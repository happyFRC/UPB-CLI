/* engine.c -- 面向用户的单星演化引擎的实现
 *
 * 对应原 Python 项目中的 engine.py。
 */

#include "engine.h"
#include "massloss.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

/* ------------------------------------------------------------------ 阶段名 */

/* 英文阶段名，对应 Python 的 _PHASE_LABELS */
static const char *PHASE_LABELS_EN[16] = {
    "MS (convective)",
    "Main sequence",
    "Hertzsprung gap",
    "First giant branch",
    "Core helium burning",
    "Early AGB",
    "Thermally pulsing AGB",
    "Naked helium star MS",
    "Naked helium star HG",
    "Naked helium star GB",
    "Helium white dwarf",
    "Carbon/oxygen white dwarf",
    "Oxygen/neon white dwarf",
    "Neutron star",
    "Black hole",
    "Massless remnant"
};

/* 中文阶段名，对应 main.py 的 PHASE_CN */
static const char *PHASE_LABELS_CN[16] = {
    "主序(全对流)",
    "主序",
    "赫氏空隙",
    "第一巨星支",
    "核心氦燃烧",
    "早期AGB",
    "热脉动AGB",
    "氦星主序",
    "氦星赫氏空隙",
    "氦星巨星支",
    "氦白矮星",
    "碳氧白矮星",
    "氧氖白矮星",
    "中子星",
    "黑洞",
    "无质量残骸"
};

const char *phase_name_en(int phase)
{
    if (phase < 0 || phase > 15) return "unknown";
    return PHASE_LABELS_EN[phase];
}

const char *phase_name_cn(int phase)
{
    if (phase < 0 || phase > 15) return "未知";
    return PHASE_LABELS_CN[phase];
}

/* 残骸阶段（不再是核燃烧恒星） */
static int is_remnant_phase(int phase)
{
    return (phase >= 10 && phase <= 15);
}

/* ------------------------------------------------------------------ 辅助 */

static double _spectral_boundary_to_value(SpectralType sp)
{
    switch (sp) {
        case SPECTRAL_O: return 30000.0;
        case SPECTRAL_B: return 10000.0;
        case SPECTRAL_A: return  7500.0;
        case SPECTRAL_F: return  6000.0;
        case SPECTRAL_G: return  5200.0;
        case SPECTRAL_K: return  3700.0;
        default:         return     0.0;
    }
}

/* 与 phases.py / engine.py 里 _spectral_type 相同的判断 */
static void _spectral_type_str(double t_eff, char out[2])
{
    if      (t_eff >= 30000.0) { out[0] = 'O'; }
    else if (t_eff >= 10000.0) { out[0] = 'B'; }
    else if (t_eff >=  7500.0) { out[0] = 'A'; }
    else if (t_eff >=  6000.0) { out[0] = 'F'; }
    else if (t_eff >=  5200.0) { out[0] = 'G'; }
    else if (t_eff >=  3700.0) { out[0] = 'K'; }
    else                       { out[0] = 'M'; }
    out[1] = '\0';
}

/* 近似中心温度 (K)。仅为诊断，不影响演化。
 *
 * 病毒定理估计 T_c ~ G M mu m_u / (k R)，无量纲系数校准到今日太阳
 * 为 1.57e7 K。退化残骸不适用病毒定理，返回 0。
 */
static double _core_temperature(double mass, double radius,
                                double x, double y, double z, int phase)
{
    double denom, mu;
    if (phase >= 7 || mass <= 0.0 || radius <= 0.0) return 0.0;
    denom = 2.0 * x + 0.75 * y + 0.5 * z;
    mu = (denom > 0.0) ? (1.0 / denom) : 0.6;
    return 1.57e7 * (mass / 1.0) * (mu / 0.62) * (1.0 / radius);
}

/* 把 TrackPoint 转成 StellarState */
static void _to_state(const TrackPoint *pt, StellarState *st)
{
    double t_eff;
    memset(st, 0, sizeof(*st));
    t_eff = pt->temperature;
    if (t_eff <= 0.0) t_eff = temperature_from_lr(pt->luminosity, pt->radius);
    st->time        = pt->time;
    st->luminosity  = pt->luminosity;
    st->radius      = pt->radius;
    st->temperature = t_eff;
    _spectral_type_str(t_eff, st->spectral_type);
    st->core_h_fraction  = 0.0;
    st->core_he_fraction = 0.0;
    st->core_temperature = 0.0;
    st->mass             = pt->mass;
    st->core_mass        = pt->core_mass;
    st->phase            = pt->phase;
    strncpy(st->phase_name, phase_name_cn(pt->phase),
            sizeof(st->phase_name) - 1);
    st->phase_name[sizeof(st->phase_name) - 1] = '\0';
    st->envelope_mass    = pt->menv;
    st->envelope_radius  = pt->renv;
}

/* 校验 (mass, metallicity)。返回 NULL 表示通过，否则返回错误信息。 */
static const char *_validate(double mass, double metallicity)
{
    static char msg[128];
    if (!isfinite(mass) || !isfinite(metallicity)) {
        return "mass and metallicity must be finite numbers";
    }
    if (mass < ENG_MASS_MIN || mass > ENG_MASS_MAX) {
        snprintf(msg, sizeof(msg),
                 "mass %g is outside the supported range [%g, %g] Msun",
                 mass, ENG_MASS_MIN, ENG_MASS_MAX);
        return msg;
    }
    if (metallicity < ENG_Z_MIN || metallicity > ENG_Z_MAX) {
        snprintf(msg, sizeof(msg),
                 "metallicity %g is outside the supported range [%g, %g]",
                 metallicity, ENG_Z_MIN, ENG_Z_MAX);
        return msg;
    }
    return NULL;
}

/* --------------------------------------------------------------- 生命周期 */

double engine_get_main_sequence_lifetime(double mass, double metallicity)
{
    Sequence sp;
    if (mass < ENG_MASS_MIN || mass > ENG_MASS_MAX) return 0.0;
    if (metallicity < ENG_Z_MIN || metallicity > ENG_Z_MAX) return 0.0;
    if (sequence_init(&sp, mass, metallicity) != 0) return 0.0;
    return sp.tm;
}

double engine_get_bgb_lifetime(double mass, double metallicity)
{
    Sequence sp;
    if (mass < ENG_MASS_MIN || mass > ENG_MASS_MAX) return 0.0;
    if (metallicity < ENG_Z_MIN || metallicity > ENG_Z_MAX) return 0.0;
    if (sequence_init(&sp, mass, metallicity) != 0) return 0.0;
    return sp.tscls[1];
}

int engine_get_phase_times(double mass, double metallicity, PhaseTimes *out)
{
    Sequence sp;
    if (out == NULL) return -1;
    memset(out, 0, sizeof(*out));
    if (mass < ENG_MASS_MIN || mass > ENG_MASS_MAX) return -1;
    if (metallicity < ENG_Z_MIN || metallicity > ENG_Z_MAX) return -1;
    if (sequence_init(&sp, mass, metallicity) != 0) return -1;

    out->t_ms       = sp.tm;
    out->t_bgb      = sp.tscls[1];
    out->t_he_i     = sp.tscls[2];
    out->t_he       = sp.tscls[3];
    out->t_bagb     = sp.t_bagb;
    out->t_tp       = sp.t_tp;
    out->t_nuclear  = sp.tn;

    out->l_zams     = sp.lums[1];
    out->l_tms      = sp.lums[2];
    out->l_bgb      = sp.lums[3];
    out->l_he_i     = sp.lums[4];
    out->l_min_he   = sp.lums[5];
    out->l_bagb     = sp.lums[7];

    out->r_zams     = seq_rzams(&sp, sp.mass);
    out->r_tms      = seq_r_tms(&sp, sp.mass);

    out->m_hook     = sp.mhook;
    out->m_he_f     = sp.mhef;
    out->m_fgb      = sp.mfgb;
    return 0;
}

/* ------------------------------------------------------------------ evolve */

/* 追加 StellarState 到结果。 */
static int result_append(EvolutionResult *r, const StellarState *st)
{
    if (r->n_track >= r->cap_track) {
        int newcap = (r->cap_track == 0) ? 1024 : r->cap_track * 2;
        StellarState *np = (StellarState *)realloc(r->evolution_track,
                                                   newcap * sizeof(StellarState));
        if (np == NULL) return -1;
        r->evolution_track = np;
        r->cap_track = newcap;
    }
    r->evolution_track[r->n_track++] = *st;
    return 0;
}

/* 时钟：返回毫秒 */
static double now_ms(void)
{
    return (double)clock() * 1000.0 / (double)CLOCKS_PER_SEC;
}

int engine_evolve(double mass, double metallicity,
                  double end_time, double step_size,
                  int track_history, int mass_loss,
                  EvolutionResult *result)
{
    EvolveOptions opt;
    Track tr;
    const char *problem;
    double t0;
    double x, y, z;
    int i, rc;

    if (result == NULL) return -1;
    memset(result, 0, sizeof(*result));

    result->mass = mass;
    result->metallicity = metallicity;
    result->end_time = end_time;
    result->step_size = step_size;

    /* 参数校验 */
    problem = _validate(mass, metallicity);
    if (problem != NULL) {
        result->success = 0;
        snprintf(result->error_message, sizeof(result->error_message),
                 "%s", problem);
        return -1;
    }
    if (end_time < 0.0) {
        result->success = 0;
        snprintf(result->error_message, sizeof(result->error_message),
                 "end_time must not be negative");
        return -1;
    }

    t0 = now_ms();

    /* 调用底层 evolve */
    opt.end_time = end_time;
    opt.step_size = step_size;
    opt.max_steps = 100000;
    opt.record = track_history ? 1 : 0;
    opt.mass_loss = mass_loss ? 1 : 0;
    opt.stop_at_remnant = 1;

    rc = evolve(mass, metallicity, &opt, &tr);
    if (rc != 0 || !tr.success) {
        result->success = 0;
        snprintf(result->error_message, sizeof(result->error_message),
                 "%s", tr.error);
        result->computation_time_ms = now_ms() - t0;
        track_free(&tr);
        return -1;
    }

    result->total_steps = tr.steps;
    result->stagnated = tr.stagnated;
    result->t_ms   = tr.t_ms;
    result->t_bgb  = tr.t_bgb;
    result->t_he_i = tr.t_he_i;
    result->t_he   = tr.t_he;
    result->t_bagb = tr.t_bagb;
    result->t_tp   = tr.t_tp;

    /* 包层成分: Y = 0.24 + 2Z, X = 1 - Y - Z */
    y = 0.24 + 2.0 * metallicity;
    x = 1.0 - y - metallicity;
    if (x < 0.0) x = 0.0;
    z = metallicity;

    /* 轨道 */
    if (track_history) {
        for (i = 0; i < tr.n_points; i++) {
            StellarState st;
            _to_state(&tr.points[i], &st);
            st.core_h_fraction  = x;
            st.core_he_fraction = y;
            st.core_temperature = _core_temperature(st.mass, st.radius,
                                                    x, y, z, st.phase);
            if (result_append(result, &st) != 0) {
                result->success = 0;
                snprintf(result->error_message, sizeof(result->error_message),
                         "out of memory building track");
                track_free(&tr);
                return -1;
            }
        }
    }

    /* 最终状态 */
    {
        const TrackPoint *fp = track_final(&tr);
        if (fp != NULL) {
            StellarState fs;
            _to_state(fp, &fs);
            fs.core_h_fraction  = x;
            fs.core_he_fraction = y;
            fs.core_temperature = _core_temperature(fs.mass, fs.radius,
                                                    x, y, z, fs.phase);
            result->final_state = fs;

            /* 保持存储的最后一个点与最终状态一致 */
            if (result->n_track > 0) {
                result->evolution_track[result->n_track - 1] = fs;
            }
        }
    }

    result->computation_time_ms = now_ms() - t0;
    result->success = 1;
    track_free(&tr);
    return 0;
}

void evolution_result_free(EvolutionResult *r)
{
    if (r == NULL) return;
    if (r->evolution_track != NULL) {
        free(r->evolution_track);
        r->evolution_track = NULL;
    }
    r->n_track = 0;
    r->cap_track = 0;
}

/* ------------------------------------------------------------- 批量/便捷 */

int engine_evolve_batch(const double *masses, const double *metallicities,
                        const double *end_times, const double *step_sizes,
                        int n, int track_history,
                        EvolutionResult *results)
{
    int i;
    if (masses == NULL || metallicities == NULL ||
        end_times == NULL || step_sizes == NULL || results == NULL)
        return 0;
    for (i = 0; i < n; i++) {
        engine_evolve(masses[i], metallicities[i],
                      end_times[i], step_sizes[i],
                      track_history, 1, &results[i]);
    }
    return n;
}

int evolve_star(double mass, double metallicity, double end_time,
                double step_size, EvolutionResult *result)
{
    return engine_evolve(mass, metallicity, end_time, step_size,
                         0, 1, result);
}

int evolve_star_detailed(double mass, double metallicity, double end_time,
                         double step_size, EvolutionResult *result)
{
    return engine_evolve(mass, metallicity, end_time, step_size,
                         1, 1, result);
}

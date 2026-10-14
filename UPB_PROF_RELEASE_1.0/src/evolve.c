/* evolve.c -- 单星演化驱动的实现
 *
 * 对应原 Python 项目中的 sse/evolve.py。
 */

#include "evolve.h"
#include "massloss.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* -------------------------------------------------------------- Track 管理 */

void track_init(Track *tr, double mass, double metallicity)
{
    if (tr == NULL) return;
    memset(tr, 0, sizeof(*tr));
    tr->mass = mass;
    tr->metallicity = metallicity;
    tr->points = NULL;
    tr->n_points = 0;
    tr->cap_points = 0;
    tr->success = 1;
}

void track_free(Track *tr)
{
    if (tr == NULL) return;
    if (tr->points != NULL) {
        free(tr->points);
        tr->points = NULL;
    }
    tr->n_points = 0;
    tr->cap_points = 0;
}

const TrackPoint *track_final(const Track *tr)
{
    if (tr == NULL) return NULL;
    if (tr->n_points > 0) return &tr->points[tr->n_points - 1];
    return &tr->last;
}

/* 追加一个记录点到 Track。cap_points 满时按 2 倍增长。 */
static int track_append(Track *tr, const TrackPoint *pt)
{
    if (tr->n_points >= tr->cap_points) {
        int newcap = (tr->cap_points == 0) ? 256 : tr->cap_points * 2;
        TrackPoint *np = (TrackPoint *)realloc(tr->points,
                                               newcap * sizeof(TrackPoint));
        if (np == NULL) return -1;
        tr->points = np;
        tr->cap_points = newcap;
    }
    tr->points[tr->n_points++] = *pt;
    return 0;
}

/* ------------------------------------------------------------------ 边界 */

/* 阶段边界时间，按升序排列。对裸氦星 (kw 7..9) 用氦星主序寿命。 */
static int boundaries(const Sequence *sp, int kw, double mt, double *out)
{
    int n = 0;
    double ts[6];
    int i, j;

    if (kw >= 7 && kw <= 9) {
        double tm_he = (mt > 0.0) ? seq_t_hems(sp, mt) : 0.0;
        if (tm_he > 0.0) { out[n++] = tm_he; }
        return n;
    }

    ts[0] = sp->tm;
    ts[1] = sp->tscls[1];
    ts[2] = sp->tscls[2];
    ts[3] = sp->tscls[2] + sp->tscls[3];
    ts[4] = sp->tscls[13];
    ts[5] = sp->tscls[14];

    /* 去掉非正值 */
    for (i = 0; i < 6; i++) {
        if (ts[i] > 0.0) {
            out[n++] = ts[i];
        }
    }
    /* 冒泡排序 */
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - 1 - i; j++) {
            if (out[j] > out[j + 1]) {
                double t = out[j];
                out[j] = out[j + 1];
                out[j + 1] = t;
            }
        }
    }
    return n;
}

/* ------------------------------------------------------------------ evolve */

int evolve(double mass, double metallicity,
           const EvolveOptions *opts, Track *out)
{
    EvolveOptions opt;
    Sequence sp, sp_cur;
    Track tr;
    double aj, mt, mc, mass0, prev_kw;
    double end_time, base, dt;
    int kw, n_steps;
    int phase_start_aj, phase_steps, prev_phase;
    int i, nb;
    double bd[6];

    /* 默认选项 */
    if (opts != NULL) {
        opt = *opts;
    } else {
        opt.end_time = 0.0;
        opt.step_size = 0.0;
        opt.max_steps = 100000;
        opt.record = 1;
        opt.mass_loss = 1;
        opt.stop_at_remnant = 1;
    }
    if (opt.max_steps <= 0) opt.max_steps = 100000;

    track_init(&tr, mass, metallicity);

    if (sequence_init(&sp, mass, metallicity) != 0) {
        tr.success = 0;
        snprintf(tr.error, sizeof(tr.error),
                 "sequence_init failed for M=%g Z=%g", mass, metallicity);
        *out = tr;
        return -1;
    }

    tr.t_ms   = sp.tm;
    tr.t_bgb  = sp.tscls[1];
    tr.t_he_i = sp.tscls[2];
    tr.t_he   = sp.tscls[3];
    tr.t_bagb = sp.t_bagb;
    tr.t_tp   = sp.t_tp;
    tr.t_nuc  = sp.tn;
    tr.l_zams = sp.lums[1];
    tr.r_zams = seq_rzams(&sp, sp.mass);

    end_time = (opt.end_time > 0.0) ? opt.end_time : 1.0e300;

    /* 初始状态 */
    aj = 0.0;
    mt = mass;
    mc = 0.0;
    kw = (mass < sp.mhook - 0.3) ? 0 : 1;
    sp_cur = sp;
    mass0 = sp.mass_in;
    prev_kw = (double)kw;

    if (opt.step_size > 0.0) {
        base = opt.step_size;
    } else {
        base = fmax(sp.tm / 1000.0, 1.0e-9);
    }

    n_steps = 0;
    phase_start_aj = (int)aj;
    phase_steps = 0;
    prev_phase = kw;

    while (n_steps < opt.max_steps) {
        Point pt;

        /* --- 有效初始质量 M0 (论文 Section 7.1) --- */
        if (kw >= 7 && kw <= 9 && prev_kw <= 6.0) {
            sequence_init(&sp_cur, mt, metallicity);
            mass0 = sp_cur.mass_in;
        } else if (kw <= 2 && fabs(mt - mass0) > 1.0e-4 * fmax(mass0, 1.0)) {
            Sequence sp_new;
            sequence_init(&sp_new, mt, metallicity);
            if (kw <= 1) {
                if (sp_cur.tm > 0.0)
                    aj *= sp_new.tm / sp_cur.tm;
            } else if (sp_cur.tscls[1] > sp_cur.tm) {
                double frac = (aj - sp_cur.tm)
                              / (sp_cur.tscls[1] - sp_cur.tm);
                if (frac < 0.0) frac = 0.0;
                if (frac > 1.0) frac = 1.0;
                aj = sp_new.tm + frac * (sp_new.tscls[1] - sp_new.tm);
            }
            sp_cur = sp_new;
            mass0 = sp_cur.mass_in;
        }
        prev_kw = (double)kw;

        /* --- 一步 hrdiag --- */
        pt = hrdiag(&sp_cur, kw, aj, mt, mc, mass0);
        kw = pt.kw;
        mc = pt.mc;
        if (pt.kw == 15 || pt.mt > 0.0) {
            mt = pt.mt;
        }

        if (kw != prev_phase) {
            prev_phase = kw;
            phase_start_aj = (int)aj;
            phase_steps = 0;
        } else {
            phase_steps++;
        }

        /* --- 记录 --- */
        {
            TrackPoint rec;
            double temp = temperature_from_lr(pt.lum, pt.r);
            rec.time = aj;
            rec.luminosity = pt.lum;
            rec.radius = pt.r;
            rec.temperature = temp;
            rec.mass = mt;
            rec.core_mass = mc;
            rec.phase = kw;
            rec.spectral_type[0] = spectral_type_from_temperature(temp)[0];
            rec.spectral_type[1] = '\0';
            rec.core_temperature = 0.0;
            rec.menv = pt.menv;
            rec.renv = pt.renv;

            if (opt.record) {
                if (track_append(&tr, &rec) != 0) {
                    tr.success = 0;
                    snprintf(tr.error, sizeof(tr.error),
                             "out of memory appending track point");
                    *out = tr;
                    return -1;
                }
            }
            tr.last = rec;
        }

        tr.steps = n_steps + 1;

        if (aj >= end_time) break;

        /* --- 到达残骸 --- */
        if (opt.stop_at_remnant
            && (kw >= EV_REMNANT_MIN && kw <= EV_REMNANT_MAX)) {
            break;
        }
        if (opt.stop_at_remnant && mt <= 0.0) break;

        /* --- 选择下一步步长 --- */
        dt = base;

        if (opt.step_size <= 0.0) {
            /* 停滞保护 */
            if (phase_steps > 60000) {
                tr.stagnated = 1;
                break;
            }

            /* 不跨过阶段边界 */
            if (kw >= 3) {
                nb = boundaries(&sp_cur, kw, mt, bd);
                for (i = 0; i < nb; i++) {
                    if (aj < bd[i] && bd[i] < aj + dt) {
                        dt = bd[i] - aj;
                        break;
                    }
                }
            }
        }

        {
            double t_next = aj + dt;
            if (end_time != 1.0e300 && t_next > end_time) {
                t_next = end_time;
            }

            /* --- 应用质量损失 --- */
            if (opt.mass_loss && mt > 0.0) {
                double rate = mass_loss_rate(kw, pt.lum, pt.r, mt, mc,
                                             metallicity, 0.0);
                mt = mt - rate * (t_next - aj);
                if (kw >= 2 && kw <= 9 && mc > 0.0) {
                    if (mt < mc) mt = mc;
                }
                if (mt < 0.0) mt = 0.0;
            }

            aj = t_next;
        }

        n_steps++;
    }

    tr.final_mass = mt;
    *out = tr;
    return 0;
}

/* massloss.c -- 质量损失与观测换算的实现
 *
 * 对应原 Python 项目中的 sse/massloss.py。
 */

#include "massloss.h"
#include <math.h>

/* ------------------------------------------------------------------ 质量损失 */

double mass_loss_rate(int kw, double lum, double r,
                      double mt, double mc, double z, double rl)
{
    double zfac, dms, dml, mew;
    double p0, dmt;

    if (mt <= 0.0 || r <= 0.0 || lum <= 0.0)
        return 0.0;

    zfac = sqrt(z / 0.02);

    /* --- Nieuwenhuijzen & de Jager (1990)，全赫罗图 --- */
    if (lum > 4000.0) {
        double x = fmin(1.0, (lum - 4000.0) / 500.0);
        dms = 9.6e-15 * x * pow(r, 0.81) * pow(lum, 1.24) * pow(mt, 0.16);
        dms *= zfac;
    } else {
        dms = 0.0;
    }

    if (kw >= 2 && kw <= 9) {
        /* --- Reimers (1975) --- */
        dml = ML_NETA * 4.0e-13 * r * lum / mt;
        if (rl > 0.0) {
            double rr = fmin(0.5, r / rl);
            dml *= 1.0 + ML_B_WIND * pow(rr, 6.0);
        }

        /* --- Vassiliadis & Wood (1993): 长周期 AGB 脉动 --- */
        if (kw == 5 || kw == 6) {
            p0 = -2.07 - 0.9 * log10(mt) + 1.94 * log10(r);
            p0 = fmin(pow(10.0, p0), 2000.0);
            dmt = -11.4 + 0.0125 * (p0 - 100.0 * fmax(mt - 2.5, 0.0));
            dmt = fmin(pow(10.0, dmt), 1.36e-9 * lum);
            dml = fmax(dml, dmt);
        }

        if (kw > 6) {
            dms = fmax(dml, 1.0e-13 * pow(lum, 1.5));
        } else {
            dms = fmax(dml, dms);
            mew = ((mt - mc) / mt)
                  * fmin(5.0, fmax(1.2, pow(lum / ML_L0, ML_KAP)));
            if (mew < 1.0) {
                /* 小氢包层的减弱 Wolf-Rayet 星风 */
                dml = 1.0e-13 * pow(lum, 1.5) * (1.0 - mew) * zfac;
                dms = fmax(dml, dms);
            }
        }
    }

    /* 从 Msun/yr 转成 Msun/Myr */
    return dms * ML_YR_PER_MYR;
}

/* -------------------------------------------------------------- 观测换算 */

double temperature_from_lr(double lum, double r)
{
    if (lum <= 0.0 || r <= 0.0)
        return 0.0;
    return 5773.15 * pow(lum / (r * r), 0.25);
}

typedef struct {
    double limit;
    const char *name;
} SpectralBound;

static const SpectralBound SPECTRAL_BOUNDS[] = {
    { 30000.0, "O" },
    { 10000.0, "B" },
    {  7500.0, "A" },
    {  6000.0, "F" },
    {  5200.0, "G" },
    {  3700.0, "K" },
};

const char *spectral_type_from_temperature(double t_eff)
{
    int i;
    int n = (int)(sizeof(SPECTRAL_BOUNDS) / sizeof(SPECTRAL_BOUNDS[0]));
    for (i = 0; i < n; i++) {
        if (t_eff >= SPECTRAL_BOUNDS[i].limit)
            return SPECTRAL_BOUNDS[i].name;
    }
    return "M";
}

double bolometric_correction(double t_eff)
{
    double x;
    if (t_eff <= 0.0)
        return 0.0;
    x = log10(t_eff);
    if (x >= 3.90)
        return -5.0 + 2.5 * (x - 3.90) * 10.0;
    if (x >= 3.70)
        return -0.10 + 4.0 * (x - 3.70);
    return 0.9 + 3.0 * (x - 3.60);
}

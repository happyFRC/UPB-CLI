/* phases.c -- 恒星属性随年龄变化的实现
 *
 * 对应原 Python 项目中的 sse/phases.py。
 */

#include "phases.h"
#include "massloss.h"

#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------ 小工具 */

/* base**exponent，base<=0 时返回 0 */
static double pow_safe(double base, double exponent)
{
    if (base <= 0.0) return 0.0;
    return pow(base, exponent);
}

/* Fortran 风格的实数幂 */
static double real_pow(double base, double exponent)
{
    if (base < 0.0) return -pow(-base, exponent);
    return pow(base, exponent);
}

/* --------------------------------------------------- MS 亮度/半径系数 (5.1) */

static double lalph(const Sequence *sp, double m)
{
    const ZCoeffs *a = &sp->c;
    double mcut = 2.0;
    double b_val;
    if (m >= mcut) {
        return (a->a45 + a->a46 * pow(m, a->a48))
               / (pow(m, 0.4) + a->a47 * pow(m, 1.9));
    }
    if (m <= 0.5) return a->a49;
    if (m <= 0.7) return a->a49 + ((0.3 - a->a49) / 0.2) * (m - 0.5);
    if (m <= a->a52)
        return 0.3 + ((a->a50 - 0.3) / (a->a52 - 0.7)) * (m - 0.7);
    if (m <= a->a53)
        return a->a50 + ((a->a51 - a->a50) / (a->a53 - a->a52))
               * (m - a->a52);
    b_val = (a->a45 + a->a46 * pow(mcut, a->a48))
            / (pow(mcut, 0.4) + a->a47 * pow(mcut, 1.9));
    return a->a51 + ((b_val - a->a51) / (mcut - a->a53)) * (m - a->a53);
}

static double lbet(const Sequence *sp, double m)
{
    const ZCoeffs *a = &sp->c;
    double beta = fmax(0.0, a->a54 - a->a55 * pow(m, a->a56));
    double a57 = a->a57;
    if (m > a57 && beta > 0.0) {
        double b_val = fmax(0.0, a->a54 - a->a55 * pow(a57, a->a56));
        beta = fmax(0.0, b_val - 10.0 * (m - a57) * b_val);
    }
    return beta;
}

static double lneta(const Sequence *sp, double m)
{
    const ZCoeffs *a = &sp->c;
    double neta;
    if (m <= 1.0)      neta = 10.0;
    else if (m >= 1.1) neta = 20.0;
    else               neta = 10.0 + 100.0 * (m - 1.0);
    return fmin(neta, a->eta_max);
}

static double ralph(const Sequence *sp, double m)
{
    const ZCoeffs *a = &sp->c;
    double b_val, c_val;
    if (m <= 0.5) return a->a62;
    if (m <= 0.65)
        return a->a62 + ((a->a63 - a->a62) / 0.15) * (m - 0.5);
    if (m <= a->a68)
        return a->a63 + ((a->a64 - a->a63) / (a->a68 - 0.65)) * (m - 0.65);
    if (m <= a->a66) {
        b_val = (a->a58 * pow(a->a66, a->a60))
                / (a->a59 * pow(a->a66, a->a61));
        return a->a64 + ((b_val - a->a64) / (a->a66 - a->a68))
               * (m - a->a68);
    }
    if (m <= a->a67)
        return (a->a58 * pow(m, a->a60)) / (a->a59 * pow(m, a->a61));
    c_val = (a->a58 * pow(a->a67, a->a60))
            / (a->a59 * pow(a->a67, a->a61));
    return c_val + a->a65 * (m - a->a67);
}

static double rbet(const Sequence *sp, double m)
{
    const ZCoeffs *a = &sp->c;
    double m2 = 2.0, m3 = 16.0;
    double v, b2, b3;
    if (m <= 1.0)        v = 1.06;
    else if (m <= a->a74) v = 1.06 + ((a->a72 - 1.06) / (a->a74 - 1.0))
                                    * (m - 1.0);
    else if (m <= m2) {
        b2 = (a->a69 * pow(m2, 3.5)) / (a->a70 + pow(m2, a->a71));
        v = a->a72 + ((b2 - a->a72) / (m2 - a->a74)) * (m - a->a74);
    } else if (m <= m3) {
        v = (a->a69 * pow(m, 3.5)) / (a->a70 + pow(m, a->a71));
    } else {
        b3 = (a->a69 * pow(m3, 3.5)) / (a->a70 + pow(m3, a->a71));
        v = b3 + a->a73 * (m - m3);
    }
    return v - 1.0;
}

static double rgamm(const Sequence *sp, double m)
{
    const ZCoeffs *a = &sp->c;
    double m1 = 1.0;
    double b1, v, c_val;
    if (m > a->a75 + 0.1) return 0.0;
    b1 = fmax(0.0, a->a76 + a->a77 * pow(fabs(m1 - a->a78), a->a79));
    if (m <= m1) {
        v = a->a76 + a->a77 * pow(fabs(m - a->a78), a->a79);
    } else if (m <= a->a75) {
        v = b1 + (a->a80 - b1)
            * pow((m - m1) / (a->a75 - m1), a->a81);
    } else {
        c_val = a->a80;
        if (a->a75 <= m1) c_val = b1;
        v = c_val - 10.0 * (m - a->a75) * c_val;
    }
    return fmax(v, 0.0);
}

static double lhook(const Sequence *sp, double m)
{
    const ZCoeffs *a = &sp->c;
    double mhook = sp->mhook;
    double b_val;
    if (m <= mhook) return 0.0;
    if (m >= a->a33)
        return fmin(a->a34 / pow(m, a->a35), a->a36 / pow(m, a->a37));
    b_val = fmin(a->a34 / pow(a->a33, a->a35),
                 a->a36 / pow(a->a33, a->a37));
    return b_val * pow((m - mhook) / (a->a33 - mhook), 0.4);
}

static double rhook(const Sequence *sp, double m)
{
    const ZCoeffs *a = &sp->c;
    double mhook = sp->mhook;
    double m2, b2;
    if (m <= mhook) return 0.0;
    if (m <= a->a42)
        return a->a43 * sqrt((m - mhook) / (a->a42 - mhook));
    if (m <= 2.0) {
        m2 = 2.0;
        b2 = (a->a38 + a->a39 * pow(m2, 3.5))
             / (a->a40 * pow(m2, 3.0) + pow(m2, a->a41)) - 1.0;
        return a->a43 + (b2 - a->a43)
               * pow((m - a->a42) / (m2 - a->a42), a->a44);
    }
    return (a->a38 + a->a39 * pow(m, 3.5))
           / (a->a40 * pow(m, 3.0) + pow(m, a->a41)) - 1.0;
}

static double mctms(const Sequence *sp, double m)
{
    double m525 = pow(m, 5.25);
    (void)sp;
    return (1.586 + m525) / (2.434 + 1.02 * m525);
}

/* ------------------------------------------------- 小包层扰动 (eq. 83-84) */

static double lpert(double m, double mew)
{
    double b = 0.002 * fmax(1.0, 2.5 / m);
    double c = 3.0;
    return ((1.0 + pow(b, c)) * pow(mew / b, c))
           / (1.0 + pow(mew / b, c));
}

static double rpert(double m, double mew, double r, double rc)
{
    double a = 0.1;
    double b = 0.006 * fmax(1.0, 2.5 / m);
    double c = 3.0;
    double q = log(r / rc);
    return ((1.0 + pow(b, c)) * pow(mew / b, c) * pow(mew, a / q))
           / (1.0 + pow(mew / b, c));
}

/* ------------------------------------------------------------ 白矮星半径 */

static double wd_radius(double mt)
{
    double arg;
    if (mt < 5.0e-6) return 0.009;
    if (mt < 5.0e-4) return 0.09;
    arg = pow(MCH_ / mt, 2.0 / 3.0) - pow(mt / MCH_, 2.0 / 3.0);
    return fmax(PH_RNS, 0.0115 * sqrt(fmax(1.48204e-6, arg)));
}

/* ------------------------------------------------------------ 裸氦星处理 */

/* 返回 r, lum, kw, mc, rc, menv, renv, k2。
 * 注意: 本函数不会把 mt 传出去（不需要——所有 kw>=10 的路径都直接
 * 在 hrdiag 里构造输出并返回）。 */
static void helium_star(const Sequence *sp, int kw, double aj, double mt,
                        double mc,
                        double *r_out, double *lum_out, int *kw_out,
                        double *mc_out, double *rc_out,
                        double *menv_out, double *renv_out, double *k2_out)
{
    double tscls[16], lums[9], gb[10];
    double tm, rzams, rx, rc, menv, renv, k2;
    double mtc, mcmax;
    int i;

    for (i = 0; i < 16; i++) tscls[i] = sp->tscls[i];
    for (i = 0; i < 9;  i++) lums[i]  = sp->lums[i];
    for (i = 0; i < 10; i++) gb[i]    = sp->gb[i];

    tm = seq_t_hems(sp, mt);
    tscls[1] = tm;
    lums[1]  = seq_l_zhe(sp, mt);
    {
        double am = fmax(0.0, 0.85 - 0.08 * mt);
        lums[2] = lums[1] * (1.0 + 0.45 + am);
    }
    gb[8] = 8.0e-5;
    gb[3] = 4.1e4;
    gb[4] = 5.5e4 / (1.0 + 0.4 * pow(mt, 4.0));
    gb[5] = 5.0;
    gb[6] = 3.0;
    gb[7] = pow(gb[3] / gb[4], 1.0 / (gb[5] - gb[6]));
    lums[6] = gb[4] * pow(gb[7], gb[5]);

#define HE_MC_FROM_L(lum_) ( \
    ((lum_) <= gb[4] * pow(pow(gb[3] / gb[4], 1.0 / (gb[5] - gb[6])), gb[5])) \
        ? (((lum_) > 0.0) ? pow((lum_) / gb[4], 1.0 / gb[5]) : 0.0) \
        : (((lum_) > 0.0) ? pow((lum_) / gb[3], 1.0 / gb[6]) : 0.0))

    {
        double mc1 = HE_MC_FROM_L(lums[2]);
        tscls[4] = tm + 1.0 / ((gb[5] - 1.0) * gb[8] * gb[4])
                   * pow(mc1, 1.0 - gb[5]);
        tscls[6] = tscls[4] - (tscls[4] - tm)
                   * pow(gb[7] / mc1, 1.0 - gb[5]);
        tscls[5] = tscls[6] + 1.0 / ((gb[6] - 1.0) * gb[8] * gb[3])
                   * pow(gb[7], 1.0 - gb[6]);
    }

    mtc = fmin(mt, 1.45 * mt - 0.31);
    if (mtc <= 0.0) mtc = mt;
    mcmax = fmin(mtc, fmax(MCH_, 0.773 * mt - 0.35));
    if (mcmax <= gb[7])
        tscls[14] = tscls[4] - 1.0 / ((gb[5] - 1.0) * gb[8] * gb[4])
                    * pow(mcmax, 1.0 - gb[5]);
    else
        tscls[14] = tscls[5] - 1.0 / ((gb[6] - 1.0) * gb[8] * gb[3])
                    * pow(mcmax, 1.0 - gb[6]);
    tscls[14] = fmax(tscls[14], tm);

    rzams = seq_r_zhe(sp, mt);
    rx = rzams;
    rc = 0.0;
    menv = 1.0e-10;
    renv = 1.0e-10;
    k2 = 0.21;

    if (aj < tm) {
        double tau = aj / tm;
        double am, l_tms_he, mc_end;
        *kw_out = 7;
        am = fmax(0.0, 0.85 - 0.08 * sp->mass);
        *lum_out = lums[1] * (1.0 + 0.45 * tau + am * tau * tau);
        am = fmax(0.0, 0.4 - 0.22 * log10(mt));
        *r_out = rx * (1.0 + am * (tau - pow(tau, 6.0)));
        l_tms_he = lums[2];
        if (l_tms_he <= gb[4] * pow(pow(gb[3] / gb[4],
                                        1.0 / (gb[5] - gb[6])), gb[5]))
            mc_end = (l_tms_he > 0.0) ? pow(l_tms_he / gb[4], 1.0 / gb[5]) : 0.0;
        else
            mc_end = (l_tms_he > 0.0) ? pow(l_tms_he / gb[3], 1.0 / gb[6]) : 0.0;
        mc = fmin(mcmax, mc_end);
        if (mt < sp->zpars[10]) *kw_out = 10;
    } else {
        double lum, rg;
        *kw_out = 8;
        lum = seq_l_from_t(sp, aj, gb[8], tscls[4], tscls[5], tscls[6], mt);
        *r_out = seq_r_he_hg(sp, mt, lum, rx, lums[2]);
        rg = seq_r_he_gb(sp, lum);
        mc = HE_MC_FROM_L(lum);
        mtc = fmin(mt, 1.45 * mt - 0.31);
        if (mtc <= 0.0) mtc = mt;
        mcmax = fmin(mtc, fmax(MCH_, 0.773 * sp->mass - 0.35));
        mc = fmax(mc, mcmax);
        if (mc >= mcmax) {
            double mcx;
            aj = 0.0;
            mc = fmax(mcmax, 1.0e-10);
            if (mc < MCH_) {
                if (sp->mass < 1.6) {
                    mt = fmax(mc, (mc + 0.31) / 1.45);
                    *kw_out = 11;
                } else {
                    mt = mc;
                    *kw_out = 12;
                }
            } else {
                if (mc < 2.5)      mcx = 0.161767 * mc + 1.067055;
                else               mcx = 0.314154 * mc + 0.686088;
                if (mc <= 5.0)      mt = mcx;
                else if (mc < 7.6)  mt = mcx + (mc - 5.0) * (mt - mcx) / 2.6;
                mc = mt;
                *kw_out = (mt <= PH_MXNS) ? 13 : 14;
            }
            *lum_out = lum;
            *mc_out = mc;
            *rc_out = *r_out;
            *menv_out = 1.0e-10;
            *renv_out = 1.0e-10;
            *k2_out = 0.21;
            return;
        }
        if (*r_out >= rg) {
            *kw_out = 9;
            *r_out = rg;
        }
    }
#undef HE_MC_FROM_L

    if (*kw_out >= 10) {
        double mcx;
        if (*kw_out >= 10 && *kw_out <= 12) {
            mc = mt;
            *lum_out = 635.0 * mt * pow(sp->z, 0.4)
                       / pow((*kw_out == 10) ? PH_AHE
                              : ((*kw_out == 11) ? PH_ACO : PH_AONE),
                             1.4);
            *r_out = wd_radius(mt);
            *mc_out = mc;
            *rc_out = *r_out;
            *menv_out = 1.0e-10;
            *renv_out = 1.0e-10;
            *k2_out = 0.21;
            return;
        }
        if (*kw_out == 13) {
            mc = mt;
            *lum_out = 0.02 * pow(mt, 2.0 / 3.0) / pow(fmax(aj, 0.1), 2.0);
            *r_out = PH_RNS;
            *mc_out = mc;
            *rc_out = *r_out;
            *menv_out = 1.0e-10;
            *renv_out = 1.0e-10;
            *k2_out = 0.21;
            return;
        }
        if (*kw_out == 14) {
            mc = mt;
            *lum_out = 1.0e-10;
            *r_out = 4.24e-6 * mt;
            *mc_out = mc;
            *rc_out = *r_out;
            *menv_out = 1.0e-10;
            *renv_out = 1.0e-10;
            *k2_out = 0.21;
            return;
        }
        mcx = 0.0; (void)mcx;
    }

    *mc_out = mc;
    *rc_out = rc;
    *menv_out = menv;
    *renv_out = renv;
    *k2_out = k2;
}

/* ------------------------------------------------------------------ mrenv */

void mrenv(int kw, double mass, double mt, double mc,
           double lum, double rad, double rc,
           double aj, double tm, double ltms, double lbgb, double lhei,
           double rzams, double rtms, double rg,
           double *menv_out, double *renv_out, double *k2_out)
{
    double logm, a, c, d, e;
    double k2z, k2bgb, k2g, menvg, renvg;
    double menv, renv, k2e;
    double teff = 0.0, tebgb = 0.0;   /* kw<=2 分支用，外层需要 */

    if (mass <= 0.0 || mt <= 0.0 || rad <= 0.0 || lum <= 0.0) {
        *menv_out = 1.0e-10;
        *renv_out = 1.0e-10;
        *k2_out   = 0.21;
        return;
    }

    logm = log10(mass);
    a = fmin(0.81, fmax(0.68, 0.68 + 0.4 * logm));
    c = fmax(-2.5, fmin(-1.5, -2.5 + 5.0 * logm));
    d = -0.1;
    e = 0.025;

    k2z = fmin(0.21, fmax(0.09 - 0.27 * logm, 0.037 + 0.033 * logm));
    if (logm > 1.3)
        k2z = k2z - 0.055 * pow(logm - 1.3, 2.0);
    k2bgb = fmin(0.15, fmin(0.147 + 0.03 * logm, 0.162 - 0.04 * logm));

    if (kw >= 3 && kw <= 6) {
        double logmt = log10(mt);
        double f = 0.208 + 0.125 * logmt - 0.035 * logmt * logmt;
        double b = 1.0e4 * pow(mt, 1.5) / (1.0 + 0.1 * pow(mt, 1.5));
        double x = pow((lum - lbgb) / b, 2.0);
        double y = (f - 0.033 * log10(lbgb)) / k2bgb - 1.0;
        k2g = (f - 0.033 * log10(lum) + 0.4 * x)
              / (1.0 + y * (lbgb / lum) + x);
    } else if (kw == 9) {
        double b = 3.0e4 * pow(mt, 1.5);
        double x = pow(fmax(0.0, lum / b - 0.5), 2.0);
        k2g = (k2bgb + 0.4 * x) / (1.0 + 0.4 * x);
    } else {
        k2g = k2bgb;
    }

    if (kw <= 2) {
        menvg = 0.5; renvg = 0.65;
    } else if (kw == 3 && lum < 3.0 * lbgb) {
        double x = fmin(3.0, lhei / lbgb);
        double tau = fmax(0.0, fmin(1.0, (x - lum / lbgb) / (x - 1.0)));
        menvg = 1.0 - 0.5 * tau * tau;
        renvg = 1.0 - 0.35 * tau * tau;
    } else {
        menvg = 1.0; renvg = 1.0;
    }

    if (rad < rg) {
        if (kw <= 6) {
            double rzr = fmax(1.0, rad / rzams);
            k2e = (k2z - e) * pow(rzr, c) + e * pow(rzr, d);
        } else if (kw == 7) {
            double tau = (tm != 0.0) ? aj / tm : 0.0;
            k2e = 0.08 - 0.03 * tau;
        } else {
            k2e = 0.08 * rzams / rad;
        }

        {
            double tauenv;
            if (kw <= 2) {
                teff = pow(lum, 0.25) / sqrt(rad);
                tebgb = pow(lbgb, 0.25) / sqrt(rg);
                tauenv = fmax(0.0, fmin(1.0, (tebgb / teff - a) / (1.0 - a)));
            } else {
                tauenv = fmax(0.0, fmin(1.0,
                              (sqrt(rad / rg) - a) / (1.0 - a)));
            }

            if (tauenv > 0.0) {
                menv = menvg * pow(tauenv, 5.0);
                renv = renvg * pow(tauenv, 1.25);
                if (kw <= 1) {
                    double x = fmax(0.0, fmin(1.0, (0.1 - logm) / 0.55));
                    double menvz = 0.18 * x + 0.82 * pow(x, 5.0);
                    double renvz = 0.4 * pow(x, 0.25) + 0.6 * pow(x, 10.0);
                    double y = 2.0 + 8.0 * x;
                    double tetms = pow(ltms, 0.25) / sqrt(rtms);
                    double tautms = fmax(0.0, fmin(1.0,
                                     (tebgb / tetms - a) / (1.0 - a)));
                    double menvt = menvg * pow(tautms, 5.0);
                    double renvt = renvg * pow(tautms, 1.25);
                    double tau = (tm != 0.0) ? aj / tm : 0.0;
                    if (tautms > 0.0 && menvt > 0.0 && renvt > 0.0) {
                        menv = menvz + pow(tau, y) * menv
                               * (menvt - menvz) / menvt;
                        renv = renvz + pow(tau, y) * renv
                               * (renvt - renvz) / renvt;
                    } else {
                        menv = 0.0;
                        renv = 0.0;
                    }
                    k2e = k2e + pow(tau, y) * pow(tauenv, 3.0)
                          * (k2g - k2e);
                } else {
                    k2e = k2e + pow(tauenv, 3.0) * (k2g - k2e);
                }
            } else {
                menv = 0.0;
                renv = 0.0;
            }
        }
    } else {
        menv = menvg;
        renv = renvg;
        k2e = k2g;
    }

    menv = menv * (mt - mc);
    renv = renv * (rad - rc);
    *menv_out = fmax(menv, 1.0e-10);
    *renv_out = fmax(renv, 1.0e-10);
    *k2_out   = k2e;
}

/* ------------------------------------------------------------------ 阈值 */

double mA_ms_threshold(const Sequence *sp)
{
    return sp->mhook - 0.3;
}

/* ------------------------------------------------------------------ hrdiag */

Point hrdiag(const Sequence *sp, int kw, double aj, double mt,
             double mc, double mass)
{
    Point out;
    double mass0, mass_f, mt0;
    double mc_prev, mcx_prev;
    double r, lum, rc, menv, renv, k2;
    double mtc;
    double *tscls, *lums, *gb;
    const double *zp;
    double tscls_local[16], lums_local[9], gb_local[10];
    double rg = 0.0;            /* 让 mrenv 调用处可见 */
    int i;

    /* 清零 */
    memset(&out, 0, sizeof(out));

    mass0 = (mass < 0.0) ? sp->mass_in : mass;
    mass_f = fmin(mass0, 100.0);
    mt0 = mt;
    (void)mt0;
    mt = fmin(mt, 100.0);
    mc_prev = mc;
    mcx_prev = (kw == 6) ? fmax(mc, sp->mc_du) : mc;

    for (i = 0; i < 16; i++) tscls_local[i] = sp->tscls[i];
    for (i = 0; i < 9;  i++) lums_local[i]  = sp->lums[i];
    for (i = 0; i < 10; i++) gb_local[i]    = sp->gb[i];
    tscls = tscls_local;
    lums  = lums_local;
    gb    = gb_local;
    zp    = sp->zpars;

    out.mass = mass_f;
    out.mt   = mt;

    r = lum = 0.0;
    rc = 0.0;
    menv = renv = 1.0e-10;
    k2 = 0.21;
    mtc = fmin(mt, 1.45 * mt - 0.31);

    if (kw <= 6) {
        double tbagb_local = tscls[2] + tscls[3];
        double thg = tscls[1] - sp->tm;
        double rzams = seq_rzams(sp, mass_f);
        double rtms = seq_r_tms(sp, mass_f);
        rg = seq_r_gb(sp, lums[3], mt);

        if (aj < tscls[1]) {
            /* --------------------------------------- MS 或 HG */
            if (aj < sp->tm) {
                /* MS */
                double tau = aj / sp->tm;
                double thook = seq_t_hook(sp, mass_f) * tscls[1];
                double zeta = 0.01;
                double tau1 = fmin(1.0, aj / thook);
                double tau2 = fmax(0.0, fmin(1.0,
                                (aj - (1.0 - zeta) * thook) / (zeta * thook)));
                double dell = lhook(sp, mass_f);
                double dtau = tau1 * tau1 - tau2 * tau2;
                double alpha = lalph(sp, mass_f);
                double beta = lbet(sp, mass_f);
                double neta = lneta(sp, mass_f);
                double lx = log10(lums[2] / lums[1]);
                double xx;
                double delr, ralpha, rbeta, rgamma, rx;
                mc = 0.0;
                if (tau > PH_TAUMIN)
                    xx = alpha * tau + beta * pow(tau, neta)
                         + (lx - alpha - beta) * tau * tau - dell * dtau;
                else
                    xx = alpha * tau + (lx - alpha) * tau * tau - dell * dtau;
                lum = lums[1] * pow(10.0, xx);

                delr = rhook(sp, mass_f);
                dtau = tau1 * tau1 * tau1 - tau2 * tau2 * tau2;
                ralpha = ralph(sp, mass_f);
                rbeta  = rbet(sp, mass_f);
                rgamma = rgamm(sp, mass_f);
                rx = log10(rtms / rzams);
                if (tau > PH_TAUMIN)
                    xx = ralpha * tau + rbeta * pow(tau, 10.0)
                         + rgamma * pow(tau, 40.0)
                         + (rx - ralpha - rbeta - rgamma)
                           * tau * tau * tau - delr * dtau;
                else
                    xx = ralpha * tau + (rx - ralpha) * tau * tau * tau
                         - delr * dtau;
                r = rzams * pow(10.0, xx);

                if (mass_f < mA_ms_threshold(sp)) {
                    double rdgen = 0.0258 * pow(1.0 + sp->x_h, 5.0 / 3.0)
                                   * pow(mass_f, -1.0 / 3.0);
                    kw = 0;
                    r = fmax(rdgen, r);
                } else {
                    kw = 1;
                }
            } else {
                /* ----------------------------------- 赫氏空隙 */
                double mcx = mc;
                double neta, tau;
                if (mass_f <= sp->mhef)
                    mc = seq_mc_from_l(sp, lums[3], mass_f);
                else if (mass_f <= sp->mfgb)
                    mc = seq_mc_he_if(sp, mass_f, zp[9]);
                else
                    mc = seq_mc_he_if(sp, mass_f, zp[10]);
                neta = mctms(sp, mass_f);
                tau = (aj - sp->tm) / thg;
                mc = ((1.0 - tau) * neta + tau) * mc;
                mc = fmax(mc, mcx);

                if (mc >= mt) {
                    aj = 0.0;
                    if (mass_f > sp->mhef) {
                        /* 零龄裸氦星 */
                        Sequence sp2;
                        mc = 0.0;
                        sequence_init(&sp2, mt, sp->z);
                        helium_star(&sp2, 7, 0.0, mt, mc,
                                    &r, &lum, &kw, &mc, &rc,
                                    &menv, &renv, &k2);
                        out.mass = mt;
                        out.mt = mt;
                        out.r = r; out.lum = lum; out.kw = kw;
                        out.mc = mc; out.rc = rc;
                        out.menv = menv; out.renv = renv; out.k2 = k2;
                        return out;
                    } else {
                        /* 零龄氦白矮星 */
                        mc = mt;
                        out.mass = mt;
                        kw = 10;
                    }
                } else {
                    double rx, tau;
                    lum = lums[2] * pow(lums[3] / lums[2], tau);
                    if (mass_f <= sp->mfgb) {
                        rx = rg;
                    } else {
                        double rmin = seq_r_min_he(sp, mass_f);
                        double ry = seq_ragbf(sp, lums[4], mt);
                        double texp, tau2;
                        rx = fmin(rmin, ry);
                        if (mass_f <= MLP_) {
                            texp = log(mass_f / MLP_) / log(sp->mfgb / MLP_);
                            rx = rg;
                            rx = rmin * pow(rx / rmin, texp);
                        }
                        tau2 = seq_tau_bl(sp, mass_f);
                        if (tau2 < TINY_) rx = ry;
                    }
                    r = rtms * pow(rx / rtms, tau);
                    kw = 2;
                }
            }
        } else if (aj < tscls[2]) {
            /* ----------------------------------- 第一巨星支 */
            kw = 3;
            lum = seq_l_from_t(sp, aj, gb[1], tscls[4], tscls[5],
                               tscls[6], mass_f);
            if (mass_f <= sp->mhef) {
                mc = seq_mc_from_l(sp, lum, mass_f);
            } else {
                double tau = (aj - tscls[1]) / (tscls[2] - tscls[1]);
                double mcx = seq_mc_he_if(sp, mass_f, zp[9]);
                double mcy = seq_mc_he_if(sp, mass_f, zp[10]);
                mc = mcx + (mcy - mcx) * tau;
            }
            r = seq_r_gb(sp, lum, mt);
            rg = r;
            if (mc >= mt) {
                aj = 0.0;
                if (mass_f > sp->mhef) {
                    Sequence sp2;
                    mc = 0.0;
                    out.mass = mt;
                    sequence_init(&sp2, mt, sp->z);
                    helium_star(&sp2, 7, 0.0, mt, mc,
                                &r, &lum, &kw, &mc, &rc,
                                &menv, &renv, &k2);
                    out.r = r; out.lum = lum; out.kw = kw;
                    out.mc = mc; out.rc = rc;
                    out.menv = menv; out.renv = renv; out.k2 = k2;
                    out.mt = mt;
                    return out;
                } else {
                    mc = mt;
                    out.mass = mt;
                    kw = 10;
                }
            }
        } else if (aj < tbagb_local) {
            /* ----------------------------------- 核心氦燃烧 */
            double mcx, tau, texp, taul, tauh, tau2, lx, ly, rx, ry, rmin;
            if (kw == 3 && mass_f <= sp->mhef) {
                /* LM 星: 氦闪瞬时，重置时钟 */
                (void)mt;
                aj = tscls[2];
            }
            if (mass_f <= sp->mhef)
                mcx = seq_mc_from_l(sp, lums[4], mass_f);
            else
                mcx = seq_mc_he_if(sp, mass_f, zp[10]);
            tau = (aj - tscls[2]) / tscls[3];
            mc = mcx + (seq_mc_bagb(sp, mass_f) - mcx) * tau;

            if (mass_f <= sp->mhef) {
                /* LM: 整个 CHeB 都是蓝的 */
                lx = lums[5];
                ly = lums[7];
                rx = seq_r_zahb(sp, mc, mt);
                rg = seq_r_gb(sp, lx, mt);
                rmin = rg * pow(sp->rmin_ratio, mass_f / sp->mhef);
                texp = fmin(fmax(0.4, rmin / rx), 2.5);
                ry = seq_ragbf(sp, ly, mt);
                if (rmin < rx)
                    taul = pow(log(rx / rmin), 1.0 / 3.0);
                else {
                    rmin = rx;
                    taul = 0.0;
                }
                tauh = pow(log(ry / rmin), 1.0 / 3.0);
                tau2 = taul * (tau - 1.0) + tauh * tau;
                r = rmin * exp(pow(fabs(tau2), 3.0));
                rg = rg + tau * (ry - rg);
                lum = lx * pow(ly / lx, pow(tau, texp));
            } else if (mass_f > sp->mfgb) {
                /* HM: 先蓝后红 (超) 巨星 */
                tau2 = seq_tau_bl(sp, mass_f);
                {
                    double tloop = tscls[2] + tau2 * tscls[3];
                    rmin = seq_r_min_he(sp, mass_f);
                    rg = seq_r_gb(sp, lums[4], mt);
                    rx = seq_ragbf(sp, lums[4], mt);
                    rmin = fmin(rmin, rx);
                    if (mass_f <= MLP_) {
                        texp = log(mass_f / MLP_) / log(sp->mfgb / MLP_);
                        rx = rg;
                        rx = rmin * pow(rx / rmin, texp);
                    } else {
                        rx = rmin;
                    }
                    texp = fmin(fmax(0.4, rmin / rx), 2.5);
                    lum = lums[4] * pow(lums[7] / lums[4], pow(tau, texp));
                    if (aj < tloop) {
                        double ly2 = lums[4] * pow(lums[7] / lums[4],
                                                    pow(tau2, texp));
                        double ry2 = seq_ragbf(sp, ly2, mt);
                        double taul2 = 0.0, tauh2 = 0.0;
                        if (fabs(rmin - rx) > TINY_)
                            taul2 = pow(log(rx / rmin), 1.0 / 3.0);
                        if (ry2 > rmin)
                            tauh2 = pow(log(ry2 / rmin), 1.0 / 3.0);
                        tau = (aj - tscls[2]) / (tau2 * tscls[3]);
                        tau2 = taul2 * (tau - 1.0) + tauh2 * tau;
                        r = rmin * exp(pow(fabs(tau2), 3.0));
                        rg = rg + tau * (ry2 - rg);
                    } else {
                        r = seq_ragbf(sp, lum, mt);
                        rg = r;
                    }
                }
            } else {
                /* IM: 先红巨星后蓝回圈 */
                tau2 = 1.0 - seq_tau_bl(sp, mass_f);
                {
                    double tloop = tscls[2] + tau2 * tscls[3];
                    if (aj < tloop) {
                        tau = (tloop - aj) / (tau2 * tscls[3]);
                        lum = lums[5] * pow(lums[4] / lums[5],
                                             pow(tau, 3.0));
                        r = seq_r_gb(sp, lum, mt);
                        rg = r;
                    } else {
                        lx = lums[5];
                        ly = lums[7];
                        rx = seq_r_gb(sp, lx, mt);
                        rmin = seq_r_min_he(sp, mt);
                        texp = fmin(fmax(0.4, rmin / rx), 2.5);
                        ry = seq_ragbf(sp, ly, mt);
                        if (rmin < rx)
                            taul = pow(log(rx / rmin), 1.0 / 3.0);
                        else {
                            rmin = rx;
                            taul = 0.0;
                        }
                        tauh = pow(log(ry / rmin), 1.0 / 3.0);
                        tau = (aj - tloop)
                              / (tscls[3] - (tloop - tscls[2]));
                        tau2 = taul * (tau - 1.0) + tauh * tau;
                        r = rmin * exp(pow(fabs(tau2), 3.0));
                        rg = rx + tau * (ry - rx);
                        lum = lx * pow(ly / lx, pow(tau, texp));
                    }
                }
            }

            if (mc >= mt) {
                /* 演化后的 MS 裸氦星 */
                double xx;
                Sequence sp2;
                kw = 7;
                xx = (aj - tscls[2]) / tscls[3];
                sequence_init(&sp2, mt, sp->z);
                out.mass = mt;
                aj = xx * seq_t_hems(&sp2, mt);
                helium_star(&sp2, kw, aj, mt, mc,
                            &r, &lum, &kw, &mc, &rc,
                            &menv, &renv, &k2);
                out.r = r; out.lum = lum; out.kw = kw;
                out.mc = mc; out.rc = rc;
                out.menv = menv; out.renv = renv; out.k2 = k2;
                out.mt = mt;
                return out;
            }
            kw = 4;
        } else {
            /* ----------------------------------- AGB */
            double mcbagb = seq_mc_bagb(sp, mass_f);
            double mcx = seq_mc_from_t(sp, tbagb_local, gb[8], tscls[7],
                                       tscls[8], tscls[9], mass_f);
            double mcmax = fmax(fmax(MCH_, 0.773 * mcbagb - 0.35),
                                1.05 * mcx);

            if (aj < tscls[13]) {
                /* 早期 AGB: H-耗尽核心质量恒定 */
                mcx = seq_mc_from_t(sp, aj, gb[8], tscls[7], tscls[8],
                                    tscls[9], mass_f);
                mc = mcbagb;
                lum = seq_l_from_mc(sp, mcx, mass_f);
                if (mt <= mc) {
                    /* 演化后的裸氦星 */
                    Sequence sp2;
                    kw = 9;
                    mt = mc;
                    out.mass = mt;
                    mc = mcx;
                    sequence_init(&sp2, mt, sp->z);
                    if (mc <= gb[7])
                        aj = tscls[4] - 1.0
                             / ((gb[5] - 1.0) * gb[8] * gb[4])
                             * pow(mc, 1.0 - gb[5]);
                    else
                        aj = tscls[5] - 1.0
                             / ((gb[6] - 1.0) * gb[8] * gb[3])
                             * pow(mc, 1.0 - gb[6]);
                    aj = fmax(aj, sp2.tm);
                    helium_star(&sp2, 9, aj, mt, mc,
                                &r, &lum, &kw, &mc, &rc,
                                &menv, &renv, &k2);
                    out.r = r; out.lum = lum; out.kw = kw;
                    out.mc = mc; out.rc = rc;
                    out.menv = menv; out.renv = renv; out.k2 = k2;
                    out.mt = mt;
                    return out;
                }
                kw = 5;
            } else {
                /* TPAGB */
                double lam, tau, mcx_du;
                kw = 6;
                mcx = seq_mc_from_t(sp, aj, gb[2], tscls[10], tscls[11],
                                    tscls[12], mass_f);
                if (mcx <= 0.0) mcx = mcmax;
                mcx = fmax(mcx, mcx_prev);
                lum = seq_l_from_mc(sp, mcx, mass_f);
                lam = fmin(0.9, 0.3 + 0.001 * pow(mass_f, 5.0));
                tau = tscls[13];
                mcx_du = seq_mc_from_t(sp, tau, gb[2], tscls[10],
                                       tscls[11], tscls[12], mass_f);
                mc = mcx - lam * (mcx - mcx_du);
                mcx_prev = mcx;
                mcmax = fmin(mt, mcmax);
            }

            r = seq_ragbf(sp, lum, mt);
            rg = r;

            if (mcmax - mcx < TINY_) {
                aj = 0.0;
                mc = mcmax;
                if (mc < MCH_) {
                    mt = mc;
                    kw = (mcbagb < 1.6) ? 11 : 12;
                    out.mass = mt;
                } else if (mcbagb < 1.6) {
                    kw = 15;
                    aj = 0.0;
                    mt = 0.0;
                    lum = 1.0e-10;
                    r = 1.0e-10;
                    out.mass = 0.0;
                } else {
                    double mcx2;
                    if (mc < 6.0)
                        mcx2 = 1.44;
                    else
                        mcx2 = 1.4512017 * mc - 6.5913737e-3 * mc * mc
                               - 6.1073371;
                    mt = mcx2;
                    mc = mt;
                    kw = (mt <= PH_MXNS) ? 13 : 14;
                    out.mass = mt;
                }
            }
        }
    }

    out.mt = mt;

    /* ------------------------------------------------ 裸氦星 */
    if (kw >= 7 && kw <= 9) {
        helium_star(sp, kw, aj, mt, mc, &r, &lum, &kw, &mc, &rc,
                    &menv, &renv, &k2);
        out.r = r; out.lum = lum; out.kw = kw;
        out.mc = mc; out.rc = rc;
        out.menv = menv; out.renv = renv; out.k2 = k2;
        return out;
    }

    /* ------------------------------------------------ 白矮星 */
    if (kw >= 10 && kw <= 12) {
        double xx;
        mc = mt;
        if (mc >= MCH_ && kw == 12) {
            kw = 13;
            mt = 1.3;
        } else if (mc >= MCH_ - 0.06 && kw <= 11) {
            kw = 15;
            aj = 0.0;
            mt = 0.0;
            lum = 1.0e-10;
            r = 1.0e-10;
        } else {
            xx = (kw == 10) ? PH_AHE : ((kw == 11) ? PH_ACO : PH_AONE);
            lum = 635.0 * mt * pow(sp->z, 0.4)
                  / pow(xx * (aj + 0.1), 1.4);
            r = wd_radius(mt);
        }
        out.mass = mt; out.mt = mt;
        out.r = r; out.lum = lum; out.kw = kw; out.mc = mc;
        return out;
    }

    if (kw == 13) {
        mc = mt;
        if (mc > PH_MXNS) {
            kw = 14;
            aj = 0.0;
        } else {
            lum = 0.02 * pow(mt, 2.0 / 3.0) / pow(fmax(aj, 0.1), 2.0);
            r = PH_RNS;
        }
        out.r = r; out.lum = lum; out.kw = kw; out.mc = mc;
        return out;
    }

    if (kw == 14) {
        mc = mt;
        lum = 1.0e-10;
        r = 4.24e-6 * mt;
        out.r = r; out.lum = lum; out.kw = kw; out.mc = mc;
        return out;
    }

    /* -------------------------------------- H 星的核心半径 */
    {
        double lx = 0.0, rx = 0.0;
        if (kw <= 1 || kw == 7) {
            rc = 0.0;
        } else if (kw <= 3) {
            if (mass_f > sp->mhef) {
                lx = seq_l_zhe(sp, mc);
                rx = seq_r_zhe(sp, mc);
                rc = rx;
            } else {
                lx = 300.0 * mc * pow(sp->z, 0.4)
                     / pow(PH_AHE * 0.1, 1.18);
                rx = 0.0115 * sqrt(fmax(1.48204e-6,
                     pow(MCH_ / mc, 2.0 / 3.0)
                     - pow(mc / MCH_, 2.0 / 3.0)));
                rc = 5.0 * rx;
            }
        } else if (kw == 4) {
            double tau = (aj - tscls[2]) / tscls[3];
            double mc_he;
            Sequence sp7;
            double am;
            if (mass_f <= sp->mhef)
                mc_he = seq_mc_from_l(sp, lums[4], mt);
            else
                mc_he = seq_mc_he_if(sp, mass_f, zp[10]);
            sequence_init(&sp7, mc, sp->z);
            am = fmax(0.0, 0.85 - 0.08 * mc);
            lx = seq_lzams(&sp7, mc) * (1.0 + 0.45 * tau + am * tau * tau);
            rx = seq_r_zhe(&sp7, mc);
            am = fmax(0.0, 0.4 - 0.22 * log10(mc));
            rx = rx * (1.0 + am * (tau - pow(tau, 6.0)));
            rc = rx;
            (void)mc_he;
        } else if (kw == 5) {
            double tbagb_v = sp->t_bagb;
            double tau = (sp->tn > tbagb_v)
                         ? 3.0 * (aj - tbagb_v) / (sp->tn - tbagb_v) : 0.0;
            double mcx = mc;
            Sequence sp9;
            sequence_init(&sp9, mc, sp->z);
            lx = seq_l_from_mc(sp, mcx, mc);
            if (tau < 1.0)
                lx = lums[2] * pow(lx / lums[2], tau);
            rx = seq_r_zhe(&sp9, mc);
            {
                double r1 = seq_r_he_hg(&sp9, mc, lx, rx, lums[2]);
                double r2 = seq_r_he_gb(&sp9, lx);
                rx = fmin(r1, r2);
            }
            rc = rx;
        } else if (kw <= 9) {
            lx = 300.0 * mc * pow(sp->z, 0.4) / pow(PH_ACO * 0.1, 1.18);
            rx = 0.0115 * sqrt(fmax(1.48204e-6,
                 pow(MCH_ / mc, 2.0 / 3.0)
                 - pow(mc / MCH_, 2.0 / 3.0)));
            rc = 5.0 * rx;
        } else {
            rc = r;
            menv = 1.0e-10;
            renv = 1.0e-10;
            k2 = 0.21;
        }

        /* 小包层扰动 (eq. 83-84) */
        if (kw >= 2 && kw <= 9 && kw != 7 && mt > 0.0) {
            double lum_pos = fmax(lum, 1.0e-30);
            double mew = ((mt - mc) / mt)
                         * fmin(5.0, fmax(1.2, pow(lum_pos / PH_LUM0,
                                                    PH_KAP)));
            if (kw >= 8)
                mew = (mtc > 0.0) ? ((mtc - mc) / mtc) * 5.0 : 5.0;
            if (mew < 1.0 && lx > 0.0 && rx > 0.0 && r > 0.0) {
                double xx = lpert(mt, mew);
                lum = lx * pow(lum_pos / lx, xx);
                xx = (r <= rx) ? 0.0 : rpert(mt, mew, r, rx);
                r = rx * pow(r / rx, xx);
            }
            rc = fmin(rc, r);
        }
    }

    /* ------------------------------------------ 对流包层 (mrenv) */
    if (kw >= 0 && kw < 10) {
        mrenv(kw, mass_f, mt, mc, lum, r, rc, aj, sp->tm,
              lums[2], lums[3], lums[4],
              seq_rzams(sp, mass_f), seq_r_tms(sp, mass_f), rg,
              &menv, &renv, &k2);
    }

    out.r = r; out.lum = lum; out.kw = kw;
    out.mc = mc; out.rc = rc;
    out.menv = menv; out.renv = renv; out.k2 = k2;
    return out;
}

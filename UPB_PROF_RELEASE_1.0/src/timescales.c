/* timescales.c -- SSE 演化时标与阶段边界量的实现
 *
 * 对应原 Python 项目中的 sse/timescales.py。
 */

#include "timescales.h"
#include "zams_coeffs.h"

#include <math.h>
#include <stddef.h>

/* ------------------------------------------------------------------ 小工具 */

/* base**exponent，base<=0 时返回 0 */
static double pow_safe(double base, double exponent)
{
    if (base <= 0.0) return 0.0;
    return pow(base, exponent);
}

/* Fortran 风格的实数幂：负底数时取实分支 (-|base|^exp) */
static double real_pow(double base, double exponent)
{
    if (base < 0.0) return -pow(-base, exponent);
    return pow(base, exponent);
}

/* ------------------------------------------------------------------ 构造 */

int sequence_init(Sequence *sp, double mass, double z)
{
    double m, lz, lzd;
    double a_h, bb, dd, p, q;
    double lums[9], tscls[16], gb[10], zp[15];
    int i;

    if (sp == NULL)      return -1;
    if (mass <= 0.0)     return -1;
    if (z <= 0.0)        return -1;

    /* 清零 */
    for (i = 0; i < 9;  i++) lums[i]  = 0.0;
    for (i = 0; i < 16; i++) tscls[i] = 0.0;
    for (i = 0; i < 10; i++) gb[i]    = 0.0;
    for (i = 0; i < 15; i++) zp[i]    = 0.0;

    sp->mass_in = mass;
    sp->mass    = (mass > 100.0) ? 100.0 : mass;
    sp->z       = z;

    if (zcoeffs(z, &sp->c) != 0)
        return -1;

    sp->mhook = sp->c.Mhook;
    sp->mhef  = sp->c.MHeF;
    sp->mfgb  = sp->c.MFGB;
    sp->zeta  = sp->c.zeta;
    lz        = log10(z);
    lzd       = sp->zeta + 1.0;
    sp->lz    = lz;
    sp->lzd   = lzd;
    sp->x_h   = 0.76 - 3.0 * z;
    sp->y_he  = 0.24 + 2.0 * z;

    m = sp->mass;

    /* --- ZAMS 与终端主序 --- */
    lums[1] = seq_lzams(sp, m);
    lums[2] = seq_l_tms(sp, m);
    tscls[1] = seq_t_bgb(sp, m);
    sp->tm   = seq_t_ms(sp, m);

    /* --- 这个质量的巨支参数 --- */
    {
        double log_d0 = 5.37 + sp->zeta * 0.135;
        a_h = pow(10.0, fmax(-4.8, fmin(-5.7 + 0.8 * m, -4.1 + 0.14 * m)));
        bb  = fmax(3.0e4, 500.0 + 1.75e4 * pow(m, 0.6));
        if (m <= 2.0) {
            dd = pow(10.0, log_d0);
            p = 6.0;
            q = 3.0;
        } else if (m < 2.5) {
            double dx = log_d0 - (0.975 * log_d0 - 0.18 * 2.5);
            double d  = log_d0 - dx * (m - 2.0) / 0.5;
            dd = pow(10.0, d);
            p = 6.0 - (m - 2.0) / 0.5;
            q = 3.0 - (m - 2.0) / 0.5;
        } else {
            double d = fmax(-1.0, 0.5 * log_d0 - 0.06 * m);
            d = fmax(d, 0.975 * log_d0 - 0.18 * m);
            dd = pow(10.0, d);
            p = 5.0;
            q = 2.0;
        }
    }
    gb[1] = a_h;
    gb[2] = 1.27e-5;
    gb[3] = bb;
    gb[4] = dd;
    gb[5] = p;
    gb[6] = q;
    gb[7] = pow(bb / dd, 1.0 / (p - q));      /* Mx */
    gb[8] = 8.0e-5;
    lums[6] = dd * pow(gb[7], p);              /* L(Mx) */

    lums[4] = seq_l_he_i(sp, m);
    lums[7] = seq_l_bagb(sp, m);

    /* 参考核心质量 (M = MHeF) */
    zp[9]  = seq_mc_he_if(sp, sp->mhef,
                          seq_mc_from_l(sp, seq_l_bgb(sp, sp->mhef), sp->mhef));
    zp[10] = seq_mc_he_if(sp, sp->mhef,
                          seq_mc_from_l(sp, lums[4], sp->mhef));

    if (m <= sp->mfgb) {
        lums[3] = seq_l_bgb(sp, m);
        /* 巨支时标 (eq. 34-36) */
        tscls[4] = tscls[1] + 1.0 / ((p - 1.0) * a_h * dd)
                   * pow(dd / lums[3], (p - 1.0) / p);
        tscls[6] = tscls[4] - (tscls[4] - tscls[1])
                   * pow(lums[3] / lums[6], (p - 1.0) / p);
        tscls[5] = tscls[6] + 1.0 / ((q - 1.0) * a_h * bb)
                   * pow(bb / lums[6], (q - 1.0) / q);
        /* 氦点火时间 (eq. 44) */
        if (lums[4] <= lums[6]) {
            tscls[2] = tscls[4] - 1.0 / ((p - 1.0) * a_h * dd)
                       * pow(dd / lums[4], (p - 1.0) / p);
        } else {
            tscls[2] = tscls[5] - 1.0 / ((q - 1.0) * a_h * bb)
                       * pow(bb / lums[4], (q - 1.0) / q);
        }
        if (m <= sp->mhef) {
            double mc1 = seq_mc_from_l(sp, lums[4], m);
            sp->mc_he_i_val = mc1;
            lums[5]  = seq_l_zahb(sp, mc1, m);
            tscls[3] = seq_t_he(sp, mc1, m);
        } else {
            sp->mc_he_i_val = seq_mc_he_if(sp, m, zp[10]);
            lums[5]  = seq_l_min_he(sp, m);
            tscls[3] = seq_t_he(sp, 0.0, m);   /* mc = -inf 时不用，见实现 */
        }
    } else {
        tscls[2] = tscls[1];
        tscls[3] = seq_t_he(sp, 0.0, m);
        lums[5]  = lums[7];
        lums[3]  = lums[4];
    }

    /* --- BGB 核心质量 --- */
    if (m <= sp->mhef)
        gb[9] = seq_mc_from_l(sp, lums[3], m);
    else if (m <= sp->mfgb)
        gb[9] = seq_mc_he_if(sp, m, zp[9]);
    else
        gb[9] = seq_mc_he_if(sp, m, zp[10]);
    sp->mc_bgb_val = gb[9];

    /* --- FAGB (EAGB) 时标 (eq. 63-70) --- */
    {
        double tbagb = tscls[2] + tscls[3];
        sp->t_bagb = tbagb;
        tscls[7] = tbagb + 1.0 / ((p - 1.0) * gb[8] * dd)
                   * pow(dd / lums[7], (p - 1.0) / p);
        tscls[9] = tscls[7] - (tscls[7] - tbagb)
                   * pow(lums[7] / lums[6], (p - 1.0) / p);
        tscls[8] = tscls[9] + 1.0 / ((q - 1.0) * gb[8] * bb)
                   * pow(bb / lums[6], (q - 1.0) / q);
    }

    /* --- 第二 dredge-up (eq. 81) --- */
    {
        double mcbagb = seq_mc_bagb(sp, m);
        double mc1    = mcbagb;
        if (mc1 >= 0.8 && mc1 < 2.25)
            mc1 = 0.44 * mc1 + 0.448;
        sp->mc_du = mc1;
        lums[8]   = seq_l_from_mc(sp, mc1, m);
        if (mc1 <= gb[7])
            tscls[13] = tscls[7] - 1.0 / ((p - 1.0) * gb[8] * dd)
                        * pow(mc1, 1.0 - p);
        else
            tscls[13] = tscls[8] - 1.0 / ((q - 1.0) * gb[8] * bb)
                        * pow(mc1, 1.0 - q);
        sp->t_tp = tscls[13];
    }

    /* --- SAGB (TPAGB) 时标 --- */
    {
        double mc1 = sp->mc_du;
        if (mc1 <= gb[7]) {
            tscls[10] = tscls[13] + 1.0 / ((p - 1.0) * gb[2] * dd)
                        * pow(dd / lums[8], (p - 1.0) / p);
            tscls[12] = tscls[10] - (tscls[10] - tscls[13])
                        * pow(lums[8] / lums[6], (p - 1.0) / p);
            tscls[11] = tscls[12] + 1.0 / ((q - 1.0) * gb[2] * bb)
                        * pow(bb / lums[6], (q - 1.0) / q);
        } else {
            tscls[10] = tscls[7];
            tscls[12] = tscls[9];
            tscls[11] = tscls[13] + 1.0 / ((q - 1.0) * gb[2] * bb)
                        * pow(bb / lums[8], (q - 1.0) / q);
        }
    }

    /* --- 达到 Mc,max 的时间 (eq. 75-76) --- */
    {
        double tau   = tscls[2] + tscls[3];
        double mc2   = seq_mc_from_t(sp, tau, gb[8], tscls[7], tscls[8], tscls[9], m);
        double mcbagb = seq_mc_bagb(sp, m);
        double mc1    = sp->mc_du;
        double mcmax  = fmax(fmax(MCH_, 0.773 * mcbagb - 0.35), 1.05 * mc2);
        sp->mc_sn = mcmax;
        if (mcmax <= mc1) {
            if (mcmax <= gb[7])
                tscls[14] = tscls[7] - 1.0 / ((p - 1.0) * gb[8] * dd)
                            * pow(mcmax, 1.0 - p);
            else
                tscls[14] = tscls[8] - 1.0 / ((q - 1.0) * gb[8] * bb)
                            * pow(mcmax, 1.0 - q);
        } else {
            double lam    = fmin(0.9, 0.3 + 0.001 * pow(m, 5.0));
            double mcmax2 = (mcmax - lam * mc1) / (1.0 - lam);
            if (mcmax2 <= gb[7])
                tscls[14] = tscls[10] - 1.0 / ((p - 1.0) * gb[2] * dd)
                            * pow(mcmax2, 1.0 - p);
            else
                tscls[14] = tscls[11] - 1.0 / ((q - 1.0) * gb[2] * bb)
                            * pow(mcmax2, 1.0 - q);
        }
        tscls[14] = fmax(sp->t_bagb, tscls[14]);
        sp->tn    = tscls[14];
    }

    /* --- LM 星 CHeB 最小半径比例 --- */
    {
        double rmin_hef = seq_r_min_he(sp, sp->mhef);
        double rgbf_zahb = seq_r_gb(sp,
            seq_l_zahb(sp, seq_mc_from_l(sp, lums[4], sp->mhef), sp->mhef),
            sp->mhef);
        sp->rmin_ratio = (rgbf_zahb != 0.0) ? rmin_hef / rgbf_zahb : 0.0;
    }

    /* 写回数组 */
    for (i = 0; i < 9;  i++) sp->lums[i]  = lums[i];
    for (i = 0; i < 16; i++) sp->tscls[i] = tscls[i];
    for (i = 0; i < 10; i++) sp->gb[i]    = gb[i];
    for (i = 0; i < 15; i++) sp->zpars[i] = zp[i];

    return 0;
}

/* ============================================================ ZAMS 辅助 */

double seq_lzams(const Sequence *sp, double m)
{
    return zams_lzams(m, sp->z);
}

double seq_rzams(const Sequence *sp, double m)
{
    return zams_rzams(m, sp->z);
}

/* ============================================================ MS/HG 边界量 */

double seq_t_bgb(const Sequence *sp, double m)
{
    const ZCoeffs *a = &sp->c;
    return (a->a1 + a->a2 * pow(m, 4.0) + a->a3 * pow(m, 5.5)
            + pow(m, 7.0))
           / (a->a4 * m * m + a->a5 * pow(m, 7.0));
}

double seq_t_hook(const Sequence *sp, double m)
{
    const ZCoeffs *a = &sp->c;
    double mu = 1.0 - 0.01 * fmax(a->a6 / pow(m, a->a7),
                                  a->a8 + a->a9 / pow(m, a->a10));
    return fmax(mu, 0.5);
}

double seq_x_factor(const Sequence *sp)
{
    double z = sp->z;
    double v1 = 0.95 - (10.0 / 3.0) * (z - 0.01);
    double v2 = fmin(0.99, 0.98 - (100.0 / 7.0) * (z - 0.001));
    return fmax(0.95, fmax(v1, v2));
}

double seq_t_ms(const Sequence *sp, double m)
{
    return fmax(seq_t_hook(sp, m), seq_x_factor(sp) * seq_t_bgb(sp, m));
}

double seq_l_tms(const Sequence *sp, double m)
{
    const ZCoeffs *a = &sp->c;
    return (a->a11 * pow(m, 3.0) + a->a12 * pow(m, 4.0)
            + a->a13 * pow(m, a->a16 + 1.8))
           / (a->a14 + a->a15 * pow(m, 5.0) + pow(m, a->a16));
}

double seq_r_tms(const Sequence *sp, double m)
{
    const ZCoeffs *a = &sp->c;
    double a17 = a->a17;
    double r_lo, r_hi, mstar;

    /* 分支 lo */
    {
        double r = (a->a18 + a->a19 * pow(m, a->a21))
                   / (a->a20 + pow(m, a->a22));
        r_lo = fmax(r, 1.5 * seq_rzams(sp, m));
    }
    /* 分支 hi */
    r_hi = (a->c1 * pow(m, 3.0) + a->a23 * pow(m, a->a26)
            + a->a24 * pow(m, a->a26 + 1.5))
           / (a->a25 + pow(m, 5.0));

    if (m <= a17) return r_lo;
    mstar = a17 + 0.1;
    if (m >= mstar) return r_hi;
    {
        double r_lo_ = (a->a18 + a->a19 * pow(a17, a->a21))
                       / (a->a20 + pow(a17, a->a22));
        r_lo_ = fmax(r_lo_, 1.5 * seq_rzams(sp, a17));
        r_hi = (a->c1 * pow(mstar, 3.0) + a->a23 * pow(mstar, a->a26)
                + a->a24 * pow(mstar, a->a26 + 1.5))
               / (a->a25 + pow(mstar, 5.0));
        return r_lo_ + (r_hi - r_lo_) * (m - a17) / 0.1;
    }
}

double seq_l_bgb(const Sequence *sp, double m)
{
    const ZCoeffs *a = &sp->c;
    return (a->a27 * pow(m, a->a31) + a->a28 * pow(m, a->c2))
           / (a->a29 + a->a30 * pow(m, a->c3) + pow(m, a->a32));
}

/* ================================================== 巨支 L-Mc 关系 (5.2) */

/* 返回 (A_H, B, D, p, q) —— 这里只用 B, D, p, q */
static void gb_params(const Sequence *sp, double m,
                      double *A_H, double *B, double *D,
                      double *p, double *q)
{
    double log_d0 = 5.37 + sp->zeta * 0.135;
    *A_H = pow(10.0, fmax(-4.8, fmin(-5.7 + 0.8 * m, -4.1 + 0.14 * m)));
    *B   = fmax(3.0e4, 500.0 + 1.75e4 * pow(m, 0.6));
    if (m <= 2.0) {
        *D = pow(10.0, log_d0);
        *p = 6.0;
        *q = 3.0;
    } else if (m < 2.5) {
        double dx = log_d0 - (0.975 * log_d0 - 0.18 * 2.5);
        double d  = log_d0 - dx * (m - 2.0) / 0.5;
        *D = pow(10.0, d);
        *p = 6.0 - (m - 2.0) / 0.5;
        *q = 3.0 - (m - 2.0) / 0.5;
    } else {
        double d = fmax(-1.0, 0.5 * log_d0 - 0.06 * m);
        d = fmax(d, 0.975 * log_d0 - 0.18 * m);
        *D = pow(10.0, d);
        *p = 5.0;
        *q = 2.0;
    }
}

double seq_l_from_mc(const Sequence *sp, double mc, double m)
{
    double A_H, B, D, p, q, mx;
    gb_params(sp, m, &A_H, &B, &D, &p, &q);
    mx = pow(B / D, 1.0 / (p - q));
    if (mc <= mx) return D * pow(mc, p);
    return B * pow(mc, q);
}

double seq_mc_from_l(const Sequence *sp, double lum, double m)
{
    double A_H, B, D, p, q, lx;
    gb_params(sp, m, &A_H, &B, &D, &p, &q);
    lx = B * pow(pow(B / D, 1.0 / (p - q)), q);
    if (lum <= lx) return pow_safe(lum / D, 1.0 / p);
    return pow_safe(lum / B, 1.0 / q);
}

/* 内部: 通用 (exp, coeff, base) 计算 */
static void gb_time_args(const Sequence *sp, double t, double a_rate,
                         double tinf1, double tinf2, double tx, double m,
                         double *exp_out, double *coeff_out, double *base_out)
{
    double A_H, B, D, p, q;
    gb_params(sp, m, &A_H, &B, &D, &p, &q);
    if (t <= tx) {
        *exp_out   = p;
        *coeff_out = D;
        *base_out  = (p - 1.0) * a_rate * D * (tinf1 - t);
    } else {
        *exp_out   = q;
        *coeff_out = B;
        *base_out  = (q - 1.0) * a_rate * B * (tinf2 - t);
    }
}

double seq_mc_from_t(const Sequence *sp, double t, double a_rate,
                     double tinf1, double tinf2, double tx, double m)
{
    double e, c, base;
    gb_time_args(sp, t, a_rate, tinf1, tinf2, tx, m, &e, &c, &base);
    if (base <= 0.0) return 0.0;
    return pow_safe(base, 1.0 / (1.0 - e));
}

double seq_l_from_t(const Sequence *sp, double t, double a_rate,
                    double tinf1, double tinf2, double tx, double m)
{
    double e, c, base;
    gb_time_args(sp, t, a_rate, tinf1, tinf2, tx, m, &e, &c, &base);
    if (base <= 0.0) return 0.0;
    return c * pow_safe(base, e / (1.0 - e));
}

/* =================================================== 巨支半径 (5.2/5.4) */

double seq_r_gb(const Sequence *sp, double lum, double m)
{
    const ZCoeffs *b = &sp->c;
    double aa = fmin(b->b4 / pow(m, b->b5), b->b6 / pow(m, b->b7));
    return aa * (pow(lum, b->b1) + b->b2 * pow(lum, b->b3));
}

double seq_ragbf(const Sequence *sp, double lum, double m)
{
    const ZCoeffs *b = &sp->c;
    double m1 = sp->mhef - 0.2;
    double xx, b50, aa;

    if (m >= sp->mhef)      xx = b->b55;
    else if (m >= m1)       xx = 1.0 + 5.0 * (b->b55 - 1.0) * (m - m1);
    else                    xx = 1.0;
    b50 = xx * b->b3;

    if (m <= m1)
        aa = b->b56 + b->b57 * m;
    else if (m >= sp->mhef)
        aa = fmin(b->b51 / pow(m, b->b52), b->b53 / pow(m, b->b54));
    else {
        double a_lo = b->b56 + b->b57 * m1;
        double a_hi = fmin(b->b51 / pow(sp->mhef, b->b52),
                           b->b53 / pow(sp->mhef, b->b54));
        aa = a_lo + 5.0 * (a_hi - a_lo) * (m - m1);
    }
    return aa * (pow(lum, b->b1) + b->b2 * pow(lum, b50));
}

double seq_r_agb(const Sequence *sp, double lum, double m)
{
    return seq_ragbf(sp, lum, m);
}

/* ==================================================== CHeB 关系 (5.3) */

double seq_l_he_i(const Sequence *sp, double m)
{
    const ZCoeffs *b = &sp->c;
    if (m < sp->mhef) {
        double l_im_hef = (b->b11 + b->b12 * pow(sp->mhef, 3.8))
                          / (b->b13 + sp->mhef * sp->mhef);
        double one_plus_alpha1 = (b->b9 * pow(sp->mhef, b->b10)) / l_im_hef;
        return b->b9 * pow(m, b->b10) / one_plus_alpha1;
    }
    return (b->b11 + b->b12 * pow(m, 3.8)) / (b->b13 + m * m);
}

double seq_l_he_ratio(const Sequence *sp, double m)
{
    const ZCoeffs *b = &sp->c;
    return (b->b14 + b->b15 * pow(m, b->b17 + 0.1))
           / (b->b16 + pow(m, b->b17));
}

double seq_l_min_he(const Sequence *sp, double m)
{
    return seq_l_he_i(sp, m) * seq_l_he_ratio(sp, m);
}

double seq_r_min_he(const Sequence *sp, double m)
{
    const ZCoeffs *b = &sp->c;
    double m28 = pow(m, b->b28);
    return (b->b24 * m + pow(b->b25 * m, b->b26) * m28)
           / (b->b27 + m28);
}

double seq_mu_zahb(const Sequence *sp, double mc, double m)
{
    if (sp->mhef == mc) return 0.0;
    return fmin(1.0, fmax(0.0, (m - mc) / (sp->mhef - mc)));
}

double seq_l_zahb(const Sequence *sp, double mc, double m)
{
    const ZCoeffs *b = &sp->c;
    double lzhe = seq_l_zhe(sp, mc);
    double l_min_hef = seq_l_min_he(sp, sp->mhef);
    double alpha2, mm;

    if (l_min_hef == lzhe) return lzhe;
    alpha2 = (b->b18 + lzhe - l_min_hef) / (l_min_hef - lzhe);
    mm = fmax(seq_mu_zahb(sp, mc, m), 1.0e-12);
    return lzhe + (1.0 + b->b20) / (1.0 + b->b20 * pow(mm, 1.6479))
           * (b->b18 * pow(mm, b->b19))
           / (1.0 + alpha2 * exp(15.0 * (m - sp->mhef)));
}

double seq_r_zahb(const Sequence *sp, double mc, double m)
{
    const ZCoeffs *b = &sp->c;
    double mm = fmax(seq_mu_zahb(sp, mc, m), 1.0e-12);
    double f = (1.0 + b->b21) * pow(mm, b->b22)
               / (1.0 + b->b21 * pow(mm, b->b23));
    return (1.0 - f) * seq_r_zhe(sp, mc)
           + f * seq_r_gb(sp, seq_l_zahb(sp, mc, m), m);
}

double seq_t_he_ms(const Sequence *sp, double mc)
{
    (void)sp;
    if (mc < 1.0e-12) mc = 1.0e-12;
    return (0.4129 + 18.81 * pow(mc, 4.0) + 1.853 * pow(mc, 6.0))
           / pow(mc, 6.5);
}

double seq_t_he_ratio(const Sequence *sp, double m, double mc)
{
    const ZCoeffs *b = &sp->c;
    if (m <= sp->mhef) {
        double mm;
        if (sp->mhef != mc)
            mm = fmax((sp->mhef - m) / (sp->mhef - mc), 1.0e-12);
        else
            mm = 1.0;
        return (b->b39 + (seq_t_he_ms(sp, mc) - b->b39) * pow(mm, b->b40))
               * (1.0 + b->b57 * exp(m * b->b56));
    }
    {
        double m5 = pow(m, 5.0);
        return (b->b41 * pow(m, b->b42) + b->b43 * m5) / (b->b44 + m5);
    }
}

double seq_thef_im(const Sequence *sp, double m)
{
    const ZCoeffs *b = &sp->c;
    double m5 = pow(m, 5.0);
    return seq_t_bgb(sp, m) * (b->b41 * pow(m, b->b42) + b->b43 * m5)
           / (b->b44 + m5);
}

double seq_thef_rel(const Sequence *sp, double m, double mc)
{
    const ZCoeffs *b = &sp->c;
    if (m <= sp->mhef) {
        double mm;
        if (mc < 0.0) {          /* mc = None 时用 0，等价于 Python */
            if (m >= sp->mhef) return seq_thef_im(sp, m);
            mc = 0.0;
        }
        if (sp->mhef != mc)
            mm = fmax((sp->mhef - m) / (sp->mhef - mc), 1.0e-12);
        else
            mm = 1.0;
        return (b->b39 + (seq_t_he_ms(sp, mc) - b->b39) * pow(mm, b->b40))
               * (1.0 + b->b57 * exp(m * b->b56));
    }
    return seq_thef_im(sp, m);
}

double seq_t_he(const Sequence *sp, double mc, double m)
{
    const ZCoeffs *b = &sp->c;
    if (m <= sp->mhef) {
        double t_he_im_hef = seq_thef_rel(sp, sp->mhef, -1.0);
        double one_plus_alpha4 = t_he_im_hef / b->b39;
        return fabs(seq_t_he_ratio(sp, m, mc) * one_plus_alpha4);
    }
    return fabs(seq_thef_rel(sp, m, -1.0));
}

double seq_tau_bl(const Sequence *sp, double m)
{
    const ZCoeffs *b = &sp->c;
    double t;

    if (m < sp->mhef) return 1.0;
    {
        double mr = sp->mhef / sp->mfgb;
        if (m <= sp->mfgb) {
            double m1 = m / sp->mfgb;
            double m2 = (m1 > 0.0) ?
                        fmax(log10(m1) / log10(mr), 1.0e-12) : 1.0e-12;
            double alpha_bl = (1.0 - b->b45 * pow(mr, 0.414))
                              * real_pow(log10(mr), -b->b46);
            t = b->b45 * pow(m1, 0.414) + alpha_bl * real_pow(m2, b->b46);
        } else {
            double r1 = 1.0 - seq_r_min_he(sp, m)
                        / seq_ragbf(sp, seq_l_he_i(sp, m), m);
            r1 = fmax(r1, 1.0e-12);
            t = (1.0 - b->b47) * pow(m, b->b48) * pow(r1, b->b49)
                / seq_f_bl(sp, sp->mfgb);
        }
    }
    t = fmin(1.0, fmax(0.0, t));
    return (t < 1.0e-10) ? 0.0 : t;
}

double seq_f_bl(const Sequence *sp, double m)
{
    const ZCoeffs *b = &sp->c;
    double rr = seq_ragbf(sp, seq_l_he_i(sp, m), m);
    if (rr <= 0.0) return 0.0;
    return pow(m, b->b48)
           * pow(fmax(0.0, 1.0 - seq_r_min_he(sp, m) / rr), b->b49);
}

/* ==================================================== 核心质量 (5.2-5.4) */

double seq_mc_bagb(const Sequence *sp, double m)
{
    const ZCoeffs *b = &sp->c;
    return pow(b->b38 + b->b36 * pow(m, b->b37), 0.25);
}

double seq_mc_he_if(const Sequence *sp, double m, double mc_ref)
{
    const ZCoeffs *a = &sp->c;
    double mcbagb = seq_mc_bagb(sp, m);
    double a3 = pow(mc_ref, 4.0)
                - a->c1_mcbgb * pow(sp->mhef, a->c2_mcbgb);
    return fmin(0.95 * mcbagb,
                pow(a3 + a->c1_mcbgb * pow(m, a->c2_mcbgb), 0.25));
}

/* ============================================ 裸氦星 (Section 6.1) */

double seq_l_zhe(const Sequence *sp, double m)
{
    (void)sp;
    /* _LZHE_A = 1.5262e4, p = 41/4, D = (0.0469, 31.18, 29.54) */
    double m15 = m * sqrt(m);
    return 1.5262e4 * pow(m, 41.0 / 4.0)
           / (0.0469 + pow(m, 6.0) * (31.18 + m15 * (29.54 + m15)));
}

double seq_r_zhe(const Sequence *sp, double m)
{
    (void)sp;
    /* 0.2391 m^4.6 / (0.0065 + (0.162 + m) m^3) */
    return 0.2391 * pow(m, 4.6)
           / (0.0065 + (0.162 + m) * pow(m, 3.0));
}

double seq_t_hems(const Sequence *sp, double m)
{
    (void)sp;
    return (0.4129 + 18.81 * pow(m, 4.0) + 1.853 * pow(m, 6.0))
           / pow(m, 6.5);
}

double seq_r_he_hg(const Sequence *sp, double m, double lum,
                   double rx, double lx)
{
    double cm;
    (void)sp;
    cm = 2.0e-3 * pow(m, 2.5) / (2.0 + pow(m, 5.0));
    return rx * pow(lum / lx, 0.2)
           + 0.02 * (exp(cm * lum) - exp(cm * lx));
}

double seq_r_he_gb(const Sequence *sp, double lum)
{
    (void)sp;
    return 0.08 * pow(lum, 0.75);
}

/* ==================================================== L_BAGB (eq. 56) */

double seq_l_bagb(const Sequence *sp, double m)
{
    const ZCoeffs *b = &sp->c;
    if (m < sp->mhef) {
        double l_im_hef = (b->b31 + b->b32 * pow(sp->mhef, b->b33 + 1.8))
                          / (b->b34 + pow(sp->mhef, b->b33));
        double one_plus_alpha3 = (b->b29 * pow(sp->mhef, b->b30)) / l_im_hef;
        return b->b29 * pow(m, b->b30) / one_plus_alpha3;
    }
    return (b->b31 + b->b32 * pow(m, b->b33 + 1.8))
           / (b->b34 + pow(m, b->b33));
}
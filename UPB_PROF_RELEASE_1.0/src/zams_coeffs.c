/* zams_coeffs.c -- Tout et al. (1996) ZAMS 拟合的实现
 *
 * 对应原 Python 项目中的 sse/zams_coeffs.py。
 *
 * 每个系数行 (alpha, beta, gamma, eta, mu) 对应 zeta 的 0..4 次多项式。
 */

#include "zams_coeffs.h"
#include <math.h>

/* --- L_ZAMS 系数: l1..l7 --- */
static const double L_ZAMS[7][5] = {
    { 3.970417e-01, -3.2913574e-01,  3.4776688e-01,  3.7470851e-01,  9.011915e-02 },
    { 8.527626e+00, -2.441225973e+01, 5.643597107e+01, 3.706152575e+01, 5.4562406e+00 },
    { 2.5546e-04,   -1.23461e-03,   -2.3246e-04,    4.5519e-04,     1.6176e-04   },
    { 5.432889e+00, -8.62157806e+00, 1.344202049e+01, 1.451584135e+01, 3.39793084e+00 },
    { 5.563579e+00, -1.032345224e+01,1.944322980e+01, 1.897361347e+01, 4.16903097e+00 },
    { 7.8866060e-01,-2.90870942e+00, 6.54713531e+00,  4.05606657e+00,  5.3287322e-01 },
    { 5.86685e-03,  -1.704237e-02,   3.872348e-02,    2.570041e-02,    3.83376e-03  },
};

/* --- R_ZAMS 系数: r1..r9 (r6 是与 zeta 无关的常数) --- */
static const double R_ZAMS[9][5] = {
    { 1.715359e+00,  6.2246212e-01, -9.2557761e-01, -1.16996966e+00, -3.0631491e-01 },
    { 6.597788e+00, -4.2450044e-01, -1.213339427e+01,-1.073509484e+01, -2.51487077e+00},
    { 1.008855000e+01,-7.11727086e+00,-3.167119479e+01,-2.424848322e+01, -5.33608972e+00},
    { 1.012495e+00,  3.2699690e-01, -9.23418e-03,   -3.876858e-02,    -4.12750e-03  },
    { 7.490166e-02,  2.410413e-02,   7.233664e-02,    3.040467e-02,     1.97741e-03  },
    { 1.077422e-02,  0.0,            0.0,             0.0,              0.0          },
    { 3.082234e+00,  9.447205e-01, -2.15200882e+00, -2.49219496e+00, -6.3848738e-01 },
    { 1.784778e+01, -7.4534569e+00, -4.896066856e+01,-4.005386135e+01, -9.09331816e+00},
    { 2.2582e-04,   -1.86899e-03,   3.88783e-03,     1.42402e-03,     -7.671e-05    },
};

/* alpha + beta*zeta + gamma*zeta^2 + eta*zeta^3 + mu*zeta^4 */
static double poly5(double zeta, const double row[5])
{
    return row[0] + zeta * (row[1] + zeta * (row[2]
                + zeta * (row[3] + zeta * row[4])));
}

void zams_lzams_coeffs(double z, double out[7])
{
    double zeta = log10(z / 0.02);
    int i;
    for (i = 0; i < 7; i++)
        out[i] = poly5(zeta, L_ZAMS[i]);
}

void zams_rzams_coeffs(double z, double out[9])
{
    double zeta = log10(z / 0.02);
    int i;
    for (i = 0; i < 9; i++)
        out[i] = poly5(zeta, R_ZAMS[i]);
}

double zams_lzams(double m, double z)
{
    double c[7];
    double mx = sqrt(m);
    double num, den;

    zams_lzams_coeffs(z, c);

    num = c[0] * pow(m, 5.0) * mx + c[1] * pow(m, 11.0);
    den = c[2] + pow(m, 3.0)
          + c[3] * pow(m, 5.0)
          + c[4] * pow(m, 7.0)
          + c[5] * pow(m, 8.0)
          + c[6] * pow(m, 9.0) * mx;
    return num / den;
}

double zams_rzams(double m, double z)
{
    double c[9];
    double mx = sqrt(m);
    double num, den;

    zams_rzams_coeffs(z, c);

    num = (c[0] * m * m + c[1] * pow(m, 6.0)) * mx
          + c[2] * pow(m, 11.0)
          + (c[3] + c[4] * mx) * pow(m, 19.0);
    den = c[5] + c[6] * m * m
          + (c[7] * pow(m, 8.0) + pow(m, 18.0) + c[8] * pow(m, 19.0)) * mx;
    return num / den;
}

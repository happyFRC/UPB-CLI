/* zcoeffs.c -- SSE 附录系数的实现
 *
 * 对应原 Python 项目中的 sse/zcoeffs.py 的 zcoeffs() 函数。
 * 逐步翻译，保持与 Python 版本一一对应。
 */

#include "zcoeffs.h"
#include <math.h>
#include <stddef.h>

/* 多项式: alpha + beta*zeta + gamma*zeta^2 + eta*zeta^3 + mu*zeta^4 */
static double poly5(double zeta, double alpha, double beta,
                    double gamma, double eta, double mu)
{
    return alpha + zeta * (beta + zeta * (gamma + zeta * (eta + zeta * mu)));
}

/* 大多数只有前几项的简化调用 */
static double poly1(double zeta, double alpha)
{
    return alpha;
}
static double poly2(double zeta, double alpha, double beta)
{
    return alpha + beta * zeta;
}
static double poly3(double zeta, double alpha, double beta, double gamma)
{
    return alpha + zeta * (beta + zeta * gamma);
}
static double poly4(double zeta, double alpha, double beta,
                    double gamma, double eta)
{
    return alpha + zeta * (beta + zeta * (gamma + zeta * eta));
}

/* 论文 eq.(18): MS 半径的 alpha_R 系数。
 * a68 必须是尚未做 a68 = min(a68, a66) 钳制之前的值。 */
static double alpha_R_from(double a58, double a59, double a60, double a61,
                           double a62, double a63, double a64, double a65,
                           double a66, double a67, double a68)
{
    double m = a66;
    double b_val, c_val;

    if (m <= 0.5)
        return a62;
    if (m <= 0.65)
        return a62 + ((a63 - a62) / 0.15) * (m - 0.5);
    if (m <= a68)
        return a63 + ((a64 - a63) / (a68 - 0.65)) * (m - 0.65);
    if (m <= a66) {
        b_val = (a58 * pow(a66, a60)) / (a59 * pow(a66, a61));
        return a64 + ((b_val - a64) / (a66 - a68)) * (m - a68);
    }
    if (m <= a67)
        return (a58 * pow(m, a60)) / (a59 * pow(m, a61));
    c_val = (a58 * pow(a67, a60)) / (a59 * pow(a67, a61));
    return c_val + a65 * (m - a67);
}

int zcoeffs(double Z, ZCoeffs *o)
{
    double zeta, sigma, rho;
    double loga17, a18p, a19p, a29p, a34_, a11p, a12p;
    double b2_, b3p;
    double Mhook, MHeF, MFGB;
    double use_alpha_R, a68_before_min;

    if (o == NULL)
        return -1;
    if (Z <= 0.0)
        return -1;

    zeta  = log10(Z / 0.02);
    sigma = log10(Z);
    rho   = zeta + 1.0;

    o->z     = Z;
    o->zeta  = zeta;
    o->sigma = sigma;
    o->rho   = rho;

    /* --- 临界质量 Mhook, MHeF, MFGB --- */
    Mhook = 1.0185 + 0.16015 * zeta + 0.0892 * zeta * zeta;
    MHeF  = 1.995  + 0.25    * zeta + 0.087  * zeta * zeta;
    MFGB  = 13.048 * pow(Z / 0.02, 0.06)
            / (1.0 + 0.0012 * pow(0.02 / Z, 1.27));
    o->Mhook = Mhook;
    o->MHeF  = MHeF;
    o->MFGB  = MFGB;

    /* --- a1..a5 (tBGB) --- */
    o->a1 = poly4(zeta, 1.593890e3, 2.053038e3, 1.231226e3, 2.327785e2);
    o->a2 = poly4(zeta, 2.706708e3, 1.483131e3, 5.772723e2, 7.411230e1);
    o->a3 = poly4(zeta, 1.466143e2, -1.048442e2, -6.795374e1, -1.391127e1);
    o->a4 = poly4(zeta, 4.141960e-2, 4.564888e-2, 2.958542e-2, 5.571483e-3);
    o->a5 = poly1(zeta, 3.426349e-1);

    /* --- a6..a10 (tMS) --- */
    o->a6  = poly4(zeta, 1.949814e1, 1.758178e0, -6.008212e0, -4.470533e0);
    o->a7  = poly1(zeta, 4.903830e0);
    o->a8  = poly4(zeta, 5.212154e-2, 3.166411e-2, -2.750074e-3, -2.271549e-3);
    o->a9  = poly4(zeta, 1.312179e0, -3.294936e-1, 9.231860e-2, 2.610989e-2);
    o->a10 = poly1(zeta, 8.073972e-1);

    /* --- a11'..a16 (L_TMS) --- */
    a11p = poly5(zeta, 1.031538e0, -2.434480e-1, 7.732821e0, 6.460705e0, 1.374484e0);
    a12p = poly5(zeta, 1.043715e0, -1.577474e0, -5.168234e0, -5.596506e0, -1.299394e0);
    o->a13 = poly4(zeta, 7.859573e2, -8.542048e0, -2.642511e1, -9.585707e0);
    o->a14 = poly5(zeta, 3.858911e3, 2.459681e3, -7.630093e1, -3.486057e2, -4.861703e1);
    o->a15 = poly4(zeta, 2.888720e2, 2.952979e2, 1.850341e2, 3.797254e1);
    o->a16 = poly4(zeta, 7.196580e0, 5.613746e-1, 3.805871e-1, 8.398728e-2);
    o->a11 = a11p * o->a14;
    o->a12 = a12p * o->a14;

    /* --- a18'..a26, a17, a18, a19 (R_TMS) --- */
    a18p = poly5(zeta, 2.187715e-1, -2.154437e0, -3.768678e0, -1.975518e0, -3.021475e-1);
    a19p = poly5(zeta, 1.466440e0, 1.839725e0, 6.442199e0, 4.023635e0, 6.957529e-1);
    o->a20 = poly5(zeta, 2.652091e1, 8.178458e1, 1.156058e2, 7.633811e1, 1.950698e1);
    o->a21 = poly4(zeta, 1.472103e0, -2.947609e0, -3.312828e0, -9.945065e-1);
    o->a22 = poly4(zeta, 3.071048e0, -5.679941e0, -9.745523e0, -3.594543e0);
    o->a23 = poly4(zeta, 2.617890e0, 1.019135e0, -3.292551e-2, -7.445123e-2);
    o->a24 = poly4(zeta, 1.075567e-2, 1.773287e-2, 9.610479e-3, 1.732469e-3);
    o->a25 = poly4(zeta, 1.476246e0, 1.899331e0, 1.195010e0, 3.035051e-1);
    o->a26 = poly4(zeta, 5.502535e0, -6.601663e-2, 9.968707e-2, 3.599801e-2);
    loga17 = fmax(0.097 - 0.1072 * (sigma + 3.0),
                  fmax(0.097, fmin(0.1461, 0.1461 + 0.1237 * (sigma + 2.0))));
    o->a17 = pow(10.0, loga17);
    o->a18 = a18p * o->a20;
    o->a19 = a19p * o->a20;

    /* --- a27..a32 (L_BGB) --- */
    o->a27 = poly4(zeta, 9.511033e1, 6.819618e1, -1.045625e1, -1.474939e1);
    o->a28 = poly4(zeta, 3.113458e1, 1.012033e1, -4.650511e0, -2.463185e0);
    a29p   = poly4(zeta, 1.413057e0, 4.578814e-1, -6.850581e-2, -5.588658e-2);
    o->a30 = poly4(zeta, 3.910862e1, 5.196646e1, 2.264970e1, 2.873680e0);
    o->a31 = poly3(zeta, 4.597479e0, -2.855179e-1, 2.709724e-1);
    o->a32 = poly3(zeta, 6.682518e0, 2.827718e-1, -7.294429e-2);
    o->a29 = pow(a29p, o->a32);

    /* --- a34..a37, a33 (Delta L hook) --- */
    o->a34 = poly4(zeta, 1.910302e-1, 1.158624e-1, 3.348990e-2, 2.599706e-3);
    o->a35 = poly4(zeta, 3.931056e-1, 7.277637e-2, -1.366593e-1, -4.508946e-2);
    o->a36 = poly4(zeta, 3.267776e-1, 1.204424e-1, 9.988332e-2, 2.455361e-2);
    o->a37 = poly4(zeta, 5.990212e-1, 5.570264e-2, 6.207626e-2, 1.777283e-2);
    o->a33 = fmin(1.4, 1.5135 + 0.3769 * zeta);
    o->a33 = fmax(0.6355 - 0.4192 * zeta, fmax(1.25, o->a33));

    /* --- a38..a44 (Delta R hook) --- */
    o->a38 = poly4(zeta, 7.330122e-1, 5.192827e-1, 2.316416e-1, 8.346941e-3);
    o->a39 = poly4(zeta, 1.172768e0, -1.209262e-1, -1.193023e-1, -2.859837e-2);
    o->a40 = poly4(zeta, 3.982622e-1, -2.296279e-1, -2.262539e-1, -5.219837e-2);
    o->a41 = poly4(zeta, 3.571038e0, -2.223625e-2, -2.611794e-2, -6.359648e-3);
    o->a42 = poly3(zeta, 1.9848e0, 1.1386e0, 3.5640e-1);
    o->a43 = poly3(zeta, 6.300e-2, 4.810e-2, 9.840e-3);
    o->a44 = poly2(zeta, 1.200e0, 2.450e0);
    o->a42 = fmin(1.25, fmax(1.1, o->a42));
    o->a44 = fmin(1.3, fmax(0.45, o->a44));

    /* --- a45..a53 (alpha_L) --- */
    o->a45 = poly4(zeta, 2.321400e-1, 1.828075e-3, -2.232007e-2, -3.378734e-3);
    o->a46 = poly4(zeta, 1.163659e-2, 3.427682e-3, 1.421393e-3, -3.710666e-3);
    o->a47 = poly4(zeta, 1.048020e-2, -1.231921e-2, -1.686860e-2, -4.234354e-3);
    o->a48 = poly4(zeta, 1.555590e0, -3.223927e-1, -5.197429e-1, -1.066441e-1);
    o->a49 = poly3(zeta, 9.7700e-2, -2.3100e-1, -7.5300e-2);
    o->a50 = poly3(zeta, 2.4000e-1, 1.8000e-1, 5.9500e-1);
    o->a51 = poly3(zeta, 3.3000e-1, 1.3200e-1, 2.1800e-1);
    o->a52 = poly3(zeta, 1.1064e0, 4.1500e-1, 1.8000e-1);
    o->a53 = poly3(zeta, 1.1900e0, 3.7700e-1, 1.7600e-1);
    o->a49 = fmax(o->a49, 0.145);
    o->a50 = fmin(o->a50, 0.306 + 0.053 * zeta);
    o->a51 = fmin(o->a51, 0.3625 + 0.062 * zeta);
    o->a52 = fmax(o->a52, 0.9);
    if (Z > 0.01) o->a52 = fmin(o->a52, 1.0);
    o->a53 = fmax(o->a53, 1.0);
    if (Z > 0.01) o->a53 = fmin(o->a53, 1.1);

    /* --- a54..a57 (beta_L) --- */
    o->a54 = poly5(zeta, 3.855707e-1, -6.104166e-1, 5.676742e0, 1.060894e1, 5.284014e0);
    o->a55 = poly5(zeta, 3.579064e-1, -6.442936e-1, 5.494644e0, 1.054952e1, 5.280991e0);
    o->a56 = poly3(zeta, 9.587587e-1, 8.777464e-1, 2.017321e-1);
    o->a57 = poly2(zeta, 1.5135e0, 3.7690e-1);
    o->a57 = fmin(1.4, o->a57);
    o->a57 = fmax(0.6355 - 0.4192 * zeta, fmax(1.25, o->a57));

    /* --- a58..a68 (R_ZAMS) --- */
    o->a58 = poly4(zeta, 4.907546e-1, -1.683928e-1, -3.108742e-1, -7.202918e-2);
    o->a59 = poly4(zeta, 4.537070e0, -4.465455e0, -1.612690e0, -1.623246e0);
    o->a60 = poly4(zeta, 1.796220e0, 2.814020e-1, 1.423325e0, 3.421036e-1);
    o->a61 = poly4(zeta, 2.256216e0, 3.773400e-1, 1.537867e0, 4.396373e-1);
    o->a62 = poly3(zeta, 8.4300e-2, -4.7500e-2, -3.5200e-2);
    o->a63 = poly3(zeta, 7.3600e-2, 7.4900e-2, 4.4260e-2);
    o->a64 = poly2(zeta, 1.3600e-1, 3.5200e-2);
    o->a65 = poly5(zeta, 1.564231e-3, 1.653042e-3, -4.439786e-3, -4.951011e-3, -1.216530e-3);
    o->a66 = poly2(zeta, 1.4770e0, 2.9600e-1);
    o->a67 = poly3(zeta, 5.210157e0, -4.143695e0, -2.120870e0);
    o->a68 = poly2(zeta, 1.1160e0, 1.6600e-1);
    o->a62 = fmax(0.065, o->a62);
    if (Z < 0.004) o->a63 = fmin(0.055, o->a63);
    o->a64 = fmax(0.091, fmin(0.121, o->a64));
    o->a66 = fmax(o->a66, fmin(1.6, -0.308 - 1.046 * zeta));
    o->a66 = fmax(0.8, fmin(0.8 - 2.0 * zeta, o->a66));
    o->a68 = fmax(0.9, fmin(o->a68, 1.0));
    use_alpha_R = (o->a68 > o->a66);
    a68_before_min = o->a68;
    o->a68 = fmin(o->a68, o->a66);

    /* --- a69..a74 (beta_R) --- */
    o->a69 = poly4(zeta, 1.071489e0, -1.164852e-1, -8.623831e-2, -1.582349e-2);
    o->a70 = poly4(zeta, 7.108492e-1, 7.935927e-1, 3.926983e-1, 3.622146e-2);
    o->a71 = poly4(zeta, 3.478514e0, -2.585474e-2, -1.512955e-2, -2.833691e-3);
    /* a72 的说明见 Python 版本注释：这里按参考实现取 zeta^2 项 */
    o->a72 = poly3(zeta, 9.132108e-1, -1.653695e-1, 3.636784e-2);
    o->a73 = poly4(zeta, 3.969331e-3, 4.539076e-3, 1.720906e-3, 1.897857e-4);
    o->a74 = poly3(zeta, 1.600e0, 7.640e-1, 3.322e-1);
    if (Z > 0.01) o->a72 = fmax(o->a72, 0.95);
    o->a74 = fmax(1.4, fmin(o->a74, 1.6));

    /* --- a75..a81 (gamma) --- */
    o->a75 = poly2(zeta, 8.109e-1, -6.282e-1);
    o->a76 = poly4(zeta, 1.192334e-2, 1.083057e-2, 1.230969e0, 1.551656e0);
    o->a77 = poly4(zeta, -1.668868e-1, 5.818123e-1, -1.105027e1, -1.668070e1);
    o->a78 = poly4(zeta, 7.615495e-1, 1.068243e-1, -2.011333e-1, -9.371415e-2);
    o->a79 = poly2(zeta, 9.409838e0, 1.522928e0);
    o->a80 = poly3(zeta, -2.7110e-1, -5.7560e-1, -8.3800e-2);
    o->a81 = poly2(zeta, 2.4930e0, 1.1475e0);
    o->a75 = fmax(1.0, fmin(o->a75, 1.27));
    o->a75 = fmax(o->a75, 0.6355 - 0.4192 * zeta);
    o->a76 = fmax(o->a76,
                  -0.1015564 - 0.2161264 * zeta - 0.05182516 * zeta * zeta);
    o->a77 = fmax(-0.3868776 - 0.5457078 * zeta - 0.1463472 * zeta * zeta,
                  fmin(0.0, o->a77));
    o->a78 = fmax(0.0, fmin(o->a78, 7.454 + 9.046 * zeta));
    o->a79 = fmin(o->a79, fmax(2.0, -13.3 - 18.6 * zeta));
    o->a80 = fmax(0.0585542, o->a80);
    o->a81 = fmin(1.5, fmax(0.4, o->a81));

    /* --- b1..b7 (R_GB, R_AGB) --- */
    o->b1 = poly3(zeta, 3.9700e-1, 2.8826e-1, 5.2930e-1);
    o->b4 = poly5(zeta, 9.960283e-1, 8.164393e-1, 2.383830e0, 2.223436e0, 8.638115e-1);
    o->b5 = poly5(zeta, 2.561062e-1, 7.072646e-2, -5.444596e-2, -5.798167e-2, -1.349129e-2);
    o->b6 = poly5(zeta, 1.157338e0, 1.467883e0, 4.299661e0, 3.130500e0, 6.992080e-1);
    o->b7 = poly5(zeta, 4.022765e-1, 3.050010e-1, 9.962137e-1, 7.914079e-1, 1.728098e-1);
    o->b1 = fmin(0.54, o->b1);
    b2_   = pow(10.0, -4.6739 - 0.9394 * sigma);
    b2_   = fmin(fmax(b2_, -0.04167 + 55.67 * Z),
                 0.4771 - 9329.21 * pow(Z, 2.94));
    o->b2 = b2_;
    b3p = fmax(-0.1451, -2.2794 - 1.5175 * sigma - 0.254 * sigma * sigma);
    o->b3 = pow(10.0, b3p);
    if (Z > 0.004)
        o->b3 = fmax(o->b3, 0.7307 + 14265.1 * pow(Z, 3.395));
    o->b4 = o->b4 + 0.1231572 * pow(zeta, 5.0);
    o->b6 = o->b6 + 0.01640687 * pow(zeta, 5.0);

    /* --- b9..b13 (L_HeI) --- */
    o->b9  = poly2(zeta, 2.751631e3, 3.557098e2);
    o->b10 = poly2(zeta, -3.820831e-2, 5.872664e-2);
    {
        double b11p = poly3(zeta, 1.071738e2, -8.970339e1, -3.949739e1);
        double b13p = poly3(zeta, 9.219293e0, -2.005865e0, -5.561309e-1);
        o->b11 = b11p * b11p;
        o->b13 = b13p * b13p;
    }
    o->b12 = poly3(zeta, 7.348793e2, -1.531020e2, -3.793700e1);

    /* --- b14'..b16', b14, b16, b17 (L_minHe) --- */
    {
        double b14p = poly3(zeta, 2.917412e0, 1.575290e0, 5.751814e-1);
        double b16p = poly3(zeta, 4.916389e0, 2.862149e0, 7.844850e-1);
        o->b15 = poly3(zeta, 3.629118e0, -9.112722e-1, 1.042291e0);
        o->b14 = pow(b14p, o->b15);
        o->b16 = pow(b16p, o->b15);
    }
    o->b17 = 1.0;
    if (zeta > -1.0)
        o->b17 = 1.0 - 0.3880523 * pow(zeta + 1.0, 2.862149);

    /* --- b18..b23 (L_ZAHB, R_ZAHB) --- */
    o->b18 = poly3(zeta, 5.496045e1, -1.289968e1, 6.385758e0);
    o->b19 = poly3(zeta, 1.832694e0, -5.766608e-2, 5.696128e-2);
    o->b20 = poly1(zeta, 1.211104e2);
    o->b21 = poly4(zeta, 2.214088e2, 2.187113e2, 1.170177e1, -2.635340e1);
    o->b22 = poly4(zeta, 2.063983e0, 7.363827e-1, 2.654323e-1, -6.140719e-2);
    o->b23 = poly4(zeta, 2.003160e0, 9.388871e-1, 9.656450e-1, 2.362266e-1);

    /* --- b24'..b28, b24, b26, b27 (R_mHe) --- */
    {
        double b24p = poly4(zeta, 1.609901e1, 7.391573e0, 2.277010e1, 8.334227e0);
        double b27p = poly4(zeta, 2.752869e0, 2.729201e-2, 4.996927e-1, 2.496551e-1);
        o->b25 = poly4(zeta, 1.747500e-1, 6.271202e-2, -2.324229e-2, -1.844559e-2);
        o->b28 = poly4(zeta, 3.518506e0, 1.112440e0, -4.556216e-1, -2.179426e-1);
        o->b24 = pow(b24p, o->b28);
        o->b26 = 5.0 - 0.09138012 * pow(Z, -0.3671407);
        o->b27 = pow(b27p, 2.0 * o->b28);
    }

    /* --- b29..b34 (L_BAGB) --- */
    o->b29 = poly3(zeta, 1.626062e2, -1.168838e1, -5.498343e0);
    o->b30 = poly3(zeta, 3.336833e-1, -1.458043e-1, -2.011751e-2);
    {
        double b31p = poly4(zeta, 7.425137e1, 1.790236e1, 3.033910e1, 1.018259e1);
        double b34p = poly4(zeta, 1.127018e1, 1.622158e0, -1.443664e0, -9.474699e-1);
        o->b32 = poly4(zeta, 9.268325e2, -9.739859e1, -7.702152e1, -3.158268e1);
        o->b33 = poly2(zeta, 2.474401e0, 3.892972e-1);
        o->b31 = pow(b31p, o->b33);
        o->b34 = pow(b34p, o->b33);
    }

    /* --- b36'..b38', b36, b37, b38 (M_c,BAGB) --- */
    {
        double b36p = poly4(zeta, 1.445216e-1, -6.180219e-2, 3.093878e-2, 1.567090e-2);
        double b37p = poly4(zeta, 1.304129e0, 1.395919e-1, 4.142455e-3, -9.732503e-3);
        double b38p = poly2(zeta, 5.114149e-1, -1.160850e-2);
        o->b36 = pow(b36p, 4.0);
        o->b37 = 4.0 * b37p;
        o->b38 = pow(b38p, 4.0);
    }

    /* --- b39..b44 (t_He) --- */
    o->b39 = poly4(zeta, 1.314955e2, 2.009258e1, -5.143082e-1, -1.379140e0);
    o->b40 = poly3(zeta, 1.823973e1, -3.074559e0, -4.307878e0);
    {
        double b41p = poly4(zeta, 2.327037e0, 2.403445e0, 1.208407e0, 2.087263e-1);
        double b44p = poly4(zeta, 2.327409e0, 6.901582e-1, -2.158431e-1, -1.084117e-1);
        o->b42 = poly2(zeta, 1.997378e0, -8.126205e-1);
        o->b43 = poly4(zeta, 1.079113e-1, 1.762409e-2, 1.096601e-2, 3.058818e-3);
        o->b40 = fmax(o->b40, 1.0);
        o->b41 = pow(b41p, o->b42);
        o->b44 = pow(b44p, 5.0);
    }

    /* --- b45..b49 (t_bl) --- */
    o->b45 = 1.0 - (2.47162 * rho - 5.401682 * rho * rho
                    + 3.247361 * rho * rho * rho);
    if (rho <= 0.0) o->b45 = 1.0;
    {
        double b46_ = poly2(zeta, 2.214315e0, -1.975747e0);
        double b48_ = poly4(zeta, 5.072525e0, 1.146189e1, 6.961724e0, 1.316965e0);
        o->b46 = -1.0 * b46_ * log10(MHeF / MFGB);
        o->b47 = 1.127733 * rho + 0.2344416 * rho * rho
                 - 0.3793726 * rho * rho * rho;
        o->b48 = -1.0 * b48_ * log10(MHeF / MFGB);
        o->b49 = poly1(zeta, 5.139740e0);
    }

    /* --- 第 5 节与 Z 无关的常数 --- */
    o->c1       = -8.672073e-2;
    o->c2       =  9.301992;
    o->c3       =  4.637345;
    o->c1_mcbgb =  9.20925e-5;
    o->c2_mcbgb =  5.402216;
    o->eta_max  = (Z > 0.0009) ? 10.0 : 20.0;

    /* --- b51'..b57', b51, b53, b55, b56, b57 (R_AGB) --- */
    {
        double b51p = poly5(zeta, 1.125124e0, 1.306486e0, 3.622359e0, 2.601976e0, 3.031270e-1);
        double b53p = poly5(zeta, 1.467794e0, 2.798142e0, 9.455580e0, 8.963904e0, 3.339719e0);
        double b56p = poly5(zeta, 1.110866e0, 9.623856e-1, 2.735487e0, 2.445602e0, 8.826352e-1);
        double b57p = poly5(zeta, -1.584333e-1, -1.728865e-1, -4.461431e-1, -3.925259e-1, -1.276203e-1);
        o->b51 = b51p - 0.1343798 * pow(zeta, 5.0);
        o->b52 = poly5(zeta, 3.349489e-1, 4.531269e-3, 1.131793e-1, 2.300156e-1, 7.632745e-2);
        o->b53 = b53p + 0.4426929 * pow(zeta, 5.0);
        o->b54 = poly5(zeta, 4.658512e-1, 2.597451e-1, 9.048179e-1, 7.394505e-1, 1.607092e-1);
        o->b55 = poly3(zeta, 1.0422e0, 1.3156e-1, 4.5000e-2);
        o->b55 = fmin(0.99164 - 743.123 * pow(Z, 2.83), o->b55);
        o->b56 = b56p + 0.1140142 * pow(zeta, 5.0);
        o->b57 = b57p - 0.01308728 * pow(zeta, 5.0);
    }

    /* --- 由 alpha_R 关系反推 a64 --- */
    if (use_alpha_R) {
        o->a64 = alpha_R_from(o->a58, o->a59, o->a60, o->a61,
                              o->a62, o->a63, o->a64, o->a65,
                              o->a66, o->a67, a68_before_min);
    }

    return 0;
}

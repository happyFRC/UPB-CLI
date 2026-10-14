/* zcoeffs.h -- SSE 附录系数（随金属度 Z 变化）
 *
 * 对应原 Python 项目中的 sse/zcoeffs.py。
 *
 * Hurley, Pols & Tout (2000), MNRAS 315, 543 的附录给出了每个系数
 *
 *     a_n, b_n = alpha + beta*zeta + gamma*zeta^2 + eta*zeta^3 + mu*zeta^4
 *
 * 其中 zeta = log10(Z / 0.02)，sigma = log10(Z)，rho = zeta + 1.0。
 *
 * 用法:
 *     ZCoeffs c;
 *     zcoeffs(0.02, &c);
 *     printf("a1 = %g\n", c.a1);
 */

#ifndef ZCOEFFS_H
#define ZCOEFFS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    /* 基本信息 */
    double z;        /* 金属度 Z */
    double zeta;     /* log10(Z / 0.02) */
    double sigma;    /* log10(Z) */
    double rho;      /* zeta + 1.0 */

    /* 临界质量 */
    double Mhook;
    double MHeF;
    double MFGB;

    /* a1..a5 (t_BGB) */
    double a1, a2, a3, a4, a5;
    /* a6..a10 (t_MS) */
    double a6, a7, a8, a9, a10;
    /* a11..a16 (L_TMS) */
    double a11, a12, a13, a14, a15, a16;
    /* a17..a26 (R_TMS) */
    double a17, a18, a19, a20, a21, a22, a23, a24, a25, a26;
    /* a27..a32 (L_BGB) */
    double a27, a28, a29, a30, a31, a32;
    /* a33..a37 (Delta L hook) */
    double a33, a34, a35, a36, a37;
    /* a38..a44 (Delta R hook) */
    double a38, a39, a40, a41, a42, a43, a44;
    /* a45..a53 (alpha_L) */
    double a45, a46, a47, a48, a49, a50, a51, a52, a53;
    /* a54..a57 (beta_L) */
    double a54, a55, a56, a57;
    /* a58..a68 (R_ZAMS) */
    double a58, a59, a60, a61, a62, a63, a64, a65, a66, a67, a68;
    /* a69..a74 (beta_R) */
    double a69, a70, a71, a72, a73, a74;
    /* a75..a81 (gamma) */
    double a75, a76, a77, a78, a79, a80, a81;

    /* b1..b7 (R_GB, R_AGB) */
    double b1, b2, b3, b4, b5, b6, b7;
    /* b9..b13 (L_HeI) */
    double b9, b10, b11, b12, b13;
    /* b14..b17 (L_minHe) */
    double b14, b15, b16, b17;
    /* b18..b23 (L_ZAHB, R_ZAHB) */
    double b18, b19, b20, b21, b22, b23;
    /* b24..b28 (R_mHe) */
    double b24, b25, b26, b27, b28;
    /* b29..b34 (L_BAGB) */
    double b29, b30, b31, b32, b33, b34;
    /* b36..b38 (M_c,BAGB) */
    double b36, b37, b38;
    /* b39..b44 (t_He) */
    double b39, b40, b41, b42, b43, b44;
    /* b45..b49 (t_bl) */
    double b45, b46, b47, b48, b49;
    /* b51..b57 (R_AGB) */
    double b51, b52, b53, b54, b55, b56, b57;

    /* 第 5 节中与 Z 无关的常数 */
    double c1;            /* R_TMS 的 c1 = -8.672073e-2 */
    double c2;            /* L_BGB 的 c2 = 9.301992 */
    double c3;            /* L_BGB 的 c3 = 4.637345 */
    double c1_mcbgb;      /* eq.46 的 c1 = 9.20925e-5 */
    double c2_mcbgb;      /* eq.46 的 c2 = 5.402216 */
    double eta_max;       /* eq.13 中 MS 光度多项式的指数上限 */
} ZCoeffs;

/* 计算 Z 对应的全部系数。Z 必须 > 0。
 * 成功返回 0，失败（Z<=0）返回 -1。 */
int zcoeffs(double Z, ZCoeffs *out);

#ifdef __cplusplus
}
#endif

#endif /* ZCOEFFS_H */

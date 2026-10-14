/* zams_coeffs.h -- Tout et al. (1996) 零龄主序拟合
 *
 * 对应原 Python 项目中的 sse/zams_coeffs.py。
 *
 * Hurley, Pols & Tout (2000) 采用 Tout et al. (1996, MNRAS 281, 257) 的
 * ZAMS 拟合作为 L 和 R 的起点。这里把系数和求值函数独立出来：
 *
 *   zeta = log10(Z / 0.02)
 *   每个多项式系数 c(zeta) = alpha + beta*zeta + gamma*zeta^2
 *                            + eta*zeta^3 + mu*zeta^4
 *
 *   L_ZAMS = (l1 m^5.5 + l2 m^11) /
 *            (l3 + m^3 + l4 m^5 + l5 m^7 + l6 m^8 + l7 m^9.5)
 *
 *   R_ZAMS = ((r1 m^2 + r2 m^6) m^0.5 + r3 m^11
 *             + (r4 + r5 m^0.5) m^19) /
 *            (r6 + r7 m^2 + (r8 m^8 + m^18 + r9 m^19) m^0.5)
 *
 * 精度: L_ZAMS 约 3%，R_ZAMS 约 1.2%（全质量范围）。
 */

#ifndef ZAMS_COEFFS_H
#define ZAMS_COEFFS_H

#ifdef __cplusplus
extern "C" {
#endif

/* 计算 ZAMS 光度 L_ZAMS (Lsun) */
double zams_lzams(double m, double z);

/* 计算 ZAMS 半径 R_ZAMS (Rsun) */
double zams_rzams(double m, double z);

/* 取 7 个 L_ZAMS 系数（按 l1..l7 顺序写入 out[0..6]） */
void zams_lzams_coeffs(double z, double out[7]);

/* 取 9 个 R_ZAMS 系数（按 r1..r9 顺序写入 out[0..8]） */
void zams_rzams_coeffs(double z, double out[9]);

#ifdef __cplusplus
}
#endif

#endif /* ZAMS_COEFFS_H */

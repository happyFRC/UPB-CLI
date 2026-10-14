/* timescales.h -- SSE 演化时标与阶段边界量
 *
 * 对应原 Python 项目中的 sse/timescales.py 的 Sequence 类。
 *
 * 给定初始质量 M 与金属度 Z，计算所有与时间 t 无关、只依赖 (M, Z) 的量：
 *   - 临界质量 Mhook / MHeF / MFGB
 *   - 特征光度 lums[0..8] = L_ZAMS, L_TMS, L_BGB, L_HeI, L_minHe, L(Mx), L_BAGB, ...
 *   - 阶段边界时间 tscls[0..15] = t_BGB, t_HeI, t_He, t_inf1, t_inf2, t_x, ...
 *   - 巨支参数 gb[0..9] = A_H, A_HHe, B, D, p, q, Mx, A_He, Mc_BGB, ...
 *   - 一些标量: tm, t_bagb, t_tp, tn, mc_du, mc_sn, rmin_ratio, mc_he_i_val, ...
 *
 * 单位: 质量 Msun, 光度 Lsun, 半径 Rsun, 时间 Myr。
 */

#ifndef TIMESCALES_H
#define TIMESCALES_H

#include "zcoeffs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* 与 Python 版一致的常量 */
#define MCH_   1.44       /* Chandrasekhar 质量 */
#define MLP_   12.0       /* R_HeI 的 HM 插值分界质量 */
#define TINY_  1.0e-14

typedef struct {
    /* --- 输入 --- */
    double mass_in;   /* 用户给的初始质量 */
    double mass;      /* 公式里用的质量，上限 100 */
    double z;

    /* --- 系数 --- */
    ZCoeffs c;

    /* --- 派生量 --- */
    double mhook;
    double mhef;
    double mfgb;
    double zeta;
    double lz;        /* log10(Z) */
    double lzd;       /* zeta + 1 */
    double x_h;       /* X = 0.76 - 3Z */
    double y_he;      /* Y = 0.24 + 2Z */

    /* --- 特征量数组 (与 Python 版下标一致) ---
     * lums[0] 恒为 0，实际用 lums[1..8]
     * tscls[0] 恒为 0，实际用 tscls[1..15]
     * gb[0]   恒为 0，实际用 gb[1..9]
     * zpars[0]恒为 0，实际用 zpars[1..14]
     */
    double lums[9];
    double tscls[16];
    double gb[10];
    double zpars[15];

    /* --- 常用标量 --- */
    double tm;          /* 主序寿命 t_MS */
    double t_bagb;      /* AGB 底端时间 */
    double t_tp;        /* TPAGB 起始时间 */
    double tn;          /* 核演化终点时间 */
    double mc_du;       /* 二 dredge-up 后的核心质量 */
    double mc_sn;       /* 超新星前最大核心质量 */
    double mc_he_i_val; /* 氦点火时的核心质量 */
    double mc_bgb_val;  /* BGB 时的核心质量 */
    double rmin_ratio;  /* 用于 LM 星 CHeB 最小半径的比例 */
} Sequence;

/* 构造。mass > 0, z > 0。
 * 成功返回 0，失败返回 -1。 */
int sequence_init(Sequence *sp, double mass, double z);

/* --- 序列相关量的函数接口 (与 Python 版 Sequence 方法一一对应) --- */

/* ZAMS 光度/半径 (Tout et al. 1996) */
double seq_lzams(const Sequence *sp, double m);
double seq_rzams(const Sequence *sp, double m);

/* 各阶段边界时间 */
double seq_t_bgb   (const Sequence *sp, double m);
double seq_t_hook  (const Sequence *sp, double m);
double seq_x_factor(const Sequence *sp);
double seq_t_ms    (const Sequence *sp, double m);

/* 特征光度 */
double seq_l_tms   (const Sequence *sp, double m);
double seq_r_tms   (const Sequence *sp, double m);
double seq_l_bgb   (const Sequence *sp, double m);

/* 巨支 L-Mc 关系 (eq. 31) 及其反函数 */
double seq_l_from_mc(const Sequence *sp, double mc, double m);
double seq_mc_from_l(const Sequence *sp, double lum, double m);

/* 巨支 L(t) / Mc(t) (eq. 33, 34) */
double seq_mc_from_t(const Sequence *sp, double t, double a_rate,
                     double tinf1, double tinf2, double tx, double m);
double seq_l_from_t (const Sequence *sp, double t, double a_rate,
                     double tinf1, double tinf2, double tx, double m);

/* 巨支半径 (eq. 37, 74) */
double seq_r_gb  (const Sequence *sp, double lum, double m);
double seq_ragbf (const Sequence *sp, double lum, double m);
double seq_r_agb (const Sequence *sp, double lum, double m);

/* CHeB 相关 */
double seq_l_he_i    (const Sequence *sp, double m);
double seq_l_he_ratio(const Sequence *sp, double m);
double seq_l_min_he  (const Sequence *sp, double m);
double seq_r_min_he  (const Sequence *sp, double m);
double seq_mu_zahb   (const Sequence *sp, double mc, double m);
double seq_l_zahb    (const Sequence *sp, double mc, double m);
double seq_r_zahb    (const Sequence *sp, double mc, double m);
double seq_t_he_ms   (const Sequence *sp, double mc);
double seq_t_he_ratio(const Sequence *sp, double m, double mc);
double seq_t_he      (const Sequence *sp, double mc, double m);
double seq_thef_rel  (const Sequence *sp, double m, double mc);
double seq_thef_im   (const Sequence *sp, double m);
double seq_tau_bl    (const Sequence *sp, double m);
double seq_f_bl      (const Sequence *sp, double m);

/* 核心质量关系 */
double seq_mc_bagb  (const Sequence *sp, double m);
double seq_mc_he_if (const Sequence *sp, double m, double mc_ref);

/* 裸氦星 (Section 6.1) */
double seq_l_zhe   (const Sequence *sp, double m);
double seq_r_zhe   (const Sequence *sp, double m);
double seq_t_hems  (const Sequence *sp, double m);
double seq_r_he_hg (const Sequence *sp, double m, double lum, double rx, double lx);
double seq_r_he_gb (const Sequence *sp, double lum);

/* L_BAGB (eq. 56) */
double seq_l_bagb(const Sequence *sp, double m);

#ifdef __cplusplus
}
#endif

#endif /* TIMESCALES_H */
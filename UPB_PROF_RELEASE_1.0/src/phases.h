/* phases.h -- 恒星属性随年龄的变化（论文第 5 节公式）
 *
 * 对应原 Python 项目中的 sse/phases.py。
 *
 * 这是 Hurley, Pols & Tout (2000) SSE 包中 hrdiag 例程的移植，加上
 * mrenv 包层估计。给定初始质量 M、金属度 Z、当前（可能被质量损失减小
 * 的）质量 mt 与年龄 aj，返回光度、半径、核心质量、核心半径、恒星类型
 * 与对流包层属性。
 *
 * 恒星类型 kw 沿用论文第 4 节：
 *   0  MS 星，M <~ 0.7（深/全对流）
 *   1  MS 星，M >~ 0.7
 *   2  赫氏空隙
 *   3  第一巨星支
 *   4  核心氦燃烧
 *   5  早期渐近巨星支
 *   6  热脉动 AGB
 *   7  裸氦星 MS
 *   8  裸氦星赫氏空隙
 *   9  裸氦星巨星支
 *   10 氦白矮星
 *   11 碳氧白矮星
 *   12 氧氖白矮星
 *   13 中子星
 *   14 黑洞
 *   15 无质量残骸
 */

#ifndef PHASES_H
#define PHASES_H

#include "timescales.h"

#ifdef __cplusplus
extern "C" {
#endif

/* 与 Python 版一致的常量 */
#define PH_LUM0     7.0e4
#define PH_KAP      (-0.5)
#define PH_AHE      4.0
#define PH_ACO      15.0
#define PH_AONE     17.0
#define PH_MXNS     3.0
#define PH_RNS      1.4e-5
#define PH_TAUMIN   5.0e-8

/* hrdiag 的输出 */
typedef struct {
    double r;         /* 半径 (Rsun) */
    double lum;       /* 光度 (Lsun) */
    int    kw;        /* 恒星类型 */
    double mc;        /* 核心质量 (Msun) */
    double rc;        /* 核心半径 (Rsun) */
    double mass;      /* 当前质量 (Msun) */
    double mt;        /* 当前总质量 (Msun) */
    double menv;      /* 对流包层质量 (Msun) */
    double renv;      /* 对流包层半径 (Rsun) */
    double k2;        /* 回转半径 */
} Point;

/* 深对流 MS 星的质量阈值 */
double mA_ms_threshold(const Sequence *sp);

/* 核心例程。
 *
 * 参数:
 *   sp   (M, Z) 对应的序列参数
 *   kw   当前恒星类型
 *   aj   当前年龄 (Myr)
 *   mt   当前总质量 (Msun)
 *   mc   上一时刻的核心质量 (用于 HG 单调性修正与核心质量测试)
 *   mass 零龄质量; 传 <0 表示用 sp->mass_in
 */
Point hrdiag(const Sequence *sp, int kw, double aj, double mt,
             double mc, double mass);

/* 对流包层 (论文第 7.2 节)。返回 menv, renv, k2。 */
void mrenv(int kw, double mass, double mt, double mc,
           double lum, double rad, double rc,
           double aj, double tm, double ltms, double lbgb, double lhei,
           double rzams, double rtms, double rg,
           double *menv_out, double *renv_out, double *k2_out);

#ifdef __cplusplus
}
#endif

#endif /* PHASES_H */

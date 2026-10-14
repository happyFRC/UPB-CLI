/* evolve.h -- 单星演化驱动
 *
 * 对应原 Python 项目中的 sse/evolve.py。
 *
 * 让一颗恒星按照 Hurley, Pols & Tout (2000) 的阶段推进：从 ZAMS 经过主序、
 * 赫氏空隙、第一巨星支、核心氦燃烧、渐近巨星支，到残骸阶段，期间按照论文
 * 第 7 节叠加质量损失。
 *
 * 默认步长为主序寿命的 1/1000，所有阶段统一使用；在需要时缩短以避免跨过
 * 阶段边界。按照论文第 7.1 节，时标、光度和核心质量用的序列参数取自有效
 * 初始质量 M0；M0 只在主序和赫氏空隙上重新锚定到当前质量，以及恒星被剥离
 * 成裸氦星时。半径公式和质量损失率使用当前质量。
 */

#ifndef EVOLVE_H
#define EVOLVE_H

#include "timescales.h"
#include "phases.h"

#ifdef __cplusplus
extern "C" {
#endif

/* 残骸类型 (不再进行核燃烧) */
#define EV_REMNANT_MIN 10
#define EV_REMNANT_MAX 15

/* 演化轨道上的一个记录点 */
typedef struct {
    double time;            /* 年龄 (Myr) */
    double luminosity;      /* 光度 (Lsun) */
    double radius;          /* 半径 (Rsun) */
    double temperature;     /* 有效温度 (K) */
    double mass;            /* 当前总质量 (Msun) */
    double core_mass;       /* 核心质量 (Msun) */
    int    phase;           /* SSE 恒星类型 */
    char   spectral_type[2];/* 光谱型单字符字符串 */
    double core_temperature;
    double menv;
    double renv;
} TrackPoint;

/* 一次演化运行的结果 */
typedef struct {
    double mass;
    double metallicity;

    TrackPoint *points;     /* 记录点数组 (若 record=0 则为 NULL) */
    int    n_points;
    int    cap_points;

    double t_ms;
    double t_bgb;
    double t_he_i;
    double t_he;
    double t_bagb;
    double t_tp;
    double t_nuc;
    double l_zams;
    double r_zams;
    double final_mass;

    int    success;
    char   error[256];
    int    steps;
    int    stagnated;       /* 核心质量/光度停止变化 */
    TrackPoint last;        /* 最后一个记录点 (即使 record=0 也有) */
} Track;

/* 演化选项 */
typedef struct {
    double end_time;        /* 停止年龄 (Myr); <=0 表示演化到残骸 */
    double step_size;       /* 固定步长 (Myr); <=0 表示默认 t_MS/1000 */
    int    max_steps;       /* 步数上限 */
    int    record;          /* 是否记录每个点 */
    int    mass_loss;       /* 是否应用质量损失 */
    int    stop_at_remnant; /* 是否在残骸处停止 */
} EvolveOptions;

/* 初始化 Track（清零、分配默认值） */
void track_init(Track *tr, double mass, double metallicity);

/* 释放 Track 内部的 points 数组 */
void track_free(Track *tr);

/* 取最后一个记录点；track 为空时返回 NULL */
const TrackPoint *track_final(const Track *tr);

/* 主演化函数。
 *
 * 返回 0 成功；-1 参数错误。
 */
int evolve(double mass, double metallicity,
           const EvolveOptions *opts, Track *out);

#ifdef __cplusplus
}
#endif

#endif /* EVOLVE_H */

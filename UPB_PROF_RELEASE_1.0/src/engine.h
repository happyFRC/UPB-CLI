/* engine.h -- 面向用户的单星演化引擎
 *
 * 对应原 Python 项目中的 engine.py。
 *
 * 这是命令行界面使用的入口。它包装 sse 包并暴露一个小巧、便于 JSON 化的 API：
 *   - StellarEvolutionEngine (evolve / get_main_sequence_lifetime / ...)
 *   - StellarState / EvolutionResult 结果容器
 *   - evolve_star / evolve_star_detailed 便捷函数
 *
 * 所有物理都在底层模块里；这里只做输入检查、单位/结构转换和计时。
 */

#ifndef ENGINE_H
#define ENGINE_H

#include "evolve.h"
#include "timescales.h"

#ifdef __cplusplus
extern "C" {
#endif

/* 支持的质量/金属度范围 */
#define ENG_MASS_MIN  0.077
#define ENG_MASS_MAX  100.0
#define ENG_Z_MIN     0.0001
#define ENG_Z_MAX     0.03

/* 光谱型：从有效温度得到的传统分类 */
typedef enum {
    SPECTRAL_O = 0,
    SPECTRAL_B,
    SPECTRAL_A,
    SPECTRAL_F,
    SPECTRAL_G,
    SPECTRAL_K,
    SPECTRAL_M
} SpectralType;

/* 演化阶段，使用论文第 4 节的整数类型 */
typedef enum {
    PHASE_MS_CONVECTIVE = 0,
    PHASE_MS            = 1,
    PHASE_HG            = 2,
    PHASE_GB            = 3,
    PHASE_CHeB          = 4,
    PHASE_EAGB          = 5,
    PHASE_TPAGB         = 6,
    PHASE_HeMS          = 7,
    PHASE_HeHG          = 8,
    PHASE_HeGB          = 9,
    PHASE_He_WD         = 10,
    PHASE_CO_WD         = 11,
    PHASE_ONe_WD        = 12,
    PHASE_NS            = 13,
    PHASE_BH            = 14,
    PHASE_MASSLESS      = 15
} EvolutionPhase;

/* 单个恒星快照 */
typedef struct {
    double time;             /* 年龄 (Myr) */
    double luminosity;       /* 光度 (Lsun) */
    double radius;           /* 半径 (Rsun) */
    double temperature;      /* 有效温度 (K) */
    char   spectral_type[2]; /* 单字符光谱型 */
    double core_h_fraction;  /* 包层氢质量分数 X */
    double core_he_fraction; /* 包层氦质量分数 Y */
    double core_temperature; /* 由核心质量估计的中心温度 (K) */
    double mass;             /* 当前总质量 (Msun) */
    double core_mass;        /* 核心质量 (Msun) */
    int    phase;            /* 恒星类型 */
    char   phase_name[40];   /* 中文阶段名（UTF-8） */
    double envelope_mass;    /* 包层质量 (Msun) */
    double envelope_radius;  /* 包层半径 (Rsun) */
} StellarState;

/* 一次演化的结果 */
typedef struct {
    /* 参数 */
    double mass;
    double metallicity;
    double end_time;
    double step_size;

    /* 结果 */
    StellarState final_state;
    StellarState *evolution_track;  /* 动态数组 */
    int    n_track;
    int    cap_track;
    int    total_steps;
    double computation_time_ms;

    int    success;
    char   error_message[256];
    int    stagnated;

    /* 阶段边界时间 (Myr) */
    double t_ms, t_bgb, t_he_i, t_he, t_bagb, t_tp;
} EvolutionResult;

/* ------------------------------------------------------------- 引擎 API */

/* 单星演化。opts 可为 NULL，使用默认选项。
 * 成功返回 0，失败返回 -1；失败时 result->error_message 有说明。 */
int engine_evolve(double mass, double metallicity,
                  double end_time, double step_size,
                  int track_history, int mass_loss,
                  EvolutionResult *result);

/* 主序寿命 t_MS (Myr)。参数越界返回 0。 */
double engine_get_main_sequence_lifetime(double mass, double metallicity);

/* BGB 寿命 t_BGB (Myr)。参数越界返回 0。 */
double engine_get_bgb_lifetime(double mass, double metallicity);

/* 所有阶段边界时间与特征光度。
 * 用一个简单的结构体带出（对应 Python 的 dict）。 */
typedef struct {
    double t_ms, t_bgb, t_he_i, t_he, t_bagb, t_tp, t_nuclear;
    double l_zams, l_tms, l_bgb, l_he_i, l_min_he, l_bagb;
    double r_zams, r_tms;
    double m_hook, m_he_f, m_fgb;
} PhaseTimes;

int engine_get_phase_times(double mass, double metallicity, PhaseTimes *out);

/* 批量演化：对 params_list 里的每一项调用 engine_evolve。
 *
 * 为了贴近 Python 的 evolve_batch，这里提供一个简化版本：
 * 调用方给出 mass/z/end/step 的数组，结果数组由调用方预先分配。
 * 成功处理 n 个返回 n（即使部分失败）。 */
int engine_evolve_batch(const double *masses, const double *metallicities,
                        const double *end_times, const double *step_sizes,
                        int n, int track_history,
                        EvolutionResult *results);

/* 释放 EvolutionResult 内部的动态数组 */
void evolution_result_free(EvolutionResult *r);

/* 便捷函数：演化为摘要结构体（不保留轨道） */
int evolve_star(double mass, double metallicity, double end_time,
                double step_size, EvolutionResult *result);

/* 便捷函数：演化为完整结果（保留轨道） */
int evolve_star_detailed(double mass, double metallicity, double end_time,
                         double step_size, EvolutionResult *result);

/* 光谱型 / 阶段名辅助 */
const char *phase_name_cn(int phase);
const char *phase_name_en(int phase);

#ifdef __cplusplus
}
#endif

#endif /* ENGINE_H */

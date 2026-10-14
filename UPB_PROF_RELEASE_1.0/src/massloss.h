/* massloss.h -- 质量损失与观测换算
 *
 * 对应原 Python 项目中的 sse/massloss.py。
 *
 * Hurley, Pols & Tout (2000) 第 7 节给出了与拟合公式配套使用的完整质量损失
 * 算法。详细模型本身是无质量损失的，所以这里在演化轨道之上叠加质量损失。
 *
 * 三个处方叠加使用：
 *   - Nieuwenhuijzen & de Jager (1990): 全赫罗图的高光度星风
 *   - Reimers (1975):                  巨星 (kw 2..9)
 *   - Vassiliadis & Wood (1993):       脉动 AGB (kw 5, 6)
 * 再加上小包层被剥离恒星的减弱 Wolf-Rayet 星风。
 *
 * 注意: 论文第 7 节的处方都以 Msun/yr 为单位，本程序内部时间用 Myr，
 *       所以最后乘 1.0e6 转换成 Msun/Myr。
 */

#ifndef MASSLOSS_H
#define MASSLOSS_H

#ifdef __cplusplus
extern "C" {
#endif

/* 与 Python 版一致的常量 */
#define ML_L0     7.0e4
#define ML_KAP    (-0.5)
#define ML_A_HE   4.0
#define ML_A_CO   15.0
#define ML_A_ONe  17.0
#define ML_B_WIND 0.0
#define ML_NETA   0.5

/* 1 Myr = 1.0e6 yr，用于把 Msun/yr 转成 Msun/Myr */
#define ML_YR_PER_MYR 1.0e6

/* 总质量损失率，单位 Msun/Myr。
 *
 * 参数:
 *   kw        恒星类型 (见 phases.h)
 *   lum, r    L (Lsun), R (Rsun)
 *   mt, mc    当前总质量、核心质量 (Msun)
 *   z         金属度
 *   rl        Roche-lobe 半径 (只在潮汐增强风里用，单星时传 0)
 */
double mass_loss_rate(int kw, double lum, double r,
                      double mt, double mc, double z, double rl);

/* 从 L (Lsun) 和 R (Rsun) 求有效温度 (K) */
double temperature_from_lr(double lum, double r);

/* 从有效温度判断光谱型 (单字符字符串，如 "G") */
const char *spectral_type_from_temperature(double t_eff);

/* 近似热改正 BC (等值)，仅用于报告，不影响演化 */
double bolometric_correction(double t_eff);

#ifdef __cplusplus
}
#endif

#endif /* MASSLOSS_H */

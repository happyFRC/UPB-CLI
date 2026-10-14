/* cli.h -- 命令行界面共用的结构与常量
 *
 * 对应原 Python 项目中的 main.py（顶部常量与单位系统部分）。
 */

#ifndef CLI_H
#define CLI_H

#include <stdio.h>

/* 版本信息 */
#define CLI_APP_NAME       "UPB-CLI"
#define CLI_VERSION_STRING "UPB-PROFESSIONAL RELEASE 1.0"
#define CLI_VERSION_CODE   "1.0"
#define CLI_RELEASE_DATE   "2026-10-14"
#define CLI_PROMPT         "UPB-CLI>"
#define CLI_MAX_BATCH      100
#define CLI_MAX_STEPS      200000

/* 单位换算常量（cgs / SI 混合，按需换算） */
#define U_MSUN_KG      1.98847e30
#define U_MJUP_KG      1.8981246e27
#define U_MEARTH_KG    5.9722e24
#define U_LSUN_W       3.828e26
#define U_RSUN_M       6.957e8
#define U_REARTH_M     6.371e6
#define U_AU_M         1.495978707e11
#define U_RSTEPHENSON  2150.0
#define U_LR136A1      8.0e6

/* --- 单位系统 --- */
typedef struct {
    int mass;          /* 1=MSun 2=MJupiter 3=MEarth 4=kg */
    int luminosity;    /* 1=LSun 2=LR136a1 3=W */
    int radius;        /* 1=RSun 2=RStephenson 3=AU 4=REarth 5=km */
    int temperature;   /* 1=K 2=C 3=F */
} Units;

/* 单位换算 */
double units_mass_from(const Units *u, double v);
double units_lum_from (const Units *u, double v);
double units_rad_from (const Units *u, double v);
double units_temp_from(const Units *u, double v);

/* 单位符号 */
const char *units_mass_symbol(const Units *u);
const char *units_lum_symbol (const Units *u);
const char *units_rad_symbol (const Units *u);
const char *units_temp_symbol(const Units *u);

/* 重置所有单位代码为 1 */
void units_reset(Units *u);

/* 把当前单位设置写入 upb_config.json；失败时静默 */
void units_save(const Units *u);

/* 从 upb_config.json 读回上次的单位设置（文件不存在或损坏则保持默认） */
void units_load(Units *u);

/* 打印一行摘要 */
void units_summary(const Units *u, char *buf, int buflen);

/* --- 路径/文件工具 --- */

/* 取 results 目录（若不存在则创建）。返回缓冲区指针。 */
const char *cli_results_dir(void);

/* 把用户给的输出文件名清洗成安全文件名（去掉路径/非法字符，确保 .txt） */
int cli_safe_filename(const char *name, char *out, int outlen);

/* 生成自动文件名: star_M{mass}_Z{z}.txt */
void cli_auto_filename(double mass, double z, char *out, int outlen);

/* --- 常用格式化 --- */
void cli_fmt(double value, int sig, char *out, int outlen);

/* 阶段名（中文） */
const char *cli_phase_cn(int kw);
const char *cli_phase_en(int kw);

/* 元素丰度: Y = 0.24 + 2Z, X = 0.76 - 3Z */
void cli_composition(double z, double *x_out, double *y_out);

#endif /* CLI_H */
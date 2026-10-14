/* cli.c -- 命令行界面共用的工具实现（单位、路径、格式化）
 *
 * 对应原 Python 项目中的 main.py 的 Units / 小工具部分。
 */

#include "cli.h"
#include "engine.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _WIN32
#  include <direct.h>
#  define CLI_MKDIR(p) _mkdir(p)
#else
#  include <unistd.h>
#  define CLI_MKDIR(p) mkdir((p), 0755)
#endif

/* ------------------------------------------------------------------ 单位 */

typedef struct {
    int code;
    const char *symbol;
    double scale;      /* 在内部单位中的大小 */
} UnitDef;

static const UnitDef MASS_UNITS[] = {
    { 0, "",       0.0 },                               /* 占位，下标 1 起 */
    { 1, "MSun",   1.0 },
    { 2, "MJupiter", U_MJUP_KG / U_MSUN_KG },
    { 3, "MEarth",  U_MEARTH_KG / U_MSUN_KG },
    { 4, "kg",      1.0 / U_MSUN_KG },
};
static const int MASS_UNITS_N = 5;

static const UnitDef LUM_UNITS[] = {
    { 0, "",     0.0 },
    { 1, "LSun", 1.0 },
    { 2, "LR136a1", U_LR136A1 },
    { 3, "W",    1.0 / U_LSUN_W },
};
static const int LUM_UNITS_N = 4;

static const UnitDef RAD_UNITS[] = {
    { 0, "",     0.0 },
    { 1, "RSun", 1.0 },
    { 2, "RStephenson_218", U_RSTEPHENSON },
    { 3, "AU",   U_AU_M / U_RSUN_M },
    { 4, "REarth", U_REARTH_M / U_RSUN_M },
    { 5, "km",   1.0 / (U_RSUN_M / 1000.0) },
};
static const int RAD_UNITS_N = 6;

static const UnitDef TEMP_UNITS[] = {
    { 0, "", 0.0 },
    { 1, "K", 1.0 },
    { 2, "C", 1.0 },   /* 仿射换算，见 units_temp_from */
    { 3, "F", 1.0 },
};
static const int TEMP_UNITS_N = 4;

static const UnitDef *find_def(const UnitDef *table, int n, int code)
{
    int i;
    for (i = 1; i < n; i++)
        if (table[i].code == code) return &table[i];
    return NULL;
}

double units_mass_from(const Units *u, double v)
{
    const UnitDef *d = find_def(MASS_UNITS, MASS_UNITS_N, u->mass);
    return (d != NULL) ? v / d->scale : v;
}
double units_lum_from(const Units *u, double v)
{
    const UnitDef *d = find_def(LUM_UNITS, LUM_UNITS_N, u->luminosity);
    return (d != NULL) ? v / d->scale : v;
}
double units_rad_from(const Units *u, double v)
{
    const UnitDef *d = find_def(RAD_UNITS, RAD_UNITS_N, u->radius);
    return (d != NULL) ? v / d->scale : v;
}
double units_temp_from(const Units *u, double v)
{
    if (u->temperature == 2) return v - 273.15;
    if (u->temperature == 3) return (v - 273.15) * 9.0 / 5.0 + 32.0;
    return v;
}

const char *units_mass_symbol(const Units *u)
{
    const UnitDef *d = find_def(MASS_UNITS, MASS_UNITS_N, u->mass);
    return (d != NULL) ? d->symbol : "?";
}
const char *units_lum_symbol(const Units *u)
{
    const UnitDef *d = find_def(LUM_UNITS, LUM_UNITS_N, u->luminosity);
    return (d != NULL) ? d->symbol : "?";
}
const char *units_rad_symbol(const Units *u)
{
    const UnitDef *d = find_def(RAD_UNITS, RAD_UNITS_N, u->radius);
    return (d != NULL) ? d->symbol : "?";
}
const char *units_temp_symbol(const Units *u)
{
    const UnitDef *d = find_def(TEMP_UNITS, TEMP_UNITS_N, u->temperature);
    return (d != NULL) ? d->symbol : "?";
}

void units_reset(Units *u)
{
    u->mass = u->luminosity = u->radius = u->temperature = 1;
}

void units_summary(const Units *u, char *buf, int buflen)
{
    snprintf(buf, buflen,
             "mass=%d(%s) luminosity=%d(%s) radius=%d(%s) temperature=%d(%s)",
             u->mass, units_mass_symbol(u),
             u->luminosity, units_lum_symbol(u),
             u->radius, units_rad_symbol(u),
             u->temperature, units_temp_symbol(u));
}

/* ------------------------------------------------------------------ 路径 */

static char g_results_dir[512] = {0};

const char *cli_results_dir(void)
{
    if (g_results_dir[0] == '\0') {
        /* 与可执行文件同目录下的 results/ */
        snprintf(g_results_dir, sizeof(g_results_dir), "results");
    }
    CLI_MKDIR(g_results_dir);
    return g_results_dir;
}

int cli_safe_filename(const char *name, char *out, int outlen)
{
    int i, j = 0;
    int len;
    if (name == NULL || out == NULL) return -1;
    /* 跳过前导空白 */
    while (*name == ' ' || *name == '\t') name++;
    len = (int)strlen(name);
    /* 去掉尾随空白 */
    while (len > 0 && (name[len - 1] == ' ' || name[len - 1] == '\t'))
        len--;
    if (len == 0) return -1;

    for (i = 0; i < len && j < outlen - 6; i++) {
        char c = name[i];
        if (c == '\\' || c == '/' || c == ':' || c == '*' ||
            c == '?' || c == '"' || c == '<' || c == '>' ||
            c == '|' || (unsigned char)c < 0x20) {
            c = '_';
        }
        out[j++] = c;
    }
    out[j] = '\0';
    /* 去掉尾随的点 */
    while (j > 0 && out[j - 1] == '.') { out[--j] = '\0'; }
    if (j == 0) return -1;

    /* 确保 .txt 后缀 */
    if (j < 4 || strcmp(out + j - 4, ".txt") != 0) {
        if (j + 4 >= outlen) return -1;
        strcpy(out + j, ".txt");
        j += 4;
    }
    if (j > 120) return -1;
    return 0;
}

void cli_auto_filename(double mass, double z, char *out, int outlen)
{
    char raw[128];
    snprintf(raw, sizeof(raw), "star_M%.4g_Z%.4g", mass, z);
    if (cli_safe_filename(raw, out, outlen) != 0) {
        snprintf(out, outlen, "star.txt");
    }
}

/* ------------------------------------------------------------------ 格式 */

void cli_fmt(double value, int sig, char *out, int outlen)
{
    if (value == 0.0) {
        snprintf(out, outlen, "0");
        return;
    }
    if (!isfinite(value)) {
        snprintf(out, outlen, "n/a");
        return;
    }
    snprintf(out, outlen, "%.*g", sig, value);
}

/* ---------------------------------------------------------------- 阶段名 */

const char *cli_phase_cn(int kw)
{
    static const char *tbl[16] = {
        "主序(全对流)", "主序", "赫氏空隙", "第一巨星支",
        "核心氦燃烧", "早期AGB", "热脉动AGB",
        "氦星主序", "氦星赫氏空隙", "氦星巨星支",
        "氦白矮星", "碳氧白矮星", "氧氖白矮星",
        "中子星", "黑洞", "无质量残骸"
    };
    if (kw < 0 || kw > 15) return "未知";
    return tbl[kw];
}

const char *cli_phase_en(int kw)
{
    return phase_name_en(kw);
}

void cli_composition(double z, double *x_out, double *y_out)
{
    double y = 0.24 + 2.0 * z;
    double x = 0.76 - 3.0 * z;
    if (x_out != NULL) *x_out = x;
    if (y_out != NULL) *y_out = y;
}
#define CLI_CONFIG_PATH "upb_config.json"

void units_save(const Units *u)
{
    FILE *fp = fopen(CLI_CONFIG_PATH, "w");
    if (fp == NULL) return;   /* 写不了就算了，不报错 */
    fprintf(fp,
            "{\n"
            "  \"mass\": %d,\n"
            "  \"luminosity\": %d,\n"
            "  \"radius\": %d,\n"
            "  \"temperature\": %d\n"
            "}\n",
            u->mass, u->luminosity, u->radius, u->temperature);
    fclose(fp);
}

void units_load(Units *u)
{
    FILE *fp = fopen(CLI_CONFIG_PATH, "r");
    int m, l, r, t;
    if (fp == NULL) return;   /* 文件不存在就保持默认 */
    /* 只抓 4 个整数，忽略所有其它字符 */
    if (fscanf(fp, "%*[^{]%*c%*[^0-9]%d%*[^0-9]%d%*[^0-9]%d%*[^0-9]%d",
               &m, &l, &r, &t) == 4) {
        if (m >= 1 && m <= 4)          u->mass        = m;
        if (l >= 1 && l <= 3)          u->luminosity  = l;
        if (r >= 1 && r <= 5)          u->radius      = r;
        if (t >= 1 && t <= 3)          u->temperature = t;
    }
    fclose(fp);
}

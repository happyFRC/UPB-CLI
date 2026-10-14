/* main.c -- UPB-CLI 命令行程序 (C 版本)
 *
 * 对应原 Python 项目中的 main.py。
 *
 * 用法:
 *   upb_cli                         进入交互模式
 *   upb_cli help                    执行一条指令后退出
 *   upb_cli evolution_single 1.0 0.02 0 0
 *   upb_cli "unit mass 2; evolution_single 1.0 0.02 0 0"
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#include "cli.h"
#include "engine.h"

/* ------------------------------------------------------------------ 封面 */

/* ------------------------------------------------------------------ 封面 */

static void print_cover(void)
{
    int width = 78;
    int i;
    char rule[128];
    for (i = 0; i < width; i++) rule[i] = '=';
    rule[width] = '\0';

    printf("%s\n", rule);
    printf("%s\n", "  ★  U P B _ C L I  恒星演化计算器  ★");
    printf("  %s\n", CLI_VERSION_STRING);
    printf("  发布日期: %s\n", CLI_RELEASE_DATE);
    printf("%s\n", rule);
    printf("  模型: Hurley, Pols & Tout (2000) MNRAS 315, 543 单星解析演化公式\n");
    printf("  指令: 输入 help 查看全部指令, exit / quit / q 退出\n");
    printf("%s\n\n", rule);
}

/* ------------------------------------------------------------ 指令解析 */

/* 按分号拆成多条指令（括号内的分号不分）。返回拆分出的条数。 */
static int split_statements(const char *line, char parts[][512], int max_parts)
{
    int n = 0, depth = 0, j = 0;
    const char *p = line;
    for (; *p; p++) {
        char c = *p;
        if (c == '(') depth++;
        else if (c == ')') { if (depth > 0) depth--; }
        if (c == ';' && depth == 0) {
            if (j > 0 && n < max_parts) {
                parts[n][j] = '\0';
                /* 去掉两端空白 */
                n++;
            }
            j = 0;
        } else {
            if (j < 511) parts[n][j++] = c;
            else if (n < max_parts - 1) { parts[n][j] = '\0'; n++; j = 0; }
        }
    }
    if (j > 0 && n < max_parts) {
        parts[n][j] = '\0';
        n++;
    }
    /* 去掉两端空白 */
    {
        int i;
        for (i = 0; i < n; i++) {
            char *s = parts[i];
            char *e;
            while (*s == ' ' || *s == '\t') s++;
            e = s + strlen(s);
            while (e > s && (e[-1] == ' ' || e[-1] == '\t')) *--e = '\0';
            if (s != parts[i]) memmove(parts[i], s, strlen(s) + 1);
        }
    }
    return n;
}

/* 把一行指令按空白和括号拆成 tokens 和 groups。
 *
 * tokens 里不含括号内容；groups 是每个 (...) 里的内容。
 * 返回: tokens 的数量。
 */
static int tokenize(const char *line,
                    char tokens[][64], int max_tokens,
                    char groups[][256], int max_groups,
                    int *n_groups_out)
{
    int n_tok = 0, n_grp = 0;
    const char *p = line;
    char buf[512];
    int blen = 0;

    /* 提取 groups: 扫描 "(...)" */
    while (*p) {
        if (*p == '(') {
            const char *q = p + 1;
            int glen = 0;
            while (*q && *q != ')' && glen < 255) {
                groups[n_grp][glen++] = *q++;
            }
            groups[n_grp][glen] = '\0';
            n_grp++;
            if (n_grp >= max_groups) break;
            if (*q == ')') q++;
            p = q;
        } else {
            p++;
        }
    }
    *n_groups_out = n_grp;

    /* 提取 tokens: 把 "(...)" 替换成空格后按空白切分 */
    p = line;
    while (*p) {
        if (*p == '(') {
            while (*p && *p != ')') p++;
            if (*p == ')') p++;
            if (blen > 0) {
                buf[blen] = '\0';
                if (n_tok < max_tokens) {
                    strncpy(tokens[n_tok], buf, 63);
                    tokens[n_tok][63] = '\0';
                    n_tok++;
                }
                blen = 0;
            }
        } else if (*p == ' ' || *p == '\t') {
            if (blen > 0) {
                buf[blen] = '\0';
                if (n_tok < max_tokens) {
                    strncpy(tokens[n_tok], buf, 63);
                    tokens[n_tok][63] = '\0';
                    n_tok++;
                }
                blen = 0;
            }
            p++;
        } else {
            if (blen < 511) buf[blen++] = *p;
            p++;
        }
    }
    if (blen > 0 && n_tok < max_tokens) {
        buf[blen] = '\0';
        strncpy(tokens[n_tok], buf, 63);
        tokens[n_tok][63] = '\0';
        n_tok++;
    }
    return n_tok;
}

/* 解析浮点数 */
static int parse_float(const char *s, double *out)
{
    char *end;
    double v;
    if (s == NULL || out == NULL) return -1;
    v = strtod(s, &end);
    if (end == s || *end != '\0') return -1;
    *out = v;
    return 0;
}

/* 解析整数 */
static int parse_int(const char *s, int *out)
{
    char *end;
    long v;
    if (s == NULL || out == NULL) return -1;
    v = strtol(s, &end, 10);
    if (end == s || *end != '\0') return -1;
    *out = (int)v;
    return 0;
}

/* ------------------------------------------------------------ 交互辅助 */

/* 询问一个问题，带默认值；返回用户输入（可能为空字符串表示默认）。 */
static void ask(const char *prompt, const char *default_value,
                char *out, int outlen)
{
    char line[256];
    printf("%s [%s]: ", prompt, default_value);
    fflush(stdout);
    if (fgets(line, sizeof(line), stdin) == NULL) {
        strncpy(out, default_value, outlen - 1);
        out[outlen - 1] = '\0';
        return;
    }
    /* 去掉换行 */
    line[strcspn(line, "\r\n")] = '\0';
    /* 去空白 */
    {
        char *s = line;
        char *e;
        while (*s == ' ' || *s == '\t') s++;
        e = s + strlen(s);
        while (e > s && (e[-1] == ' ' || e[-1] == '\t')) *--e = '\0';
        if (*s == '\0')
            strncpy(out, default_value, outlen - 1);
        else
            strncpy(out, s, outlen - 1);
        out[outlen - 1] = '\0';
    }
}

/* ============================================================ 指令实现 */

/* --- help --- */
static int cmd_help(void)
{
    printf(
"恒星演化命令行程序 (C 版本) 指令列表\n"
"==============================================================================\n"
"  evolution_single [mass] [metallicity] [step_size] [end_time] [filename]\n"
"      跑单个恒星演化\n"
"      参数: 质量(%.3f-%.0f) 金属度(%.4f-%.2f) 步长(Myr) 终点(Myr) 文件名(可选)\n"
"      示例: evolution_single 1.0 0.02 0 0\n"
"      示例: evolution_single 5.0 0.02 0.1 120 my_star\n"
"\n"
"  evolution_multiple (m,z,step,end[,name]) ...\n"
"      批量演化 (最多 %d 个)\n"
"      示例: evolution_multiple (1.0,0.02,0,0) (2.0,0.02,0,0)\n"
"\n"
"  phases [mass] [metallicity]\n"
"      显示各演化阶段的时间与特征光度\n"
"      示例: phases 1.0 0.02\n"
"\n"
"  results                列出结果文件\n"
"  results open           打开结果文件夹\n"
"  remove [filename] [--yes]   删除某个结果文件\n"
"\n"
"  unit [quantity] [code] 查看或设置输出单位\n"
"      mass        1=MSun 2=MJupiter 3=MEarth 4=kg\n"
"      luminosity  1=LSun 2=LR136a1 3=W\n"
"      radius      1=RSun 2=RStephenson_218 3=AU 4=REarth 5=km\n"
"      temperature 1=K 2=C 3=F\n"
"      unit reset 恢复所有单位代码为 1\n"
"\n"
"  help / version / about / exit / quit / q\n"
"==============================================================================\n",
        ENG_MASS_MIN, ENG_MASS_MAX, ENG_Z_MIN, ENG_Z_MAX, CLI_MAX_BATCH);
    return 0;
}

/* --- version --- */
static int cmd_version(void)
{
    printf("%s  %s\n", CLI_APP_NAME, CLI_VERSION_STRING);
    printf("  版本号    : %s\n", CLI_VERSION_CODE);
    printf("  发布日期  : %s\n", CLI_RELEASE_DATE);
    printf("  物理模型  : Hurley, Pols & Tout (2000) MNRAS 315, 543 (SSE)\n");
    printf("  语言/标准 : C89/C99 (MinGW)\n");
    printf("  结果目录  : %s\n", cli_results_dir());
    return 0;
}

/* --- about --- */
static int cmd_about(void)
{
    printf("关于演化模型\n"
"------------------------------------------------------------------------------\n"
"物理模型\n"
"  单星演化采用 Hurley, Pols & Tout (2000, MNRAS 315, 543) 的解析拟合公式。\n"
"  零龄主序 (ZAMS) 的 L、R 取自 Tout et al. (1996) 的拟合。\n"
"\n"
"演化阶段 (论文第 4 节的恒星类型)\n"
"  0/1 主序 -> 2 赫氏空隙 -> 3 第一巨星支 ->\n"
"  4 核心氦燃烧 -> 5 早期渐近巨星支 -> 6 热脉动 AGB\n"
"  剥离后成为裸氦星 7/8/9; 最终形成残骸:\n"
"  10 氦白矮星 / 11 碳氧白矮星 / 12 氧氖白矮星 / 13 中子星 /\n"
"  14 黑洞 / 15 无质量残骸\n"
"\n"
"质量损失 (论文第 7 节)\n"
"  Nieuwenhuijzen & de Jager (1990) + Reimers (1975) +\n"
"  Vassiliadis & Wood (1993) + 小包层类 Wolf-Rayet 星风\n"
"\n"
"适用范围\n"
"  质量 %.3f - %.0f Msun, 金属度 Z = %.4f - %.2f\n"
"  元素丰度按 Y = 0.24 + 2Z, X = 0.76 - 3Z\n"
"------------------------------------------------------------------------------\n",
        ENG_MASS_MIN, ENG_MASS_MAX, ENG_Z_MIN, ENG_Z_MAX);
    return 0;
}

/* --- unit --- */
static int cmd_unit(Units *u, char tokens[][64], int n_tok)
{
    int code;
    if (n_tok < 2) {
        char sum[256];
        units_summary(u, sum, sizeof(sum));
        printf("当前输出单位: %s\n", sum);
        return 0;
    }
    if (strcmp(tokens[1], "reset") == 0 || strcmp(tokens[1], "clear") == 0) {
        units_reset(u);
        units_save(u);
        printf("[ok] 所有单位代码已恢复为 1\n");
        return 0;
    }
    if (n_tok < 3) {
        printf("[错误] 用法: unit [quantity] [code]\n");
        return 1;
    }
    if (parse_int(tokens[2], &code) != 0) {
        printf("[错误] 单位代码必须是整数\n");
        return 1;
    }
    if (strcmp(tokens[1], "mass") == 0 && code >= 1 && code <= 4) {
        u->mass = code;
    } else if ((strcmp(tokens[1], "luminosity") == 0 ||
                strcmp(tokens[1], "lum") == 0) && code >= 1 && code <= 3) {
        u->luminosity = code;
    } else if (strcmp(tokens[1], "radius") == 0 && code >= 1 && code <= 5) {
        u->radius = code;
    } else if ((strcmp(tokens[1], "temperature") == 0 ||
                strcmp(tokens[1], "temp") == 0) && code >= 1 && code <= 3) {
        u->temperature = code;
    } else {
        printf("[错误] 未知单位类别或代码: %s %d\n", tokens[1], code);
        return 1;
    }
    units_save(u);
    printf("[ok] 已更新 %s = %d\n", tokens[1], code);
    return 0;
}

/* --- phases --- */
static int cmd_phases(double mass, double z)
{
    PhaseTimes ph;
    double x, y;
    if (engine_get_phase_times(mass, z, &ph) != 0) {
        printf("[错误] 参数越界或计算失败\n");
        return 1;
    }
    cli_composition(z, &x, &y);
    printf("恒星参数: M = %g Msun, Z = %g  (X=%.4f, Y=%.4f)\n",
           mass, z, x, y);
    printf("主序寿命 t_MS = %.2f Myr\n\n", ph.t_ms);
    printf("演化阶段 (时间单位: Myr)\n");
    printf("--------------------------------------------------------------------------\n");
    printf("%-22s%14s%14s%14s\n", "阶段", "开始", "结束", "持续时间");
    printf("--------------------------------------------------------------------------\n");
    printf("%-22s%14.2f%14.2f%14.2f\n", "主序 (MS)", 0.0, ph.t_ms, ph.t_ms);
    printf("%-22s%14.2f%14.2f%14.2f\n", "赫氏空隙 (HG)", ph.t_ms, ph.t_bgb,
           ph.t_bgb - ph.t_ms);
    printf("%-22s%14.2f%14.2f%14.2f\n", "第一巨星支 (GB)", ph.t_bgb, ph.t_he_i,
           ph.t_he_i - ph.t_bgb);
    printf("%-22s%14.2f%14.2f%14.2f\n", "核心氦燃烧 (CHeB)", ph.t_he_i,
           ph.t_he_i + ph.t_he, ph.t_he);
    printf("%-22s%14.2f%14.2f%14.2f\n", "早期 AGB (EAGB)",
           ph.t_he_i + ph.t_he, ph.t_tp, ph.t_tp - ph.t_he_i - ph.t_he);
    printf("%-22s%14.2f%14.2f%14.2f\n", "热脉动 AGB (TPAGB)",
           ph.t_tp, ph.t_nuclear, ph.t_nuclear - ph.t_tp);
    printf("--------------------------------------------------------------------------\n");
    printf("核演化终点 t_nuc = %.2f Myr\n", ph.t_nuclear);
    printf("\n特征光度 (单位: Lsun)\n");
    printf("  L_ZAMS  = %.6g\n", ph.l_zams);
    printf("  L_TMS   = %.6g\n", ph.l_tms);
    printf("  L_BGB   = %.6g\n", ph.l_bgb);
    printf("  L_HeI   = %.6g\n", ph.l_he_i);
    printf("  L_BAGB  = %.6g\n", ph.l_bagb);
    printf("  临界质量: Mhook=%.4f  MHeF=%.4f  MFGB=%.4f Msun\n",
           ph.m_hook, ph.m_he_f, ph.m_fgb);
    return 0;
}

/* --- results --- */

/* 简单统计文件行数（不含注释与表头） */
static int count_rows(const char *path)
{
    FILE *fp = fopen(path, "r");
    int rows = 0;
    char line[2048];
    if (fp == NULL) return -1;
    while (fgets(line, sizeof(line), fp) != NULL) {
        char *s = line;
        while (*s == ' ' || *s == '\t') s++;
        if (*s == '\0' || *s == '\n' || *s == '#' || *s == '\r') continue;
        rows++;
    }
    fclose(fp);
    return (rows > 0) ? (rows - 1) : 0;
}

static int cmd_results(int open_dir)
{
    const char *dir = cli_results_dir();
    if (open_dir) {
#ifdef _WIN32
        char cmd[600];
        snprintf(cmd, sizeof(cmd), "start \"\" \"%s\"", dir);
        system(cmd);
#else
        char cmd[600];
        snprintf(cmd, sizeof(cmd), "xdg-open \"%s\" &", dir);
        system(cmd);
#endif
        printf("[ok] 已尝试打开结果文件夹: %s\n", dir);
        return 0;
    }

    printf("结果文件夹: %s\n", dir);
#ifdef _WIN32
    {
        char cmd[600];
        snprintf(cmd, sizeof(cmd), "dir /B /O-D \"%s\\*.txt\" 2>NUL", dir);
        system(cmd);
    }
#else
    {
        char cmd[600];
        snprintf(cmd, sizeof(cmd), "ls -t \"%s\"/*.txt 2>/dev/null", dir);
        system(cmd);
    }
#endif
    return 0;
}

/* --- remove --- */
static int cmd_remove(const char *filename, int skip_confirm)
{
    char safe[256];
    char path[600];
    if (cli_safe_filename(filename, safe, sizeof(safe)) != 0) {
        printf("[错误] 文件名非法\n");
        return 1;
    }
    snprintf(path, sizeof(path), "%s/%s", cli_results_dir(), safe);
    {
        FILE *fp = fopen(path, "r");
        if (fp == NULL) {
            printf("[错误] 文件不存在: %s\n", path);
            return 1;
        }
        fclose(fp);
    }
    if (!skip_confirm) {
        char line[16];
        printf("确认删除 %s ? (y/N): ", safe);
        fflush(stdout);
        if (fgets(line, sizeof(line), stdin) == NULL) return 0;
        if (line[0] != 'y' && line[0] != 'Y') {
            printf("已取消\n");
            return 0;
        }
    }
    if (remove(path) == 0) {
        printf("[ok] 已删除 %s\n", safe);
        return 0;
    }
    printf("[错误] 删除失败\n");
    return 1;
}

/* ---------------------------------------------------------- 主演化例程 */

/* 写结果文件。 */
static int write_track_file(const char *path, const EvolutionResult *r,
                            double mass, double z, double dt, double end,
                            const Units *u)
{
    FILE *fp = fopen(path, "w");
    int i;
    char tbuf[64];
    double x, y;
    time_t now = time(NULL);
    struct tm *tm_info = localtime(&now);

    if (fp == NULL) return -1;

    cli_composition(z, &x, &y);
    strftime(tbuf, sizeof(tbuf), "%Y-%m-%d %H:%M:%S", tm_info);

    fprintf(fp, "# ============================================================================\n");
    fprintf(fp, "# %s 恒星演化输出文件\n", CLI_APP_NAME);
    fprintf(fp, "# 程序版本   : %s  发布日期: %s\n",
            CLI_VERSION_STRING, CLI_RELEASE_DATE);
    fprintf(fp, "# 物理模型   : Hurley, Pols & Tout (2000) MNRAS 315, 543 SSE\n");
    fprintf(fp, "# ----------------------------------------------------------------------------\n");
    fprintf(fp, "# 初始质量   : %.6g Msun\n", mass);
    fprintf(fp, "# 金属度 Z   : %.6g   (X = %.4f, Y = %.4f)\n", z, x, y);
    fprintf(fp, "# 步长       : %.6g Myr\n", dt);
    fprintf(fp, "# 主序寿命   : %.4f Myr\n", r->t_ms);
    fprintf(fp, "# 红巨星底端 : %.4f Myr\n", r->t_bgb);
    fprintf(fp, "# 氦点火     : %.4f Myr\n", r->t_he_i);
    fprintf(fp, "# 演化终点   : %s\n",
            (end <= 0.0) ? "残骸阶段" : "用户指定");
    fprintf(fp, "# 演化点数   : %d   步数: %d\n", r->n_track, r->total_steps);
    fprintf(fp, "# 最终状态   : t = %.4f Myr, 阶段 = %s (type %d)\n",
            r->final_state.time, r->final_state.phase_name,
            r->final_state.phase);
    fprintf(fp, "# 最终质量   : %.6g Msun   核心质量: %.6g Msun\n",
            r->final_state.mass, r->final_state.core_mass);
    fprintf(fp, "# 生成时间   : %s\n", tbuf);
    fprintf(fp, "# 单位代码   : mass=%d lum=%d radius=%d temp=%d\n",
            u->mass, u->luminosity, u->radius, u->temperature);
    fprintf(fp, "# ----------------------------------------------------------------------------\n");
    fprintf(fp, "# 数据列: 1=年龄(Myr) 2=亮度 3=半径 4=温度 5=质量\n");
    fprintf(fp, "# ============================================================================\n");
    fprintf(fp, "%16s %20s %20s %16s %20s\n",
            "年龄(Myr)",
            "亮度",
            "半径",
            "温度",
            "质量");

    for (i = 0; i < r->n_track; i++) {
        const StellarState *s = &r->evolution_track[i];
        fprintf(fp, "%16.4f %20.10g %20.10g %16.2f %20.10g\n",
                s->time,
                units_lum_from(u, s->luminosity),
                units_rad_from(u, s->radius),
                units_temp_from(u, s->temperature),
                units_mass_from(u, s->mass));
    }
    fprintf(fp, "\n");
    fclose(fp);
    return 0;
}

/* 单次演化 + 输出文件。返回 0 成功。 */
static int run_evolution(const Units *u, double mass, double z,
                         double step, double end,
                         const char *filename, int banner)
{
    EvolutionResult r;
    double t_ms;
    double dt;
    char safe[256];
    char path[600];
    double x, y;

    t_ms = engine_get_main_sequence_lifetime(mass, z);
    if (t_ms <= 0.0) {
        printf("[错误] 参数越界\n");
        return 1;
    }

    /* 规则: step > 0 时 end_time 不能超过主序寿命 */
    if (step > 0.0 && end > 0.0 && end > t_ms) {
        printf("  [注意] end_time = %g Myr 超过 t_MS = %.2f Myr; 已限制\n",
               end, t_ms);
        end = t_ms;
    }

    /* 默认步长 */
    if (step <= 0.0) {
        dt = t_ms / 1000.0;
    } else {
        dt = step;
    }

    cli_composition(z, &x, &y);
    if (banner) {
        printf("恒星演化: M = %g Msun, Z = %g  (X=%.4f, Y=%.4f)\n",
               mass, z, x, y);
        printf("  主序寿命 t_MS : %.2f Myr\n", t_ms);
        printf("  步长          : %s\n",
               (step <= 0.0) ? "默认 t_MS/1000" : "固定");
        printf("  演化终点      : %s\n",
               (end <= 0.0) ? "残骸阶段" : "用户指定");
        printf("--------------------------------------------------------------------------\n");
    }

    /* 调用引擎 */
    if (engine_evolve(mass, z, end, step, 1, 1, &r) != 0) {
        printf("[错误] 演化失败: %s\n", r.error_message);
        evolution_result_free(&r);
        return 1;
    }

    /* 输出文件 */
    if (filename != NULL && filename[0] != '\0') {
        if (cli_safe_filename(filename, safe, sizeof(safe)) != 0) {
            snprintf(safe, sizeof(safe), "star.txt");
        }
    } else {
        cli_auto_filename(mass, z, safe, sizeof(safe));
    }
    snprintf(path, sizeof(path), "%s/%s", cli_results_dir(), safe);
    write_track_file(path, &r, mass, z, dt, end, u);

    /* 摘要 */
    {
        const StellarState *fs = &r.final_state;
        printf("演化完成  步数 %d, 记录点数 %d, 耗时 %.1f ms\n",
               r.total_steps, r.n_track, r.computation_time_ms);
        printf("  最终年龄  : %.2f Myr\n", fs->time);
        printf("  最终阶段  : %s (type %d, %s)\n",
               fs->phase_name, fs->phase, phase_name_en(fs->phase));
        printf("  最终质量  : %.6g %s\n",
               units_mass_from(u, fs->mass), units_mass_symbol(u));
        printf("  最终光度  : %.6g %s\n",
               units_lum_from(u, fs->luminosity), units_lum_symbol(u));
        printf("  最终半径  : %.6g %s\n",
               units_rad_from(u, fs->radius), units_rad_symbol(u));
        printf("  最终温度  : %.6g %s (%s 型)\n",
               units_temp_from(u, fs->temperature),
               units_temp_symbol(u), fs->spectral_type);
        if (r.stagnated) {
            printf("  [注意] 演化被截断 (未到残骸)\n");
        }
        printf("  输出文件  : %s\n", path);
    }
    evolution_result_free(&r);
    return 0;
}

/* --- evolution_single --- */
static int cmd_evolution_single(Units *u, char tokens[][64], int n_tok,
                                int interactive)
{
    double mass = 1.0, z = 0.02, step = 0.0, end = 0.0;
    char filename[256] = "";

    if (n_tok < 2) {
        if (!interactive) {
            printf("[错误] 非交互模式下必须给出参数\n");
            return 1;
        }
        {
            char buf[128];
            printf("交互式输入恒星参数 (直接回车使用默认值)\n");
            ask("  质量 (Msun)", "1.0", buf, sizeof(buf));
            parse_float(buf, &mass);
            ask("  金属度 Z", "0.02", buf, sizeof(buf));
            parse_float(buf, &z);
            ask("  步长 (Myr, 0=默认 t_MS/1000)", "0", buf, sizeof(buf));
            parse_float(buf, &step);
            ask("  演化终点 (Myr, 0=演化到残骸)", "0", buf, sizeof(buf));
            parse_float(buf, &end);
            ask("  输出文件名 (回车=自动命名)", "", buf, sizeof(buf));
            strncpy(filename, buf, sizeof(filename) - 1);
        }
    } else {
        if (n_tok > 6) {
            printf("[错误] 参数过多\n");
            return 1;
        }
        if (parse_float(tokens[1], &mass) != 0) {
            printf("[错误] 质量必须是数字\n");
            return 1;
        }
        if (n_tok >= 3) parse_float(tokens[2], &z);
        if (n_tok >= 4) parse_float(tokens[3], &step);
        if (n_tok >= 5) parse_float(tokens[4], &end);
        if (n_tok >= 6) strncpy(filename, tokens[5], sizeof(filename) - 1);
    }

    if (mass < ENG_MASS_MIN || mass > ENG_MASS_MAX) {
        printf("[错误] 质量 %g 超出范围 [%.3f, %.0f]\n",
               mass, ENG_MASS_MIN, ENG_MASS_MAX);
        return 1;
    }
    if (z < ENG_Z_MIN || z > ENG_Z_MAX) {
        printf("[错误] 金属度 Z=%g 超出范围 [%.4f, %.2f]\n",
               z, ENG_Z_MIN, ENG_Z_MAX);
        return 1;
    }
    if (step < 0.0 || end < 0.0) {
        printf("[错误] 步长和终点不能为负\n");
        return 1;
    }
    return run_evolution(u, mass, z, step, end,
                         (filename[0] != '\0') ? filename : NULL, 1);
}

/* --- evolution_multiple --- */
static int cmd_evolution_multiple(Units *u,
                                  char groups[][256], int n_grp)
{
    int i;
    printf("批量演化: 共 %d 个恒星\n", n_grp);
    printf("==============================================================================\n");
    for (i = 0; i < n_grp; i++) {
        char fields[5][64];
        int nf = 0;
        char *p = groups[i];
        char *tok;
        double mass = 1.0, z = 0.02, step = 0.0, end = 0.0;
        char filename[256] = "";
        char *saveptr = NULL;

        /* 按逗号拆分 */
        tok = strtok_r(p, ",", &saveptr);
        while (tok != NULL && nf < 5) {
            char *s = tok;
            char *e;
            while (*s == ' ' || *s == '\t') s++;
            e = s + strlen(s);
            while (e > s && (e[-1] == ' ' || e[-1] == '\t')) *--e = '\0';
            strncpy(fields[nf], s, 63);
            fields[nf][63] = '\0';
            nf++;
            tok = strtok_r(NULL, ",", &saveptr);
        }
        if (nf < 4) {
            printf("[跳过] 参数组 (%s) 字段数应为 4 或 5\n", groups[i]);
            continue;
        }
        parse_float(fields[0], &mass);
        parse_float(fields[1], &z);
        parse_float(fields[2], &step);
        parse_float(fields[3], &end);
        if (nf >= 5) strncpy(filename, fields[4], sizeof(filename) - 1);

        printf("\n[%d/%d] M = %g Msun, Z = %g\n", i + 1, n_grp, mass, z);
        run_evolution(u, mass, z, step, end,
                      (filename[0] != '\0') ? filename : NULL, 0);
    }
    return 0;
}

/* ============================================================ 主执行入口 */

/* 返回值: -1 表示退出 */
static int execute(Units *u, const char *line, int interactive)
{
    char tokens[16][64];
    char groups[CLI_MAX_BATCH][256];
    int n_tok, n_grp;
    const char *cmd;

    n_tok = tokenize(line, tokens, 16, groups, CLI_MAX_BATCH, &n_grp);
    if (n_tok == 0) return 0;
    cmd = tokens[0];

    if (strcmp(cmd, "exit") == 0 || strcmp(cmd, "quit") == 0 ||
        strcmp(cmd, "q") == 0) {
        return -1;
    }
    if (strcmp(cmd, "help") == 0 || strcmp(cmd, "?") == 0) {
        cmd_help(); return 0;
    }
    if (strcmp(cmd, "version") == 0) { cmd_version(); return 0; }
    if (strcmp(cmd, "about") == 0)   { cmd_about();   return 0; }
    if (strcmp(cmd, "unit") == 0) {
        return cmd_unit(u, tokens, n_tok);
    }
    if (strcmp(cmd, "phases") == 0) {
        double mass = (n_tok > 1) ? atof(tokens[1]) : 1.0;
        double z    = (n_tok > 2) ? atof(tokens[2]) : 0.02;
        return cmd_phases(mass, z);
    }
    if (strcmp(cmd, "results") == 0) {
        int open_dir = (n_tok > 1 && strcmp(tokens[1], "open") == 0);
        return cmd_results(open_dir);
    }
    if (strcmp(cmd, "remove") == 0) {
        if (n_tok < 2) { printf("用法: remove [filename] [--yes]\n"); return 1; }
        {
            int skip = (n_tok > 2 && strcmp(tokens[2], "--yes") == 0);
            return cmd_remove(tokens[1], skip);
        }
    }
    if (strcmp(cmd, "evolution_single") == 0) {
        return cmd_evolution_single(u, tokens, n_tok, interactive);
    }
    if (strcmp(cmd, "evolution_multiple") == 0) {
        if (n_grp == 0) {
            printf("[错误] evolution_multiple 需要参数组, 例如 (1.0,0.02,0,0)\n");
            return 1;
        }
        return cmd_evolution_multiple(u, groups, n_grp);
    }
    printf("[错误] 未知指令: %s (输入 help 查看)\n", cmd);
    return 1;
}

/* ============================================================ 入口/REPL */

static int repl(Units *u)
{
    char line[1024];
    printf("输入 help 查看指令列表; exit / quit / q 退出程序\n\n");
    while (1) {
        char parts[8][512];
        int n, i;
        printf("%s ", CLI_PROMPT);
        fflush(stdout);
        if (fgets(line, sizeof(line), stdin) == NULL) {
            printf("\n");
            break;
        }
        line[strcspn(line, "\r\n")] = '\0';
        if (line[0] == '\0') continue;
        n = split_statements(line, parts, 8);
        for (i = 0; i < n; i++) {
            if (execute(u, parts[i], 1) < 0)
                goto done;
        }
    }
done:
    printf("再见!\n");
    return 0;
}

int main(int argc, char **argv)
{
    Units u;
    int i;
    int show_banner = 1;
    int interactive = 0;
    int n_commands = 0;
    char *commands[32];
    int rc = 0;

    units_reset(&u);
    units_load(&u);

    for (i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (strcmp(a, "--no-banner") == 0 || strcmp(a, "-n") == 0 ||
            strcmp(a, "--quiet") == 0) {
            show_banner = 0;
        } else if (strcmp(a, "--interactive") == 0 || strcmp(a, "-i") == 0) {
            interactive = 1;
        } else {
            if (n_commands < 32) commands[n_commands++] = argv[i];
        }
    }

    if (show_banner) print_cover();

    if (interactive || n_commands == 0) {
        return repl(&u);
    }

    /* 一次性执行: 把多条命令行参数拼成分号分隔的字符串 */
    {
        char joined[4096] = "";
        char parts[32][512];
        int n, j;
        for (i = 0; i < n_commands; i++) {
            if (i > 0) strcat(joined, " ");
            strncat(joined, commands[i],
                    4095 - strlen(joined));
        }
        n = split_statements(joined, parts, 32);
        for (j = 0; j < n; j++) {
            int ret = execute(&u, parts[j], 0);
            if (ret < 0) break;
            if (ret != 0) rc = ret;
        }
    }
    return rc;
}

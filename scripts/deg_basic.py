#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LKFlow2 Free-Basic · 差异表达分析 DEGs（经典转录组模块 2.1 的免费基础版）
==========================================================================
零代码差异表达分析（纯 Python，无需 R）：
  - 输入：表达矩阵 CSV（行=基因，列=样本）+ 分组表 CSV（sample, group）
        或内置 --demo 合成数据（自动生成"病例 vs 对照"表达谱）
  - 统计：Welch t 检验 + Benjamini-Hochberg FDR + log2FC
  - 输出：火山图 / 显著基因 Top 热图 / DEG 结果表(CSV) / 汇总报告(MD)

说明: 免费基础版提供"单分组对比"差异分析（支持配对/非配对 t 检验）。
      多分组两两对比、limma 经验贝叶斯、DESeq2 负二项模型、批次协变量校正、
      富集分析串联等进阶功能请使用 LKFlow2 完整软件（模块 2.1 DEGs）。
"""
import argparse
import os

import numpy as np
import pandas as pd

FIG_DPI = 150


def _log(msg):
    print(f"[LKFlow-Free/deg] {msg}", flush=True)


# --------------------------------------------------------------------------
# 演示数据：合成"病例 vs 对照"表达矩阵
# --------------------------------------------------------------------------
def make_demo_deg(out_dir, n_genes=12000, n_case=30, n_ctrl=30, seed=42):
    """生成演示表达矩阵 + 分组表。返回 (expr_path, group_path)。"""
    os.makedirs(out_dir, exist_ok=True)
    rng = np.random.default_rng(seed)
    n = n_case + n_ctrl
    genes = [f"GENE{i:05d}" for i in range(1, n_genes + 1)]
    samples = [f"case_{i+1}" for i in range(n_case)] + [f"ctrl_{i+1}" for i in range(n_ctrl)]
    group = ["case"] * n_case + ["ctrl"] * n_ctrl

    # 基线表达（log2 尺度）+ 样本级个体漂移（更接近真实数据）
    mu = rng.normal(4.0, 2.0, n_genes)
    expr = rng.normal(loc=mu[:, None], scale=0.75, size=(n_genes, n))
    expr += rng.normal(0, 0.4, size=n)[None, :]
    expr = np.clip(expr, 0, 20)
    # 注入 ~300 个差异基因（150 up / 150 down），中等幅度 → 火山图 y 轴均衡不空顶
    n_diff = 300
    idx = rng.choice(n_genes, size=n_diff, replace=False)
    fc = rng.choice([-1, 1], size=n_diff) * rng.uniform(0.9, 2.2, n_diff)
    case_mask = np.array([g == "case" for g in group])
    case_idx = np.where(case_mask)[0]
    expr[np.ix_(idx, case_idx)] += fc[:, None]

    df = pd.DataFrame(expr, index=genes, columns=samples)
    expr_path = os.path.join(out_dir, "expression_matrix.csv")
    df.round(4).to_csv(expr_path)

    grp = pd.DataFrame({"sample": samples, "group": group})
    grp_path = os.path.join(out_dir, "group_info.csv")
    grp.to_csv(grp_path, index=False)
    _log(f"✔ 演示数据: {expr_path} / {grp_path}  (差异基因 {n_diff} 个，幅度 |log2FC| 1~3.5)")
    return expr_path, grp_path


# --------------------------------------------------------------------------
# 输入读取
# --------------------------------------------------------------------------
def read_inputs(expr_csv, group_csv, case="case", ctrl="ctrl"):
    expr = pd.read_csv(expr_csv, index_col=0)
    grp = pd.read_csv(group_csv)
    # 兼容单列 / 双列分组表
    if "group" in grp.columns and "sample" in grp.columns:
        pass
    elif grp.shape[1] >= 2:
        grp.columns = ["sample", "group"] + list(grp.columns[2:])
    elif grp.shape[1] == 1:
        grp = grp.rename(columns={grp.columns[0]: "sample"})
        grp["group"] = grp["sample"].map(lambda s: case if s.startswith(case) else ctrl)
    common = [s for s in grp["sample"] if s in expr.columns]
    if not common:
        raise ValueError("分组表样本与表达矩阵列名无交集，请检查")
    grp = grp[grp["sample"].isin(common)].reset_index(drop=True)
    expr = expr[common]
    return expr, grp


# --------------------------------------------------------------------------
# 差异检验
# --------------------------------------------------------------------------
def run_degs(expr, grp, case="case", ctrl="ctrl", paired=False):
    """Welch t 检验（配对时用配对 t 检验）+ BH-FDR。返回 DataFrame。"""
    from scipy import stats
    from scipy.stats import false_discovery_control

    case_samp = grp.loc[grp["group"] == case, "sample"].tolist()
    ctrl_samp = grp.loc[grp["group"] == ctrl, "sample"].tolist()
    a = expr[case_samp].values.T          # (n_case, n_genes)
    b = expr[ctrl_samp].values.T          # (n_ctrl, n_genes)
    if a.shape[0] < 2 or b.shape[0] < 2:
        raise ValueError("每组样本数需 ≥2")

    if paired:
        if a.shape[0] != b.shape[0]:
            raise ValueError("配对检验要求两组样本数一致")
        # 按原始样本顺序配对（case_i ↔ ctrl_i）
        n = min(a.shape[0], b.shape[0])
        a, b = a[:n], b[:n]
        t, p = stats.ttest_rel(a, b, axis=0)
    else:
        t, p = stats.ttest_ind(a, b, axis=0, equal_var=False)

    with np.errstate(divide="ignore", invalid="ignore"):
        mean_a = np.nanmean(a, axis=0)
        mean_b = np.nanmean(b, axis=0)
        log2fc = mean_a - mean_b          # 输入为 log2 尺度 → 均值差即 log2FC
    padj = false_discovery_control(np.nan_to_num(p, nan=1.0))
    padj = np.clip(padj, 1e-300, 1.0)

    res = pd.DataFrame({
        "gene": expr.index,
        "log2FC": log2fc,
        "mean_case": mean_a,
        "mean_ctrl": mean_b,
        "stat": t,
        "pvalue": np.clip(p, 1e-300, 1.0),
        "padj": padj,
    })
    res["direction"] = np.select([res.log2FC > 0, res.log2FC < 0], ["up", "down"], "ns")
    return res


# --------------------------------------------------------------------------
# 出图
# --------------------------------------------------------------------------
def plot_volcano(res, out_path, fdr=0.05, lfc=1.0):
    import matplotlib.pyplot as plt

    from plot_style import DOWN, NS, UP, apply_style, save_fig

    apply_style()
    sig = (res.padj < fdr) & (res.log2FC.abs() >= lfc)
    up = sig & (res.log2FC > 0)
    dn = sig & (res.log2FC < 0)
    nl = -np.log10(res.pvalue.clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    # NS 底层（小、淡） → 显著点上层（大、带白边）
    ax.scatter(res.log2FC[~sig], nl[~sig], s=3, c=NS, alpha=0.45,
               linewidths=0, label=f"NS (n={int((~sig).sum())})", zorder=1)
    ax.scatter(res.log2FC[dn], nl[dn], s=9, c=DOWN, alpha=0.9,
               linewidths=0.4, edgecolors="white", label=f"Down (n={int(dn.sum())})", zorder=3)
    ax.scatter(res.log2FC[up], nl[up], s=9, c=UP, alpha=0.9,
               linewidths=0.4, edgecolors="white", label=f"Up (n={int(up.sum())})", zorder=3)
    ax.axvline(lfc, color="#999", lw=0.8, ls="--", zorder=2)
    ax.axvline(-lfc, color="#999", lw=0.8, ls="--", zorder=2)
    ax.axhline(-np.log10(fdr), color="#999", lw=0.8, ls="--", zorder=2)
    ax.grid(True, ls=":", lw=0.5, color="#E3E6EC", zorder=0)

    # 标注 top 基因（up 前 3 + down 前 3，按 padj 取，文字贴点放置不出界）
    label_df = pd.concat([res[up].nsmallest(3, "padj"), res[dn].nsmallest(3, "padj")])
    ymax = nl.max() * 1.12
    xmax = max(2.6, res.log2FC.abs().max() + 0.4)
    for i, (_, r) in enumerate(label_df.iterrows()):
        y_pt = -np.log10(max(r.pvalue, 1e-300))
        side = 1 if r.log2FC > 0 else -1
        tx = r.log2FC + side * 0.22
        ha = "left" if side > 0 else "right"
        ty = min(y_pt + 1.0 + (i % 2) * 1.4, ymax * 0.94)
        if tx > xmax * 0.97:
            tx, ha = r.log2FC - 0.22, "right"
        elif tx < -xmax * 0.97:
            tx, ha = r.log2FC + 0.22, "left"
        ax.annotate(r.gene, (r.log2FC, y_pt), xytext=(tx, ty),
                    fontsize=7, color="#37474F", ha=ha,
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="#90A4AE"))
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("log2 Fold Change (case vs ctrl)")
    ax.set_ylabel("-log10 P-value")
    ax.set_title("Differential Expression Volcano")
    ax.legend(loc="lower right")
    save_fig(fig, out_path)
    _log(f"✔ 火山图: {out_path}")


def plot_deg_heatmap(expr, res, out_path, top_n=50):
    """Top 差异基因 z-score 热图（按 padj 排序），蓝-白-橙发散配色。"""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    from plot_style import DIVERGING, apply_style, save_fig

    apply_style()
    top = res.sort_values("padj").head(top_n).gene.values
    sub = expr.loc[top].copy()
    # 行 z-score
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, 1), axis=0)
    z = z.clip(-2.5, 2.5)

    fig, ax = plt.subplots(figsize=(9, 10.5))
    cmap = LinearSegmentedColormap.from_list("bwo", DIVERGING)
    im = ax.imshow(z.values, aspect="auto", cmap=cmap, vmin=-2.5, vmax=2.5,
                   interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("Samples")
    ax.set_ylabel("Genes")
    ax.set_title(f"Top {top_n} DEGs (by padj) — z-score")
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(0.8)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("z-score", fontsize=9)
    cb.outline.set_linewidth(0.6)
    save_fig(fig, out_path)
    _log(f"✔ Top DEG 热图: {out_path}")


# --------------------------------------------------------------------------
# 主流程
# --------------------------------------------------------------------------
def run_deg(demo=False, expr_csv=None, group_csv=None, out_dir=None,
            case="case", ctrl="ctrl", fdr=0.05, lfc=1.0, paired=False, top_n=50,
            min_mean=0.5):
    import textwrap

    if not out_dir:
        out_dir = os.path.join(os.getcwd(), "lkflow_out", "deg")
    os.makedirs(out_dir, exist_ok=True)

    if demo:
        expr_csv, group_csv = make_demo_deg(out_dir)
        if case == "case" and ctrl == "ctrl":
            pass
        _log(f"演示模式: 病例={case} vs 对照={ctrl}")

    _log(f"读取表达矩阵: {expr_csv}")
    expr, grp = read_inputs(expr_csv, group_csv, case=case, ctrl=ctrl)
    # 低表达过滤（标准做法）：均值过低的基因方差趋零，会产生病态 p 值
    n0 = expr.shape[0]
    expr = expr[expr.mean(axis=1) >= min_mean]
    _log(f"基因 {n0} → {expr.shape[0]}（低表达过滤 mean ≥ {min_mean}）")
    _log(f"样本 {expr.shape[1]}  |  {case}: {int((grp.group == case).sum())}   "
         f"{ctrl}: {int((grp.group == ctrl).sum())}")

    res = run_degs(expr, grp, case=case, ctrl=ctrl, paired=paired)
    res = res.sort_values(["padj", "pvalue"])
    res_csv = os.path.join(out_dir, "deg_results_all_genes.csv")
    res.round(6).to_csv(res_csv, index=False)

    sig = res[(res.padj < fdr) & (res.log2FC.abs() >= lfc)]
    sig_up, sig_dn = int((sig.direction == "up").sum()), int((sig.direction == "down").sum())
    _log(f"显著基因 (padj<{fdr} & |log2FC|≥{lfc}): {len(sig)}  (up {sig_up} / down {sig_dn})")

    plot_volcano(res, os.path.join(out_dir, "deg_volcano.png"), fdr=fdr, lfc=lfc)
    if len(sig) >= 5:
        plot_deg_heatmap(expr, sig, os.path.join(out_dir, "deg_heatmap_top50.png"),
                         top_n=min(top_n, len(sig)))
    elif len(sig) > 0:
        plot_deg_heatmap(expr, res, os.path.join(out_dir, "deg_heatmap_top50.png"),
                         top_n=min(top_n, len(sig)))

    # 报告
    top = sig.head(30)
    lines = [
        "# DEG 差异表达分析报告（LKFlow2 Free-Basic · 经典转录组 deg）",
        "",
        f"- 输入矩阵: {expr_csv}",
        f"- 分组: {case} (n={int((grp.group == case).sum())}) vs {ctrl} (n={int((grp.group == ctrl).sum())})",
        f"- 检验: {'配对' if paired else 'Welch 非配对'} t 检验 + BH-FDR",
        f"- 阈值: padj < {fdr} 且 |log2FC| ≥ {lfc}",
        "",
        "## 统计结果",
        "",
        f"- 差异基因总数: **{len(sig)}**（上调 {sig_up} / 下调 {sig_dn}）",
        f"- Top DEG 热图基因数: {min(top_n, len(sig))}",
        "",
        "## Top 15 显著基因",
        "",
        "| gene | log2FC | mean_case | mean_ctrl | pvalue | padj | direction |",
        "|---|---|---|---|---|---|---|",
    ]
    for _, r in top.head(15).iterrows():
        lines.append(f"| {r.gene} | {r.log2FC:.2f} | {r.mean_case:.2f} | {r.mean_ctrl:.2f} "
                     f"| {r.pvalue:.2e} | {r.padj:.2e} | {r.direction} |")
    lines += [
        "",
        "## 文件清单",
        "- 图: deg_volcano.png / deg_heatmap_top50.png",
        "- 表: deg_results_all_genes.csv（全部基因检验结果）",
        "",
        "## 下一步（进阶功能请在 LKFlow2 完整软件中体验）",
        "- 多分组两两对比、limma 经验贝叶斯、DESeq2、批次/协变量校正",
        "- 差异基因 → 韦恩图交集(venn) → 机器学习筛选(ml) → 富集/免疫浸润/生存验证 一键串联",
        "- 网络药理学(3.x)、PPI(5.1)、GO+KEGG(5.2)、GSEA/GSVA/免疫浸润/ROC/KM 生存(7.x)",
    ]
    with open(os.path.join(out_dir, "summary_report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    _log(f"✔ 完成！结果目录: {os.path.abspath(out_dir)}")
    return out_dir


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", action="store_true", help="使用内置演示数据")
    ap.add_argument("--expr", default=None, help="表达矩阵 CSV（行=基因 列=样本）")
    ap.add_argument("--group", default=None, help="分组表 CSV（sample, group）")
    ap.add_argument("--case", default="case", help="病例组标签")
    ap.add_argument("--ctrl", default="ctrl", help="对照组标签")
    ap.add_argument("--out", default=None)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--lfc", type=float, default=1.0)
    ap.add_argument("--paired", action="store_true")
    a = ap.parse_args()
    if not a.demo and not (a.expr and a.group):
        ap.error("请用 --demo 或同时提供 --expr 与 --group")
    run_deg(demo=a.demo, expr_csv=a.expr, group_csv=a.group, out_dir=a.out,
            case=a.case, ctrl=a.ctrl, fdr=a.fdr, lfc=a.lfc, paired=a.paired)

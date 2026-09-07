#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LKFlow2 Free-Basic · 孟德尔随机化（MR）基础因果分析
====================================================
纯 Python 实现 Two-sample MR 核心方法，无需 R 环境：
  - 数据对齐 harmonise（等位基因方向校正）
  - 5 种因果估计：IVW（固定/随机效应）、MR-Egger、加权中位数、加权众数(直方图近似)
  - 敏感性检验：Cochran's Q 异质性、Egger 截距多效性、留一法 leave-one-out
  - 输出：散点图 / 森林图 / 漏斗图 / 留一森林图 + 5 张结果表 + 汇总报告

输入格式（CSV，逗号分隔，含表头）:
  SNP, effect_allele, other_allele, eaf, beta, se, pval
  （两表需含相同 SNP；未匹配等位基因将按互补方向自动翻转校正）

说明: 免费基础版支持"单暴露 vs 单结局"摘要数据 MR。
      多基因批量扫描(全基因组)、多变量 MR、双向 MR、共定位 coloc、
      药物靶点(pQTL)MR、多效性分层等进阶功能请使用 LKFlow2 完整软件。
"""
import argparse
import os
import sys
import textwrap

import numpy as np
import pandas as pd

FIG_DPI = 150
METHODS = ["IVW", "MR-Egger", "Weighted median", "Weighted mode", "Simple mode"]


def _log(msg):
    print(f"[LKFlow-Free/mr] {msg}", flush=True)


# --------------------------------------------------------------------------
# 工具函数：加权回归
# --------------------------------------------------------------------------
def _wlm(x, y, w, intercept=True):
    """加权最小二乘。返回 (coef, se)。"""
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    w = np.asarray(w, float).ravel()
    if intercept:
        X = np.column_stack([np.ones_like(x), x])
    else:
        X = x[:, None]
    sw = np.sqrt(w)
    Xw, yw = X * sw[:, None], y * sw
    coef, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    resid = yw - Xw @ coef
    dof = len(y) - X.shape[1]
    sigma2 = (resid @ resid) / dof if dof > 0 else 0.0
    cov = sigma2 * np.linalg.inv(Xw.T @ Xw)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    return coef, se


def _weighted_median(theta, w):
    """加权中位数（Bowden 2016）。"""
    order = np.argsort(theta)
    theta_s, w_s = theta[order], w[order]
    cw = np.cumsum(w_s) / w_s.sum()
    idx = np.searchsorted(cw, 0.5)
    if cw[idx] == 0.5:
        return (theta_s[idx] + theta_s[idx + 1]) / 2.0
    return theta_s[idx]


def _weighted_mode(theta, w, bins=100):
    """加权众数（直方图近似，取权重最大的 bin 质心）。"""
    lo, hi = np.quantile(theta, [0.01, 0.99])
    edges = np.linspace(lo, hi, bins + 1)
    tot = np.zeros(bins)
    for i in range(bins):
        m = (theta >= edges[i]) & (theta < edges[i + 1])
        if m.any():
            tot[i] = w[m].sum()
    idx = int(np.argmax(tot))
    return (edges[idx] + edges[idx + 1]) / 2.0


# --------------------------------------------------------------------------
# 数据对齐
# --------------------------------------------------------------------------
def harmonise(exposure, outcome):
    """按 SNP 内连接并对齐效应等位基因方向。返回合并后数据。"""
    req = ["SNP", "effect_allele", "other_allele", "beta", "se", "pval"]
    for name, df in (("exposure", exposure), ("outcome", outcome)):
        miss = [c for c in req if c not in df.columns]
        if miss:
            raise ValueError(f"{name} 缺少列: {miss}（需要 {req}；eaf 可选）")

    m = pd.merge(exposure, outcome, on="SNP", suffixes=("_exp", "_out"), how="inner")
    _log(f"合并后共 {len(m)} 个匹配 SNP")
    if len(m) < 3:
        raise ValueError("匹配 SNP 过少(<3)，请检查两表 SNP 是否一致")

    ea, oa = m.effect_allele_exp.values, m.other_allele_exp.values
    eao, oao = m.effect_allele_out.values, m.other_allele_out.values
    same = (ea == eao)                    # 效应等位基因一致 → 无需处理
    comp = (ea == oao) & (oa == eao)      # 互补颠倒 → 翻转 outcome
    ok = same | comp

    n_bad = int((~ok).sum())
    if n_bad:
        _log(f"⚠ {n_bad} 个 SNP 等位基因无法对齐，已剔除")
    m = m[ok].copy().reset_index(drop=True)

    idx = m.index[(m.effect_allele_exp.values != m.effect_allele_out.values)]
    if len(idx):
        m.loc[idx, "beta_out"] = -m.loc[idx, "beta_out"]
        if "eaf_out" in m.columns:
            m.loc[idx, "eaf_out"] = 1 - m.loc[idx, "eaf_out"]
        ea_tmp = m.loc[idx, "effect_allele_out"].copy()
        m.loc[idx, "effect_allele_out"] = m.loc[idx, "other_allele_out"]
        m.loc[idx, "other_allele_out"] = ea_tmp
    if len(m) < 3:
        raise ValueError("对齐后 SNP 过少(<3)")
    return m.reset_index(drop=True)


# --------------------------------------------------------------------------
# MR 方法
# --------------------------------------------------------------------------
def mr_estimates(m):
    """5 法估计。返回方法表。"""
    bx, by = m.beta_exp.values.astype(float), m.beta_out.values.astype(float)
    sey = m.se_out.values.astype(float)
    # 以 outcome 的 se 为权重（两样本 MR 常用）
    theta = by / bx
    w = 1.0 / sey**2

    res = []
    # IVW（固定效应：无截距加权回归 by ~ -1 + bx）
    coef_ivw, se_ivw = _wlm(bx, by, w, intercept=False)
    b_ivw, s_ivw = float(coef_ivw[0]), float(se_ivw[0])
    # 异质性 Cochran's Q
    q = float(np.sum(w * (by - b_ivw * bx) ** 2))
    q_p = float(1.0 if len(bx) <= 1 else _chi2_p(q, len(bx) - 1))
    # 随机效应校正
    se_ivw_re = s_ivw * np.sqrt(max(1.0, q / max(1, len(bx) - 1))) if len(bx) > 1 else s_ivw
    res.append(("IVW", b_ivw, s_ivw, q_p))

    # Egger
    coef_eg, se_eg = _wlm(bx, by, w, intercept=True)
    b_eg, s_eg = float(coef_eg[1]), float(se_eg[1])
    intercept_eg, se_int = float(coef_eg[0]), float(se_eg[0])
    z_int = intercept_eg / se_int
    p_int = 2 * (1 - _norm_cdf(abs(z_int)))
    res.append(("MR-Egger", b_eg, s_eg, p_int))

    # 加权中位数
    b_wm = _weighted_median(theta, w)
    se_wm = _median_se(theta, w)
    res.append(("Weighted median", b_wm, se_wm, float("nan")))

    # 加权众数 / 简单众数
    b_wmode = _weighted_mode(theta, w)
    b_smode = _weighted_mode(theta, np.ones_like(w))
    res.append(("Weighted mode", b_wmode, _mode_se(theta, w), float("nan")))
    res.append(("Simple mode", b_smode, _mode_se(theta, np.ones_like(w)), float("nan")))

    df = pd.DataFrame(res, columns=["method", "beta", "se", "sensitivity_p"])
    df["OR"] = np.exp(df.beta)
    df["OR_lci"] = np.exp(df.beta - 1.96 * df.se)
    df["OR_uci"] = np.exp(df.beta + 1.96 * df.se)
    df["pval"] = [2 * (1 - _norm_cdf(abs(b / s))) if s > 0 else np.nan
                  for b, s in zip(df.beta, df.se)]
    df["se"] = df.se.where(df.se > 0, np.nan)
    return df, dict(q=q, q_df=len(bx) - 1, q_p=q_p,
                    egger_intercept=intercept_eg, egger_intercept_se=se_int,
                    egger_intercept_p=p_int, n_snp=len(bx))


def _median_se(theta, w):
    """加权中位数法 SE（简化：以 θ 的加权 MAD 近似，保留展示用途）。"""
    order = np.argsort(theta)
    tw, t = theta[order], w[order]
    cw = np.cumsum(tw) / tw.sum()
    lo, hi = tw[np.searchsorted(cw, 0.16)], tw[np.searchsorted(cw, 0.84)]
    return max(abs(hi - lo) / 1.349, 1e-9)


def _mode_se(theta, w):
    s = np.std(theta, ddof=1) if len(theta) > 1 else 1e-9
    return max(s / np.sqrt(len(theta)), 1e-9)


def _norm_cdf(x):
    from scipy.stats import norm
    return norm.cdf(x)


def _chi2_p(x, df):
    from scipy.stats import chi2
    return float(chi2.sf(x, df))


# --------------------------------------------------------------------------
# 留一法
# --------------------------------------------------------------------------
def leave_one_out(m):
    rows = []
    bx, by = m.beta_exp.values.astype(float), m.beta_out.values.astype(float)
    sey = m.se_out.values.astype(float)
    for i in range(len(m)):
        keep = np.ones(len(m), bool)
        keep[i] = False
        coef, se = _wlm(bx[keep], by[keep], 1.0 / sey[keep] ** 2, intercept=False)
        rows.append((m.SNP[i], float(coef[0]), float(se[0])))
    df = pd.DataFrame(rows, columns=["excluded_SNP", "beta", "se"])
    df["OR"] = np.exp(df.beta)
    df["OR_lci"] = np.exp(df.beta - 1.96 * df.se)
    df["OR_uci"] = np.exp(df.beta + 1.96 * df.se)
    return df


# --------------------------------------------------------------------------
# 绘图
# --------------------------------------------------------------------------
def plot_scatter(m, ax, ivw_b, egger_b, title):
    bx, by = m.beta_exp.values.astype(float), m.beta_out.values.astype(float)
    sex, sey = m.se_exp.values.astype(float), m.se_out.values.astype(float)
    xs = np.linspace(bx.min() * 0.8, bx.max() * 1.2, 100)
    ax.errorbar(bx, by, xerr=1.96 * sex, yerr=1.96 * sey, fmt="o", ms=4,
                color="#5B7DB1", ecolor="#B9C6DA", elinewidth=1, capsize=0, zorder=2)
    ax.plot(xs, ivw_b * xs, color="#C0392B", lw=2, label=f"IVW slope={ivw_b:.3f}")
    ax.plot(xs, egger_b * xs, color="#27AE60", lw=2, ls="--", label=f"Egger slope={egger_b:.3f}")
    ax.set_xlabel("SNP effect on exposure (β)")
    ax.set_ylabel("SNP effect on outcome (β)")
    ax.set_title(title)
    ax.legend(fontsize=8, frameon=False)


def _fmt_p(p):
    try:
        if pd.isna(p):
            return "n/a"
        return f"{p:.3g}"
    except Exception:
        return "n/a"


def plot_forest(mr_df, ax, title):
    rows = mr_df.copy()
    rows = rows.sort_values("beta")
    y = np.arange(len(rows))
    ax.errorbar(rows.OR, y, xerr=[rows.OR - rows.OR_lci, rows.OR_uci - rows.OR],
                fmt="o", ms=6, color="#1B4F72", ecolor="#1B4F72", elinewidth=1.5,
                capsize=3, zorder=3)
    ax.axvline(1.0, color="gray", lw=1, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.method} (p={_fmt_p(r.pval)})" for _, r in rows.iterrows()],
                       fontsize=8)
    ax.set_xlabel("OR (95% CI)")
    ax.set_title(title)


def plot_funnel(m, ax, ivw_b):
    bx, by = m.beta_exp.values.astype(float), m.beta_out.values.astype(float)
    sey = m.se_out.values.astype(float)
    theta, inv = by / bx, bx / sey
    ax.scatter(theta, inv, s=18, color="#5B7DB1", zorder=2)
    ax.axvline(ivw_b, color="#C0392B", lw=1.5)
    ax.set_xlabel("θ (β_outcome / β_exposure)")
    ax.set_ylabel("Instrument strength (β_exp / se_out)")
    ax.set_title("Funnel plot")


def plot_loo(loo, ax, ivw_b, ivw_se):
    rows = loo.sort_values("beta")
    y = np.arange(len(rows))
    ax.errorbar(rows.OR, y, xerr=[rows.OR - rows.OR_lci, rows.OR_uci - rows.OR],
                fmt="s", ms=4, color="#7D6608", ecolor="#B7950B", elinewidth=1, capsize=2, zorder=3)
    ax.axvline(np.exp(ivw_b), color="#C0392B", lw=1.5)
    ax.axvline(np.exp(ivw_b - 1.96 * ivw_se), color="#C0392B", lw=0.8, ls=":")
    ax.axvline(np.exp(ivw_b + 1.96 * ivw_se), color="#C0392B", lw=0.8, ls=":")
    ax.set_yticks(y)
    ax.set_yticklabels(rows.excluded_SNP, fontsize=6)
    ax.set_xlabel("OR (leave-one-out IVW)")
    ax.set_title("Leave-one-out analysis")


# --------------------------------------------------------------------------
# 演示数据
# --------------------------------------------------------------------------
def make_demo_mr_data(out_dir="./demo_mr", n_snp=50, seed=42):
    """生成自洽的演示 GWAS 摘要数据（带部分多效性 SNP）。"""
    rng = np.random.default_rng(seed)
    os.makedirs(out_dir, exist_ok=True)

    bases = ["A", "C", "G", "T"]
    snps = [f"rs{1000000 + i:07d}" for i in range(n_snp)]
    ea = rng.choice(bases, n_snp)
    oa = np.array([rng.choice([x for x in bases if x != e]) for e in ea])

    beta_exp = rng.normal(0.05, 0.02, n_snp)
    se_exp = rng.uniform(0.006, 0.012, n_snp)
    se_out = rng.uniform(0.008, 0.014, n_snp)

    causal_b = 0.22  # 潜在因果效应（正相关）
    noise = rng.normal(0, se_out * 0.6, n_snp)
    beta_out = causal_b * beta_exp + noise
    # 制造 6 个多效性 SNP（对 Egger 截距有贡献）
    pleio_idx = rng.choice(n_snp, size=6, replace=False)
    beta_out[pleio_idx] += rng.normal(0.03, 0.01, 6)

    def _table(beta, se, eaf_allele):
        df = pd.DataFrame({
            "SNP": snps,
            "effect_allele": eaf_allele[0],
            "other_allele": eaf_allele[1],
            "eaf": rng.uniform(0.1, 0.9, n_snp).round(5),
            "beta": np.round(beta, 5),
            "se": np.round(se, 5),
            "pval": np.round(2 * (1 - _norm_cdf(np.abs(beta / se))), 6),
        })
        return df

    exp = _table(beta_exp, se_exp, (ea, oa))
    # 结局表：部分 SNP 故意让等位基因互补颠倒，演示 harmonise 自动翻转
    flip = rng.choice([True, False], n_snp)
    out_ea = np.where(flip, oa, ea)
    out_oa = np.where(flip, ea, oa)
    out = _table(np.where(flip, -beta_out, beta_out), se_out, (out_ea, out_oa))
    out["eaf"] = np.where(flip, 1 - out.eaf, out.eaf).round(5)

    exp.to_csv(os.path.join(out_dir, "exposure.csv"), index=False)
    out.to_csv(os.path.join(out_dir, "outcome.csv"), index=False)
    _log(f"✔ 演示数据已生成: {os.path.abspath(out_dir)}/exposure.csv, outcome.csv")
    return out_dir


# --------------------------------------------------------------------------
# 主流程
# --------------------------------------------------------------------------
def run_mr(exposure_csv, outcome_csv, exposure_name="exposure", outcome_name="outcome",
           out_dir=None):
    _log("LKFlow2 Free-Basic · 孟德尔随机化因果分析启动")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    exp = pd.read_csv(exposure_csv)
    out = pd.read_csv(outcome_csv)
    _log(f"暴露 {len(exp)} SNP / 结局 {len(out)} SNP")

    m = harmonise(exp, out)

    mr_df, extra = mr_estimates(m)
    loo = leave_one_out(m)

    base = out_dir or os.path.join("lkflow_out", "mr")
    os.makedirs(base, exist_ok=True)

    mr_df.to_csv(os.path.join(base, "mr_results.csv"), index=False)
    m[["SNP", "beta_exp", "se_exp", "beta_out", "se_out", "effect_allele_exp", "effect_allele_out"]].to_csv(
        os.path.join(base, "snps_used.csv"), index=False)
    pd.DataFrame([{"Q_stat": extra["q"], "df": extra["q_df"], "Q_pval": extra["q_p"],
                   "egger_intercept": extra["egger_intercept"],
                   "egger_intercept_se": extra["egger_intercept_se"],
                   "egger_intercept_p": extra["egger_intercept_p"],
                   "n_snp": extra["n_snp"]}]).to_csv(os.path.join(base, "sensitivity.csv"), index=False)
    loo.to_csv(os.path.join(base, "leave_one_out.csv"), index=False)

    ivw = mr_df[mr_df.method == "IVW"].iloc[0]
    egger = mr_df[mr_df.method == "MR-Egger"].iloc[0]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    plot_scatter(m, axes[0, 0], ivw.beta, egger.beta,
                 f"MR scatter: {exposure_name} → {outcome_name}")
    plot_forest(mr_df, axes[0, 1], "Causal estimates (OR)")
    plot_funnel(m, axes[1, 0], ivw.beta)
    plot_loo(loo, axes[1, 1], ivw.beta, ivw.se)
    fig.suptitle(f"MR analysis · {exposure_name} → {outcome_name}", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(base, f"mr_figures.{ext}"), dpi=FIG_DPI)
    plt.close(fig)

    # 单图（便于直接插入文档）
    figs = [("mr_scatter", lambda: plot_scatter(m, plt.gca(), ivw.beta, egger.beta,
                                                 f"{exposure_name} → {outcome_name}")),
            ("mr_forest", lambda: plot_forest(mr_df, plt.gca(), "Causal estimates (OR)")),
            ("mr_funnel", lambda: plot_funnel(m, plt.gca(), ivw.beta)),
            ("mr_leave_one_out", lambda: plot_loo(loo, plt.gca(), ivw.beta, ivw.se))]
    for name, fn in figs:
        f2, ax2 = plt.subplots(figsize=(6, 5))
        fn()
        f2.tight_layout()
        f2.savefig(os.path.join(base, f"{name}.png"), dpi=FIG_DPI)
        plt.close(f2)

    _write_mr_report(base, exposure_name, outcome_name, mr_df, extra, m, loo)
    _log(f"✔ 完成！结果目录: {os.path.abspath(base)}")


def _write_mr_report(base, exp_name, out_name, mr_df, extra, m, loo):
    ivw = mr_df[mr_df.method == "IVW"].iloc[0]
    dir_sign = "增加" if ivw.beta > 0 else "降低"
    lines = [f"# MR 分析报告（LKFlow2 Free-Basic）",
             "",
             f"- 暴露: {exp_name}  |  结局: {out_name}",
             f"- 匹配工具变量 SNP 数: {extra['n_snp']}",
             f"- IVW 因果估计 β={ivw.beta:.4f}，OR={ivw.OR:.3f} "
             f"(95% CI {ivw.OR_lci:.3f}-{ivw.OR_uci:.3f})，p={ivw.pval:.3g}",
             f"- 方向解读: 暴露每变化一个单位，结局风险{'增加' if ivw.OR >= 1 else '降低'} "
             f"{(max(ivw.OR, 1 / ivw.OR) - 1) * 100:.1f}%",
             f"- 异质性: Cochran's Q={extra['q']:.2f} (p={extra['q_p']:.3g})",
             f"- 多效性: Egger 截距={extra['egger_intercept']:.4f} "
             f"(p={extra['egger_intercept_p']:.3g})",
             "",
             "## MR 方法结果", ""]
    lines.append("| method | beta | se | OR | OR_lci | OR_uci | pval |")
    lines.append("|---|---|---|---|---|---|---|")
    for _, r in mr_df.iterrows():
        lines.append(f"| {r.method} | {r.beta:.4f} | {r.se:.4f} | {r.OR:.3f} "
                     f"| {r.OR_lci:.3f} | {r.OR_uci:.3f} | {r.pval:.3g} |")
    lines += ["",
              "## 文件清单",
              "- 图: mr_figures.png/pdf, mr_scatter.png, mr_forest.png, mr_funnel.png, mr_leave_one_out.png",
              "- 表: mr_results.csv, sensitivity.csv, leave_one_out.csv, snps_used.csv",
              "",
              "## 下一步（进阶功能请在 LKFlow2 完整软件中体验）",
              "- 暴露-结局全基因组批量筛选、差异基因∩MR 交集、森林图汇总",
              "- 多变量 MR(MVMR)、双向 MR、共定位(coloc)、药物靶点(pQTL)MR、多效性分层",
              "- STROBE-MR 20 条投稿清单自动生成 + AI 结果解读"]
    with open(os.path.join(base, "summary_report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p_data = sub.add_parser("demo-data")
    p_data.add_argument("--out", default="./demo_mr")
    p_data.add_argument("--n-snp", type=int, default=50)
    p_mr = sub.add_parser("mr")
    p_mr.add_argument("--exposure", required=True)
    p_mr.add_argument("--outcome", required=True)
    p_mr.add_argument("--out", default=None)
    args = ap.parse_args()
    if args.cmd == "demo-data":
        make_demo_mr_data(out_dir=args.out, n_snp=args.n_snp)
    else:
        run_mr(exposure_csv=args.exposure, outcome_csv=args.outcome, out_dir=args.out)

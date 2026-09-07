#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LKFlow2 Free-Basic · 机器学习特征基因筛选（经典转录组模块 6.x 的免费基础版）
=============================================================================
零代码机器学习建模与候选基因筛选（纯 Python / scikit-learn）：
  - 输入：表达矩阵 CSV（行=基因，列=样本）+ 分组表 CSV（sample, group）
        或内置 --demo 合成数据（病例 vs 对照，含 20 个真实信号基因）
  - 流程：单变量特征预筛(ANOVA F) → 8 种分类器 5 折交叉验证 → 平均 ROC-AUC 排序
         → ROC 曲线 / AUC 条形图 / 随机森林特征重要性 Top 基因
  - 输出：模型对比图 / ROC 图 / 特征重要性图 / 模型评分表 / Top 基因表 / 汇总报告

说明: 免费基础版提供"分类建模 + 8 种算法对比 + 特征重要性"基础闭环。
      机器学习-10/127/8/14/207 大规模算法组合寻优、SHAP 解释、
      Lasso 特征筛选、预后生存模型(6.7)、外部验证请在 LKFlow2 完整软件中体验。
"""
import argparse
import os

import numpy as np
import pandas as pd

FIG_DPI = 150


def _log(msg):
    print(f"[LKFlow-Free/ml] {msg}", flush=True)


# --------------------------------------------------------------------------
# 演示数据
# --------------------------------------------------------------------------
def make_demo_ml(out_dir, n_genes=1500, n_case=30, n_ctrl=30, n_signal=20, seed=1):
    os.makedirs(out_dir, exist_ok=True)
    rng = np.random.default_rng(seed)
    n = n_case + n_ctrl
    genes = [f"GENE{i:05d}" for i in range(1, n_genes + 1)]
    samples = [f"case_{i+1}" for i in range(n_case)] + [f"ctrl_{i+1}" for i in range(n_ctrl)]
    group = ["case"] * n_case + ["ctrl"] * n_ctrl

    mu = rng.normal(4.0, 2.0, n_genes)
    expr = rng.normal(loc=mu[:, None], scale=0.9, size=(n_genes, n))
    expr = np.clip(expr, 0, 20)
    sig = rng.choice(n_genes, size=n_signal, replace=False)
    # 信号较弱 + 样本个体漂移噪声 → 模型 AUC 约 0.8~0.99 有区分度
    fc = rng.choice([-1, 1], size=n_signal) * rng.uniform(0.7, 1.5, n_signal)
    case_idx = np.where(np.array([g == "case" for g in group]))[0]
    expr[np.ix_(sig, case_idx)] += fc[:, None]
    expr += rng.normal(0, 0.7, size=n)[None, :]   # 样本级批次/个体噪声

    df = pd.DataFrame(expr, index=genes, columns=samples)
    ep = os.path.join(out_dir, "expression_matrix.csv")
    df.round(4).to_csv(ep)
    grp = pd.DataFrame({"sample": samples, "group": group})
    gp = os.path.join(out_dir, "group_info.csv")
    grp.to_csv(gp, index=False)
    _log(f"✔ 演示数据: {ep}  (信号基因 {n_signal} 个，|log2FC| 0.6~1.3)")
    return ep, gp


def read_inputs(expr_csv, group_csv, case="case", ctrl="ctrl"):
    expr = pd.read_csv(expr_csv, index_col=0)
    grp = pd.read_csv(group_csv)
    if grp.shape[1] >= 2 and "group" not in grp.columns:
        grp.columns = ["sample", "group"] + list(grp.columns[2:])
    common = [s for s in grp["sample"] if s in expr.columns]
    grp = grp[grp["sample"].isin(common)].reset_index(drop=True)
    expr = expr[common]
    y = np.array([1 if g == case else 0 for g in grp["group"]])
    # 仅保留 case/ctrl 两组的样本
    keep = grp["group"].isin([case, ctrl])
    return expr.loc[:, grp.loc[keep, "sample"]], np.array([1 if g == case else 0
                                                           for g in grp.loc[keep, "group"]])


# --------------------------------------------------------------------------
# 模型定义
# --------------------------------------------------------------------------
def get_models():
    from sklearn.ensemble import (AdaBoostClassifier, ExtraTreesClassifier,
                                  GradientBoostingClassifier, RandomForestClassifier)
    from sklearn.linear_model import LogisticRegression
    from sklearn.naive_bayes import GaussianNB
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.svm import SVC

    return {
        "Logistic": LogisticRegression(max_iter=2000, C=1.0),
        "RandomForest": RandomForestClassifier(n_estimators=300, random_state=0, n_jobs=-1),
        "ExtraTrees": ExtraTreesClassifier(n_estimators=300, random_state=0, n_jobs=-1),
        "SVM": SVC(kernel="rbf", probability=True, random_state=0),
        "KNN": KNeighborsClassifier(n_neighbors=5),
        "GradientBoosting": GradientBoostingClassifier(n_estimators=200, random_state=0),
        "AdaBoost": AdaBoostClassifier(n_estimators=200, random_state=0),
        "NaiveBayes": GaussianNB(),
    }


# --------------------------------------------------------------------------
# 主流程
# --------------------------------------------------------------------------
def run_ml(demo=False, expr_csv=None, group_csv=None, out_dir=None,
           case="case", ctrl="ctrl", n_keep=150, cv=5, seed=0):
    import warnings
    warnings.filterwarnings("ignore", message=".*probability.*deprecated.*")
    warnings.filterwarnings("ignore", message=".*SVC.*probability.*")
    try:
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.feature_selection import SelectKBest, f_classif
        from sklearn.metrics import auc, roc_auc_score, roc_curve
        from sklearn.model_selection import StratifiedKFold, cross_val_predict
    except Exception as e:
        raise RuntimeError("缺少 scikit-learn，请先运行: python3 scripts/lkflow_free.py install") from e

    if not out_dir:
        out_dir = os.path.join(os.getcwd(), "lkflow_out", "ml")
    os.makedirs(out_dir, exist_ok=True)

    if demo:
        expr_csv, group_csv = make_demo_ml(out_dir)
    _log(f"读取表达矩阵: {expr_csv}")
    expr, y = read_inputs(expr_csv, group_csv, case=case, ctrl=ctrl)
    X_all = expr.values.T.astype(float)          # (n_samples, n_genes)
    _log(f"样本 {X_all.shape[0]} × 基因 {X_all.shape[1]}  |  病例 {int(y.sum())} vs 对照 {int((y == 0).sum())}")

    # 特征预筛：ANOVA F 取 top n_keep
    selector = SelectKBest(f_classif, k=min(n_keep, X_all.shape[1]))
    X = selector.fit_transform(X_all, y)
    sel_genes = expr.index[selector.get_support()]
    _log(f"单变量预筛后保留 {X.shape[1]} 个基因（ANOVA F top {min(n_keep, X_all.shape[1])}）")

    # 8 种模型 5 折交叉验证 OOF 概率
    skf = StratifiedKFold(n_splits=cv, shuffle=True, random_state=seed)
    rows, rocs = [], {}
    for name, clf in get_models().items():
        try:
            proba = cross_val_predict(clf, X, y, cv=skf, method="predict_proba", n_jobs=-1)[:, 1]
            auroc = roc_auc_score(y, proba)
            rows.append({"model": name, "mean_AUC": round(auroc, 4), "n_features": X.shape[1]})
            fpr, tpr, _ = roc_curve(y, proba)
            rocs[name] = (fpr, tpr, auroc)
            _log(f"  {name:<16} AUC = {auroc:.4f}")
        except Exception as e:
            _log(f"  {name} 失败: {e}")
    scores = pd.DataFrame(rows).sort_values("mean_AUC", ascending=False).reset_index(drop=True)
    scores.to_csv(os.path.join(out_dir, "ml_model_scores.csv"), index=False)

    # ---- ROC 图（Top 6 模型，按 AUC 降序；最佳模型加粗高亮） ----
    import matplotlib.pyplot as plt

    from plot_style import DOWN, UP, apply_style, save_fig

    apply_style()
    top = scores.head(6)
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    cmap = plt.get_cmap("YlGnBu")
    for rank, name in enumerate(top.model):
        fpr, tpr, a = rocs[name]
        if rank == 0:                       # 最佳模型：橙色加粗高亮
            color, lw, alpha = UP, 2.8, 1.0
        else:
            color, lw, alpha = cmap(0.32 + 0.12 * rank), 1.6, 0.85
        ax.plot(fpr, tpr, lw=lw, alpha=alpha, color=color,
                label=f"{name} (AUC={a:.3f})")
    ax.plot([0, 1], [0, 1], ls="--", color="#9AA1AC", lw=1)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("Machine Learning — ROC (5-fold CV, out-of-fold)")
    ax.legend(loc="lower right")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.05)
    save_fig(fig, os.path.join(out_dir, "ml_roc.png"))

    # ---- AUC 条形图（最优高亮橙，其余蓝；数值标签） ----
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    x = np.arange(len(scores))
    colors = [UP if i == 0 else DOWN for i in range(len(scores))]
    ax.bar(x, scores.mean_AUC, color=colors, alpha=0.92, width=0.62)
    ax.set_xticks(x)
    ax.set_xticklabels(scores.model, rotation=28, ha="right")
    ax.set_ylabel("Mean ROC-AUC (5-fold CV)")
    ax.set_ylim(0.4, 1.06)
    ax.set_title("Model comparison — 8 classifiers")
    ax.axhline(0.5, color="#9AA1AC", lw=0.8, ls="--")
    for xi, v in zip(x, scores.mean_AUC):
        ax.text(xi, v + 0.012, f"{v:.3f}", ha="center", fontsize=8.5)
    save_fig(fig, os.path.join(out_dir, "ml_auc_bar.png"))

    # ---- 特征重要性（RandomForest 在预筛基因上） ----
    rf = RandomForestClassifier(n_estimators=500, random_state=seed, n_jobs=-1)
    rf.fit(X, y)
    imp = pd.DataFrame({"gene": sel_genes, "importance": rf.feature_importances_})
    imp = imp.sort_values("importance", ascending=False).reset_index(drop=True)
    imp.to_csv(os.path.join(out_dir, "ml_feature_importance.csv"), index=False)

    top20 = imp.head(20)
    fig, ax = plt.subplots(figsize=(7.2, 7.6))
    ypos = range(len(top20))[::-1]
    ax.barh(ypos, top20.importance, color=UP, alpha=0.9, height=0.66)
    ax.set_yticks(list(ypos))
    ax.set_yticklabels(top20.gene)
    ax.set_xlabel("Random Forest feature importance")
    ax.set_title("Top 20 candidate genes")
    ax.grid(axis="x", ls=":", lw=0.5, color="#E3E6EC", zorder=0)
    ax.set_axisbelow(True)
    save_fig(fig, os.path.join(out_dir, "ml_feature_importance.png"))

    # ---- 报告 ----
    top_genes = top20.gene.head(20).tolist()
    lines = [
        "# 机器学习特征基因筛选报告（LKFlow2 Free-Basic · 经典转录组 ml）",
        "",
        f"- 输入: {expr_csv}",
        f"- 分组: {case} vs {ctrl}（{int(y.sum())} vs {int((y == 0).sum())} 样本）",
        f"- 特征预筛: ANOVA F top {X.shape[1]}",
        f"- 建模: 8 种分类器 × {cv} 折交叉验证（out-of-fold ROC-AUC）",
        "",
        "## 模型表现（按平均 AUC 排序）",
        "",
        "| model | mean_AUC |",
        "|---|---|",
    ]
    for _, r in scores.iterrows():
        lines.append(f"| {r.model} | {r.mean_AUC:.4f} |")
    lines += [
        "",
        f"## Top 20 候选基因（RandomForest 特征重要性）",
        "",
    ]
    for i, g in enumerate(top_genes, 1):
        lines.append(f"{i}. {g}")
    lines += [
        "",
        "## 文件清单",
        "- 图: ml_roc.png / ml_auc_bar.png / ml_feature_importance.png",
        "- 表: ml_model_scores.csv / ml_feature_importance.csv",
        "",
        "## 下一步（进阶功能请在 LKFlow2 完整软件中体验）",
        "- 机器学习-10/127/8/14/207 大规模算法组合自动寻优、Lasso 特征筛选、SHAP 解释图",
        "- 预后生存模型(6.7 机器学习-100)、外部数据集独立验证",
        "- 与 deg(差异) → venn(交集) 串联做「DEG∩ML 特征」核心基因",
    ]
    with open(os.path.join(out_dir, "summary_report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    _log(f"✔ 完成！结果目录: {os.path.abspath(out_dir)}")
    return out_dir


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--expr", default=None)
    ap.add_argument("--group", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--case", default="case")
    ap.add_argument("--ctrl", default="ctrl")
    ap.add_argument("--n-keep", type=int, default=150)
    a = ap.parse_args()
    if not a.demo and not (a.expr and a.group):
        ap.error("请用 --demo 或同时提供 --expr 与 --group")
    run_ml(demo=a.demo, expr_csv=a.expr, group_csv=a.group, out_dir=a.out,
           case=a.case, ctrl=a.ctrl, n_keep=a.n_keep)

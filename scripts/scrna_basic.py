#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LKFlow2 Free-Basic · 单细胞 scRNA-seq 全自动基础流程
====================================================
流程: 读取(demo/自定义) → QC 质控 → 归一化 → 高变基因 → PCA → 聚类 → UMAP/tSNE
      → 每群 marker 基因 → 出图出表 → 一键汇总报告

输入支持: 10x 目录 / 10x h5 / h5ad / mtx 目录 / csv,txt 表达矩阵（基因×细胞 或 细胞×基因）
输出: PNG/PDF 图 + CSV 表 + summary_report.md
"""
import argparse
import os
import sys
import textwrap

import numpy as np
import pandas as pd

OUT_SUBDIR = "scrna"
FIG_DPI = 150


def _log(msg):
    print(f"[LKFlow-Free/scrna] {msg}", flush=True)


def read_adata(input_path=None, demo=False, sample_name="sample"):
    """读取数据：--demo 下载 PBMC3K；否则按扩展名/目录结构自动识别。"""
    import scanpy as sc

    if demo:
        _log("下载演示数据 PBMC 3K（约 2700 细胞，需联网，仅首次）…")
        adata = sc.datasets.pbmc3k()
        adata.obs["sample"] = "pbmc3k"
        return adata, "PBMC3K-demo"

    if not input_path:
        raise ValueError("请提供 --input 数据路径，或使用 --demo 体验演示数据")

    p = os.path.abspath(input_path)
    if not os.path.exists(p):
        raise FileNotFoundError(f"输入路径不存在: {p}")
    _log(f"读取输入: {p}")
    if os.path.isdir(p):
        # 10x 目录
        if os.path.exists(os.path.join(p, "matrix.mtx")) or os.path.exists(os.path.join(p, "matrix.mtx.gz")):
            adata = sc.read_10x_mtx(p, var_names="gene_symbols", make_unique=True)
        elif any(f.endswith(".h5") for f in os.listdir(p)):
            h5 = [f for f in os.listdir(p) if f.endswith(".h5")][0]
            adata = sc.read_10x_h5(os.path.join(p, h5))
        else:
            raise ValueError("无法识别的目录结构：需要 10x mtx 目录（matrix.mtx）或目录内 .h5 文件")
    else:
        ext = os.path.splitext(p)[1].lower()
        if ext == ".h5":
            adata = sc.read_10x_h5(p)
        elif ext == ".h5ad":
            adata = sc.read_h5ad(p)
        elif ext in (".csv", ".txt", ".tsv", ".gz"):
            df = pd.read_csv(p, sep=None, engine="python", index_col=0)
            adata = sc.AnnData(df.T.astype(np.float32))  # 默认 基因×细胞 则转置
        else:
            raise ValueError(f"不支持的扩展名: {ext}（支持目录/h5/h5ad/csv/txt/mtx）")

    if adata.n_vars == 0 or adata.n_obs == 0:
        raise ValueError("数据为空，请检查输入文件格式")
    adata.obs["sample"] = sample_name
    _log(f"载入完成: {adata.n_obs} 细胞 × {adata.n_vars} 基因")
    return adata, sample_name


def run_qc(adata, min_genes=200, max_mt=20.0):
    """过滤 + 质控评分。返回 QC 汇总表。"""
    import scanpy as sc

    _log("计算质控指标（基因数/UMI/线粒体占比）…")
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    before = adata.n_obs
    keep = (adata.obs["n_genes_by_counts"] >= min_genes) & (adata.obs["pct_counts_mt"] <= max_mt)
    adata = adata[keep.values, :].copy()
    after = adata.n_obs
    _log(f"QC 过滤: {before} → {after} 细胞（基因数≥{min_genes}，线粒体≤{max_mt}%）")
    if after < 50:
        raise RuntimeError(f"过滤后仅剩 {after} 细胞，请调低 --min-genes 或调高 --max-mt 重试")
    return adata, pd.DataFrame({"stage": ["before", "after"], "n_cells": [before, after]})


def run_pipeline(adata, qc_stats, out_dir, sample_name, n_top_markers=20):
    import scanpy as sc

    sc.settings.figdir = out_dir
    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=FIG_DPI, dpi_save=300, frameon=False, figsize=(6, 5))

    # --- 质控小提琴图 ---
    _log("出 QC 图…")
    for metric, ylab in [("n_genes_by_counts", "基因数"), ("total_counts", "UMI 总数"),
                         ("pct_counts_mt", "线粒体 %")]:
        try:
            sc.pl.violin(adata, metric, groupby="sample", rotation=90, show=False,
                         save=f"_{metric}_qc.png")
        except Exception as e:
            _log(f"QC 图 {metric} 跳过: {e}")

    # --- 标准流程 ---
    _log("归一化 + log1p + 高变基因…")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", subset=False)

    _log("PCA 降维…")
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

    _log("UMAP + tSNE…")
    sc.tl.umap(adata)
    try:
        sc.tl.tsne(adata, n_pcs=30, perplexity=min(30, adata.n_obs - 1))
        have_tsne = True
    except Exception as e:
        have_tsne = False
        _log(f"tSNE 跳过（可忽略）: {e}")

    _log("Leiden 聚类（分辨率 1.0）…")
    _cluster(adata)

    n_clusters = adata.obs["leiden"].nunique()
    _log(f"得到 {n_clusters} 个细胞群")

    # --- 图 ---
    color_use = None
    _log("出 UMAP / tSNE 聚类图…")
    try:
        sc.pl.umap(adata, color="leiden", palette="tab20", show=False,
                   save="_leiden.png", title=f"{sample_name} Leiden clusters")
        sc.pl.umap(adata, color="leiden", palette="tab20", show=False,
                   save="_leiden.pdf", title=f"{sample_name} Leiden clusters")
    except Exception as e:
        _log(f"UMAP 图失败: {e}")
    if have_tsne:
        try:
            sc.pl.tsne(adata, color="leiden", palette="tab20", show=False, save="_leiden.png")
        except Exception:
            pass

    # --- marker 基因 ---
    _log("每群 marker 基因（Wilcoxon 秩和，耗时与细胞数相关）…")
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", n_genes=n_top_markers)
    de = sc.get.rank_genes_groups_df(adata, group=None)
    de.to_csv(os.path.join(out_dir, "cluster_markers_all.csv"), index=False)

    # 每群 top marker
    top = de.sort_values("pvals_adj").groupby("group").head(20).copy()
    top.to_csv(os.path.join(out_dir, "cluster_markers_top20.csv"), index=False)

    # DotPlot（每群 top3）
    _log("出 marker DotPlot / 热图…")
    n_show = min(8, n_clusters)
    groups_show = list(adata.obs["leiden"].cat.categories)[:n_show]
    markers = {}
    for g in groups_show:
        sub = de[de["group"] == g].sort_values("pvals_adj")
        markers[g] = sub["names"].head(3).tolist()
    genes = list(dict.fromkeys([g for v in markers.values() for g in v]))
    if genes:
        try:
            sc.pl.dotplot(adata, var_names=genes, groupby="leiden", show=False,
                          save="_top_marker_dotplot.png")
        except Exception as e:
            _log(f"DotPlot 跳过: {e}")
        try:
            sc.pl.heatmap(adata, var_names=genes, groupby="leiden", swap_axes=True,
                          show=False, save="_top_marker_heatmap.png")
        except Exception as e:
            _log(f"Heatmap 跳过: {e}")

    # 群占比堆叠条
    _log("出细胞群占比图…")
    _stacked_bar(adata, out_dir)

    # --- 表 ---
    qc_stats.to_csv(os.path.join(out_dir, "qc_filter_stats.csv"), index=False)
    membership = adata.obs[["leiden"]].copy()
    membership.to_csv(os.path.join(out_dir, "cluster_membership.csv"))
    sc.write(os.path.join(out_dir, "result.h5ad"), adata)

    _write_summary(out_dir, sample_name, adata, n_top_markers)
    return out_dir


def _cluster(adata):
    """聚类：优先 leidenalg，退回 louvain。"""
    try:
        import scanpy as sc
        sc.tl.leiden(adata, resolution=1.0, key_added="leiden", flavor="leidenalg")
    except Exception as e1:
        try:
            import scanpy as sc
            sc.tl.leiden(adata, resolution=1.0, key_added="leiden", flavor="igraph")
        except Exception as e2:
            import scanpy as sc
            sc.tl.louvain(adata, resolution=1.0, key_added="leiden")
    adata.obs["leiden"] = adata.obs["leiden"].astype(str).astype("category")


def _stacked_bar(adata, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # 中文字体回退，避免标签显示为方框
    plt.rcParams["font.sans-serif"] = ["PingFang SC", "Hiragino Sans GB",
                                       "Heiti SC", "Arial Unicode MS", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    counts = adata.obs["leiden"].value_counts().sort_index()
    frac = counts / counts.sum()
    fig, ax = plt.subplots(figsize=(max(4, 1 + 0.35 * len(counts)), 4))
    bars = ax.bar(range(len(frac)), frac.values * 100, color=plt.cm.tab20.colors[: len(frac)])
    ax.set_xticks(range(len(frac)))
    ax.set_xticklabels(frac.index, rotation=45, ha="right")
    ax.set_ylabel("细胞占比 (%)")
    ax.set_title("Clusters composition")
    for b, f in zip(bars, frac.values):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.5, f"{f * 100:.1f}%",
                ha="center", fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out_dir, f"cluster_composition.{ext}"), dpi=FIG_DPI)
    plt.close(fig)


def _write_summary(out_dir, sample_name, adata, n_top_markers):
    n_clusters = adata.obs["leiden"].nunique()
    files = sorted(os.listdir(out_dir))
    imgs = [f for f in files if f.endswith((".png", ".pdf"))]
    tables = [f for f in files if f.endswith(".csv")]
    report = textwrap.dedent(f"""\
        # {sample_name} · 单细胞基础分析报告（LKFlow2 Free-Basic）

        - 细胞数: {adata.n_obs}
        - 基因数: {adata.n_vars}
        - 聚类数: {n_clusters}
        - 聚类方法: Leiden (resolution=1.0)

        ## 结果文件
        - 图: {', '.join(imgs)}
        - 表: {', '.join(tables)}
        - 完整对象: result.h5ad（可用 scanpy 继续分析）

        ## 下一步
        - 本结果为免费基础版输出。自动细胞注释(SingleR/Azimuth/LLM 六法)、细胞通讯
          (CellChat/NicheNet)、轨迹(Monocle)、CNV 恶性细胞、虚拟敲除/扰动、药物预测
          等进阶分析，请使用 **LKFlow2 完整软件**（支持一键全自动流水线）。
        """)
    with open(os.path.join(out_dir, "summary_report.md"), "w", encoding="utf-8") as fh:
        fh.write(report)


def run_scrna(demo=False, input_path=None, sample_name="sample", out_dir=None,
              min_genes=200, max_mt=20.0, n_top_markers=20):
    _log("LKFlow2 Free-Basic · 单细胞 scRNA-seq 全自动分析启动")
    if not demo and not input_path:
        print("用法: python3 scripts/lkflow_free.py scrna --demo    # 演示数据")
        print("      python3 scripts/lkflow_free.py scrna --input <数据路径> --sample-name 名字")
        sys.exit(2)

    base = out_dir or os.path.join("lkflow_out", f"{OUT_SUBDIR}_{sample_name}")
    os.makedirs(base, exist_ok=True)

    adata, real_name = read_adata(input_path=input_path, demo=demo, sample_name=sample_name)
    adata, qc_stats = run_qc(adata, min_genes=min_genes, max_mt=max_mt)
    run_pipeline(adata, qc_stats, base, real_name, n_top_markers=n_top_markers)
    _log(f"✔ 完成！结果目录: {os.path.abspath(base)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--input", default=None)
    ap.add_argument("--sample-name", default="sample")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    run_scrna(demo=args.demo, input_path=args.input, sample_name=args.sample_name, out_dir=args.out)

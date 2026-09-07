#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LKFlow2 Free-Basic · 空间转录组基础分析
========================================
流程: 加载(Visium 10x 目录 / h5ad / 合成演示) → QC → 归一化 → PCA → 聚类
      → UMAP + 空间定位图(Spatial) → marker → 汇总报告

说明: 免费基础版覆盖"读取 → 质控 → 聚类 → 空间可视化"。
      SpaGCN/BANKSY 空间域、RCTD/CARD/Cell2location 等 10 种反卷积、MISTy、
      空间通讯与统计检验等进阶功能请使用 LKFlow2 完整软件。
"""
import argparse
import os
import sys
import textwrap

import numpy as np
import pandas as pd

FIG_DPI = 150


def _log(msg):
    print(f"[LKFlow-Free/spatial] {msg}", flush=True)


def _synth_demo(rows=60, cols=40, n_domains=6, n_genes=900):
    """合成演示：60×40 网格按行切成 6 个空间域，每域约 30 个基因特异上调。"""
    rng = np.random.default_rng(2026)
    xs, ys = np.meshgrid(np.arange(cols), np.arange(rows))
    xs, ys = xs.ravel(), ys.ravel()
    domain = np.minimum((ys // (rows // n_domains)).astype(int), n_domains - 1)

    gene_names = [f"GENE{i}" for i in range(1, n_genes + 1)]
    base = rng.gamma(shape=2.0, scale=1.0, size=n_genes)  # 基因基线表达
    up_genes = {}
    for d in range(n_domains):
        up_genes[d] = rng.choice(np.arange(n_genes), size=30, replace=False)

    counts = np.zeros((len(xs), n_genes), dtype=np.float32)
    for i in range(len(xs)):
        lam = base.copy()
        lam[up_genes[domain[i]]] *= 6.0
        counts[i] = rng.poisson(lam).astype(np.float32)

    import scanpy as sc
    adata = sc.AnnData(counts)
    adata.var_names = gene_names
    adata.obs["domain"] = [f"D{d}" for d in domain]
    adata.obs["x"] = xs
    adata.obs["y"] = ys
    # spatial 坐标（像素）
    adata.obsm["spatial"] = np.column_stack([xs * 50.0, ys * 50.0])
    adata.uns["spatial"] = {
        "synth_demo": {"scalefactors": {"tissue_hires_scalef": 1.0, "spot_diameter_fullres": 40.0}}
    }
    return adata


def _read_spatial(input_path=None, demo=False, sample_name="spatial_sample"):
    import scanpy as sc

    if demo:
        _log("生成合成演示数据（60×40 网格 · 6 个空间域）…")
        return _synth_demo(), "spatial-demo"

    p = os.path.abspath(input_path)
    if not os.path.exists(p):
        raise FileNotFoundError(f"输入路径不存在: {p}")
    _log(f"读取空间数据: {p}")
    if os.path.isdir(p):
        # Visium spaceranger outs 目录
        h5_candidates = [f for f in os.listdir(p) if f.endswith(".h5") and "filtered" in f]
        if not h5_candidates:
            h5_candidates = [f for f in os.listdir(p) if f.endswith(".h5")]
        if not h5_candidates:
            raise ValueError("Visium 目录中未找到 filtered_feature_bc_matrix.h5")
        adata = sc.read_visium(p, count_file=h5_candidates[0])
    else:
        ext = os.path.splitext(p)[1].lower()
        if ext == ".h5ad":
            adata = sc.read_h5ad(p)
        elif ext in (".csv", ".txt", ".tsv"):
            raise ValueError("CSV 空间数据缺少空间坐标，请提供 10x Visium 目录或带 spatial 信息的 .h5ad")
        else:
            raise ValueError(f"不支持的扩展名: {ext}（支持 Visium 10x 目录 / .h5ad）")

    if adata.n_obs == 0:
        raise ValueError("数据为空")
    _log(f"载入完成: {adata.n_obs} spots × {adata.n_vars} 基因")
    return adata, sample_name


def run_spatial(demo=False, input_path=None, sample_name="spatial_sample", out_dir=None):
    _log("LKFlow2 Free-Basic · 空间转录组基础分析启动")
    if not demo and not input_path:
        print("用法: python3 scripts/lkflow_free.py spatial --demo    # 合成演示数据")
        print("      python3 scripts/lkflow_free.py spatial --input <Visium 10x 目录或.h5ad> --sample-name 名字")
        sys.exit(2)

    import scanpy as sc
    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=FIG_DPI, dpi_save=300, frameon=False, figsize=(6, 5))

    base = out_dir or os.path.join("lkflow_out", f"spatial_{sample_name}")
    os.makedirs(base, exist_ok=True)
    sc.settings.figdir = base

    adata, real_name = _read_spatial(input_path=input_path, demo=demo, sample_name=sample_name)
    adata.obs["sample"] = real_name

    # QC
    _log("质控（基因数/UMI）…")
    if "total_counts" not in adata.obs:
        sc.pp.calculate_qc_metrics(adata, log1p=False, inplace=True)
    before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=200)
    qc = pd.DataFrame({"stage": ["before", "after"], "n_spots": [before, adata.n_obs]})
    qc.to_csv(os.path.join(base, "qc_stats.csv"), index=False)
    _log(f"QC 过滤: {before} → {adata.n_obs} spots")
    if adata.n_obs < 30:
        raise RuntimeError("过滤后 spot 过少，请检查输入数据")

    # 标准流程
    _log("归一化 + 高变基因 + PCA + 聚类…")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", subset=False)
    sc.tl.pca(adata, n_comps=30, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata)
    _cluster(adata)
    n_clusters = adata.obs["leiden"].nunique()
    _log(f"得到 {n_clusters} 个空间聚类")

    # 图
    _log("出 UMAP 与 空间定位图…")
    try:
        sc.pl.umap(adata, color="leiden", palette="tab20", show=False, save="_leiden.png")
    except Exception as e:
        _log(f"UMAP 图跳过: {e}")
    spatial_ok = False
    try:
        sc.pl.spatial(adata, color="leiden", palette="tab20", show=False,
                      save="_spatial_leiden.png", spot_size=1.2)
        spatial_ok = True
    except Exception as e:
        _log(f"sc.pl.spatial 不可用（缺少组织图像?）→ 改用坐标散点兜底: {e}")
    if not spatial_ok:
        # 兜底：无组织图像时，直接用 spot 坐标画聚类散点图（同样能展示空间定位）
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            plt.rcParams["font.sans-serif"] = ["PingFang SC", "Hiragino Sans GB",
                                               "Heiti SC", "Arial Unicode MS", "DejaVu Sans"]
            plt.rcParams["axes.unicode_minus"] = False
            if "spatial" not in adata.obsm:
                raise RuntimeError("obsm 中无 spatial 坐标")
            coords = np.asarray(adata.obsm["spatial"], dtype=float)
            cats = list(adata.obs["leiden"].cat.categories)
            cmap = plt.cm.tab20.colors
            fig, ax = plt.subplots(figsize=(6.5, 5.2))
            for i, c in enumerate(cats):
                sel = adata.obs["leiden"] == c
                ax.scatter(coords[sel.values, 0], coords[sel.values, 1],
                           s=6, color=cmap[i % len(cmap)], label=str(c), linewidths=0)
            ax.set_aspect("equal", adjustable="datalim")
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title("Spatial clusters (demo)", fontsize=12)
            ax.legend(title="cluster", fontsize=8, markerscale=2, loc="center left",
                      bbox_to_anchor=(1.0, 0.5), frameon=False)
            fig.tight_layout()
            for ext in ("png", "pdf"):
                fig.savefig(os.path.join(base, f"spatial_clusters.{ext}"), dpi=FIG_DPI)
            plt.close(fig)
            _log("✔ 已生成空间定位兜底图 spatial_clusters.png")
        except Exception as e2:
            _log(f"空间定位图兜底失败: {e2}")

    # marker（简化为每群 top 基因，便于免费版快速体验）
    _log("每群 marker 基因…")
    try:
        sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", n_genes=20)
        de = sc.get.rank_genes_groups_df(adata, group=None)
        de.to_csv(os.path.join(base, "cluster_markers_top20.csv"), index=False)
        n_show = min(8, n_clusters)
        markers = {}
        for g in list(adata.obs["leiden"].cat.categories)[:n_show]:
            sub = de[de["group"] == g].sort_values("pvals_adj")
            markers[g] = sub["names"].head(3).tolist()
        genes = list(dict.fromkeys([g for v in markers.values() for g in v]))
        if genes:
            sc.pl.dotplot(adata, var_names=genes, groupby="leiden", show=False,
                          save="_top_marker_dotplot.png")
    except Exception as e:
        _log(f"marker 分析跳过: {e}")

    # 表
    adata.obs[["leiden"]].to_csv(os.path.join(base, "cluster_membership.csv"))
    try:
        adata.write(os.path.join(base, "result.h5ad"))
    except Exception:
        pass

    files = sorted(os.listdir(base))
    report = textwrap.dedent(f"""\
        # {real_name} · 空间转录组基础分析报告（LKFlow2 Free-Basic）

        - Spots: {adata.n_obs}
        - 基因数: {adata.n_vars}
        - 聚类数: {n_clusters}

        ## 结果文件
        - 图: {', '.join(f for f in files if f.endswith(('.png', '.pdf')))}
        - 表: {', '.join(f for f in files if f.endswith('.csv'))}

        ## 下一步（进阶功能请在 LKFlow2 完整软件中体验）
        - 空间域识别（SpaGCN/BANKSY）、SVG 空间可变基因
        - 细胞反卷积（RCTD/CARD/Cell2location 等 10 法，可用单细胞注释结果做参考定位细胞类型）
        - 区域功能与共表达、Moran's I/LISA 空间统计、MISTy/空间通讯/邻域分析
        """)
    with open(os.path.join(base, "summary_report.md"), "w", encoding="utf-8") as fh:
        fh.write(report)

    _log(f"✔ 完成！结果目录: {os.path.abspath(base)}")


def _cluster(adata):
    import scanpy as sc
    try:
        sc.tl.leiden(adata, resolution=0.5, key_added="leiden", flavor="leidenalg")
    except Exception:
        try:
            sc.tl.leiden(adata, resolution=0.5, key_added="leiden", flavor="igraph")
        except Exception:
            sc.tl.louvain(adata, resolution=0.5, key_added="leiden")
    adata.obs["leiden"] = adata.obs["leiden"].astype(str).astype("category")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--input", default=None)
    ap.add_argument("--sample-name", default="spatial_sample")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    run_spatial(demo=args.demo, input_path=args.input, sample_name=args.sample_name, out_dir=args.out)

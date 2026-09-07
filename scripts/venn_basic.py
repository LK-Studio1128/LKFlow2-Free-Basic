#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LKFlow2 Free-Basic · 多集合交集 VENN（经典转录组模块 4.x 的免费基础版）
======================================================================
零代码韦恩图与交集提取（纯 Python）：
  - 输入：2~3 个基因列表文件（txt/csv，一行一个基因；或首列含表头亦可）
        或内置 --demo 合成三组列表（DEG / 网药靶点 / 机器学习特征 风格）
  - 输出：韦恩图(venn2/venn3) / 各区域基因明细 CSV / 核心交集 CSV / 汇总报告
  - 支持任意大小写基因名自动去重、交集自动导出 interGenes.List.txt（与完整软件 4.x 同名兼容）

说明: 免费基础版支持 2~3 集合韦恩图。
      VENN+1 二次追加交集、VENN-Auto 批量自动交集、4 集合以上欧拉图、
      与单细胞/MR 模块的跨组学自动交集请在 LKFlow2 完整软件中体验。
"""
import argparse
import os

import numpy as np
import pandas as pd

FIG_DPI = 150


def _log(msg):
    print(f"[LKFlow-Free/venn] {msg}", flush=True)


# --------------------------------------------------------------------------
# 演示数据：三组基因列表（DEG up / 药物靶点 / ML 特征）
# --------------------------------------------------------------------------
def make_demo_venn(out_dir, seed=7):
    """生成 3 个演示基因列表文件，返回路径列表。"""
    os.makedirs(out_dir, exist_ok=True)
    rng = np.random.default_rng(seed)
    base = np.array([f"GENE{i:05d}" for i in range(1, 4001)])
    # DEG up 550 个；药物靶点 620 个；ML 特征 300 个，构造真实重叠
    deg = set(rng.choice(base, size=550, replace=False))
    drug = set(rng.choice(base, size=620, replace=False))
    # 让 ml 与 deg 有强重叠（真实流程 ML 用 DEG 输入）
    shared_de = list(deg)[:240]
    ml = set(rng.choice(np.setdiff1d(base, list(deg | drug)), size=80, replace=False))
    ml |= set(rng.choice(shared_de, size=180, replace=False))
    ml |= set(rng.choice(list(drug), size=40, replace=False))

    paths = []
    for name, genes in [("A_DEG_up.txt", deg), ("B_drug_target.txt", drug), ("C_ml_features.txt", ml)]:
        p = os.path.join(out_dir, name)
        pd.Series(sorted(genes)).to_csv(p, index=False, header=False)
        paths.append(p)
    _log(f"✔ 演示列表: {[os.path.basename(p) for p in paths]}")
    return paths


def read_gene_list(path):
    """读取基因列表（每行一个；容忍表头/空行/大小写）。返回大写规范化去重 set。"""
    raw = pd.read_csv(path, header=None, dtype=str)
    col = raw.iloc[:, 0].dropna().astype(str).str.strip()
    col = col[col.str.len() > 0]
    # 若首行是常见表头则剔除
    if len(col) and col.iloc[0].upper() in ("GENE", "GENES", "SYMBOL", "GENEID"):
        col = col.iloc[1:]
    return set(g.upper() for g in col)


# --------------------------------------------------------------------------
# 出图
# --------------------------------------------------------------------------
def plot_venn(sets, names, out_path):
    """2 集合用 venn2，3 集合用 venn3。Okabe-Ito 色盲友好配色。"""
    import matplotlib.pyplot as plt

    from plot_style import CAT, apply_style

    apply_style()
    try:
        from matplotlib_venn import venn2, venn3
    except Exception:
        raise RuntimeError("缺少 matplotlib-venn，请先运行: python3 scripts/lkflow_free.py install "
                           "（或 pip install matplotlib-venn）")

    palette = [CAT[5], CAT[4], CAT[0]]  # 橙 / 蓝 / 橙黄（Okabe-Ito）
    fig, ax = plt.subplots(figsize=(8, 7))

    if len(sets) == 2:
        a, b = sets[0], sets[1]
        v = venn2(subsets=(len(a - b), len(b - a), len(a & b)), set_labels=names, ax=ax)
        for i, c in enumerate(["10", "01", "11"]):
            if v.get_patch_by_id(c) is not None:
                v.get_patch_by_id(c).set_alpha(0.5)
                v.get_patch_by_id(c).set_color(palette[i if i < 2 else 1])
    elif len(sets) == 3:
        a, b, c = sets[0], sets[1], sets[2]
        only = [len(a - b - c), len(b - a - c), len(c - a - b),
                len((a & b) - c), len((a & c) - b), len((b & c) - a), len(a & b & c)]
        v = venn3(subsets=tuple(only), set_labels=names, ax=ax)
        for i, pid in enumerate(["100", "010", "001", "110", "101", "011", "111"]):
            if v.get_patch_by_id(pid) is not None:
                v.get_patch_by_id(pid).set_alpha(0.55)
                v.get_patch_by_id(pid).set_color(palette[i % 3])
    else:
        raise ValueError("免费版韦恩图支持 2~3 个集合；4+ 集合请用 LKFlow2 完整软件 VENN-Auto")

    # 圆内计数放大加粗、区域标签加深
    for txt in v.set_labels:
        if txt is not None:
            txt.set_fontsize(12)
            txt.set_fontweight("medium")
    for txt in v.subset_labels:
        if txt is not None:
            txt.set_fontsize(11)
            txt.set_fontweight("bold")
            txt.set_color("#263238")

    ax.set_title("Venn Diagram — overlapping genes")
    from plot_style import save_fig
    save_fig(fig, out_path)
    _log(f"✔ 韦恩图: {out_path}")


# --------------------------------------------------------------------------
# 主流程
# --------------------------------------------------------------------------
def run_venn(files=None, names=None, demo=False, out_dir=None):
    os.makedirs(out_dir or "./lkflow_out/venn", exist_ok=True)
    if not out_dir:
        out_dir = os.path.join(os.getcwd(), "lkflow_out", "venn")
    os.makedirs(out_dir, exist_ok=True)

    if demo:
        files = make_demo_venn(out_dir)
    if files is None or len(files) < 2:
        raise ValueError("需提供 ≥2 个基因列表文件（--demo 或 --files a.txt b.txt [c.txt]）")
    if len(files) > 3:
        raise ValueError("免费版韦恩图支持 2~3 个集合；4+ 集合请用 LKFlow2 完整软件 VENN-Auto")
    if names is None:
        names = [os.path.splitext(os.path.basename(f))[0] for f in files]
    names = [str(n).strip() for n in names]

    sets = [read_gene_list(f) for f in files]
    for i, (n, s) in enumerate(zip(names, sets)):
        _log(f"集合 {i+1}「{n}」: {len(s)} 个基因  <- {files[i]}")

    # 交集表
    rows = []
    core = None
    import itertools
    for k in range(2, len(sets) + 1):
        for comb in itertools.combinations(range(len(sets)), k):
            inter = set.intersection(*[sets[j] for j in comb])
            region = "∩".join(names[j] for j in comb)
            rows.append({"region": region, "n_genes": len(inter), "genes": ";".join(sorted(inter))})
            if k == len(sets):
                core = inter
    region_df = pd.DataFrame(rows)
    region_csv = os.path.join(out_dir, "venn_regions.csv")
    region_df.to_csv(region_csv, index=False)
    _log(f"✔ 区域交集表: {region_csv}")

    # 核心交集（全集合交集），与完整软件 4.x 同名
    if core is not None and len(core) > 0:
        core_list = os.path.join(out_dir, "interGenes.List.txt")
        pd.Series(sorted(core)).to_csv(core_list, index=False, header=False)
        _log(f"✔ 核心交集 {len(core)} 个基因: {core_list}")

    # 图
    png = os.path.join(out_dir, "venn.png")
    plot_venn(sets, names, png)

    # 报告
    lines = [
        "# VENN 交集分析报告（LKFlow2 Free-Basic · 经典转录组 venn）",
        "",
        "## 输入集合",
        "",
        "| # | 名称 | 基因数 | 来源文件 |",
        "|---|---|---|---|",
    ]
    for i, (n, f) in enumerate(zip(names, files)):
        lines.append(f"| {i+1} | {n} | {len(sets[i])} | {f} |")
    lines += ["", "## 交集区域", "", "| 区域 | 基因数 | 基因示例 |", "|---|---|---|"]
    for _, r in region_df.iterrows():
        example = ";".join(r.genes.split(";")[:6]) + (" …" if r.n_genes > 6 else "")
        lines.append(f"| {r.region} | {r.n_genes} | {example} |")
    lines += ["", "## 文件清单",
              f"- 图: venn.png",
              f"- 表: venn_regions.csv（全部区域基因明细）"]
    if core is not None and len(core) > 0:
        lines.append(f"- 核心交集: interGenes.List.txt（{len(core)} 个基因，可直接喂给下游模块/完整软件）")
    lines += ["",
              "## 下一步（进阶功能请在 LKFlow2 完整软件中体验）",
              "- VENN+1 二次追加、VENN-Auto 批量自动交集、4+ 集合欧拉图",
              "- DEG×网药×机器学习×MR 跨组学自动交集与韦恩图联动"]
    with open(os.path.join(out_dir, "summary_report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    _log(f"✔ 完成！结果目录: {os.path.abspath(out_dir)}")
    return out_dir


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--files", nargs="*", default=None, help="2~3 个基因列表文件")
    ap.add_argument("--names", nargs="*", default=None, help="集合显示名（与文件顺序对应）")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run_venn(files=a.files, names=a.names, demo=a.demo, out_dir=a.out)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LKFlow2 Free-Basic · 统一出版级绘图风格
========================================
基于 Nature/Science 系期刊通用图表规范与 Okabe-Ito 色盲友好调色板：
  - 去顶/右 边框，刻度朝内，legend 无框
  - Okabe-Ito 主色（上/下调、分类色），发散数据用 蓝-白-红
  - 200 DPI 保存、tight bbox
所有出图脚本在绘图前调用 apply_style() 即可统一风格。
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# ---- Okabe-Ito 色盲友好调色板 ----
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
             "#0072B2", "#D55E00", "#CC79A7", "#000000"]
UP = "#D55E00"      # 上调 / 显著+
DOWN = "#0072B2"    # 下调 / 显著-
NS = "#B9C0CC"      # 不显著
ACCENT = "#009E73"  # 通用强调（绿）
CAT = OKABE_ITO[:7]  # 分类图默认循环色

# 蓝-白-红 发散 colormap（z-score / 热图）
DIVERGING = ["#0072B2", "#FFFFFF", "#D55E00"]


def apply_style():
    """全局 rcParams：出版级基础风格。"""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.titleweight": "medium",
        "axes.labelsize": 10.5,
        "axes.linewidth": 0.9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.prop_cycle": plt.cycler(color=CAT),
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "legend.handlelength": 1.6,
        "lines.linewidth": 1.8,
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.08,
    })


def save_fig(fig, path):
    """统一保存：200 DPI + tight。"""
    fig.savefig(path, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)

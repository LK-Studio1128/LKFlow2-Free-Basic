#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LKFlow2 Free-Basic · 统一命令行入口
===================================
零代码生信基础分析（免费开源体验版）
  经典转录组线（对应完整软件模块 1–7 的基础版）
  - deg       差异表达 DEGs（2.1 基础版）：火山图 / Top 热图 / 差异表
  - venn      多集合交集（4.x 基础版）：2~3 集合韦恩图 / 交集基因导出
  - ml        机器学习特征筛选（6.x 基础版）：8 分类器对比 / ROC / Top 基因
  单细胞 / 空间 / MR（对应完整软件模块 8–10 的基础版）
  - scrna     单细胞 scRNA-seq 全自动流程（demo / 自定义输入）
  - spatial   空间转录组基础分析（demo / 自定义 Visium 输入）
  - mr        孟德尔随机化 MR 因果分析（暴露/结局 GWAS 摘要）
  - install   自动创建虚拟环境并安装依赖
  - all       一条命令跑完全流程（--demo）

用法示例:
  python3 scripts/lkflow_free.py install
  python3 scripts/lkflow_free.py all --demo
  python3 scripts/lkflow_free.py deg   --demo
  python3 scripts/lkflow_free.py venn  --demo
  python3 scripts/lkflow_free.py ml    --demo
  python3 scripts/lkflow_free.py scrna --demo
  python3 scripts/lkflow_free.py spatial --demo
  python3 scripts/lkflow_free.py mr-demo-data --out ./demo_mr
  python3 scripts/lkflow_free.py mr --exposure ./demo_mr/exposure.csv --outcome ./demo_mr/outcome.csv
"""
import argparse
import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS_DIR = os.path.join(REPO_ROOT, "scripts")
ENV_NAME = "lkflow_env"

# 分析依赖（MR 只需要前 4 个；单细胞/空间额外需要 scanpy 系列；
# deg/venn/ml 需要 scikit-learn / matplotlib-venn）
REQUIRED_CORE = ["numpy", "pandas", "scipy", "matplotlib"]
REQUIRED_EXTRA = ["scanpy", "h5py", "scikit-learn", "matplotlib-venn"]
# 单细胞聚类后端：优先 leidenalg（提供二进制轮子），失败则退回 python-louvain，仍失败则仅警告
CLUSTER_BACKENDS = [["leidenalg"], ["python-louvain"]]


def _find_env(venv=None):
    """在显式路径、仓库内 lkflow_env、用户目录 ~/.lkflow_env 中查找虚拟环境。"""
    candidates = []
    if venv:
        candidates.append(venv)
    candidates += [os.path.join(REPO_ROOT, ENV_NAME), os.path.join(os.path.expanduser("~"), "." + ENV_NAME)]
    for c in candidates:
        py = os.path.join(c, "bin", "python")
        if sys.platform.startswith("win"):
            py = os.path.join(c, "Scripts", "python.exe")
        if os.path.isfile(py):
            return c, py
    return None, None


def _venv_python(venv_dir):
    if sys.platform.startswith("win"):
        return os.path.join(venv_dir, "Scripts", "python.exe")
    return os.path.join(venv_dir, "bin", "python")


def cmd_install(args):
    """创建虚拟环境并安装全部依赖。"""
    venv_dir = args.venv or os.path.join(REPO_ROOT, ENV_NAME)
    print(f"[1/4] 准备虚拟环境: {venv_dir}")
    if not os.path.exists(_venv_python(venv_dir)):
        subprocess.check_call([sys.executable, "-m", "venv", venv_dir])
    py = _venv_python(venv_dir)

    pkgs = list(REQUIRED_CORE) + (list(REQUIRED_EXTRA) if not args.core_only else [])
    print(f"[2/4] 安装依赖: {' '.join(pkgs)}  （首次约需几分钟）")
    pip = [py, "-m", "pip", "install", "--upgrade"]
    subprocess.check_call(pip + ["pip"])
    if pkgs:
        subprocess.check_call(pip + pkgs)
    if not args.core_only:
        for backend in CLUSTER_BACKENDS:
            try:
                subprocess.check_call(pip + backend)
                break
            except subprocess.CalledProcessError:
                print(f"⚠ 聚类后端 {backend[0]} 安装失败，尝试下一个…")
        else:
            print("⚠ 未成功安装任何聚类后端（leidenalg/python-louvain），"
                  "单细胞与空间聚类步骤将不可用，可稍后重跑 install。")
    print(f"[3/4] 校验依赖…")
    check = subprocess.run([py, os.path.join(SCRIPTS_DIR, "deps_check.py")], capture_output=True, text=True)
    print(check.stdout)
    if check.returncode != 0:
        print(check.stderr)
        print("⚠ 部分依赖校验失败，见上方信息。")
        sys.exit(1)
    print(f"[4/4] ✔ 依赖就绪。虚拟环境: {venv_dir}")
    print("下一步示例:\n"
          "  python3 scripts/lkflow_free.py all --demo\n"
          "  python3 scripts/lkflow_free.py deg --demo     # 差异表达\n"
          "  python3 scripts/lkflow_free.py venn --demo    # 韦恩交集\n"
          "  python3 scripts/lkflow_free.py ml --demo      # 机器学习筛选\n"
          "  python3 scripts/lkflow_free.py scrna --demo\n"
          "  python3 scripts/lkflow_free.py spatial --demo\n"
          "  python3 scripts/lkflow_free.py mr-demo-data --out ./demo_mr && python3 scripts/lkflow_free.py mr --exposure ./demo_mr/exposure.csv --outcome ./demo_mr/outcome.csv")


def _need_rerun_in_env(argv, cmd):
    """若当前解释器缺少分析依赖且仓库/用户目录存在虚拟环境，则用虚拟环境重执行。"""
    if cmd == "install":
        return None
    try:
        import numpy  # noqa
        return None  # 当前环境已有依赖，直接运行
    except Exception:
        pass
    env_dir, env_py = _find_env(None)
    if env_dir and os.path.isfile(env_py):
        print(f"检测到虚拟环境 {env_dir}，自动切换到其中的 Python 执行…\n")
        return env_py
    return None


def main(argv=None):
    parser = argparse.ArgumentParser(prog="lkflow_free", description="LKFlow2 Free-Basic 生信基础分析（免费开源版）")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_install = sub.add_parser("install", help="创建虚拟环境并安装依赖")
    p_install.add_argument("--venv", default=None, help="虚拟环境目录（默认仓库内 lkflow_env 或 ~/.lkflow_env）")
    p_install.add_argument("--core-only", action="store_true", help="只装 MR 所需核心依赖（不装 scanpy）")

    p_deg = sub.add_parser("deg", help="差异表达 DEGs（经典转录组模块 2.1 基础版）")
    p_deg.add_argument("--demo", action="store_true", help="使用内置演示数据（病例 vs 对照）")
    p_deg.add_argument("--expr", default=None, help="表达矩阵 CSV（行=基因 列=样本）")
    p_deg.add_argument("--group", default=None, help="分组表 CSV（sample, group）")
    p_deg.add_argument("--case", default="case", help="病例组标签(默认 case)")
    p_deg.add_argument("--ctrl", default="ctrl", help="对照组标签(默认 ctrl)")
    p_deg.add_argument("--fdr", type=float, default=0.05, help="显著阈值 padj(默认0.05)")
    p_deg.add_argument("--lfc", type=float, default=1.0, help="显著阈值 |log2FC|(默认1)")
    p_deg.add_argument("--paired", action="store_true", help="配对 t 检验(默认 Welch 非配对)")
    p_deg.add_argument("--out", default=None, help="输出目录（默认 ./lkflow_out/deg）")

    p_venn = sub.add_parser("venn", help="多集合交集韦恩图（经典转录组模块 4.x 基础版）")
    p_venn.add_argument("--demo", action="store_true", help="使用内置三组演示列表")
    p_venn.add_argument("--files", nargs="*", default=None, help="2~3 个基因列表文件（一行一个基因）")
    p_venn.add_argument("--names", nargs="*", default=None, help="集合显示名（可选）")
    p_venn.add_argument("--out", default=None, help="输出目录（默认 ./lkflow_out/venn）")

    p_ml = sub.add_parser("ml", help="机器学习特征基因筛选（经典转录组模块 6.x 基础版）")
    p_ml.add_argument("--demo", action="store_true", help="使用内置演示数据（含 20 个真实信号基因）")
    p_ml.add_argument("--expr", default=None, help="表达矩阵 CSV（行=基因 列=样本）")
    p_ml.add_argument("--group", default=None, help="分组表 CSV（sample, group）")
    p_ml.add_argument("--case", default="case", help="病例组标签(默认 case)")
    p_ml.add_argument("--ctrl", default="ctrl", help="对照组标签(默认 ctrl)")
    p_ml.add_argument("--n-keep", type=int, default=150, help="ANOVA 预筛保留基因数(默认150)")
    p_ml.add_argument("--out", default=None, help="输出目录（默认 ./lkflow_out/ml）")

    p_scrna = sub.add_parser("scrna", help="单细胞 scRNA-seq 全自动基础分析")
    p_scrna.add_argument("--demo", action="store_true", help="使用内置演示数据（自动下载 PBMC 3K）")
    p_scrna.add_argument("--input", default=None, help="输入数据：10x 目录 / .h5 / .h5ad / mtx 目录 / csv,txt 矩阵")
    p_scrna.add_argument("--sample-name", default="sample", help="样本名（用于输出命名与报告）")
    p_scrna.add_argument("--out", default=None, help="输出目录（默认 ./lkflow_out/scrna_<样本名>）")
    p_scrna.add_argument("--min-genes", type=int, default=200, help="细胞最少基因数阈值(默认200)")
    p_scrna.add_argument("--max-mt", type=float, default=20.0, help="线粒体占比上限 %(默认20)")
    p_scrna.add_argument("--n-top-markers", type=int, default=20, help="每群保留 marker 数(默认20)")

    p_spatial = sub.add_parser("spatial", help="空间转录组基础分析")
    p_spatial.add_argument("--demo", action="store_true", help="使用内置合成演示数据（6 空间域网格）")
    p_spatial.add_argument("--input", default=None, help="10x Visium spaceranger 输出目录（含 filtered_feature_bc_matrix.h5 与 spatial/）或 .h5ad")
    p_spatial.add_argument("--sample-name", default="spatial_sample", help="样本名")
    p_spatial.add_argument("--out", default=None, help="输出目录（默认 ./lkflow_out/spatial_<样本名>）")

    p_mrdata = sub.add_parser("mr-demo-data", help="生成 MR 演示/模板 GWAS 摘要数据")
    p_mrdata.add_argument("--out", default="./demo_mr", help="输出目录")
    p_mrdata.add_argument("--n-snp", type=int, default=50, help="工具变量 SNP 数(默认50)")

    p_mr = sub.add_parser("mr", help="孟德尔随机化因果分析")
    p_mr.add_argument("--demo", action="store_true", help="使用内置演示 GWAS 数据一键体验")
    p_mr.add_argument("--exposure", default=None, help="暴露 GWAS 摘要 CSV")
    p_mr.add_argument("--outcome", default=None, help="结局 GWAS 摘要 CSV")
    p_mr.add_argument("--exposure-name", default="exposure", help="暴露名（报告/图标题）")
    p_mr.add_argument("--outcome-name", default="outcome", help="结局名（报告/图标题）")
    p_mr.add_argument("--out", default=None, help="输出目录（默认 ./lkflow_out/mr）")

    p_all = sub.add_parser("all", help="全流程体验（--demo 演示数据）")
    p_all.add_argument("--demo", action="store_true", help="使用演示数据跑六大流程")
    p_all.add_argument("--out", default=None, help="输出根目录（默认 ./lkflow_out）")

    args = parser.parse_args(argv)
    sys.path.insert(0, SCRIPTS_DIR)

    if args.cmd == "install":
        cmd_install(args)
        return

    env_py = _need_rerun_in_env(argv, args.cmd)
    if env_py:
        subprocess.check_call([env_py, os.path.abspath(__file__)] + argv)
        return

    if args.cmd == "deg":
        from deg_basic import run_deg
        run_deg(demo=args.demo, expr_csv=args.expr, group_csv=args.group, out_dir=args.out,
                case=args.case, ctrl=args.ctrl, fdr=args.fdr, lfc=args.lfc, paired=args.paired)
    elif args.cmd == "venn":
        from venn_basic import run_venn
        run_venn(files=args.files, names=args.names, demo=args.demo, out_dir=args.out)
    elif args.cmd == "ml":
        from ml_basic import run_ml
        run_ml(demo=args.demo, expr_csv=args.expr, group_csv=args.group, out_dir=args.out,
               case=args.case, ctrl=args.ctrl, n_keep=args.n_keep)
    elif args.cmd == "scrna":
        from scrna_basic import run_scrna
        run_scrna(demo=args.demo, input_path=args.input, sample_name=args.sample_name,
                  out_dir=args.out, min_genes=args.min_genes, max_mt=args.max_mt,
                  n_top_markers=args.n_top_markers)
    elif args.cmd == "spatial":
        from spatial_basic import run_spatial
        run_spatial(demo=args.demo, input_path=args.input, sample_name=args.sample_name, out_dir=args.out)
    elif args.cmd == "mr-demo-data":
        from mr_basic import make_demo_mr_data
        out = make_demo_mr_data(out_dir=args.out, n_snp=args.n_snp)
        print("\n下一步: python3 scripts/lkflow_free.py mr --exposure "
              f"{os.path.join(out, 'exposure.csv')} --outcome {os.path.join(out, 'outcome.csv')}")
    elif args.cmd == "mr":
        from mr_basic import run_mr
        if args.demo:
            from mr_basic import make_demo_mr_data
            demo_dir = make_demo_mr_data(out_dir="./demo_mr", n_snp=50)
            args.exposure = os.path.join(demo_dir, "exposure.csv")
            args.outcome = os.path.join(demo_dir, "outcome.csv")
        if not args.exposure or not args.outcome:
            parser.error("请用 --demo 或同时提供 --exposure 与 --outcome")
        run_mr(exposure_csv=args.exposure, outcome_csv=args.outcome,
               exposure_name=args.exposure_name, outcome_name=args.outcome_name, out_dir=args.out)
    elif args.cmd == "all":
        run_all(demo=args.demo, out_root=args.out)


def run_all(demo=True, out_root=None):
    base = out_root or os.path.join(REPO_ROOT, "lkflow_out")
    os.makedirs(base, exist_ok=True)

    if demo:
        print("\n" + "=" * 66 + "\n 第 1/6 步 · 差异表达 DEGs（经典转录组 2.1）\n" + "=" * 66)
        from deg_basic import run_deg
        run_deg(demo=True, out_dir=os.path.join(base, "deg"))

        print("\n" + "=" * 66 + "\n 第 2/6 步 · 多集合交集 VENN（经典转录组 4.x）\n" + "=" * 66)
        from venn_basic import run_venn
        run_venn(demo=True, out_dir=os.path.join(base, "venn"))

        print("\n" + "=" * 66 + "\n 第 3/6 步 · 机器学习筛选（经典转录组 6.x）\n" + "=" * 66)
        from ml_basic import run_ml
        run_ml(demo=True, out_dir=os.path.join(base, "ml"))

        print("\n" + "=" * 66 + "\n 第 4/6 步 · 单细胞 scRNA-seq 基础流程\n" + "=" * 66)
        from scrna_basic import run_scrna
        run_scrna(demo=True, input_path=None, sample_name="scrna_demo",
                  out_dir=os.path.join(base, "scrna"))

        print("\n" + "=" * 66 + "\n 第 5/6 步 · 空间转录组基础流程\n" + "=" * 66)
        from spatial_basic import run_spatial
        run_spatial(demo=True, input_path=None, sample_name="spatial_demo",
                    out_dir=os.path.join(base, "spatial"))

        print("\n" + "=" * 66 + "\n 第 6/6 步 · 孟德尔随机化 MR\n" + "=" * 66)
        from mr_basic import make_demo_mr_data, run_mr
        demo_dir = make_demo_mr_data(out_dir=os.path.join(base, "demo_mr_data"), n_snp=50)
        run_mr(exposure_csv=os.path.join(demo_dir, "exposure.csv"),
               outcome_csv=os.path.join(demo_dir, "outcome.csv"),
               exposure_name="DemoExposure", outcome_name="DemoOutcome",
               out_dir=os.path.join(base, "mr"))
    else:
        print("自定义数据模式：请分别运行 deg/venn/ml/scrna/spatial/mr 子命令（见 --help）")

    print("\n✔ 全部完成！结果目录: " + base)


if __name__ == "__main__":
    main(sys.argv[1:])

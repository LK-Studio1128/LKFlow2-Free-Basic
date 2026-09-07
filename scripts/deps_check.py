#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LKFlow2 Free-Basic · 依赖自检（供 install 子命令调用）"""
import importlib
import sys

CORE = ["numpy", "pandas", "scipy", "matplotlib"]
CLASSIC = ["sklearn", "matplotlib_venn"]        # deg / venn / ml
EXTRA = ["scanpy", "h5py"]                      # scrna / spatial


def _ver(m):
    return getattr(m, "__version__", "?")


def main():
    ok = True
    print("── 核心依赖 ──")
    for name in CORE:
        try:
            m = importlib.import_module(name)
            print(f"  ✔ {name} {_ver(m)}")
        except Exception as e:
            ok = False
            print(f"  ✘ {name} 缺失: {e}")
    print("── 经典转录组扩展(deg/venn/ml) ──")
    for name in CLASSIC:
        try:
            m = importlib.import_module(name)
            print(f"  ✔ {name} {_ver(m)}")
        except Exception as e:
            print(f"  ✘ {name} 缺失: {e}")
    print("── 单细胞/空间扩展 ──")
    have_scanpy = False
    for name in EXTRA:
        try:
            m = importlib.import_module(name)
            print(f"  ✔ {name} {_ver(m)}")
            if name == "scanpy":
                have_scanpy = True
        except Exception as e:
            print(f"  ✘ {name} 缺失: {e}")
    if have_scanpy:
        try:
            import leidenalg
            print(f"  ✔ 聚类后端 leidenalg {_ver(leidenalg)}")
        except Exception:
            try:
                import community  # python-louvain
                print(f"  ✔ 聚类后端 python-louvain {_ver(community)}")
            except Exception:
                print("  ⚠ 未找到聚类后端（leidenalg / python-louvain），单细胞聚类将不可用")
    print("✔ 依赖就绪" if ok else "✘ 存在缺失依赖")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

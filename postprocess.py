"""
Pintle Injector 後処理スクリプト
=================================
計算済みケースファイルから各種データを抽出・可視化する。
Matplotlib + PyFluent を組み合わせた後処理例。
"""

import ansys.fluent.core as pyfluent
from ansys.fluent.core import launch_fluent
import matplotlib.pyplot as plt
import numpy as np
import os

CASE_FILE = os.path.join(os.path.dirname(__file__), "output", "pintle_injector.cas.h5")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output", "postprocess")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# Fluent 起動 (結果読み込みのみ)
# ============================================================

print("=== Fluent 起動・ケース読み込み ===")
solver = launch_fluent(
    product_version="25.2.0",
    mode="solver",
    dimension=2,
    ui_mode="no_gui_or_graphics",
    processor_count=4,
    insecure_mode=True,
)
solver.settings.file.read_case_data(file_name=CASE_FILE)

# ============================================================
# 1. 軸上 (r=0) の変数プロファイル
# ============================================================

print("=== 軸上プロファイル抽出 ===")

# ライン (直線) サンプル定義: x=0 から x=0.2m (軸方向), r=0
# 座標系: x=軸方向, y=半径方向 (generate_mesh.py 準拠)
_surf = solver.settings.results.surfaces
_surf.line_surface.create("axis_line")
_surf.line_surface["axis_line"].p0 = [0.0, 0.0, 0.0]   # (x_start, y, z)
_surf.line_surface["axis_line"].p1 = [0.20, 0.0, 0.0]   # (x_end, y, z)

# XY プロットデータ出力
_plot = solver.settings.results.plot
_plot.xy_plot.create("axis_temperature")
_plot.xy_plot["axis_temperature"].y_axis_function = "total-temperature"
_plot.xy_plot["axis_temperature"].surfaces_list = ["axis_line"]
_plot.xy_plot["axis_temperature"].x_axis_function = "Direction Vector"
_plot.xy_plot["axis_temperature"].plot_direction.direction_vector.x_component = 1.0
_plot.xy_plot["axis_temperature"].plot_direction.direction_vector.y_component = 0.0

axis_temp_file = os.path.join(OUTPUT_DIR, "axis_temperature.xy")
_plot.xy_plot["axis_temperature"].write_to_file(file_name=axis_temp_file)

# ============================================================
# 2. 横断面 (各 z 位置) での半径方向プロファイル
# ============================================================

print("=== 断面プロファイル抽出 ===")

z_positions = [0.02, 0.05, 0.10, 0.15, 0.20]  # [m] サンプル z 位置

for z in z_positions:
    surface_name = f"cross_section_z{int(z*1000):03d}mm"
    _surf.iso_surface.create(surface_name)
    _surf.iso_surface[surface_name].field = "x-coordinate"   # 軸方向 = x
    _surf.iso_surface[surface_name].iso_values = [z]

    # 混合分率プロット (PDF テーブル未設定時はスキップ)
    try:
        xy_plot_name = f"radial_profile_z{int(z*1000):03d}mm"
        _plot.xy_plot.create(xy_plot_name)
        _plot.xy_plot[xy_plot_name].y_axis_function = "mean-mixture-fraction"
        _plot.xy_plot[xy_plot_name].surfaces_list = [surface_name]
        _plot.xy_plot[xy_plot_name].x_axis_function = "Direction Vector"
        _plot.xy_plot[xy_plot_name].plot_direction.direction_vector.x_component = 0.0
        _plot.xy_plot[xy_plot_name].plot_direction.direction_vector.y_component = 1.0
        out_file = os.path.join(OUTPUT_DIR, f"{xy_plot_name}.xy")
        _plot.xy_plot[xy_plot_name].write_to_file(file_name=out_file)
    except Exception as e:
        print(f"  混合分率プロットスキップ z={z:.2f}m: {e}")

# ============================================================
# 3. 混合効率計算
# ============================================================

print("=== 混合効率計算 ===")

# 化学量論混合分率 (LOX/RP-1 の場合)
# RP-1: C12H23.4, 理論空燃比 (O/F) ≈ 3.4 → f_stoich = 1/(1+3.4) ≈ 0.227
F_STOICH = 0.227

try:
    _rd = solver.settings.solution.report_definitions
    _rd.surface.create("outlet_mf_avg")
    _rd.surface["outlet_mf_avg"].report_type = "facet-average"
    _rd.surface["outlet_mf_avg"].field = "mean-mixture-fraction"
    _rd.surface["outlet_mf_avg"].surface_names = ["chamber_outlet"]
    f_avg = _rd.compute(report_defs=["outlet_mf_avg"])
    print(f"  出口平均混合分率: {f_avg}")
    print(f"  化学量論混合分率: {F_STOICH:.4f}")
    print(f"  当量比 (φ): {float(f_avg) / F_STOICH:.3f}")
except Exception as e:
    print(f"  混合効率計算スキップ (PDF テーブル未設定): {e}")
    f_avg = None

# ============================================================
# 4. Matplotlib によるプロファイル可視化
# ============================================================

print("=== Matplotlib 可視化 ===")

def read_fluent_xy(filepath):
    """Fluent XY プロットファイルを読み込む (コメント行スキップ)"""
    x_vals, y_vals = [], []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("(") or line.startswith(")") or not line:
                continue
            try:
                parts = line.split()
                x_vals.append(float(parts[0]))
                y_vals.append(float(parts[1]))
            except (ValueError, IndexError):
                continue
    return np.array(x_vals), np.array(y_vals)

# 軸上温度プロファイル
if os.path.exists(axis_temp_file):
    x, T = read_fluent_xy(axis_temp_file)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x * 1000, T, "r-", linewidth=2, label="Axis Total Temperature")
    ax.set_xlabel("z [mm]", fontsize=12)
    ax.set_ylabel("Total Temperature [K]", fontsize=12)
    ax.set_title("Pintle Injector — Axial Temperature Profile", fontsize=13)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "axis_temperature_profile.png"), dpi=150)
    plt.close(fig)
    print("  軸上温度プロファイル保存完了")

# 断面 混合分率プロファイル
fig, ax = plt.subplots(figsize=(8, 5))
colors = plt.cm.viridis(np.linspace(0, 1, len(z_positions)))

for z, color in zip(z_positions, colors):
    xy_file = os.path.join(OUTPUT_DIR, f"radial_profile_z{int(z*1000):03d}mm.xy")
    if os.path.exists(xy_file):
        r, mf = read_fluent_xy(xy_file)
        ax.plot(r * 1000, mf, color=color, linewidth=1.8, label=f"z = {z*1000:.0f} mm")

ax.axhline(F_STOICH, color="k", linestyle="--", linewidth=1.2, label=f"Stoich. (f={F_STOICH:.3f})")
ax.set_xlabel("r [mm]", fontsize=12)
ax.set_ylabel("Mean Mixture Fraction [-]", fontsize=12)
ax.set_title("Pintle Injector — Radial Mixture Fraction Profiles", fontsize=13)
ax.grid(True, linestyle="--", alpha=0.5)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "radial_mixture_fraction.png"), dpi=150)
plt.close(fig)
print("  断面混合分率プロファイル保存完了")

# ============================================================
# Fluent 終了
# ============================================================

solver.exit()
print(f"=== 後処理完了: 結果は {OUTPUT_DIR} に保存 ===")

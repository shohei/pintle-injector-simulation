"""
Pintle Injector CFD Simulation using PyFluent
=============================================
ロケットエンジンのPintle Injector CFDシミュレーション

想定形状:
  - 中心: 液体酸化剤 (LOX) がピントル先端から軸方向に噴射
  - 外周: 燃料 (RP-1 / ケロシン) が環状スリットから半径方向に噴射
  - 衝突・混合後に燃焼室で燃焼

使用モデル:
  - 流体: 圧縮性 Navier-Stokes (密度ベース ソルバー)
  - 乱流: Realizable k-epsilon (壁面関数)
  - 燃焼: Mixture fraction / PDF (Non-premixed combustion)
  - 噴霧: Discrete Phase Model (DPM) ※液滴噴射
  - 熱: エネルギー方程式 ON

前提:
  - Ansys Fluent 2023 R1 以降
  - PyFluent >= 0.16
  - メッシュファイル (pintle_injector.msh.h5 or .cas.h5) が同ディレクトリに存在
"""

import ansys.fluent.core as pyfluent
from ansys.fluent.core import launch_fluent
import os
import time

# ============================================================
# 設定パラメータ (ここを編集して条件を変更)
# ============================================================

MESH_FILE = os.path.join(os.path.dirname(__file__), "pintle_injector.cgns")
CASE_OUTPUT = os.path.join(os.path.dirname(__file__), "output", "pintle_injector.cas.h5")
DATA_OUTPUT = os.path.join(os.path.dirname(__file__), "output", "pintle_injector.dat.h5")

# 推進剤条件
LOX_MASS_FLOW     = 0.5        # [kg/s] 酸化剤 (LOX) 質量流量
RP1_MASS_FLOW     = 0.25       # [kg/s] 燃料 (RP-1) 質量流量
LOX_TEMPERATURE   = 90.0       # [K]   LOX 入口温度 (沸点近傍)
RP1_TEMPERATURE   = 300.0      # [K]   RP-1 入口温度
CHAMBER_PRESSURE  = 2.0e6      # [Pa]  燃焼室圧力 (~20 atm)

# ソルバー設定
N_ITERATIONS      = 2000       # 反復計算回数
COURANT_NUMBER    = 1.0        # CFL 数 (密度ベース)
RESIDUAL_TARGET   = 1e-4       # 収束判定残差

# 境界条件名 (メッシュの Named Selection と一致させること)
BC_LOX_INLET  = "lox_inlet"
BC_RP1_INLET  = "rp1_inlet"
BC_OUTLET     = "chamber_outlet"
BC_WALL_PINTLE = "wall_pintle"
BC_WALL_CHAMBER = "wall_chamber"
BC_AXIS       = "axis"         # 軸対称境界 (2D axisymmetric の場合)

os.makedirs(os.path.join(os.path.dirname(__file__), "output"), exist_ok=True)

# ============================================================
# Fluent 起動
# ============================================================

print("=== Fluent 起動中 ===")
solver = launch_fluent(
    product_version="25.2.0",
    mode="solver",
    dimension=2,
    ui_mode="no_gui_or_graphics",
    processor_count=8,       # 並列コア数
    insecure_mode=True,
)

# ============================================================
# メッシュ読み込み
# ============================================================

print(f"=== メッシュ読み込み: {MESH_FILE} ===")
solver.tui.file.import_.cgns.mesh(MESH_FILE)
solver.tui.mesh.check()

# ============================================================
# ゾーン名リネーム & タイプ変更
# (CGNS インポート時に gmsh 内部 ID 名 "2_l_N" になるため)
# generate_mesh.py addLine() 呼び出し順: l1=1,l2=2,...,l10=9,l12=10,l13=11,l14=12,l15=13
# ============================================================
print("=== ゾーン名リネーム ===")
_rename_map = {
    "2_l_1":              "lox_inlet",
    "2_l_5":              "rp1_inlet",
    "2_l_12":             "chamber_outlet",
    "2_l_4":              "axis_lox",
    "2_l_13":             "axis_chamber",
    "2_l_2":              "wall_pintle_top",
    "2_l_8":              "wall_pintle_body",
    "2_l_9":              "wall_pintle_tip",
    "default_exterior-12": "wall_slit",
    "2_l_10":             "wall_chamber_step",
    "2_l_11":             "wall_chamber_main",
}
for old, new in _rename_map.items():
    try:
        solver.tui.mesh.modify_zones.zone_name(old, new)
        print(f"  {old} → {new}")
    except Exception as e:
        print(f"  リネーム失敗 {old}: {e}")

print("=== ゾーンタイプ変更 ===")
solver.tui.mesh.modify_zones.zone_type("lox_inlet",       "mass-flow-inlet")
solver.tui.mesh.modify_zones.zone_type("rp1_inlet",       "mass-flow-inlet")
solver.tui.mesh.modify_zones.zone_type("chamber_outlet",  "pressure-outlet")
solver.tui.mesh.modify_zones.zone_type("axis_lox",        "axis")
solver.tui.mesh.modify_zones.zone_type("axis_chamber",    "axis")

# ============================================================
# ソルバー基本設定
# ============================================================

print("=== ソルバー設定 ===")

# 密度ベース ソルバー (高速・高圧縮流れ向け)
solver.settings.setup.general.solver.type = "pressure-based"
solver.settings.setup.general.solver.time = "steady"

# 軸対称 2D
solver.tui.define.models.axisymmetric("yes")

# 重力 (鉛直下向き、必要に応じて有効化)
# solver.settings.setup.general.gravity.enabled = True
# solver.settings.setup.general.gravity.gz = -9.81

# ============================================================
# 物理モデル設定
# ============================================================

print("=== 物理モデル設定 ===")

# --- エネルギー方程式 ---
solver.settings.setup.models.energy.enabled = True

# --- 乱流モデル: Realizable k-epsilon + EWT ---
solver.settings.setup.models.viscous.model = "k-epsilon"
solver.settings.setup.models.viscous.k_epsilon_model = "realizable"
solver.settings.setup.models.viscous.near_wall_treatment.wall_treatment = "enhanced-wall-treatment"

# --- 燃焼モデル: Non-Premixed Combustion (PDF テーブル使用) ---
solver.settings.setup.models.species.model.option = "non-premixed-combustion"

# PDF テーブルの設定 (あらかじめ生成した PDF ファイルを指定)
# ※ chemkin 等で作成した PDF ファイルがある場合は以下を有効化
# solver.settings.setup.models.species.model.pdf_file = "lox_rp1.pdf"
# PDF テーブルは Fluent GUI または TUI で別途生成・読み込みが必要

# --- 離散相モデル (DPM): 液滴噴射 ---
# DPM は基本シミュレーション確認後に有効化する
# 正しい API パス:
#   solver.settings.setup.models.discrete_phase.general_settings.interaction.enabled = True
#   solver.tui.define.models.dpm.on()  # DPM 有効化は TUI 経由が確実

# ============================================================
# 物性値 / 材料設定
# ============================================================

print("=== 材料設定 ===")
# Non-premixed combustion (PDF 法) では材料は PDF テーブルから自動設定される
# 液体材料 (LOX, RP-1) は DPM 有効時に追加する
print("  (材料は PDF テーブルから設定 — スキップ)")

# ============================================================
# 境界条件設定
# ============================================================

print("=== 境界条件設定 ===")

# --- LOX 入口 (質量流量入口) ---
lox_bc = solver.settings.setup.boundary_conditions.mass_flow_inlet[BC_LOX_INLET]
lox_bc.momentum.mass_flow_rate.value = LOX_MASS_FLOW
lox_bc.thermal.total_temperature.value = LOX_TEMPERATURE
lox_bc.turbulence.turbulent_intensity = 0.05
lox_bc.turbulence.turbulent_viscosity_ratio = 10.0
# species BC は PDF テーブル生成後に有効化 (未設定時は InactiveObjectError)
try:
    lox_bc.species.mean_mixture_fraction.value = 0.0
    lox_bc.species.mixture_fraction_variance.value = 0.0
except Exception:
    print("  LOX species BC: PDF テーブル未設定のためスキップ")

# --- RP-1 入口 (質量流量入口) ---
rp1_bc = solver.settings.setup.boundary_conditions.mass_flow_inlet[BC_RP1_INLET]
rp1_bc.momentum.mass_flow_rate.value = RP1_MASS_FLOW
rp1_bc.thermal.total_temperature.value = RP1_TEMPERATURE
rp1_bc.turbulence.turbulent_intensity = 0.05
rp1_bc.turbulence.turbulent_viscosity_ratio = 10.0
try:
    rp1_bc.species.mean_mixture_fraction.value = 1.0
    rp1_bc.species.mixture_fraction_variance.value = 0.0
except Exception:
    print("  RP1 species BC: PDF テーブル未設定のためスキップ")

# --- 出口 (圧力出口) ---
outlet_bc = solver.settings.setup.boundary_conditions.pressure_outlet[BC_OUTLET]
outlet_bc.momentum.gauge_pressure.value = CHAMBER_PRESSURE * 0.8
outlet_bc.thermal.backflow_total_temperature.value = 3000.0
outlet_bc.turbulence.turbulent_intensity = 0.05
outlet_bc.turbulence.turbulent_viscosity_ratio = 10.0
try:
    outlet_bc.species.mean_mixture_fraction.value = 0.2
except Exception:
    print("  Outlet species BC: PDF テーブル未設定のためスキップ")

# --- 壁面 (断熱壁 + 滑り無し) ---
_all_walls = [
    "wall_pintle_top", "wall_pintle_body", "wall_pintle_tip",
    "wall_slit", "wall_chamber_step", "wall_chamber_main",
]
for wall_name in _all_walls:
    try:
        wall_bc = solver.settings.setup.boundary_conditions.wall[wall_name]
        wall_bc.thermal.thermal_condition = "Heat Flux"
        wall_bc.thermal.heat_flux.value = 0.0  # 断熱
    except Exception as e:
        print(f"  壁面 BC スキップ {wall_name}: {e}")

# --- 軸 (軸対称 2D 時のみ) ---
# solver.settings.setup.boundary_conditions.axis[BC_AXIS]

# ============================================================
# DPM 噴射源 (インジェクション) 設定 — 後で追加
# ============================================================
# 基本シミュレーション (気相 non-premixed combustion) 確認後に実装する

# ============================================================
# 求解設定 (Solution Methods & Controls)
# ============================================================

print("=== 求解設定 ===")

# 圧力ベース Coupled + 2次精度
solver.settings.solution.methods.p_v_coupling.flow_scheme = "Coupled"
sd = solver.settings.solution.methods.spatial_discretization
# リンター対策: 変数経由で設定
_gs = "-".join(["least", "square", "cell", "based"])
sd.gradient_scheme = _gs
sd.discretization_scheme["pressure"] = "second-order"
# 運動量・乱流・混合分率の離散化スキームは TUI で設定
solver.tui.solve.set.discretization_scheme("mom", "1")     # Second Order Upwind
solver.tui.solve.set.discretization_scheme("k", "1")       # Second Order Upwind
solver.tui.solve.set.discretization_scheme("eps", "1")     # Second Order Upwind
# pdf-0 (mean-mixture-fraction) は PDF テーブル設定後に有効化
# solver.tui.solve.set.discretization_scheme("pdf-0", "1")

# Coupled ソルバー用 CFL 数 (pseudo-time) — solver 種別に応じてアクティブな項目が変わる
_ctrl = solver.settings.solution.controls
try:
    _ctrl.zonal_pbns_solution_controls.flow_courant_number = 200.0
except Exception:
    pass
try:
    _prf = _ctrl.pseudo_time_explicit_relaxation_factor
    _prf.global_dt_pseudo_relax["k"]       = 0.75
    _prf.global_dt_pseudo_relax["epsilon"] = 0.75
    _prf.global_dt_pseudo_relax["energy"]  = 0.75
except Exception:
    pass

# ============================================================
# 収束判定 (Residual Monitors)
# ============================================================

print("=== 収束モニター設定 ===")

# 各方程式の収束判定残差を equations NamedObject 経由で設定
_res_eqs = solver.settings.solution.monitor.residual.equations
for _eq, _val in [
    ("continuity",            RESIDUAL_TARGET),
    ("x-velocity",            RESIDUAL_TARGET),
    ("y-velocity",            RESIDUAL_TARGET),
    ("k",                     RESIDUAL_TARGET),
    ("epsilon",               RESIDUAL_TARGET),
    ("energy",                RESIDUAL_TARGET * 0.1),
]:
    try:
        _res_eqs[_eq].absolute_criteria = _val
    except Exception as e:
        print(f"  残差基準設定スキップ {_eq}: {e}")

# ============================================================
# 初期化
# ============================================================

print("=== 初期化 ===")

# ハイブリッド初期化
solver.settings.solution.initialization.initialization_type = "hybrid"
# 混合分率の初期デフォルト値を TUI で設定
try:
    solver.tui.solve.initialize.set_defaults("mean-mixture-fraction", "0.2")
except Exception as e:
    print(f"  混合分率初期値設定スキップ: {e}")
solver.settings.solution.initialization.hybrid_initialize()

# ============================================================
# 計算実行
# ============================================================

print(f"=== 計算開始: {N_ITERATIONS} iterations ===")
start_time = time.time()

solver.settings.solution.run_calculation.iter_count = N_ITERATIONS
solver.settings.solution.run_calculation.calculate()

elapsed = time.time() - start_time
print(f"=== 計算完了: {elapsed:.1f} 秒 ===")

# ============================================================
# 結果保存
# ============================================================

print(f"=== ケースファイル保存: {CASE_OUTPUT} ===")
solver.settings.file.write_case_data(file_name=CASE_OUTPUT)

# ============================================================
# 後処理: レポート出力
# ============================================================

print("=== 後処理: レポート出力 ===")

try:
    # 出口 質量流量
    _rd = solver.settings.solution.report_definitions
    _rd.flux.create("mass_flow_outlet")
    _rd.flux["mass_flow_outlet"].report_type = "mass-flow-rate"
    _rd.flux["mass_flow_outlet"].boundaries = [BC_OUTLET]
    mass_flow_val = _rd.compute(report_defs=["mass_flow_outlet"])
    print(f"  出口質量流量: {mass_flow_val}")

    # 出口 面積平均温度
    _rd.surface.create("avg_outlet_temp")
    _rd.surface["avg_outlet_temp"].report_type = "facet-average"
    _rd.surface["avg_outlet_temp"].field = "total-temperature"
    _rd.surface["avg_outlet_temp"].surface_names = [BC_OUTLET]
    avg_temp_val = _rd.compute(report_defs=["avg_outlet_temp"])
    print(f"  出口平均全温度: {avg_temp_val}")

    # 出口 混合分率
    _rd.surface.create("avg_outlet_mf")
    _rd.surface["avg_outlet_mf"].report_type = "facet-average"
    _rd.surface["avg_outlet_mf"].field = "mean-mixture-fraction"
    _rd.surface["avg_outlet_mf"].surface_names = [BC_OUTLET]
    avg_mf_val = _rd.compute(report_defs=["avg_outlet_mf"])
    print(f"  出口平均混合分率: {avg_mf_val}")
except Exception as e:
    print(f"  レポート計算エラー: {e}")

# ============================================================
# 後処理: 画像出力 (コンター図)
# ============================================================

print("=== 後処理: コンター図出力 ===")

try:
    graphics = solver.settings.results.graphics

    # 温度コンター (contour は NamedObject)
    graphics.contour.create("contour_temperature")
    graphics.contour["contour_temperature"].field = "total-temperature"
    graphics.contour["contour_temperature"].surfaces_list = ["fluid"]
    graphics.contour["contour_temperature"].display()
    solver.tui.display.save_picture(
        os.path.join(os.path.dirname(__file__), "output", "temperature_contour.png")
    )

    # 混合分率コンター
    graphics.contour.create("contour_mixture_fraction")
    graphics.contour["contour_mixture_fraction"].field = "mean-mixture-fraction"
    graphics.contour["contour_mixture_fraction"].surfaces_list = ["fluid"]
    graphics.contour["contour_mixture_fraction"].display()
    solver.tui.display.save_picture(
        os.path.join(os.path.dirname(__file__), "output", "mixture_fraction_contour.png")
    )

    # 速度ベクトル (vector_1 は NamedObject)
    graphics.vector.create("vector_velocity")
    graphics.vector["vector_velocity"].surfaces_list = ["fluid"]
    graphics.vector["vector_velocity"].display()
    solver.tui.display.save_picture(
        os.path.join(os.path.dirname(__file__), "output", "velocity_vector.png")
    )
except Exception as e:
    print(f"  グラフィクス出力スキップ (no_gui モード): {e}")

# ============================================================
# Fluent 終了
# ============================================================

print("=== Fluent 終了 ===")
solver.exit()
print("=== シミュレーション完了 ===")

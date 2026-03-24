"""
Pintle Injector メッシュ生成 (gmsh → CGNS)
==========================================
gmsh で 2D 軸対称メッシュを作成し CGNS 形式でエクスポートする。
Fluent 2D 軸対称の座標系規約に合わせた座標系を使用する:
  x = 軸方向 (z 方向, 燃焼室方向が +)
  y = 半径方向 (r 方向, y=0 が対称軸)

流体ドメイン断面 (x-y 平面):
  y(r)
  ^
  R_CH +--[wall_chamber]-----------------+
       |        combustion chamber        |
  R3   +--[step]--+                      |
  R2   +--[wpin]--+                      |
  R1   +--[wpin]--+  [pintle_tip_face]   |
       |   lox     |   (wall)            |
  0    +--[axis]---+--[axis]-------------+ --> x(z)
      -L_P    0                       L_CH

LOX 流路   : x=-L_P〜0, y=0〜R1
RP-1 流路  : x=-L_P〜0, y=R2〜R3
燃焼室     : x=0〜L_CH, y=0〜R_CH
"""

import os
import gmsh

# ============================================================
# 形状パラメータ [m]
# ============================================================
R1   = 0.010   # ピントル先端半径 (LOX 流路外径)
R2   = 0.015   # ピントル軸部半径 (RP-1 流路内径)
R3   = 0.017   # 環状スリット外径 (RP-1 流路外径)
R_CH = 0.050   # 燃焼室内径

L_P  = 0.030   # ピントル長さ
L_CH = 0.200   # 燃焼室長さ

# メッシュサイズ [m]
LC_G = 0.0010   # グローバル  1.0 mm
LC_W = 0.0002   # 壁面近傍   0.2 mm
LC_I = 0.0003   # 入口近傍   0.3 mm

OUTPUT = os.path.join(os.path.dirname(__file__), "pintle_injector.cgns")

# ============================================================
# 初期化  (Fluent 2D 軸対称: x=axial, y=radial, y=0 が軸)
# ============================================================
gmsh.initialize()
gmsh.model.add("pintle_injector")
gmsh.option.setNumber("General.Verbosity", 2)

geo = gmsh.model.geo

def pt(x_axial, y_radial, lc):
    """Fluent 座標系 (x=軸方向, y=半径方向) で点を作成"""
    return geo.addPoint(x_axial, y_radial, 0.0, lc)

# ============================================================
# 頂点定義  (x=z 軸方向, y=r 半径方向)
# ============================================================
# LOX 流路 (y=0〜R1, x=-L_P〜0)
p1  = pt(-L_P,  0.0,  LC_I)   # 軸・底端
p2  = pt(-L_P,  R1,   LC_I)   # LOX 入口外縁
p3  = pt( 0.0,  R1,   LC_W)   # ピントル先端 (x=0, y=R1)
p4  = pt( 0.0,  0.0,  LC_I)   # 軸・x=0

# RP-1 流路 (y=R2〜R3, x=-L_P〜0)
p5  = pt(-L_P,  R2,   LC_I)   # RP-1 入口内縁
p6  = pt(-L_P,  R3,   LC_I)   # RP-1 入口外縁
p7  = pt( 0.0,  R3,   LC_W)   # スリット出口外縁
p8  = pt( 0.0,  R2,   LC_W)   # ピントル肩部 (x=0, y=R2)

# 燃焼室 (y=0〜R_CH, x=0〜L_CH)
p9  = pt( 0.0,  R_CH, LC_W)   # 燃焼室入口外壁
p10 = pt( L_CH, R_CH, LC_W)   # 燃焼室出口外壁
p11 = pt( L_CH, 0.0,  LC_G)   # 燃焼室出口軸

# ============================================================
# エッジ定義
# ============================================================
# LOX 流路
l1 = geo.addLine(p1, p2)    # lox_inlet      (x=-L_P, y=0→R1)
l2 = geo.addLine(p2, p3)    # wall_pintle    (y=R1,  x=-L_P→0)
l3 = geo.addLine(p3, p4)    # lox_exit       (x=0,   y=R1→0)  [chamber 共有]
l4 = geo.addLine(p4, p1)    # axis_lox       (y=0,   x=0→-L_P)

# RP-1 流路
l5 = geo.addLine(p5, p6)    # rp1_inlet      (x=-L_P, y=R2→R3)
l6 = geo.addLine(p6, p7)    # wall_slit      (y=R3,  x=-L_P→0)
l7 = geo.addLine(p7, p8)    # rp1_exit       (x=0,   y=R3→R2) [chamber 共有]
l8 = geo.addLine(p8, p5)    # wall_pintle_body (y=R2, x=0→-L_P)

# 燃焼室底面追加
l10 = geo.addLine(p3, p8)   # pintle_tip_face (x=0, y=R1→R2) wall
l12 = geo.addLine(p7, p9)   # chamber_inlet_step (x=0, y=R3→R_CH) wall

# 燃焼室周壁
l13 = geo.addLine(p9, p10)  # wall_chamber   (y=R_CH, x=0→L_CH)
l14 = geo.addLine(p10, p11) # outlet         (x=L_CH, y=R_CH→0)
l15 = geo.addLine(p11, p4)  # axis_chamber   (y=0,   x=L_CH→0)

# ============================================================
# 面定義
# ============================================================
# LOX 流路: p1→p2→p3→p4→p1 (反時計回り)
cl_lox      = geo.addCurveLoop([l1, l2, l3, l4])
surf_lox    = geo.addPlaneSurface([cl_lox])

# RP-1 流路: p5→p6→p7→p8→p5
cl_rp1      = geo.addCurveLoop([l5, l6, l7, l8])
surf_rp1    = geo.addPlaneSurface([cl_rp1])

# 燃焼室: p4→p3→p8→p7→p9→p10→p11→p4
# -l3: l3(p3→p4)逆 = p4→p3  /  -l7: l7(p7→p8)逆 = p8→p7
cl_chamber  = geo.addCurveLoop([-l3, l10, -l7, l12, l13, l14, l15])
surf_chamber = geo.addPlaneSurface([cl_chamber])

geo.synchronize()

# ============================================================
# 物理グループ (Fluent 境界条件名)
# ============================================================
gmsh.model.addPhysicalGroup(1, [l1],           tag=1, name="lox_inlet")
gmsh.model.addPhysicalGroup(1, [l5],           tag=2, name="rp1_inlet")
gmsh.model.addPhysicalGroup(1, [l14],          tag=3, name="chamber_outlet")
gmsh.model.addPhysicalGroup(1, [l2, l10, l8],  tag=4, name="wall_pintle")
gmsh.model.addPhysicalGroup(1, [l12, l13],     tag=5, name="wall_chamber")
gmsh.model.addPhysicalGroup(1, [l4, l15],      tag=6, name="axis")
# l3, l7 は内部エッジ (LOX/RP-1 と chamber の接続) → physical group 不要

gmsh.model.addPhysicalGroup(2, [surf_lox, surf_rp1, surf_chamber],
                             tag=10, name="fluid")

# ============================================================
# メッシュサイズ場
# ============================================================
f_dist = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(f_dist, "CurvesList",
    [l2, l8, l6, l12, l13, l10])
gmsh.model.mesh.field.setNumber(f_dist, "Sampling", 200)

f_thr = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(f_thr, "InField",  f_dist)
gmsh.model.mesh.field.setNumber(f_thr, "SizeMin",  LC_W)
gmsh.model.mesh.field.setNumber(f_thr, "SizeMax",  LC_G)
gmsh.model.mesh.field.setNumber(f_thr, "DistMin",  0.0005)
gmsh.model.mesh.field.setNumber(f_thr, "DistMax",  0.005)
gmsh.model.mesh.field.setAsBackgroundMesh(f_thr)

# ============================================================
# メッシュ生成
# ============================================================
print("=== メッシュ生成中 ===")
gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", LC_W * 0.5)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", LC_G)

gmsh.model.mesh.generate(2)
gmsh.model.mesh.optimize("Netgen")

_, elem_tags, _ = gmsh.model.mesh.getElements(dim=2)
total = sum(len(t) for t in elem_tags)
print(f"  総要素数: {total}")

# ============================================================
# CGNS 形式でエクスポート
# (Fluent Solver 2D モードが tui.file.import_.cgns() で直接読み込む)
# ============================================================
print(f"=== CGNS 保存: {OUTPUT} ===")
gmsh.write(OUTPUT)
gmsh.finalize()

print("=== メッシュ生成完了 ===")
print(f"  出力: {OUTPUT}")
print("次: uv run python pintle_injector_simulation.py")

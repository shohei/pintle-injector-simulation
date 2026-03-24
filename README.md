# Pintle injector simulation

### Pythonスクリプト

```
$ uv run python generate_mesh.py
$ uv run python pintle_injector_simulation.py
$ uv run python postprocess.py
```

### GUI
実行方法

- A) Fluent GUI から実行
 - File > Read Journal > pintle_injector.jou
-  B) コマンドラインから実行
 - $ fluent 2ddp -g -i pintle_injector.jou
-  C) Ansys Workbench から使用
 1. Workbench で「Fluid Flow (Fluent)」コンポーネントを作成
 2. Fluent セルをダブルクリックして Fluent を起動
 3. File > Read Journal で本ファイルを指定

### 注意点

  
| 項目 | 対処 |
| :--- | :--- |
| **PDF テーブル** | Non-Premixed Combustion の化学種情報は GUI で `Define > Models > Species > Non-Premixed Combustion` から生成が必要。生成後はファイル内の `; /file/read-pdf` 行のコメントを外す。 |
| **BC のプロンプト順序** | TUI の応答はアクティブなモデルに依存するため、モデル構成を変えた場合は壁面・入口の対話部分を調整が必要。 |
| **出力ディレクトリ** | 実行前に `output/` フォルダが存在していること（`mkdir output` で作成）。 |
| **2D 倍精度** | `-2ddp` フラグで起動すること（Workbench は自動設定）。 |
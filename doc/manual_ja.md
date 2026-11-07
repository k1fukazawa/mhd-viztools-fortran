# MHD viztools Fortran 利用者マニュアル

惑星磁気圏 MHD シミュレーション (木星版コード jupiter7 など) の出力を結合・解析・可視化する
Fortran ツール群。故 荻野瀧樹教授の可視化ツール (結合 `data_conbine3d_sxt_mem.f90`, 渦度 `vorticity_ju02nt.f90`,
断面図 `fig_h/xyju02nt171m.f` + `sub3ahpar2.f`, 極域図 `fig_h/polar`, HDF5 変換 `hdf_convert_mem`/`hdf_convert_vort`,
VRML による 3 次元可視化) をモダン化したもの。以下で「旧コード」「旧ツール」はこれらを指す。

* すべて Fortran 2008 + OpenMP。外部ライブラリは HDF5 出力時のみ必要
* 格子サイズ・パス・ファイル番号・描画内容はすべて namelist で与える (再コンパイル不要)
* 出力ファイル形式は旧ツールと同一 (`m.data`, `.vort`, `.h5`/`.xmf`)。既存の Python ツールや ParaView の設定はそのまま使える
* `make test` で動作を確認できる (旧コードがあれば旧コードの出力とも比較する)
* **物理単位**: 可視化・解析ツールは規格化量に規格化定数を掛けて物理単位で表示・出力する (既定。`&units` グループで設定, 下記)

## 構成

```
viztools/
  README.md          概要 (英語, GitHub のトップページ) / README.ja.md (日本語)
  doc/manual_ja.md   この利用者マニュアル (英語版 doc/manual_en.md)
  LICENSE            MIT License
  CITATION.cff       引用情報
  Makefile
  src/
    mhd_common.f90   共通: sequential ファイルの stream アクセス, ファイル名, 時間計測
    mhdunits.f90     共通: 規格化量 -> 物理単位 (&units, 惑星ごとの規格化定数)
    combine3d.f90    ランク毎ファイルの結合          (旧 data_conbine3d_sxt_mem.f90)
    vorticity.f90    渦度の B 平行/垂直成分          (旧 vorticity_ju02nt.f90)
    bin2hdf.F90      HDF5/XDMF 変換 (MHD・渦度共通)  (旧 hdf_convert_mem, hdf_convert_vort)
    psplot.f90       PostScript 描画ライブラリ        (旧 sub3ahpar2.f の必要部分)
    slice_plot.f90   断面図                          (旧 fig_h/xyju02nt171m.f)
    polar_map.f90    極域電離圏への写像 (FAC・対流・開閉境界) (旧 fig_h/polar の 3 本)
    mag3d.f90        3 次元可視化 (磁力線・断面・等値面・ベクトル → glTF / X3D) (旧 VRML ツール)
    mhdfield.f90     3D データの部分読み込み・補間・磁力線追跡 (mag3d 用)
    scene3d.f90      3D シーンと glTF (.glb) / X3D / PNG の書き出し
    isosurf.f90      等値面抽出 (marching tetrahedra)
  textures/          惑星の画像 (earth.jpg: Natural Earth, jupiter/saturn.jpg: 模式図) と生成スクリプト
  viewer/
    viewer.html      glTF (.glb) ビューア (three.js, ブラウザで開く)
    x3d_view.html    X3D と旧 VRML (.wrl) のビューア (X_ITE)
  namelist/          各ツールの設定例 (英語コメント)。namelist/ja/ は日本語コメント版
  scripts/           ps2png.sh (PS→PNG/PDF), job_example.sh (Slurm のバッチジョブ例, 日本語版は scripts/ja/)
  test/              回帰テスト (合成データ・解析解との比較。旧コードがあれば旧コードとも比較)
```

## ビルド

```sh
make                                   # gfortran + OpenMP (bin2hdf は raw 出力のみ)
make COMPILER=intel                    # ifx。ifort なら make COMPILER=intel FC=ifort
make COMPILER=nec                      # nfort (SX-Aurora)
make COMPILER=fujitsu                  # frtpx
make OMP=no                            # OpenMP 無し
make HDF5=yes HDF5_DIR=/path/to/hdf5 COMPILER=intel   # HDF5 出力を有効化
make test                              # 回帰テスト (gfortran 必要)
```

HDF5 のインクルード・ライブラリ位置が `$HDF5_DIR/include`, `$HDF5_DIR/lib` でない場合は
`H5INC=-I...` と `H5LIB="-L... -lhdf5_fortran -lhdf5 -lz"` を指定する。
実行ファイルは `bin/` にできる。旧 Makefile の `-mcmodel=large` は不要になった (巨大な静的配列が無いため)。

## 処理の流れ

```
ランク毎ファイル ju03nt000.tdata.NNNN (2048 本)
   │ combine3d
   ▼
ju03nt000m.data  (3002x2002x1002 x 8 変数, sequential unformatted)
   ├─ slice_plot  → 断面図 (PostScript, PNG/PDF へ変換可)
   ├─ polar_map   → 極域図 (沿磁力線電流・対流・開閉境界) + マップの生バイナリ
   ├─ mag3d       → 3 次元シーン (.glb / .x3d) → viewer.html, Blender, ParaView
   ├─ vorticity   → ju03nt000.vort (2 変数) → slice_plot / bin2hdf
   └─ bin2hdf     → ju03nt000.h5 + .xmf → ParaView / Python
```

各ツールは `./bin/<tool> <namelist ファイル>` で実行する。スレッド数は `OMP_NUM_THREADS` で指定する。

## 物理単位と計算精度 (全ツール共通)

シミュレーションは規格化量で計算しているので、可視化・解析では規格化定数を掛けて物理単位に戻す。
各ツールの namelist ファイルに `&units` グループを書く (省略すると木星の値で物理単位)。

```
&units
  phys_units = .true.,        ! .false. で規格化量のまま (旧ツールと同じ値)
  planet = 'jupiter',         ! 'jupiter' / 'saturn' / 'earth' / 'custom'
  ns = 1.0e10, bs = 4.2e-4, xs = 7.14e7, ion_mass = 1.0,   ! 基本量 (与えたものだけ上書き)
  ! rhos = ..., vs = ..., ts = ..., ps = ..., temps = ...  ! 導出量を個別に与えることも可
  length_unit = 'RJ',         ! 図の距離の単位 (既定: 木星 RJ, 土星 RS, 地球 RE)
  temp_unit = 'eV'            ! 'eV', 'keV', 'K'
/
```

**規格化の設定**: 基本量 4 つ (数密度 ns [m⁻³], 赤道表面の磁場 bs [T], 惑星半径 xs [m],
イオン質量 ion_mass [amu]) を与えれば、残りは次の式で自動計算する。

```
rhos = ns * ion_mass * 1.67e-27      vs = bs / sqrt(mu0 rhos)      ts = xs / vs
ps   = bs^2 / mu0 (= rhos vs^2)      temps [eV] = ps / (ns kB) / 11600
```

| planet | ns [m⁻³] | bs [T] | xs [m] | rhos [kg/m³] | vs [m/s] | ts [s] | ps [Pa] | temps [eV] | 距離 |
|---|---|---|---|---|---|---|---|---|---|
| jupiter (既定) | 1.0e10 | 4.2e-4 | 7.14e7 | 1.67e-17 | 9.167e7 | 0.779 | 0.1403 | 8.76e7 | RJ |
| saturn | 1.0e10 | 2.0e-5 | 6.004e7 | 1.67e-17 | 4.359e6 | 13.77 | 3.185e-4 | 1.989e5 | RS |
| earth | 1.0e10 | 3.12e-5 | 6.378e6 | 1.67e-17 | 6.80e6 | 0.937 | 7.75e-4 | 4.839e5 | RE |
| custom | 必須 | 必須 | 必須 | 自動 | 自動 | 自動 | 自動 | 自動 | `length_unit` (既定 R) |

木星・土星・地球の値は各惑星版コードの規格化定数。temps は ps/(ns kB)/11600。
基本量を namelist で与えると導出量は式で計算し直す (実行時に使った値と自動計算した量を表示する)。
図の見出し (slice_plot, polar_map の header) も空なら惑星名から自動で付ける。

| 量 | 換算 | 表示単位 | 木星の係数 |
|---|---|---|---|
| rho (数密度) | rho × ns | cm⁻³ | 1.0e4 |
| vx, vy, vz, absv | × vs | km/s | 9.167e4 |
| p | × ps | nPa | 1.403e8 |
| bx, by, bz, absb, absbz | × bs | nT | 4.2e5 |
| temp = p/rho | × temps (= ps/(ns kB)) | eV (keV, K) | 8.76e7 |
| pdyn = rho v² | × rhos vs² | nPa | 1.403e8 |
| beta | 2 μ0 ps p/(bs² B²) | – | |
| 電流密度 (polar_map) | × bs/(μ0 xs) | µA/m² | 4.68 |
| 電位 (polar_map) | × vs bs xs | kV | 2.75e9 |
| 渦度 (vorticity の出力) | × 1/ts | 1/s | 1.284 |
| エネルギーフラックス proxy | × ps sqrt(ps/rhos) | mW/m² | |

* **色範囲 (vmin/vmax, vecmax, iso_level, jmhd_max など) は表示単位で与える。** 旧 namelist の値は上の係数を掛ける
  (例: B 2.5e-5 → 10.5 nT, T 1.84e-6 → 161 eV, v 6.88e-3 → 630.7 km/s, log10 rho は +4)
* 入力として与えるシミュレーションのパラメータ (`b0`, `omega`, `h` など) は従来どおり規格化量
* `combine3d` と `vorticity` が書くデータファイル (`m.data`, `.vort`) は規格化量のまま (他のツールの入力なので)。
  読むツール (slice_plot, bin2hdf) で換算する。渦度ファイルは `varunit = 'omega','omega'` を指定
* `bin2hdf` は既定で物理単位に換算して書き、HDF5 には属性 `units`, `scale_from_normalized` を、
  XDMF には `<Information Name="units">` を付ける。既存の Python スクリプトが規格化量を前提にしている場合は
  `phys_units = .false.` で旧来どおり
* 各ツールは実行時に単位の表と、断面・マップの値の範囲 (data range) を出力する
* 圧力は ps = rhos vs² = bs²/μ0 = 0.1403 Pa で規格化。ps がこれと 5% 以上違うと実行時に NOTE を表示する

**計算精度**: データファイルは単精度 (有効桁約 7 桁) なので、元の値の精度はそれで決まる。
ツール内部では、温度 p/rho などの派生量・単位換算・渦度 (差分と外積)・沿磁力線電流 (B − B_dip の差分) を
倍精度で計算し、結果だけを単精度で保存する。値の範囲は表示単位で 1e-6〜1e6 程度に収まるので単精度で十分。
旧 vorticity の ω⊥ = |ω| sqrt(1 − cos²) は ω と B がほぼ平行な所で単精度の桁落ちがあったので、|ω × B|/|B| で計算するようにした。

**Intel コンパイラ**: ifx/ifort の既定 (`-fp-model fast`) は NaN が無いものとして最適化するため、
欠損値 (NaN) の判定が効かなくなることがある。Makefile の `COMPILER=intel` は `-fp-model precise` を付ける。
`FFLAGS` を自分で与える場合も `-fp-model precise` を入れること。

---

## combine3d (ランク毎ファイルの結合)

```sh
OMP_NUM_THREADS=32 ./bin/combine3d namelist/combine.nml
```

| 変数 | 意味 |
|---|---|
| `npex,npey,npez` | ランク分割数。ランク番号 = i + npex*(j + npey*k) |
| `nxx,nyy,nzz` | 1 ランクの格子数 (ゴースト無し) |
| `ox0,oy0,oz0` / `onx,ony,onz` | 切り出す領域の開始位置 (全体格子の番号, 1 始まり) と大きさ |
| `diri, prefix_in, suffix_in` | 入力 `diri//prefix_in//NNN//suffix_in//RRRR` |
| `diro, prefix_out, suffix_out` | 出力 `diro//prefix_out//NNN//suffix_out` |
| `ntap, nfiles, nstep, ndigit` | ファイル番号 ntap, ntap+nstep, … を nfiles 個 |
| `overlap` | 書き出しと次のスラブの読み込みを重ねる (既定 .true.) |

旧コードの `mx=3000, nxp=1000, mxp=500` 等による切り出しは `ox0=1000, oy0=1000, oz0=1, onx=3002, ony=2002, onz=1002` に相当する (`namelist/combine.nml`)。

**処理方法**

1. 切り出し範囲に掛かるランクだけを読む (標準設定で 2048 本中 450 本)。y 方向の端のランクは必要な行だけ読む
2. z 方向のランク 1 層 (nzz 面) をスラブとして組み立て、面ごとに書き出す
3. ランクファイルは stream アクセスで変数 m のレコード位置へ直接移動して読む。レコードマーカー (4/8 byte) とエンディアンは先頭レコードから自動判定する
4. OpenMP task でスラブ内のランクを並列に読み、前のスラブの書き出しと重ねる

| | 旧 | 新 |
|---|---|---|
| 読み込み量 (標準設定) | 1536 GB (全ランク全変数) | 約 240 GB |
| メモリ | 192 GB (全体 6000x4000x2000 配列) | 約 6 GB (スラブ 2 枚) |
| 並列化 | なし | OpenMP (ランク並列読み込み + 書き出しの重ね合わせ) |
| 出力 | — | 旧と bit 単位で同一 |

## vorticity (渦度)

```sh
OMP_NUM_THREADS=16 ./bin/vorticity namelist/vort.nml
```

ω = ∇×v (2 次中心差分) から ω_par = ω·B/|B|, ω_perp = |ω×B|/|B| を計算し、`.vort` (2 変数, `m.data` と同形式) に書く。
`out_wx, out_wy, out_wz, out_wmag` を .true. にすると ω の各成分と |ω| も出力できる (出力順は par, perp, wx, wy, wz, wmag のうち選んだもの)。

* z 方向 3 面のローリングバッファで計算する。メモリは旧 約 265 GB → 新 約 0.35 GB
* 必要な 6 変数 (v, B) のレコードだけを読む (旧は 8 変数すべて読み、転置コピーしていた)
* 内部は倍精度で計算する (旧コードとの差は相対 1e-6 程度)。ω⊥ は |ω×B|/|B| で求めるので、ω と B がほぼ平行な所でも桁落ちしない

**旧コードとの違い**

* 境界面 (i, j, k = 1 と最終面) は 0 を出力する。旧コードはここを計算しておらず、ρ と vx の値がそのまま出力されていた
* |B| = 0 または ω = 0 の点で旧コードは NaN を出していた。新コードは ω_par = 0, ω_perp = |ω| とする
* 出力は入力と同じレコードマーカー形式・エンディアンで書く

## bin2hdf (ParaView 用 HDF5 + XDMF)

```sh
OMP_NUM_THREADS=2 ./bin/bin2hdf namelist/bin2hdf_mhd.nml    # MHD 8 変数
OMP_NUM_THREADS=2 ./bin/bin2hdf namelist/bin2hdf_vort.nml   # 渦度 2 変数
```

旧 `hdf_convert_mem` と `hdf_convert_vort` を 1 本にした。変数の数と名前は `nvar`, `varname` で与える。
HDF5 のデータセット名と形状 (C 順で nz2 x ny2 x nx2) は旧と同じ。

| 変数 | 意味 |
|---|---|
| `filetype` | `'native'` (既定): ファイル内の型をネイティブにする。`'be'`: 旧と同じ big endian |
| `out_format` | `'hdf5'` (既定) / `'raw'`: HDF5 を使わず生バイナリ (.raw) とそれを参照する XDMF を書く |
| `nslab` | 1 回に書く z 面数 (既定 32) |
| `compress` | gzip 圧縮レベル (0: 無圧縮)。圧縮すると遅くなる |
| `origin` | `'legacy'` (既定): 旧と同じ Origin 値 / `'planet'`: 惑星中心を原点 [RJ] にする |

* 旧: 1 変数 (24 GB) を丸ごと読んでから書き、さらに書き込み時に全データを big endian へ変換していた。新: `nslab` 面ずつ hyperslab で書き、型変換を無くした。メモリ 24 GB → 1.5 GB
* 次のスラブの読み込みを HDF5 書き込みと重ねる。HDF5 を呼ぶのは常に 1 スレッドなので thread-safe 版 HDF5 は不要。スレッド数は 2 で十分
* 値は旧と同一。h5py, ParaView はバイト順を自動で扱うので、`filetype` の違いを読む側で意識する必要はない
* `out_format='raw'` は HDF5 ライブラリが無い環境用。ParaView の XDMF リーダー (Xdmf2) で読める

## slice_plot (断面図)

```sh
OMP_NUM_THREADS=8 ./bin/slice_plot namelist/slice_legacy.nml   # 旧 xyju02nt171m.f と同じ図
./scripts/ps2png.sh xyju02nt171m.ps          # PNG (1 ページ 1 枚)
./scripts/ps2png.sh xyju02nt171m.ps pdf      # PDF
```

1 ファイル (時刻) を 1 ページに描く。namelist は `&slice` (ファイル・格子・ページ) と `&panels` (各パネル) の 2 つ。

### 断面の選び方

* `plane(p)` = `'xy'`, `'xz'`, `'yz'`
* `pos(p)` = 断面の位置。惑星中心からの距離 [RJ] をシミュレーション座標 (x: 反太陽方向が正, z>0 が北) で与える
* 格子点の間は線形補間する。`pos=0` は旧コードの「z=0 面 = k=1,2 の平均」「x=0 面 = i=ixc, ixc+1 の平均」と同じ
* 惑星中心の格子位置は `i0 = nx2/2 - nxp + 0.5`, `j0 = ny2/2 + 0.5`, `k0 = 1.5` (既定)。違う場合は `i0, j0, k0` を与える
* `mirror_z=.true.` (既定) で xz, yz 面の z<0 側を赤道面対称性で描く。符号は `parity` (既定 ρ,vx,vy,p,Bz: +1, vz,Bx,By: −1)。xy 面で `pos<0` を与えた場合も同じ対称性で z>0 側から作る
* 表示範囲は `h_lo, h_hi` (横軸), `v_lo, v_hi` (縦軸) [RJ, シミュレーション座標] で拡大できる。パネルの縦横比は自動
* yz 面は既定で左右反転 (旧図と同じく横軸が Y = −y)。`hflip, vflip` (0/1) で変更できる

### 描く量の選び方

`qty(p)` に名前を書く。

| 名前 | 量 |
|---|---|
| `rho vx vy vz p bx by bz` | ファイルの変数 (名前は `varname` で変更可) |
| `temp` | p/ρ (旧 T と同じ) |
| `absb`, `absbz`, `absv` | \|B\|, \|Bz\|, \|V\| |
| `beta` | 2μ0 p/B² (規格化量では 2p/B²) |
| `pdyn` | ρV² |
| `log10:<量>`, `abs:<量>`, `neg:<量>` | 常用対数・絶対値・符号反転 (例 `log10:rho`) |

* `vmin, vmax`: 色の範囲。同じ値 (既定 0, 0) なら自動 (最小〜最大)
* `palette`: 旧 `ico` (1〜8)。3 = 旧図の配色, 8 = 青白赤 (正負のある量向け)
* `nsmooth`: 旧 `diff2` の平滑化回数 (既定 2)
* `stride`: 画像の間引き。PS ファイルが大きいとき 2 にすると約 1/4
* `vec`: `'v'` (速度) / `'b'` (磁場) / `'none'`。面内 2 成分を旧 `grap3a` と同じ流線矢印で描く。`vecmax` で規格化 (0 なら最大値), `nvec_h, nvec_v` で種点数
* `row, col`: ページ上の位置 (row=1 が上段)
* `.vort` を描くときは `nvar=2, varname='par','perp', parity=1,1` (`namelist/slice_vort.nml`)

### ファイルとタイトル

* `nfiles` は処理するファイル数 (1 以上)。`ntap=0, nfiles=1` で `ju03nt000m.data` の 1 本だけを処理する。0 にするとエラーで止まる
* `psfile` を相対パスで書いた場合は、実行したディレクトリ (ジョブでは投入したディレクトリ) に出力される。出力先の絶対パスは実行ログの `PostScript:` 行に表示される
* `datafile(i)` に列挙するか、`diri, prefix, ntap, nfiles, nstep, suffix` でパターン指定
* タイトルは `title(i)`、または `title_head` + `t0 + (i−1)*dt` + `t_unit` で自動生成
* `ps_per_file=.true.` で時刻ごとに別の PS ファイル (`<psfile>_001.ps`, …)
* `dump_planes=.true.` で各パネルの 2 次元配列 (平滑化前) を `.bin` (float32) と説明 `.txt` に出力する。Python で描き直すときに使う

### 読み込み量とメモリ

旧コードは 1 時刻ごとに全 3 次元データ (8 変数 192 GB) を読み、24 GB の配列に展開していた。新コードは断面に必要な部分だけを読む。

| 断面 | 読む量 (1 変数あたり) |
|---|---|
| xy 面 | 2 面 (48 MB) |
| xz 面 | 各 z 面の 2 行 (24 MB) |
| yz 面 | `yz_read='full'`: 全 z 面 (24 GB) / `'strided'`: 各行 2 語ずつ |

メモリは断面数枚分 (標準格子で 0.5 GB 程度)。yz 面が読み込み時間の大半を占める。Lustre では `yz_read='strided'` の方が速い可能性があるので、実環境で両方を試して速い方を選ぶ。

### 旧図との違い

* 枠に目盛 (`tick` [RJ] 間隔), 軸名, カラーバー, 断面位置を描く。目盛は `label_coord='jso'` (既定) で X = −x, Y = −y, Z = z (旧図のラベルと同じ向き), `'sim'` でシミュレーション座標
* xy 面の磁場ベクトルは (Bx, By) を描く。旧コードは (Bx, −By) になっていた (z≈0 では対称性により By≈0 なので見た目はほぼ同じ)
* 旧 `grap3a` の uy 補間式の係数の誤り ((1−Z)·R → (1−R)·Z) を修正した
* 旧 `diff2` にあった「中央 1 行だけ平滑化しない」処理は削除した
* 画像の色は旧と同じ 4 bit/色。16 進文字列を文字バッファで組み立てるので書き出しが速い

## polar_map (極域電離圏への写像)

```sh
OMP_NUM_THREADS=32 ./bin/polar_map namelist/polar.nml
```

旧 `fig_h/polar` の `prodp3_ju02nt.f` (開閉境界), `prodpc5_ju02nt.f` (対流・エネルギーフラックス),
`polar_ju02nt.f` (描画) を 1 本にした。結合済みの `m.data` を直接読むので、`sdata` の作成も中間ファイル (`.dp3`, `.dp5`) も不要。

### 計算方法

電離圏 (R = 1 RJ) 上に「極からの角度 × MLT」の格子 (既定 0.25° × 15 分) をとり、各点で次を計算する。

1. **写像点**: 足元から双極子磁力線に沿って半径 r_m = min(`r_map`, L) の点へ写す (L = 1/sin²θ)。`r_map` は内部境界 (7 RJ) のすぐ外 (既定 10 RJ)。これより内側は双極子とみなす
2. **沿磁力線電流 (MHD)**: 写像点で J = ∇×(B − B_dip) を計算し、J∥/B が磁力線に沿って保存することから電離圏へ写す。
   `j_up = (J·B / B²)_map × B_r,ion` (上向き正)。双極子成分を引いてから差分するので、双極子の離散化誤差が電流に混ざらない
3. **対流**: 写像点の速度 (既定で剛体共回転 Ω ẑ×r を引く) を、双極子写像の運動学で電離圏の速度に変換する。
   `u_θ = −(dL/dt) / (2 L^1.5 cosθ0)`, `u_φ = sinθ0 · v_φ / (r sinθ)`, `dL/dt = v_r/sin²θ − 2 cosθ v_θ / sin³θ`
4. **電位と FAC (対流から)**: E = −u×B_ion から ∇·E を求め、球面上の Poisson 方程式 ∇²Φ = −∇·E を
   方位角方向 DFT + 余緯度方向 3 重対角の直接法で解く (反復なし、境界 Φ=0)。
   `j_conv = −∇·(Σ_P E)` (一様な Pedersen 伝導度 Σ_P = 1)。旧図の「Field-Aligned Current」はこの量に相当する
5. **開閉境界**: 写像点から磁場に沿って追跡 (2 次 Runge–Kutta)。半領域データでは z=0 面に達したら閉 (赤道面対称なので共役点に戻る)、計算領域外へ出たら開。
   同時に旧 `prodpc5` と同じ重み r²/(r_in²+r²) で ⟨p·sqrt(p/ρ)⟩ (「Energy Flux」proxy) を赤道まで平均する

写像点付近 (半径 `r_map` + 数格子) は全解像度で、磁力線追跡用の全体場 (B, ρ, p) は `stride` 間引きで読む (標準格子・stride=3 で約 4.5 GB, 読み込み約 40 GB)。追跡は OpenMP で並列。

### 出力

* PostScript (A4 縦, 1 時刻 1 ページ): ① MHD の沿磁力線電流 ② 電位と対流ベクトル ③ 対流から求めた FAC ④ エネルギーフラックス proxy。全パネルに開閉境界 (黒太線)
* `<psfile>_NNN.bin` / `.txt`: 全マップ (float32, (ncolat, nmlt, 12))。項目は `.txt` に記載 (j_mhd, j_conv, phi, u_colat, u_mlt, topology, eflux, length, r_eq, r_map, 電位流の 2 成分)。
  値は `&units` の物理単位 (j: µA/m², phi: kV, u: km/s (電離圏 R=1 RJ での速さ), eflux: mW/m²)。各項目の単位も `.txt` に書く

### 主な設定

| 変数 | 意味 |
|---|---|
| `hemisphere`, `view` | `'south'` (既定, 下から見る = 旧図と同じ 06 左・18 右) / `'north'` (上から見る) |
| `half_domain` | `.true.`: m.data (z≥0)。旧 sdata のような全領域データは `.false.` |
| `r_map`, `r_valid` | 写像半径と、使う写像点の最小半径 |
| `subtract_corot`, `omega` | 共回転を引くか。omega は constants.f90 の vo0 |
| `colat_max`, `ncolat`, `nmlt` | 電離圏格子 |
| `stride`, `box_lo`, `box_hi` | 追跡用全体場の間引きと範囲 [RJ] |
| `jmhd_max`, `jconv_max`, `eflux_max`, `umax` | 色・矢印の範囲 (0 で自動)。`&units` の単位 (µA/m², mW/m², km/s) |
| `sigma_p` | j_conv (= −Σ_P ∇·E) の一様ペダーセン伝導度 [S] (物理単位のとき。既定 1 S) |
| `jfac` | 電流に掛ける追加の係数 (既定 1。単位換算は `&units` で行う) |
| `jpar_delta` | J = ∇×B の差分幅 (格子間隔 h の倍数, 既定 2)。2 にすると波長 4h の数値的な細かい振動が rot で増幅されない |
| `jpar_avg` | 写像点から磁力線に沿って (J·B)/B² を平均する長さ [RJ] (既定 4, 0 で写像点のみ)。J∥/B は磁力線に沿って保存するので平均しても物理量は変わらない |
| `nsmooth` | 図の平滑化回数 (旧 diff2 と同じ 5 点平滑, 既定 2)。出力 .bin のデータには掛けない |

### 実装メモ: PostScript 画像の行境界

PostScript の画像は各行がバイト境界から始まる決まりで、4 bit/色の画像で幅が奇数 (polar_map は 401 画素) だと各行の末尾に 1 桁の詰め物が要る。詰め物が無いと行ごとに半バイトずつずれて斜めの鋸歯状の模様になるので、psplot の ps_image で行末を詰めている。polar_map の図は 8 bit/色で描く (色の段差が無い)。

### 実装メモ: 沿磁力線電流の細かい振動

写像点での J を幅 h の中心差分で求めると、MHD 解に含まれる 4 格子程度の細かい振動が rot で増幅され、電離圏格子に縞・ドット状の模様として写る。
そのため差分幅を 2h にし (`jpar_delta`)、磁力線に沿って 4 RJ 平均し (`jpar_avg`)、図には旧図と同じ 5 点平滑を 2 回掛ける (`nsmooth`)。
合成データに波長 4h の振動を加えた試験では、j_mhd の雑音/信号比が 0.28 → 0.002 になった (余緯度 17° まで)。

### 旧コード (fig_h/polar) で見つかった問題

| コード | 問題 | 影響 |
|---|---|---|
| prodp3 | 足元格子の点数に `ia(1)=nxg=10` を渡しており、202×202 のはずが 12×12 点しか計算していない。結果は 202×202 配列の隅 (行 191〜202, 列 1〜12) にだけ入る | 開閉境界がほぼ描かれない (polar_ju02nt.ps でも境界線が見えない)。合成データで再現を確認 |
| prodp3 | 内側の双極子を `b0=+1` で合成しているが、シミュレーションは `b0=−1`。r ≈ 4.3 RJ で両者が打ち消し合い、追跡がそこで往復する | 開閉の判定が成り立たない |
| prodp3 | `ardr` が未初期化 | 判定が実行環境に依存 |
| prodpc5 | 写像点 (r=10) での速度の x, y 成分をそのまま電離圏の水平成分として使っている | z 成分が捨てられ、子午面方向の写像が正しくない。写像点の余緯度が 54.7° より大きい (足元が約 15° より低緯度) と、外向きの運動が極向きになるのに赤道向きとして描かれる (符号が逆) |
| prodpc5 | 電場の写像を等方的な sqrt(B) の比だけで扱っている | 双極子写像は子午面方向と方位角方向で縮尺が約 2 倍違うので、渦度 (= FAC) が歪む |
| polar | 速度に人工的な窓関数 `1/(1+R^8)` と `(1−r²)^0.25` を掛けてから rot をとっている | 窓関数の勾配で、外縁 (70° 付近) に偽の FAC のリングができる (旧図の外周の赤い帯) |
| polar | Poisson 方程式を ADI 法で解いているが、反復パラメータ固定・最大 100 回で打ち切り。正方形全体 (極冠の外の角も含む) を解いている | 最も大きな構造 (最低次モード) は真の解の 8% にしかならない (旧ルーチンで確認)。電位図の正方形の角の等値線や外縁に偏った等値線はこのため |
| polar | 発散なし流速を `max|u|/max|v_df|` で再規格化 | 振幅の誤りが図から分からなくなる |
| polar | FAC を「流速の渦度 × 一様伝導度」で求めている | MHD の電流そのものではない。新コードでは MHD の J∥ を写像した値を主に描き、対流からの値は比較用に残した |

### 検証 (`test/polar_test.py`, `make test` に含む)

双極子 + 一様場 + 局在したねじれ電流 + (共回転 + dL/dt = 一定 の流れ) の合成データで解析解と比較した。

| 項目 | 結果 |
|---|---|
| j_mhd (解析的 rot B から計算した値) | 相対誤差 2.7×10⁻³ |
| 写像した対流 u_colat (解析解 −c/(2L√(L−1))) | 相対誤差 2.4×10⁻³ |
| 開閉境界 (同じ場を python で追跡した境界 13.95°) | 格子 1 個 (0.25°) 以内 |
| 共回転の写像 u_mlt = Ω sinθ, 電位 Ω(sin²θ − sin²θmax) | 相対誤差 5×10⁻⁷, 3×10⁻⁵ |
| 北半球と南半球 (鏡映) | 完全一致 |

## mag3d (3 次元可視化: 旧 VRML ツールの後継)

```sh
OMP_NUM_THREADS=16 ./bin/mag3d namelist/mag3d.nml      # → ju02nt3d_001.glb, ...
```

荻野教授の VRML サブルーチン群 (`zvrsub*.f`) と、それを使った `zvrmag*.f` (磁力線), `zvrmpx*.f` (多重断面), `mcub*.f` (等値面), `zvrveld.f` (速度ベクトル) を 1 本にまとめた。結合済みの `m.data` を直接読み、1 時刻を 1 つの 3 次元シーンとして書き出す。

### 出力形式の選択

| 形式 | 長所 | 短所 | 本ツールでの扱い |
|---|---|---|---|
| VRML 2.0 (.wrl) | 旧資産。テキストで生成が簡単 | 規格・ビューアとも開発終了。ASCII なので大きく遅い | 出力しない。旧ファイルは `x3d_view.html` (X_ITE) でそのまま表示できる |
| **X3D (.x3d)** | VRML の ISO 後継。ノード構成がほぼ同じで移行しやすい。文字・視点・アニメーションを規格内で持てる | XML で大きい (VRML と同程度)。対応ソフトが少ない | `out_x3d=.true.` で出力。画像は外部 PNG |
| **glTF 2.0 (.glb)** | バイナリで小さく速い (X3D の約半分)。three.js / Blender / ParaView / Windows 3D ビューアなど対応が広い。WebGL にそのまま渡せる | 文字 (ラベル) は規格外 (extras に入れてビューアで表示) | **既定**。`viewer/viewer.html` で表示 |
| ParaView (HDF5/XDMF) | 体積レンダリング・対話的な流線など本格的な解析 | 結果の共有にはソフトが必要 | 既存の `bin2hdf` で対応 |

大規模データでは「重い処理 (磁力線追跡・等値面抽出) はスパコン上の Fortran で行い、軽いジオメトリ (.glb) だけを手元のブラウザで見る」という分担にした。3 次元データそのものを手元に持ってくる必要はない。

### 描けるもの

| 要素 | namelist | 旧ルーチン |
|---|---|---|
| 磁力線 (極域の足元格子から, 緯度 × MLT) | `fl_polar`, `fl_hemi`, `fl_lat_min/max`, `fl_nlat`, `fl_nmlt`, `fl_r_start` | ainte1 / ainte1d |
| 磁力線 (赤道面の格子から, 両方向) | `fl_equator`, `eq_xmin` … `eq_ny` | ainte21 / ainte21d |
| 磁力線 (任意の点から) | `fl_npt`, `fl_px/py/pz` | — |
| 磁力線 (3 次元の格子点から。惑星に繋がらない線用) | `fl_box`, `fl_box_lo/hi`, `fl_box_n`, `fl_keep` | — |
| 開閉境界付近の磁力線だけ選ぶ | `fl_select='ocb'` | ainte1a + zsub33d |
| 磁力線の色 | `fl_color='topology'` (閉: 緑 / 開: 青 / 太陽風: 紫) または物理量名 | 固定色 |
| 断面 (カラー画像, 半透明) | `pl_plane`, `pl_pos`, `pl_qty`, `pl_vmin/vmax`, `pl_alpha`, `pl_alpha_mode` | pixel1 / pixeco / pixect |
| 多重断面 (擬似ボリューム, 3 方向まで同時) | `stk_plane(1:3)`, `stk_n`, `stk_lo/hi`, `stk_qty`, `stk_vmin/vmax`, `stk_tra`, `stk_alpha` (最大不透明度, 既定 0.4) (例: `stk_plane = 'yz','xz','xy'`)。ビューアの「多重断面」スライダーで全組の透明度を一度に変えられる (「断面」スライダーも掛かる) | zvrmpx |
| 選択断面 (ビューアで 1 枚を選んで表示, 3 組まで) | `sel_plane`, `sel_n`, `sel_lo/hi`, `sel_qty`, `sel_vmin/vmax`, `sel_alpha`, `sel_palette` | — |
| 等値面 | `iso_qty`, `iso_level`, `iso_color`, `iso_alpha`, `iso_stride`, `iso_max_tri` (既定 200 万。超えると自動で間引く), `iso_smooth` (取る前の平滑化回数, 既定 1) | csuba |
| ベクトル矢印 (速度・磁場) | `arr_plane`, `arr_pos`, `arr_vec`, `arr_n` | arrow / zvrveld |
| 枠・座標軸・目盛・タイトル・惑星 | `draw_frame`, `draw_planet`, `draw_inner`, `label_size` | linebe/lineen, symblv2, sphere1 |
| 惑星の画像・土星の環 | `planet_texture` ('auto' / ファイル名 / 'none'), `planet_lon0`, `draw_rings`, `ring_r1/r2` | — |
| カラーバー | ビューアが自動表示。`colorbar_3d=.true.` でシーン内にも描く | — |

物理量の名前は slice_plot と同じ (rho, vx … bz, temp, absb, absbz, absv, beta, pdyn, log10:<量> など)。半領域データ (z≥0) は、磁力線・断面・等値面とも parity を使って z<0 側まで描く。
`pl_plane(n)` を与えた番号には `pl_qty(n)` も必ず与える (無いとエラー。`pl_plane(2)= ..., pl_qty(1)= ...` のような添字の書き間違いを検出する)。
値は `&units` の物理単位 (`pl_vmin/vmax`, `stk_vmin/vmax`, `iso_level`, `fl_vmin/vmax`, `arr_vmax` もこの単位)。
各断面の上に「量 [単位] 色範囲」のラベルを付ける (`plane_labels=.false.` で消す)。
実行ログには断面ごとに色範囲とデータの範囲 (data range) が出る。データが無い・一定・色範囲外のときは WARNING が出るので、
断面が真っ黒・透明になったときはまずここを見る。

### 惑星に繋がらない磁力線

赤道面・任意点・3 次元格子点 (`fl_box`) から両方向に追跡した磁力線は、両端の行き先で分類する。

| 層の名前 | 分類 | 色 (`fl_color='topology'`) |
|---|---|---|
| `fieldlines_closed` | 両端とも惑星 | 緑 |
| `fieldlines_open` | 片端が惑星, 片端が領域外 | 青 |
| `fieldlines_solarwind` | 両端とも領域外 (太陽風・IMF) | 紫 |
| `fieldlines_detached` | 惑星にも領域外にも出ない (閉じたループ, 最大長に達した線) | 灰 |

`fl_keep` で残す線を選ぶ (`'unconnected'` = solarwind + detached)。極域の足元から描く線 (`fl_polar`) は常に描く。
出発点に戻ってくる閉じたループは 1 周で止める。例:

```
fl_box = .true., fl_box_lo = -60.0, -60.0, -30.0, fl_box_hi = 150.0, 60.0, 30.0, fl_box_n = 8, 6, 3,
fl_keep = 'unconnected',
```

### 選択断面

多重断面と同じく同じ向きの断面を `sel_n` 枚作っておき、ビューアの「選択断面」スライダーで表示する 1 枚を選ぶ
(「全部表示」で全部)。色範囲を与えなければ全体場の最小〜最大で全ての面に共通。断面の位置は glb の node extras
(`group`, `index`, `pos`) に入っている。

```
sel_plane(1) = 'yz', sel_n(1) = 21, sel_lo(1) = -50.0, sel_hi(1) = 250.0, sel_qty(1) = 'temp',
sel_vmin(1) = 0.0, sel_vmax(1) = 300.0,
```

### カラーバー

断面・多重断面・選択断面・物理量で着色した磁力線・矢印の色範囲は glb の `scenes[0].extras.legends` に入れ、
viewer.html が画面右下に表示する (見えている層の分だけ。「カラーバー」チェックで消せる)。「画像を保存」にも入る。
Blender や X3D で見るときは `colorbar_3d = .true.` でシーンの中 (枠の下) にもカラーバーを描く。

### 惑星の画像

`planet_texture = 'auto'` (既定) は `textures/<planet>.jpg` (または .png) を球に貼る (北極が画像の上端、
画像の中央 = 経度 0 が太陽方向。`planet_lon0` [deg] で回す)。土星は環 (C 環〜A 環, 1.24〜2.27 RS) も描く。

| ファイル | 内容 |
|---|---|
| `textures/earth.jpg` | Natural Earth I (パブリックドメイン, naturalearthdata.com) を 1024×512 にしたもの |
| `textures/jupiter.jpg` | 帯構造と大赤斑を模した**模式図** (`textures/make_textures.py` で生成。実画像ではない) |
| `textures/saturn.jpg` | 帯構造を模した**模式図** (同上) |

実際の画像に替えるときは、同じ名前の正距円筒図法 (横:縦 = 2:1) の JPEG/PNG で置き換えるか、
`planet_texture = '/path/to/map.jpg'` で指定する。パブリックドメインの画像として NASA/JPL Photojournal の
Cassini による木星の全球図 (PIA07782) などがある。画像は glb に埋め込まれる (X3D では `<outbase>_img1.jpg` を書く)。

### 処理方法

* 全体場は `stride` 間引きで、惑星から `r_fine` 以内は全解像度で、必要な変数だけを読む (行単位で必要な行だけ読むので、読み込み量は約 1/stride²)
* 磁力線は 2 次 Runge–Kutta で追跡し、曲率と長さに応じて点を間引いて出力する (旧 VRML と同じ考え方)。足元から `fl_r_start` までは双極子磁力線で描く
* 等値面は marching tetrahedra。格子の辺ごとに頂点を共有させ、法線をスカラー場の勾配から求めるので、滑らかに陰影がつき、ファイルも小さい
* 磁力線追跡・断面のサンプリング・等値面抽出は OpenMP 並列
* 座標は表示座標 (X 太陽向き, Z 北) から Y-up の glTF/X3D 座標へ回転で写す。旧 VRML は (x, z, y) と書いていたため鏡像になっていた

### 見方

* **glb**: `viewer/viewer.html` をブラウザで開き、「ファイルを開く」で .glb を選ぶ (複数選ぶと時系列アニメーション)。ドラッグ & ドロップも可。
  層ごとの表示切替、断面・等値面の不透明度、視点 (太陽側・尾部・夕方側・北・斜め)、白背景、PNG 保存ができる。
  http サーバ経由なら `viewer.html?files=a_001.glb,a_002.glb` で直接開ける (`python3 -m http.server`)。
  Blender (File → Import → glTF) で読めば高品質な静止画・動画を作れる。ParaView 5.10 以降でも開ける
* **x3d / 旧 wrl**: `viewer/x3d_view.html?file=demo.x3d` (http サーバ経由。PNG を同じ場所に置く)。旧 VRML ツールの出力 (.wrl) もこのページで表示できる (時系列表示に使う Inline + TimeSensor の仕組みも動作を確認)
* どちらのビューアも three.js / X_ITE を CDN から読むのでインターネット接続が必要 (オフラインでは各ファイル冒頭のコメント参照)

### 性能 (合成データでの実測, 2 コアの計算機)

| データ | 設定 | 時間 | メモリ | 出力 |
|---|---|---|---|---|
| 722×482×242 × 8 変数 (2.7 GB) | stride 2, 磁力線 1920 本, 断面 2, 等値面 1 (38 万三角形), 矢印 | 1 スレッド 9.7 s / 2 スレッド 6.1 s | 0.44 GB | .glb 13 MB, .x3d 23 MB |

旧 VRML ツールの出力例は磁力線約 1000 本と断面 3 枚で 4.4 MB の ASCII だった。標準格子 (3002×2002×1002) では stride 4 で全体場 1 変数 0.38 GB になる。

### 検証

`make test` で、合成データから .glb と .x3d を書き、glTF の構造を公式 validator (gltf-validator, node が使える場合) で確認する (エラー・警告 0)。ブラウザ表示は Chromium (three.js 0.170 / X_ITE 16.5) で確認した。

---

## ファイル形式

* `m.data` / `.vort`: Fortran sequential unformatted。`do m=1,nvar; do k=1,nz2; write(u) f(1:nx2,1:ny2,k)`。1 レコード = 1 つの z 面 (nx2*ny2 個の float32)
* 新ツールはこれを stream アクセスで読み、(m, k) のレコードへ直接移動する。位置は `(irec-1)*(4*nx2*ny2 + 2*marker) + marker + 1` (irec = (m−1)*nz2 + k)
* レコードマーカーの大きさとエンディアンは先頭・末尾マーカーから自動判定する (gfortran の `-frecord-marker=8`, `-fconvert=big-endian` で作ったファイルでも読める)

## テスト

```sh
make && cd test && ./run_tests.sh [旧コードのあるディレクトリ]
```

小さい合成データで各ツールを実行し、解析解や独立な計算と比較する (python3 + numpy が必要)。
combine3d (通常 / big endian / 8 byte マーカー × 1, 2, 4 スレッド) は生成式の値と、vorticity は numpy で計算した値と、
polar_map は解析解と比較する。bin2hdf (raw, HDF5) は入力と bit 単位で一致することを確認する。
旧コード (`data_conbine3d_sxt_mem.f90`, `vorticity_ju02nt.f90`, `fig_h/`) のあるディレクトリを引数に与えると、
旧コードを小さい格子に書き換えてビルド・実行し、その出力とも比較する
(combine3d は bit 一致、vorticity は倍精度化のため相対 1e-4 以内)。
物理単位への換算は bin2hdf (raw/XDMF/HDF5 属性) と slice_plot (規格化版との比) で確認する。
旧コードがあれば slice_plot は旧・新の PostScript と PNG を `test/work/s/` に出力するので目で比較する。

## 性能 (小規模テストでの実測)

2 コアの計算機でページキャッシュを毎回破棄して測定。実機 (Lustre, 多コア) での値は別途測る必要がある。

| ツール | テスト規模 | 旧 | 新 (1 thread) | 新 (2 threads) |
|---|---|---|---|---|
| combine3d | 256 ランク, 入力 4.1 GB → 出力 545 MB | 6.6〜8.1 s, 549 MB | 0.8〜1.7 s, 69 MB | 1.1〜1.4 s |
| vorticity | 512x256x130 | 1.24 s, 718 MB | 0.71 s, 10 MB | 0.77 s |
| bin2hdf | 512x256x130 x 8 変数 | 1.0〜1.2 s | 0.64〜0.99 s | 0.79〜1.07 s |
| slice_plot | 302x202x102, 4 パネル | — (全データ読み込み) | 0.2 s | — |

combine3d の主な効果は読み込み量の削減 (標準設定で 1/6) とメモリ削減なので、実機ではテストより差が大きくなると見込まれる。
OpenMP の効果は並列ファイルシステムで同時ストリーム数が効く場合に出る。`OMP_NUM_THREADS` を 8, 16, 32 … と変えて計測し、最も速い値を使う。

## 保守メモ

* **格子・パスの変更**: namelist だけを書き換える
* **描く量の追加**: `slice_plot.f90` の `mark_needs` (必要な変数) と `eval_qty` (計算式) に case を 1 つずつ足す。mag3d 用には `mhdfield.f90` の `qty_parse` / `qty_needs` / `qty_eval` に足す
* **カラーマップの追加**: `psplot.f90` の `palette_rgb` に case を足す
* **並列化の注意**: 並列領域内で可変長文字列を返す関数や実行時書式の内部 write を使わない (gfortran で競合する)。`make_fname`, `zpad` は固定長を返すので `trim()` して使う
* **古い Intel コンパイラ**: 内部手続きからホスト関連付けの配列を書き換える書き方を避ける (ifx 2023 で正しく動かなかった)。mag3d の必要変数の収集はこのためメインプログラムで直接行っている
* **NEC/Fujitsu のレコードマーカー**: 自動判定するので設定不要。判定できない場合は「record length mismatch」で止まる (格子サイズの設定誤りもこのエラーになる)

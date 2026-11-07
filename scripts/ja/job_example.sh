#!/bin/bash
#==============================================================================
# viztools のバッチジョブ例 (一般的な Slurm 用)
#
#   * パーティション名・時間・メモリ・コア数は各計算機の利用案内に合わせて書き換える
#   * データファイル (m.data など) のあるディレクトリで投入する:
#         cd /path/to/data && sbatch /path/to/viztools/scripts/job_example.sh
#     namelist の入力・出力ディレクトリは既定で './' (投入したディレクトリ)
#   * すべてのツールは 1 ノードの OpenMP 並列。スレッド数は I/O の同時ストリーム数にもなる
#==============================================================================
#SBATCH --job-name=viztools
#SBATCH --partition=test             # 例: 各計算機のパーティション名に変更
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32           # OpenMP スレッド数
#SBATCH --mem=64G                    # 標準格子 (3002x2002x1002) の mag3d / polar_map は 8〜16 GB 程度
#SBATCH --time=02:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

VIZ=${VIZ:-$HOME/viztools}           # viztools を置いた場所
BIN=$VIZ/bin
NML=$VIZ/namelist/ja                  # 自分用に写した namelist のディレクトリに変えてもよい

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OMP_STACKSIZE=64M
export OMP_PROC_BIND=close OMP_PLACES=cores

# 使うものだけ残す
srun $BIN/combine3d  $NML/combine.nml          # ランク毎ファイル -> 結合した m.data
srun $BIN/vorticity  $NML/vort.nml             # 渦度 (.vort)
srun $BIN/slice_plot $NML/slice_legacy.nml     # 断面図 (PostScript)
srun $BIN/polar_map  $NML/polar.nml            # 極域電離圏図 (PostScript)
srun $BIN/mag3d      $NML/mag3d.nml            # 3 次元シーン (.glb)
srun $BIN/bin2hdf    $NML/bin2hdf_mhd.nml      # ParaView 用 HDF5 + XDMF

#!/bin/bash
#==============================================================================
# example batch job for viztools (generic Slurm)
#
#   * change the partition, time, memory and cores for your system
#   * submit in the directory of the data files (m.data etc.):
#         cd /path/to/data && sbatch /path/to/viztools/scripts/job_example.sh
#     input/output directories in the namelists are './' (the submit directory) by default
#   * all tools are OpenMP parallel on one node; the threads are also the concurrent I/O streams
#==============================================================================
#SBATCH --job-name=viztools
#SBATCH --partition=test             # e.g.: change to a partition of your system
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32           # OpenMP threads
#SBATCH --mem=64G                    # mag3d / polar_map on the standard grid (3002x2002x1002) need about 8-16 GB
#SBATCH --time=02:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

VIZ=${VIZ:-$HOME/viztools}           # where viztools is installed
BIN=$VIZ/bin
NML=$VIZ/namelist                    # or a directory with your own copies of the namelists

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OMP_STACKSIZE=64M
export OMP_PROC_BIND=close OMP_PLACES=cores

# keep only the tools you need
srun $BIN/combine3d  $NML/combine.nml          # per-rank files -> combined m.data
srun $BIN/vorticity  $NML/vort.nml             # vorticity (.vort)
srun $BIN/slice_plot $NML/slice_legacy.nml     # slice plots (PostScript)
srun $BIN/polar_map  $NML/polar.nml            # polar ionosphere maps (PostScript)
srun $BIN/mag3d      $NML/mag3d.nml            # 3-D scene (.glb)
srun $BIN/bin2hdf    $NML/bin2hdf_mhd.nml      # HDF5 + XDMF for ParaView

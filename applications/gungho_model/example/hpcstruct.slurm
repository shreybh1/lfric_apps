#!/bin/bash
#SBATCH --job-name=lfric
#SBATCH --time=0:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --account=n02-NGARCH
#SBATCH --account=pawsey0835
#SBATCH --partition=debug

export OMP_NUM_THREADS=1

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

module use /software/projects/pawsey0835/lparisi/setonix/2025.08/modules/zen3/gcc/14.2.0
module use /software/projects/pawsey0835/ddeeptimahanti/setonix/2025.08/modules/zen3/gcc/14.2.0
module load boost/1.86.0-b7uxjf5 \
  dyninst/13.0.0-3y45tzu \
  intel-tbb/2021.12.0-pkpjkld \
  xerces-c/3.3.0-transcoder-gnuiconv \
  elfutils/0.190-xlg3dux \
  yaml-cpp/0.8.0-tzx3y2d \
  hpctoolkit/2024.01.1-l44fgik

export MPICH_OFI_RMA_STARTUP_CONNECT=1

echo "Job started " $(date +"%T") 
hpctoolkit_dir=$(find . -name "hpctoolkit-gungho_model*")
hpcstruct ${hpctoolkit_dir} 
hpcprof ${hpctoolkit_dir}
echo "Job ended " $(date +"%T") # start time

sacct -j $SLURM_JOB_ID --format=JobID,User,NNodes,Elapsed,MaxRSS

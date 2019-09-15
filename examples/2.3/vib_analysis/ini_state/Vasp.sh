#!/bin/sh
 #PBS -N 107V1S_iV
 #PBS -l nodes=1:ppn=8
 #PBS -q batch
 #PBS -l walltime=240:00:00

ulimit -s unlimited

export OMP_NUM_THREADS=1
export I_MPI_COMPATIBILITY=4
source /opt/intel/composer_xe_2015.6.233/bin/compilervars.sh intel64
source /opt/intel/impi/5.0.3.049/intel64/bin/mpivars.sh

cd $PBS_O_WORKDIR

export PAR_RUN="/opt/intel/impi/5.0.3.049/intel64/bin/mpirun"
#export VASP="/home/apps/vasp541_gam"
#export VASP="/home/apps/vasp541_gam_wannier"
#export VASP="/home/apps/vasp541_ncl"
#export VASP="/home/apps/vasp541_ncl_wannier"
#export VASP="/home/apps/vasp541_std"
#export VASP="/home/apps/vasp541_std_wannier"
#export VASP="/home/apps/vasp544_gam"
#export VASP="/home/apps/vasp544_ncl"
#export VASP="/home/apps/vasp544_std"
#export VASP="/home/swj/bin/GF_surface_mpi_v5.0"
#export VASP="/home/swj/bin/TBA_read_wannier90_only_v5.0"
export VASP="/home/zhangzd/software/vasp.5.4.1/bin/vasp_std"

cp $PBS_NODEFILE node
NCORE=`cat node | wc -l`

date > output.$PBS_JOBID

$PAR_RUN -machinefile node -np $NCORE $VASP >log 2>error

date >> output.$PBS_JOBID

date > ~/job_detail/$PBS_JOBNAME.$PBS_JOBID
pwd >> ~/job_detail/$PBS_JOBNAME.$PBS_JOBID
tail -2 output.$PBS_JOBID >> ~/job_detail/$PBS_JOBNAME.$PBS_JOBID







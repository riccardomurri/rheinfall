#! /bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -pe mpich2_rsh 8
#$ -j y

# load modules
source /etc/profile
source /panfs/panfs0.ften.es.hpcn.uzh.ch/share/software/Modules/default/init/sh

# load MPI - must match what the binary was compiled with!
module load gcc/4.4.3

# MVAPICH is not supported by modules, apparently
export PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/gcc-4.4.1/mvapich2-1.2p1/bin:$PATH
export LD_LIBRARY_PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/gcc-4.4.1/mvapich2-1.2p1/lib:$LD_LIBRARY_PATH
export LD_RUN_PATH=$LD_LIBRARY_PATH


echo === environment info ===
echo node: $(hostname)
echo hostfile: $PE_HOSTFILE
echo === hostfile contents ===
test -n "$PE_HOSTFILE" && cat $PE_HOSTFILE
#echo === environment ===
#env | sort 
echo =========================
echo


# main
set -e #-x

LD_LIBRARY_PATH=$HOME/sw/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

PATH=$HOME/sw/bin:$PATH
export PATH

# use Google's tcmalloc
LD_PRELOAD=$HOME/sw/lib/libtcmalloc.so
export LD_PRELOAD


# create hostfile with 1 MPI rank per host
uniq < $TMPDIR/machines > $TMPDIR/hostfile
NHOSTS=$(wc -l $TMPDIR/hostfile | cut -d' ' -f1)


echo
echo === pure OpenMP ===
ulimit -c 0
ulimit -v 24000000
./src/rank_omp $@


echo
echo === pure MPI ===
mpirun_rsh -rsh -hostfile $TMPDIR/machines -np $NSLOTS \
    ./src/rank_mpi $@


# hybrid mode measurements
for np in 1 2 4 8; do
    export OMP_NUM_THREADS=$(expr $NSLOTS / $np)
    echo
    echo === OpenMP:$OMP_NUM_THREADS MPI:$np ===
    mpirun_rsh -rsh -hostfile $TMPDIR/hostfile -np $np \
        ./src/rank_mpi_omp $@

echo === Done. ===

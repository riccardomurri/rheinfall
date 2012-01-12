#! /bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -pe orte 16
#$ -j y


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
#set -e #-x

LD_LIBRARY_PATH=$HOME/sw/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

PATH=$HOME/sw/bin:$PATH
export PATH

# use Google's tcmalloc
LD_PRELOAD=$HOME/sw/lib/libtcmalloc.so
export LD_PRELOAD


# hybrid mode measurements
export OMP_NUM_THREADS=$(expr $NSLOTS / 2)
echo
echo === OpenMP:$OMP_NUM_THREADS MPI:$(expr $NSLOTS / $OMP_NUM_THREADS) ===
mpirun --pernode ./src/rank_mpi_omp $@


echo
echo === pure MPI ===
mpirun --bind-to-core -np $NSLOTS ./src/rank_mpi $@


echo === Done. ===

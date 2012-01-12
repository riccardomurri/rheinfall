#! /bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -pe orte 8
#$ -j y


echo === environment info ===
echo node: $(hostname)
echo hostfile: $PE_HOSTFILE
echo === hostfile contents ===
test -n "$PE_HOSTFILE" && cat $PE_HOSTFILE
#echo === environment ===
#env | sort 
echo === ulimit -a ===
ulimit -a
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


echo
echo === pure OpenMP ===
ulimit -c 0
ulimit -v 15000000
ulimit -s 10240
./src/rank_omp $@


# hybrid mode measurements
for np in 1 2 4 8; do
    export OMP_NUM_THREADS=$(expr $NSLOTS / $np)
    echo
    echo === OpenMP:$OMP_NUM_THREADS MPI:$np ===
    mpirun --bind-to-core -np $np ./src/rank_mpi_omp $@
done


echo
echo === pure MPI ===
mpirun --bind-to-core -np $NSLOTS ./src/rank_mpi $@


echo === Done. ===

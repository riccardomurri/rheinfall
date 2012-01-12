#! /bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -q mpi.q
#$ -pe orte 8
#$ -j y

set -e

PROG=$(basename "$0"  .sh)

export slots_per_node=${slots_per_node:-`grep -c ^processor /proc/cpuinfo`}

LD_LIBRARY_PATH=$HOME/sw/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

PATH=$HOME/sw/bin:$PATH
export PATH

#OMP_NUM_THREADS=$slots_per_node
#export OMP_NUM_THREADS

ulimit -c 0
exec $HOME/sw/bin/mpirun \
        --pernode -np $(expr $NSLOTS / $slots_per_node) \
        -x PATH \
        -x LD_LIBRARY_PATH \
    ./src/rank_mpi_omp $@

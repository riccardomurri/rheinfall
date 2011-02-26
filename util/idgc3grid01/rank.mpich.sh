#! /bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -q mpi.q
#$ -pe mpich 16
#$ -j n

set -e

PROG=$(basename "$0"  .sh)

export slots_per_node=${slots_per_node:-`grep -c ^processor /proc/cpuinfo`}


LD_LIBRARY_PATH=$HOME/sw/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

PATH=$HOME/sw/bin:/opt/gridengine/bin/lx26-amd64:$PATH
export PATH


#OMP_NUM_THREADS=$(expr $slots_per_node)
#export OMP_NUM_THREADS


# create hostfile, with one rank per node
hostfile=`mktemp -t ${PROG}.XXXXXX`
trap "{ rm -f $hostfile; }" EXIT INT ABRT QUIT
hostfile=`mktemp`
cat $PE_HOSTFILE | \
    (while read host rest; do
      echo "${host}:1" >> $hostfile
    done)

ulimit -c 0
ulimit -v 15000000
exec $HOME/sw/bin/mpiexec \
         -print-rank-map \
        -genvlist HOME,PATH,LD_LIBRARY_PATH \
        -f $hostfile \
        -np $(expr ${NSLOTS:-8} / $slots_per_node ) \
    $HOME/rheinfall/src/rank_mpi $@

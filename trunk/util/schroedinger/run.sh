#! /bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -pe openmpi2 8
#$ -j y
#$ -R y

exec="$1"; shift

flavor="$(basename "$(dirname "$exec")" )"
compiler="$(echo $flavor | cut -d- -f1)"
if [ -z "$compiler" ]; then
    compiler=gcc443
fi
mpi="$(echo $flavor | cut -d- -f2)"
if [ -z "$mpi" ]; then
    mpi=ompi
fi

# load modules
source /panfs/panfs0.ften.es.hpcn.uzh.ch/share/software/Modules/default/init/sh

# load MPI - must match what the binary was compiled with!
case "$mpi" in
    ompi) 
        case "$compiler" in
            gcc*) module load mpi/openmpi/gcc ;;
            icc|intel) module load mpi/openmpi/intel ;;
        esac
        ;;
    mvapich) # MVAPICH is not supported by modules, apparently
        export PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/gcc-4.4.1/mvapich2-1.4rc2/bin:$PATH
        export LD_LIBRARY_PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/gcc-4.4.1/mvapich2-1.4rc2/lib:$LD_LIBRARY_PATH
        export LD_RUN_PATH=$LD_LIBRARY_PATH
        ;;
    impi|intel) # Intel MPI is not supported by modules, apparently
        export PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/intel_mpi/3.2.1.009/bin64:$PATH
        export LD_LIBRARY_PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/intel_mpi/3.2.1.009/lib64:$LD_LIBRARY_PATH
        export LD_RUN_PATH=$LD_LIBRARY_PATH
        ;;
esac

# load the compiler and libs
case "$compiler" in
    gcc443)
        module load gcc/4.4.3
        ;;
    gcc450)
        module load gcc/4.5.0
        ;;
    icc|intel)
        module load intel/comp/11.1.064
        ;;
esac

cpus_per_node=$(grep -i -c '^processor' /proc/cpuinfo)

echo === running info ===
echo exec: $exec
echo compiler: $compiler
echo mpi: $mpi
echo === environment info ===
echo node: $(hostname)
echo NSLOTS: $NSLOTS
echo CPUs per node: $cpus_per_node
echo hostfile: $PE_HOSTFILE
echo === hostfile contents ===
test -n "$PE_HOSTFILE" && cat $PE_HOSTFILE
#echo === environment ===
#env | sort 
echo === ulimit -a ===
ulimit -a
echo =========================
echo


## main
#set -ex

sw="$HOME/sw/${compiler}-${mpi}"

LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

PATH=${sw}/bin:$PATH
export PATH

# use Google's tcmalloc
LD_PRELOAD=${sw}/lib/libtcmalloc.so
export LD_PRELOAD

# no core dumps; especially with MPI, they just fill the home directory...
ulimit -c 0

case "$exec" in
    *rank|*rank-omp) # no MPI... just run it!
        $exec "$@"
        ;;
    *rank-mpi) # pure MPI
        case "$mpi" in 
            none) $exec "$@" ;;
            ompi) mpirun -np $NSLOTS --bynode --nooversubscribe --bind-to-core \
                -x LD_PRELOAD -x LD_LIBRARY_PATH -x PATH \
                $exec "$@" ;;
            mvapich) mpirun_rsh -rsh -hostfile $TMPDIR/machines -np $NSLOTS $exec "$@" ;;
            intel|impi) mpiexec -n $NSLOTS \
                -genvlist LD_PRELOAD,LD_LIBRARY_PATH,PATH \
                -env I_MPI_DEVICE rdssm \
                $exec "$@" ;;
        esac
        ;;
    *rank-mpi-omp) # hybrid MPI+OpenMP
        np=$(expr $NSLOTS / $cpus_per_node)
        case "$mpi" in 
            none) $exec "$@" ;;
            ompi) mpirun -np $np \
                --bynode --nooversubscribe \
                -x LD_PRELOAD -x LD_LIBRARY_PATH -x PATH \
                $exec "$@" ;;
            mvapich) 
                # filter hostfile so that we only have one occurrence of each host
                uniq < $TMPDIR/machines > $TMPDIR/machines.uniq
                mpirun_rsh -rsh -hostfile $TMPDIR/machines.uniq -np $np \
                $exec "$@" ;;
            intel|impi) mpiexec  -perhost 1 \
                -genvlist LD_PRELOAD,LD_LIBRARY_PATH,PATH \
                -env I_MPI_DEVICE rdssm \
                -n $np \
                $exec "$@" ;;
        esac
        ;;
esac
rc=$?

echo === Done: exit code $rc ===
exit $rc

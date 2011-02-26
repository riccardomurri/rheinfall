#! /bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -pe openmpi2 8
#$ -j y
#$ -R y

PROG="$(basename $0)"

usage () {
cat <<EOF
Usage: $PROG [COMPILER]-[MPILIB]-[REVNO] PROG [ARGS ...]

Setup the correct modules and paths for supporting binaries compiled
by COMPILER with MPILIB, and finally run
'run/COMPILER-MPILIB-REVNO/PROG' with the given ARGS.

The COMPILER tag selects the compiler type (gcc412, gcc441, gcc443,
gcc450, intel); second tag MPILIB selects which MPI library will be
used (ompi, mvapich, intel, or none).

EOF
}

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

## helper functions
die () {
  rc="$1"
  shift
  (echo -n "$PROG: ERROR: ";
      if [ $# -gt 0 ]; then echo "$@"; else cat; fi) 1>&2
  exit $rc
}

have_command () {
  type "$1" >/dev/null 2>/dev/null
}

require_command () {
  if ! have_command "$1"; then
    die 1 "Could not find required command '$1' in system PATH. Aborting."
  fi
}

_ () {
    echo
    echo ==== "$@" ...;
}


## load modules
source /panfs/panfs0.ften.es.hpcn.uzh.ch/share/software/Modules/default/init/sh

supported_compilers='gcc412 gcc441 gcc443 gcc450 icc'
supported_mpilibs='openmpi mvapich intel none'

# load MPI - must match what the binary was compiled with!
case "$mpi" in
    ompi|openmpi) 
        case "$compiler" in
            gcc450)    module load mpi/openmpi/gcc-4.5.0 ;;
            gcc*) module load mpi/openmpi/gcc ;;
            icc|intel) module load mpi/openmpi/intel ;;
        esac
        ;;
    mvapich) # MVAPICH is not supported by modules, apparently
        case "$compiler" in
            gcc412)
                export PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/gcc/mvapich2-1.4rc2/bin:$PATH
                export LD_LIBRARY_PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/gcc/mvapich2-1.4rc2/lib:$LD_LIBRARY_PATH
                export LD_RUN_PATH=$LD_LIBRARY_PATH
                ;;
            gcc*)
                export PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/gcc-4.4.1/mvapich2-1.4rc2/bin:$PATH
                export LD_LIBRARY_PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/gcc-4.4.1/mvapich2-1.4rc2/lib:$LD_LIBRARY_PATH 
                export LD_RUN_PATH=$LD_LIBRARY_PATH
                ;;
        esac
        ;;
    impi|intel) # Intel MPI is not supported by modules, apparently
        export PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/intel_mpi/3.2.1.009/bin64:$PATH
        export LD_LIBRARY_PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/mpi-libs/intel_mpi/3.2.1.009/lib64:$LD_LIBRARY_PATH
        export LD_RUN_PATH=$LD_LIBRARY_PATH
        ;;
    none)
        : # nothing to do
        ;;
    *) 
        die 1 "Unknown MPI library '${mpi}' - please choose one of: $supported_mpilibs"
        ;;
esac

# load the compiler and libs
case "$compiler" in
    gcc412)
        module load gcc/4.1.2
        ;;
    gcc441)
        module load gcc/4.4.1
        ;;
    gcc443)
        module load gcc/4.4.3
        ;;
    gcc450)
        module load gcc/4.5.0
        ;;
    icc|intel)
        module load intel/comp/11.1.064
        ;;
    *)
        die 1 "Unknown compiler flavor '${compiler}' - please choose one of: $supported_compilers"
        ;;
esac


## environment information

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

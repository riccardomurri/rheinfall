#! /bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -pe openmpi2 8
#$ -j y
#$ -R y

PROG="$(basename $0)"

usage () {
cat <<EOF
Usage: $PROG [COMPILER]-[MPILIB]-[REVNO]/PROG [ARGS ...]

Setup the correct modules and paths for supporting binaries compiled
by COMPILER with MPILIB, and finally run
'run/COMPILER-MPILIB-REVNO/PROG' with the given ARGS.

The COMPILER tag selects the compiler type (gcc412, gcc441, gcc443,
gcc450, intel); second tag MPILIB selects which MPI library will be
used (ompi, mvapich, intel, or none).

EOF
}

exec="$1"; shift

source $HOME/rheinfall/util/schroedinger/functions.sh \
    || { echo 1>&2 "Cannot load 'functions.sh' - aborting."; exit 1; }

set_mpi_and_compiler_flavor "$(basename "$(dirname "$exec")" )"


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

case "$mpi" in
    ompi) # system's OpenMPI, all setup should be done by modules
        ;;
    ompi*) # my own OpenMPI
        prepare_openmpi_environment
        ;;
    *) # anything else, ignore
        ;;
esac

# use Google's tcmalloc
LD_PRELOAD=${sw}/lib/libtcmalloc.so
export LD_PRELOAD


case "$exec" in
    *rank*-mpi-omp) # hybrid MPI+OpenMP
        np=$(expr $NSLOTS / $cpus_per_node)
        case "$mpi" in 
            none)
                echo 1>&2 "Cannot run MPI binary '$exec' with MPI flavor 'none'. Aborting."
                exit 2
                ;;
            ompi*) 
                mpirun -np $np \
                    --bynode --nooversubscribe $gmca $hostfile \
                    -x LD_PRELOAD -x LD_LIBRARY_PATH -x PATH \
                    $HOME/bin/_exec.sh $exec "$@" 
                ;;
            mvapich*) 
                # filter hostfile so that we only have one occurrence of each host
                uniq < $TMPDIR/machines > $TMPDIR/machines.uniq
                mpirun_rsh -rsh -hostfile $TMPDIR/machines.uniq -np $np \
                    $HOME/bin/_exec.sh $exec "$@" 
                ;;
            intel|impi) 
                mpiexec  -perhost 1 \
                    -genvlist LD_PRELOAD,LD_LIBRARY_PATH,PATH \
                    -env I_MPI_DEVICE rdssm \
                    -n $np \
                    $HOME/bin/_exec.sh $exec "$@"
                ;;
        esac
        ;;
    *rank*-mpi) # pure MPI
        case "$mpi" in 
            none) 
                echo 1>&2 "Cannot run MPI binary '$exec' with MPI flavor 'none'. Aborting."
                exit 2
                ;;
            ompi*) 
                mpirun -np $NSLOTS --gmca btl tcp,sm,self \
                    --bynode --nooversubscribe --bind-to-core $gmca $hostfile \
                    -x LD_PRELOAD -x LD_LIBRARY_PATH -x PATH \
                    $HOME/bin/_exec.sh $exec "$@" 
                ;;
            mvapich*) 
                mpirun_rsh -rsh -hostfile $TMPDIR/machines -np $NSLOTS $HOME/bin/_exec.sh $exec "$@" 
                ;;
            intel|impi) 
                mpiexec -n $NSLOTS \
                    -genvlist LD_PRELOAD,LD_LIBRARY_PATH,PATH \
                    -env I_MPI_DEVICE rdssm \
                    $HOME/bin/_exec.sh $exec "$@" 
                ;;
        esac
        ;;
    *rank*) # no MPI... just run it!
        $HOME/bin/_exec.sh $exec "$@"
        ;;
esac
rc=$?

echo === Done: exit code $rc ===
exit $rc

#! /bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -R y
#$ -notify 
#$ -j y
#$ -pe openmpi2 8

PROG="$(basename $0)"

usage () {
cat <<EOF
Usage: $PROG [REVNO]/COMPILER-MPILIB/PROG [ARGS ...]

Setup the correct modules and paths for supporting binaries compiled
by COMPILER with MPILIB, and finally run
'run/REVNO/COMPILER-MPILIB/PROG' with the given ARGS.

The COMPILER tag selects the compiler type (gcc412, gcc434, gcc443,
gcc450, intel); second tag MPILIB selects which MPI library will be
used (ompi, parastation, mvapich, intel, or none).

EOF
}

exec="$1"; shift

libdir=$(pwd)/util/schroedinger
source $libdir/functions.sh \
    || { echo 1>&2 "Cannot load 'functions.sh' - aborting."; exit 1; }

spec="$(dirname "$exec")"
if expr "$spec" : 'run/' >/dev/null 2>&1; then
    spec="$(echo $spec | cut -c5-)"
fi
set_mpi_and_compiler_flavor "$spec"

# ensure $exec is an absolute path
if ! expr match "$exec" '/' >/dev/null 2>&1; then
  exec="$(pwd)/$exec"
fi

## detect array jobs and act accordingly

# XXX: option parsing is rheinfall-specific
opts=''
files=''
while [ $# -gt 0 ]; do
    case "$1" in
        -t|--transpose) opts="$opts -t" ;;
        -v|--verbose)   opts="$opts -v" ;;
        -m|--memory)    opts="$opts -m $2"; shift ;;
        -p|--modulus)   opts="$opts -p $2"; shift ;;
        --) shift; break ;;
        *) files="$files $1" ;;
    esac
    shift
done
files="$files $*"

nth () { declare -i n=$1; while [ $n -gt 0 ]; do shift; let n--; done; echo "$1"; }

if [ "_$SGE_TASK_ID" != '_undefined' ]; then
    args_from="$(nth $SGE_TASK_ID $files)"
    files="$(cat $args_from)"
    echo SGE Task $SGE_TASK_ID: replaced argument list with contents of file "'$args_from': $files"
fi


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
echo === start date ===
date
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
                    $libdir/_exec.sh $exec $opts $files 
                ;;
            para*) 
                # filter hostfile so that we only have one occurrence of each host
                uniq < $TMPDIR/machines > $TMPDIR/machines.uniq
                mpiexec -np $np --hostfile $TMPDIR/machines.uniq \
                    -e PATH,LD_LIBRARY_PATH,LD_PRELOAD --loopnodesfirst \
                    $libdir/_exec.sh $exec $opts $files 
                ;;
            mvapich*) 
                # filter hostfile so that we only have one occurrence of each host
                uniq < $TMPDIR/machines > $TMPDIR/machines.uniq
                mpirun_rsh -rsh -hostfile $TMPDIR/machines.uniq -np $np \
                    $libdir/_exec.sh $exec $opts $files 
                ;;
            intel|impi) 
                mpiexec  -perhost 1 \
                    -genvlist LD_PRELOAD,LD_LIBRARY_PATH,PATH \
                    -env I_MPI_DEVICE rdssm \
                    -n $np \
                    $libdir/_exec.sh $exec $opts $files
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
                    $libdir/_exec.sh $exec $opts $files 
                ;;
            para*) 
                mpiexec -np $NSLOTS --hostfile $TMPDIR/machines \
                    -e PATH,LD_LIBRARY_PATH,LD_PRELOAD --loopnodesfirst \
                    $libdir/_exec.sh $exec $opts $files 
                ;;
            mvapich*) 
                mpirun_rsh -rsh -hostfile $TMPDIR/machines -np $NSLOTS $libdir/_exec.sh $exec $opts $files 
                ;;
            intel|impi) 
                mpiexec -n $NSLOTS \
                    -genvlist LD_PRELOAD,LD_LIBRARY_PATH,PATH \
                    -env I_MPI_DEVICE rdssm \
                    $libdir/_exec.sh $exec $opts $files 
                ;;
        esac
        ;;
    *rank*) # no MPI... 
        if type -a parallel >/dev/null 2>&1; then
            # use GNU parallel to run jobs
            parallel --jobs -1 $exec $opts ::: $files
        else
            # no 'parallel' available, execute jobs sequentially
            $exec $opts $files
        fi
esac
rc=$?

echo === end date ===
date
echo === Done: exit code $rc ===
exit $rc

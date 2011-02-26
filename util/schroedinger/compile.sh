#! /bin/bash
#$ -l s_cpu=1800
#$ -cwd
#$ -S /bin/bash
#$ -v PATH,LD_LIBRARY_PATH
#$ -j y
#$ -N compile

PROG="$(basename $0)"

usage () {
cat <<EOF
Usage: $PROG [COMPILER]-[MPILIB]-[REVNO] [CXXFLAGS ...]

Setup the correct modules and paths for compiling 'Rheinfall' with the
given COMPILER and MPILIB.  Compile the code and store the results in
'~/sw/COMPILER-MPILIB-REVNO'.

The COMPILER tag selects the compiler type (gcc412, gcc441, gcc443,
gcc450, intel; default: gcc443); second tag MPILIB selects which MPI
library will be used (ompi, mvapich, intel, or none; default:
openmpi).  Third argument REVNO defaults to the current Bazaar
repository revno, but can be overridden to tag variants of the code.

EOF
}
case `hostname` in
        login*) usage; exit 1 ;;
        *) : ;;
esac

source $HOME/rheinfall/util/schroedinger/functions.sh \
    || { echo 1>&2 "Cannot load 'functions.sh' - aborting."; exit 1; }


## parse command-line 
set_mpi_and_compiler_flavor "$@"; shift
rev="$(echo $flavor | cut -d- -f3)"
if [ -z "$rev" ]; then
    rev="$( cat .bzr/branch/last-revision | (read revno rest; echo r$revno) )"
fi

CXXFLAGS="$*"
if [ -z "$CXXFLAGS" ]; then
  CXXFLAGS='-O3 -DNDEBUG'
fi
CFLAGS="$CXXFLAGS"



## environment information

show_build_information


echo === Compiling Rheinfall ... ===

top_src_dir="$HOME/rheinfall"
#build_dir="/lustre/ESPFS/scratch/oci/murri/rheinfall.${flavor}"
build_dir="$HOME/data/tmp/rheinfall.${flavor}"

set -e
rm -rf "$build_dir"
mkdir -p "$build_dir"
cd "$build_dir"
$top_src_dir/configure --with-ge=yes --with-boost=$sw --with-gmp=$sw \
    CXXFLAGS="$CXXFLAGS $cxxflags" CFLAGS="$CFLAGS $cflags"
set +e
make
rc=$?


echo === Saving executables in "${top_src_dir}/run/${compiler}-${mpi}-${rev}/" ===

mkdir -pv "${top_src_dir}/run/${compiler}-${mpi}-${rev}/"
cp -av \
    $(find src.c++ src.c -type f -perm /u+x) \
    "${top_src_dir}/run/${compiler}-${mpi}-${rev}/"


echo === Run unit tests ===

ln -s ${top_src_dir}/src.c++/test src.c++/test
make check


echo === Done: exit code $rc ===

cd "${top_src_dir}" \
    && rm -rf "$build_dir"

exit $rc

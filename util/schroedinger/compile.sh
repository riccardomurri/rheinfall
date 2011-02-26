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
Usage: $PROG [REVNO]/[COMPILER]-[MPILIB] [CXXFLAGS ...]

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

CXXFLAGS="$*"
if [ -z "$CXXFLAGS" ]; then
  CXXFLAGS='-O3 -DNDEBUG'
fi
CFLAGS="$CXXFLAGS"



## environment information

show_build_information

# figure out which source revision to check out
revno=$(echo $rev | cut -d+ -f1 | cut -c2-)
extra=$(echo $rev | cut -d+ -f2)
if [ -n "$extra" -a "$revno" -ne "$(bzr revno)" ]; then
    die 1 "Asked to tag revision '$rev' but sources are at BZR revno $(bzr revno)"
fi

echo === Compiling Rheinfall r$revno ... ===

top_src_dir="$HOME/rheinfall"
#build_dir="/lustre/ESPFS/scratch/oci/murri/rheinfall.${flavor}"
build_dir="$HOME/data/tmp/rheinfall.${flavor}"

set -e
rm -rf "$build_dir"
# check out BZR sources or compile from working directory
if [ -n "$extra" ]; then
    mkdir -p "$build_dir"
    src_dir="$top_src_dir"
else
    bzr co ./ "$build_dir" -r $revno
    src_dir="$build_dir"
    cd "$build_dir"
    autoreconf -i -s
fi
cd "$build_dir"
${src_dir}/configure --with-ge=yes --with-mpi=yes --with-boost=$sw --with-gmp=$sw \
    CXXFLAGS="$CXXFLAGS $cxxflags" CFLAGS="$CFLAGS $cflags"
set +e
make
rc=$?


echo === Saving executables in "${top_src_dir}/run/${rev}/${compiler}-${mpi}/" ===

mkdir -pv "${top_src_dir}/run/${rev}/${compiler}-${mpi}/"
cp -av \
    $(find src.c++ src.c -type f -perm /u+x) \
    "${top_src_dir}/run/${rev}/${compiler}-${mpi}/"


echo === Run unit tests ===

ln -s ${top_src_dir}/src.c++/test src.c++/test
make check


echo === Done: exit code $rc ===

cd "${top_src_dir}" \
    && rm -rf "$build_dir"

exit $rc

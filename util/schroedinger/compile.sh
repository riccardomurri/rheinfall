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
'run/REVNO/COMPILER-MPILIB'.

The COMPILER tag selects the compiler type (gcc412, gcc434, gcc443,
gcc450, intel; default: gcc443); second tag MPILIB selects which MPI
library will be used (ompi, parastation, mvapich, intel, or none;
default: openmpi).  Third argument REVNO defaults to the current
Bazaar repository revno, but can be overridden to tag variants of the
code.

EOF
}
case `hostname` in
        login*) usage; exit 1 ;;
        *) : ;;
esac

libdir=$(pwd)/util/schroedinger
source $libdir/functions.sh \
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
# if there's no field-separator, `cut` outputs the same string
# regardless of `-f`, so delete $extra here if it's the case
if [ "_r$revno" = "_$extra" ]; then
    extra=''
fi
if [ -n "$extra" ]; then
    if [  "$current_rev" = 'UNKNOWN' \
        -o "$revno" -ne "$(echo $current_rev | cut -c2-)" ]; 
    then
        die 1 "Cannot tag non-current revision: request for '$rev' but sources are at revision $current_rev"
    fi
fi

echo === Compiling Rheinfall r$revno ... ===

top_src_dir="$(pwd)"
#build_dir="/lustre/ESPFS/scratch/oci/murri/rheinfall.${flavor}"
build_dir="$HOME/data/tmp/rheinfall.${flavor}"

set -e
rm -rf "$build_dir"
# check out sources or compile from working directory
if [ -n "$extra" ]; then
    mkdir -p "$build_dir"
    src_dir="$top_src_dir"
else
    src_dir="$build_dir"
    if test -d .bzr; then
        bzr checkout -r $revno ./ "$build_dir"
    elif test -d .svn; then
        svn export -r $revno ./ "$build_dir"
    else
        die 1 "Cannot determine SCM used in directory '$(pwd)'"
    fi
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

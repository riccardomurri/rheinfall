#! /bin/sh
#$ -l s_cpu=1500
#$ -cwd
#$ -S /bin/bash
#$ -v PATH,LD_LIBRARY_PATH
#$ -j y
#$ -N compile

flavor="$1"; shift
compiler="$(echo $flavor | cut -d- -f1)"
if [ -z "$compiler" ]; then
    compiler=gcc443
fi
mpi="$(echo $flavor | cut -d- -f2)"
if [ -z "$mpi" ]; then
    mpi=ompi
fi
rev="$(echo $flavor | cut -d- -f3)"
if [ -z "$rev" ]; then
    rev="$( cat .bzr/branch/last-revision | (read revno rest; echo r$revno) )"
fi

# load modules
source /panfs/panfs0.ften.es.hpcn.uzh.ch/share/software/Modules/default/init/sh

# load MPI - must match what the binary was compiled with!
case "$mpi" in
    ompi|openmpi) 
        case "$compiler" in
            gcc450)    module load mpi/openmpi/gcc-4.5.0 ;;
            gcc*)      module load mpi/openmpi/gcc ;;
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
    none) # no MPI
        while type -a mpicc >/dev/null 2>&1; do
            p=$(dirname $(which mpicc) )
            echo 1>&2 "MPI compiler 'mpicc' found in PATH at '$p'; removing since no MPI was requested ..."
            export PATH=$(echo $PATH | sed -e "s|:${p}||;s|${p}:||;")
        done
        ;;
    *) 
        die 1 "Unknown MPI library '${mpi}' - please choose one of 'ompi', 'mvapich', 'intel' or 'none'"
        ;;
esac

# load the compiler and libs
module load binutils/2.20.1
case "$compiler" in
    gcc412)
        module load gcc/4.1.2
        CXX=g++
        toolset=gcc
        cxxflags='-march=nocona'
        ;;
    gcc441)
        module load gcc/4.4.1
        CXX=g++
        toolset=gcc
        cxxflags='-march=nocona'
        ;;
    gcc443)
        module load gcc/4.4.3
        CXX=g++
        toolset=gcc
        cxxflags='-march=native'
        ;;
    gcc450)
        module load gcc/4.5.0
        CXX=g++
        toolset=gcc
        cxxflags='-march=native'
        ;;
    icc|intel)
        module load intel/comp/11.1.064
        CXX=icc
        toolset=intel
        cxxflags='-O3 -xHOST'
        ;;
    *) 
        die 1 "Unknown compiler flavor '${compiler} - please choose one of 'gcc443', 'gcc450', 'icc'"
        ;;
esac


echo === running info ===
echo flavor: $flavor
echo compiler: $compiler
echo mpi: $mpi
echo node: $(hostname)
echo === compiler information ===

which ${CXX}
${CXX} --version

which mpicxx
echo

which as
as --version
echo

echo === compile log ===
# if [ "$mpi" != "none" ]; then
#     echo == Compiling Boost libraries ...
#     #./src/install.sh
#     rm -rf ~/sw/lib/libboost* #~/sw/src/boost_*/bin.v2/libs/mpi
#     pushd ~/sw/src/boost_*
#     ./bjam toolset=${toolset} link=static threading=multi install
#     popd
# fi

echo == Compiling Rheinfall ...

top_src_dir="$HOME/rheinfall"
sw="$HOME/sw/${compiler}-${mpi}"
build_dir="/lustre/ESPFS/scratch/oci/murri/rheinfall.${flavor}"

set -e
mkdir -p "$build_dir"
cd "$build_dir"
$top_src_dir/configure CXXFLAGS="-O3 -DNDEBUG $cxxflags" --with-boost=$sw --with-gmp=$sw
make

echo == Saving executables in "${top_src_dir}/run/${compiler}-${mpi}-${rev}/"
mkdir -pv "${top_src_dir}/run/${compiler}-${mpi}-${rev}/"
cp -av \
    $(find src.c++ src.c -type f -perm /u+x) \
    "${top_src_dir}/run/${compiler}-${mpi}-${rev}/"
cd "${top_src_dir}" \
    && rm -rf "$build_dir"

echo === Done. ===

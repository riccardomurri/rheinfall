#! /bin/sh
#$ -l h_cpu=40000 
#$ -cwd 
#$ -S /bin/bash
#$ -j y
#$ -N prereq


## Customization section

# comment any of the following lines if you already have
# the required version (or a newer one) installed in the
# standard PATH/LD_LIBRARY_PATH
#
#CYTHON=0.11.1
#PYTHON=2.5.4
#SWIG=1.3.35
LINBOX=1.1.7rc0
GMP=5.0.1
GIVARO=3.3.2
ATLAS=3.8.3
BOOST=1.43.0 # http://surfnet.dl.sourceforge.net/project/boost/boost/1.43.0/boost_1_43_0.tar.gz
TCMALLOC=1.6 # http://google-perftools.googlecode.com/files/google-perftools-1.6.tar.gz


## No customization should be necessary further down here

PROG="$(basename $0)"

usage () {
cat <<EOF
Usage: $PROG [COMPILER[-MPILIB]]

Download and install all required software for 'rheinfall' to
work into the "`pwd`/sw-COMPILER-MPILIB" directory.

First argument FLAVOR selects the compiler type (gcc443, gcc450, intel);
second argument MPILIB selects which MPI library will be used (ompi,
mvapich, intel, or none).

EOF
}


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


## parse command-line 
flavor="$1"; shift
compiler="$(echo $flavor | cut -d- -f1)"
if [ -z "$compiler" ]; then
    compiler=gcc443
fi
mpi="$(echo $flavor | cut -d- -f2)"
if [ -z "$mpi" ]; then
    mpi=ompi
fi

src_home="${3:-$HOME/sw/src}"
case "${src_home}" in
    /*) : ;;
    *) src_home="`pwd`/$src_home" ;;
esac


# load modules
source /panfs/panfs0.ften.es.hpcn.uzh.ch/share/software/Modules/default/init/sh

supported_compilers='gcc412 gcc441 gcc443 gcc450 icc'
supported_mpilibs='openmpi mvapich intel none'

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
    none) # no MPI
        while type -a mpicc >/dev/null 2>&1; do
            p=$(dirname $(which mpicc) )
            echo 1>&2 "MPI compiler 'mpicc' found in PATH at '$p'; removing since no MPI was requested ..."
            export PATH=$(echo $PATH | sed -e "s|:${p}||;s|${p}:||;")
        done
        ;;
    *) 
        die 1 "Unknown MPI library '${mpi}' - please choose one of: $supported_mpilibs"
        ;;
esac

# load the compiler and libs
module load binutils/2.20.1
case "$compiler" in
    gcc412)
        module load gcc/4.1.2
        export CC=gcc
        export CXX=g++
        cflags='-O3 -march=nocona'
        toolset=gcc
        ;;
    gcc441)
        module load gcc/4.4.1
        export CC=gcc
        export CXX=g++
        cflags='-O3 -march=nocona'
        toolset=gcc
        ;;
    gcc443)
        module load gcc/4.4.3
        export CC=gcc
        export CXX=g++
        cflags='-O3 -march=native'
        toolset=gcc
        ;;
    gcc450)
        module load gcc/4.5.0
        export CC=gcc
        export CXX=g++
        cflags='-O3 -march=native'
        toolset=gcc
        ;;
    icc|intel)
        module load intel/comp/11.1.064
        export CC=icc
        export CXX=icpc
        cflags='-O3 -xHOST'
        toolset=intel
        ;;
    *) 
        die 1 "Unknown compiler flavor '${compiler}' - please choose one of: $supported_compilers"
        ;;
esac


echo === running info ===
echo flavor: $flavor
echo compiler: $compiler
echo mpi: $mpi
echo node: $(hostname)
echo === compiler information ===

which ${CXX}
set -e # exit on error here, in case the compiler does not have enough licences
${CXX} --version
set +e

if [ "x${mpi}" != 'xnone' ]; then
    which mpicxx
    mpicxx --version
    echo
fi

which as
as --version
echo


echo === compile log ===

## main

require_command bzip2
require_command gzip
require_command make
require_command tar
#require_command wget

set -e

sw="$HOME/sw/${compiler}-${mpi}"
mkdir -p "$sw"
cd "$sw"
mkdir -p build


# paths
PATH=${sw}/bin:$PATH
LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH
PYTHONPATH=${sw}/lib/python

# misc
concurrent_make="make -j $(grep -c '^processor' /proc/cpuinfo)"

if grep -q '^cpu MHz' /proc/cpuinfo; then
    mhz=`grep '^cpu MHz' /proc/cpuinfo | head -1 | cut -d: -f2 | cut -d. -f1 | tr -d ' '`
else
    bogomips=`fgrep bogomips /proc/cpuinfo | head -1 | cut -d: -f2 | cut -d. -f1`
    mhz=`expr $bogomips / 2`
fi

case `uname -m` in
        i?86) bits=32 ;;
        x86_64) bits=64 ;;
        *) die 1 "Unknown architecture `uname -m`: is it 32-bit or 64-bit?" ;;
esac

# Python
if [ -n "${PYTHON}" ]; then
    _ Installing Python ${PYTHON}
    cd ${sw}/build/
    #wget -N http://www.python.org/ftp/python/${PYTHON}/Python-${PYTHON}.tar.bz2
    set -x
    mkdir -p Python-${PYTHON}
    cd Python-${PYTHON}
    ${src_home}/Python-${PYTHON}/configure --prefix=${sw} \
        CC=${CC} CXX=${CXX} CFLAGS="${cflags}"
    $concurrent_make
    make install
    set +x
fi # PYTHON

# SWIG
if [ -n "${SWIG}" ]; then
    _ Downloading SWIG ${SWIG}
    cd ${sw}/build/
    #wget -N http://mirror.switch.ch/ftp/ubuntu/pool/main/s/swig1.3/swig1.3_${SWIG}.orig.tar.gz
    set -x
    mkdir -p swig-${SWIG}
    cd swig-${SWIG}
    ${src_home}/swig-${SWIG}/configure --prefix=${sw} \
        --with-python=${sw}/bin \
        --without-allegrocl \
        --without-chicken \
        --without-clisp \
        --without-csharp \
        --without-gcj \
        --without-guile \
        --without-java \
        --without-lua \
        --without-mzscheme \
        --without-ocaml \
        --without-octave \
        --without-perl5 \
        --without-php4 \
        --without-pike \
        --without-r \
        --without-ruby \
        --without-rxspencer \
        --without-tcl \
         CC=${CC} CXX=${CXX} CFLAGS="${cflags}";
    $concurrent_make
    make install
    set +x
fi # SWIG

# Cython
if [ -n "$CYTHON" ]; then
    _ Installing Cython-${CYTHON}
    cd ${sw}/build/
    #wget -N http://www.cython.org/Cython-${CYTHON}.tar.gz
    set -x
    tar -xzf ${src_home}/Cython-${CYTHON}.tar.gz
    cd Cython-${CYTHON}
    python setup.py build
    python setup.py install --home=${sw}
    
    PYTHONPATH=$PYTHONPATH:`pwd`; export PYTHONPATH
    PATH=$PATH:`pwd`/bin; export PATH
    set +x
fi # CYTHON

# GMP
if [ -n "$GMP" ]; then
    _ Installing GMP ${GMP}
    cd ${sw}/build
    #wget -N ftp://sunsite.cnlab-switch.ch/mirror/gnu/gmp/gmp-${GMP}.tar.bz2
    set -x
    mkdir -p gmp-${GMP}
    cd gmp-${GMP}
    ${src_home}/gmp-${GMP}/configure --prefix=${sw} \
        --enable-cxx \
        CC=${CC} CXX=${CXX} CFLAGS="${cflags}"
    $concurrent_make
    make install
    set +x
fi # GMP

# GIVARO (cfr. http://groups.google.com/group/linbox-use/browse_thread/thread/82673844f6921271)
if [ -n "$GIVARO" ]; then
    _ Installing Givaro ${GIVARO}
    cd ${sw}/build
    #wget -N http://www-lmc.imag.fr/CASYS/LOGICIELS/givaro/Downloads/givaro-${GIVARO}.tar.gz
    tar -xzf ${src_home}/givaro-${GIVARO}.tar.gz
    givaro_vers=${GIVARO%%rc[0-9]}
    set -x
    #mkdir -p givaro-$givaro_vers
    cd givaro-$givaro_vers
    # work around bug in ./configure: the test for GMP cannot find it
    # unless it's in the LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH
    # a bug in GIVARO's `configure` script prevents out-of-src compilation:
    #${src_home}/givaro-${givaro_vers}/configure --prefix=${sw} \
    ./configure --prefix=${sw} \
        --enable-shared ${GMP:+"--with-gmp=${sw}"} \
        CC=${CC} CXX=${CXX} CFLAGS="${cflags}"
    $concurrent_make
    make install
    set +x
fi # GIVARO

# ATLAS
if [ -n "$ATLAS" ]; then
    cd ${sw}/build
    #wget "http://switch.dl.sourceforge.net/sourceforge/math-atlas/atlas${ATLAS}.tar.bz2" \
    #     -O atlas${ATLAS}.tar.bz2
    set -x
    tar -xjf ${src_home}/atlas${ATLAS}.tar.bz2
    cd ATLAS
    mkdir -p BLDdir
    cd BLDdir
    # ATLAS documentation advices against changing compilers or compilation flags
    # (see: http://math-atlas.sourceforge.net/atlas_install/atlas_install.html#SECTION00042000000000000000)
    ../configure -v 2 \
        -b ${bits} \
        -m ${mhz} \
        -D c -DPentiumCPS=${mhz} \
        -Si cputhrchk 0 \
        -Ss pmake "$concurrent_make" \
        --prefix=${sw}
    make build
    (cd lib; make cshared cptshared && cp -a *.so ${sw}/lib)
    #make check
    #make time
    make install
    set +x
fi # ATLAS

# LinBox
if [ -n "$LINBOX" ]; then
    _ Installing LinBox ${LINBOX}
    cd ${sw}/build
    #wget -N http://linalg.org/linbox-${LINBOX}.tar.gz
    set -x
    # heck, neither LinBox is not capable of out-of-srcdir compilation...
    #mkdir -p linbox-${LINBOX}
    tar -xzf ${src_home}/linbox-${LINBOX}.tar.gz
    cd linbox-${LINBOX}
    ./configure --prefix=${sw} \
        ${ATLAS:+"--with-blas=${sw}"} \
        ${GMP:+"--with-gmp=${sw}"} \
        ${GIVARO:+"--with-givaro=${sw}"} \
        CC=${CC} CXX=${CXX} CFLAGS="${cflags}";
    $concurrent_make
    make install
    set +x
fi # LINBOX

# Boost
if [ -n "$BOOST" ]; then
    _ Installing Boost ${BOOST}
    cd ${sw}/build
    boost_file=$(echo boost_$BOOST | tr . _)
    #wget "http://surfnet.dl.sourceforge.net/project/boost/boost/${BOOST}/${boost_file}.tar.gz" \
    #     -O "${boost_file}.tar.gz"
    set -x 
    tar -xzf  "${src_home}/${boost_file}.tar.gz"
    cd ${boost_file}
    ./bootstrap.sh --prefix=${sw} --with-libraries=mpi,serialization,test link=static threading=multi
    cat >> project-config.jam <<EOF
# Boost will not build Boost.MPI unless it is explicitly 
# told to by the following line:
using mpi : mpicxx ;
EOF
    # first, need to build `bjam`
    (cd tools/jam/src; sh ./build.sh)
    # then, build Boost with the new `bjam`
    PATH=$(pwd)/tools/jam/src/bin.$(uname -s | tr A-Z a-z)$(uname -m):$PATH
    export PATH
    bjam --prefix=${sw} toolset=${toolset} link=static threading=multi install
    set +x
fi # BOOST


# Google perftools
if [ -n "$TCMALLOC" ]; then
    _ Installing Google PerfTools $TCMALLOC ...
    cd ${sw}/build
    set -x
    #wget -N http://google-perftools.googlecode.com/files/google-perftools-${TCMALLOC}.tar.gz
    mkdir -p google-perftools-${TCMALLOC}
    cd google-perftools-${TCMALLOC}
    ( 
      # Google perftools error out when compiling with ICC, so revert to GCC 4.1 in this case
      if test "x${CC}" = "xicc"; then
          module load gcc/4.1.2
          export CC=gcc
          export CXX=g++
          cflags='-O3 -march=nocona'
      fi
      ${src_home}/google-perftools-${TCMALLOC}/configure --prefix=${sw} \
          --enable-frame-pointers \
          CC=${CC} CXX=${CXX} CFLAGS="${cflags}" CXXFLAGS="${cflags}"
    )
    $concurrent_make
    $concurrent_make install
    set +x
fi


_ All done.
#! /bin/bash
#$ -l s_rt=40000 
#$ -cwd 
#$ -S /bin/bash
#$ -j y
#$ -N prereq


## Customization section

# comment any of the following lines if you already have
# the required version (or a newer one) installed in the
# standard PATH/LD_LIBRARY_PATH
#
#GMP=5.0.1
#GIVARO=3.3.2
#ATLAS=3.8.3
#LINBOX=1.1.7rc0
BOOST=1.43.0 # http://surfnet.dl.sourceforge.net/project/boost/boost/1.43.0/boost_1_43_0.tar.gz
#TCMALLOC=1.6 # http://google-perftools.googlecode.com/files/google-perftools-1.6.tar.gz


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
case `hostname` in
        login*) usage; exit 1 ;;
        *) : ;;
esac

source $HOME/rheinfall/util/schroedinger/functions.sh \
    || { echo 1>&2 "Cannot load 'functions.sh' - aborting."; exit 1; }


## parse command-line 
set_mpi_and_compiler_flavor "$@"; shift

src_home="${2:-$HOME/sw/src}"
case "${src_home}" in
    /*) : ;;
    *) src_home="`pwd`/$src_home" ;;
esac


## run compilation

show_build_information

echo === compile log ===

## main

require_command bzip2
require_command gzip
require_command make
require_command tar
#require_command wget

set -e

# paths
mkdir -p "$sw"
cd "$sw"
mkdir -p build
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
    set -x
    #wget -N http://www.python.org/ftp/python/${PYTHON}/Python-${PYTHON}.tar.bz2
    mkdir -p Python-${PYTHON}
    cd Python-${PYTHON}
    ${src_home}/Python-${PYTHON}/configure --prefix=${sw} \
        CC=${CC} CXX=${CXX} CFLAGS="${cflags} ${std_cflags}"
    $concurrent_make
    make install
    set +x
fi # PYTHON

# SWIG
if [ -n "${SWIG}" ]; then
    _ Downloading SWIG ${SWIG}
    cd ${sw}/build/
    set -x
    #wget -N http://mirror.switch.ch/ftp/ubuntu/pool/main/s/swig1.3/swig1.3_${SWIG}.orig.tar.gz
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
         CC=${CC} CXX=${CXX} CFLAGS="${cflags} ${std_cflags}";
    $concurrent_make
    make install
    set +x
fi # SWIG

# Cython
if [ -n "$CYTHON" ]; then
    _ Installing Cython-${CYTHON}
    cd ${sw}/build/
    set -x
    #wget -N http://www.cython.org/Cython-${CYTHON}.tar.gz
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
    rm -rf gmp-${GMP}
    mkdir -p gmp-${GMP}
    set -x
    #wget -N ftp://sunsite.cnlab-switch.ch/mirror/gnu/gmp/gmp-${GMP}.tar.bz2    
    cd gmp-${GMP}
    ${src_home}/gmp-${GMP}/configure --prefix=${sw} \
        --enable-cxx \
        CC=${CC} CFLAGS="${cflags} ${std_cflags}" \
        CXX=${CXX} CXXFLAGS="${cxxflags} ${std_cxxflags}"
    rm -rf ${sw}/include/gmp* ${sw}/lib/libgmp* # avoid `libtool` errors
    $concurrent_make
    make install
    set +x
fi # GMP

# GIVARO (cfr. http://groups.google.com/group/linbox-use/browse_thread/thread/82673844f6921271)
if [ -n "$GIVARO" ]; then
    _ Installing Givaro ${GIVARO}
    givaro_vers=${GIVARO%%rc[0-9]}
    cd ${sw}/build
    rm -rf givaro-${givaro_vers}
    #wget -N http://www-lmc.imag.fr/CASYS/LOGICIELS/givaro/Downloads/givaro-${GIVARO}.tar.gz
    set -x
    tar -xzf ${src_home}/givaro-${GIVARO}.tar.gz
    cd givaro-$givaro_vers
    # work around bug in ./configure: the test for GMP cannot find it
    # unless it's in the LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH
    # a bug in GIVARO's `configure` script prevents out-of-src compilation:
    #${src_home}/givaro-${givaro_vers}/configure --prefix=${sw} \
    ./configure --prefix=${sw} \
        --enable-shared ${GMP:+"--with-gmp=${sw}"} \
        CC=${CC} CFLAGS="${cflags} ${std_cflags}" \
        CXX=${CXX} CXXFLAGS="${cxxflags} ${std_cxxflags}"
    rm -rf ${sw}/include/givaro* ${sw}/lib/libgivaro* # avoid `libtool` errors
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
    tar -xjf ${src_home}/atlas${ATLAS}.tar.bz2 --keep-newer-files
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
    rm -rf linbox-${LINBOX}
    #wget -N http://linalg.org/linbox-${LINBOX}.tar.gz
    set -x
    # heck, neither LinBox is not capable of out-of-srcdir compilation...
    tar -xzf ${src_home}/linbox-${LINBOX}.tar.gz
    mkdir -p linbox-${LINBOX}
    cd linbox-${LINBOX}
    ./configure --prefix=${sw} \
        --with-blas=${sw} \
        --with-gmp=${sw} \
        --with-givaro=${sw} \
        CC=${CC} CFLAGS="${cflags} ${std_cflags}" \
        CXX=${CXX} CXXFLAGS="${cxxflags} ${std_cxxflags}";
    rm -rf ${sw}/include/linbox* ${sw}/lib/liblinbox* # avoid `libtool` errors
    $concurrent_make
    make install
    set +x
fi # LINBOX

# Boost
if [ -n "$BOOST" ]; then
    _ Installing Boost ${BOOST}
    # clean up old installation
    rm -rf ${sw}/include/boost ${sw}/lib/libboost_*
    # now build new one
    cd ${sw}/build
    boost_file=$(echo boost_$BOOST | tr . _)
    #wget "http://surfnet.dl.sourceforge.net/project/boost/boost/${BOOST}/${boost_file}.tar.gz" \
    #     -O "${boost_file}.tar.gz"
    set -x 
    tar -xzf  "${src_home}/${boost_file}.tar.gz"
    cd ${boost_file}
    # build Boost.MPI for homogeneous clusters (same arch, so avoid pack/unpack)
    #sed -e 's|^//#define BOOST_MPI_HOMOGENEOUS|#define BOOST_MPI_HOMOGENEOUS|' \
    #    -i boost/mpi/config.hpp
    ./bootstrap.sh --prefix=${sw} --with-libraries=mpi,serialization,test \
        toolset=${toolset} variant=release threading=multi
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
    bjam --prefix=${sw} toolset=${toolset} threading=multi variant=release install
    set +x
fi # BOOST


# Google perftools
if [ -n "$TCMALLOC" ]; then
    _ Installing Google PerfTools $TCMALLOC ...
    cd ${sw}/build
    rm -rf google-perftools-${TCMALLOC}
    mkdir -p google-perftools-${TCMALLOC}
    set -x
    #wget -N http://google-perftools.googlecode.com/files/google-perftools-${TCMALLOC}.tar.gz
    cd google-perftools-${TCMALLOC}
    if [ "_$CC" = '_icc' ]; then
        # need to downgrade CFLAGS on icc -- "-fast" breaks the build
        cflags='-xHOST -O3'
        cxxflags='-xHOST -O3'
    fi
    ${src_home}/google-perftools-${TCMALLOC}/configure --prefix=${sw} \
        --enable-frame-pointers --disable-debugalloc \
        CC=${CC} CFLAGS="${cflags} ${std_cflags}" \
        CXX=${CXX} CXXFLAGS="${cxxflags} ${std_cxxflags}"
    $concurrent_make clean all
    make install
    set +x
fi


_ All done.

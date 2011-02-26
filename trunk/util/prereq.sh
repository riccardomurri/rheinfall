#! /bin/sh
#
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
#TBB=20100915oss # http://www.threadingbuildingblocks.org/uploads/77/161/3.0%20update%203/tbb30_20100915oss_lin.tgz

# undefine when running on a non-homogeneous cluster
# (i.e., nodes may have different arch and/or memory size)
BOOST_MPI_HOMOGENEOUS=yes


## No customization should be necessary further down here

PROG="$(basename $0)"

usage () {
cat <<EOF
Usage: $PROG

Download and install all required software for 'mgn.sh' to
work into the "`pwd`/sw" directory.

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

if [ $# -ne 0 ]; then
  usage
  exit 0
fi


## main

require_command bzip2
require_command gzip
require_command make
require_command tar
require_command wget

set -e

mkdir -p sw/src
cd sw

# paths
sw=`pwd`
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
    cd ${sw}/src/
    wget -N http://www.python.org/ftp/python/${PYTHON}/Python-${PYTHON}.tar.bz2
    set -x
    tar -xjf Python-${PYTHON}.tar.bz2
    cd Python-${PYTHON}
    ./configure --prefix=${sw}
    $concurrent_make
    make install
    set +x
fi # PYTHON

# SWIG
if [ -n "${SWIG}" ]; then
    _ Downloading SWIG ${SWIG}
    cd ${sw}/src/
    wget -N http://mirror.switch.ch/ftp/ubuntu/pool/main/s/swig1.3/swig1.3_${SWIG}.orig.tar.gz
    set -x
    tar -xzf swig1.3_${SWIG}.orig.tar.gz
    cd swig-${SWIG}
    ./configure --prefix=${sw} \
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
        ;
    $concurrent_make
    make install
    set +x
fi # SWIG

# Cython
if [ -n "$CYTHON" ]; then
    _ Installing Cython-${CYTHON}
    cd ${sw}/src/
    wget -N http://www.cython.org/Cython-${CYTHON}.tar.gz
    set -x
    tar -xzf Cython-${CYTHON}.tar.gz
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
    cd ${sw}/src
    wget -N ftp://sunsite.cnlab-switch.ch/mirror/gnu/gmp/gmp-${GMP}.tar.bz2
    set -x
    tar -xjf gmp-${GMP}.tar.bz2
    cd gmp-${GMP}
    ./configure --prefix=${sw} --enable-cxx
    $concurrent_make
    make install
    set +x
fi # GMP

# GIVARO (cfr. http://groups.google.com/group/linbox-use/browse_thread/thread/82673844f6921271)
if [ -n "$GIVARO" ]; then
    _ Installing Givaro ${GIVARO}
    cd ${sw}/src
    wget -N http://www-lmc.imag.fr/CASYS/LOGICIELS/givaro/Downloads/givaro-${GIVARO}.tar.gz
    set -x
    tar -xzf givaro-${GIVARO}.tar.gz
    cd givaro-${GIVARO%%rc[0-9]}
    # work around bug in ./configure: the test for GMP cannot find it
    # unless it's in the LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH
    ./configure  --prefix=${sw} --enable-shared ${GMP:+"--with-gmp=${sw}"}
    $concurrent_make
    make install
    set +x
fi # GIVARO

# ATLAS
if [ -n "$ATLAS" ]; then
    cd ${sw}/src
    wget "http://switch.dl.sourceforge.net/sourceforge/math-atlas/atlas${ATLAS}.tar.bz2" \
         -O atlas${ATLAS}.tar.bz2
    set -x
    tar -xjf atlas${ATLAS}.tar.bz2
    cd ATLAS
    mkdir -p BLDdir
    cd BLDdir
    ../configure -v 2 -b ${bits} -m ${mhz} -D c -DPentiumCPS=${mhz} -Si cputhrchk 0 --prefix=${sw}
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
    cd ${sw}/src
    wget -N http://linalg.org/linbox-${LINBOX}.tar.gz
    set -x
    tar -xzf linbox-${LINBOX}.tar.gz
    cd linbox-${LINBOX}
    ./configure --prefix=${sw} \
        ${ATLAS:+"--with-blas=${sw}/lib"} \
        ${GMP:+"--with-gmp=${sw}"} \
        ${GIVARO:+"--with-givaro=${sw}"} \
        ;
    $concurrent_make
    make install
    set +x
fi # LINBOX

# Boost
if [ -n "$BOOST" ]; then
    _ Installing Boost ${BOOST}
    cd ${sw}/src
    boost_file=$(echo boost_$BOOST | tr . _)
    wget "http://surfnet.dl.sourceforge.net/project/boost/boost/${BOOST}/${boost_file}.tar.gz" \
        -O "${boost_file}.tar.gz"
    set -x 
    tar -xzf  "${boost_file}.tar.gz"
    cd ${boost_file}
    # build Boost.MPI for homogeneous clusters (same arch, so avoid pack/unpack)
    if [ "x$BOOST_MPI_HOMOGENEOUS" = "xyes" ]; then
        sed -e 's|^//#define BOOST_MPI_HOMOGENEOUS|#define BOOST_MPI_HOMOGENEOUS|' \
            -i boost/mpi/config.hpp
    fi
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
    ./bjam --prefix=${sw} link=static threading=multi install
    set +x
fi # BOOST


# Google perftools
if [ -n "$TCMALLOC" ]; then
    _ Installing Google PerfTools $TCMALLOC ...
    cd ${sw}/src
    set -x
    wget -N http://google-perftools.googlecode.com/files/google-perftools-${TCMALLOC}.tar.gz
    tar -xzf "google-perftools-${TCMALLOC}.tar.gz"
    cd google-perftools-${TCMALLOC}
    ./configure --prefix=${sw} --enable-frame-pointers
    $concurrent_make
    $concurrent_make install
    set +x
fi


# Intel TBB
if [ -n "$TBB" ]; then
    _ Installing Intel TBB $TBB ...
    cd ${sw}/src/
    set -x 
    case "$TBB" in
        20100915oss) # TBB download URL changes with release ...
            wget -N "http://www.threadingbuildingblocks.org/uploads/77/161/3.0%20update%203/tbb30_20100915oss_lin.tgz"
            wget -N "http://www.threadingbuildingblocks.org/uploads/77/161/3.0%20update%203/tbb30_20100915oss_src.tgz"
            ;;
        *)
            die 1 "Unknown download URL for Intel TBB ${TBB}"
            ;;
    esac
    mkdir -p ${sw}/opt
    for tgz in tbb*.tgz; do
        tar -xzf "$tgz" -C "${sw}/opt"
    done
    set +x
fi


_ All done.

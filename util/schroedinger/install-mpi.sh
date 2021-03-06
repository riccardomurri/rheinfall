#! /bin/bash
#$ -l h_cpu=1800 
#$ -cwd 
#$ -S /bin/bash
#$ -j y
#$ -N install-mpi


## Customization section

# comment any of the following lines if you already have
# the required version (or a newer one) installed in the
# standard PATH/LD_LIBRARY_PATH
#
OPENMPI=1.4.3   # http://www.open-mpi.org/software/ompi/v1.4/downloads/openmpi-1.4.3.tar.bz2
MVAPICH2=1.5    # svn co https://mvapich.cse.ohio-state.edu/svn/mpi/mvapich2/branches/1.5/ mvapich2 
MPIP=3.3        # wget http://switch.dl.sourceforge.net/project/mpip/mpiP/mpiP-${MPIP}/mpiP-${MPIP}.tar.gz
LIBUNWIND=1.0.1 # wget http://download.savannah.gnu.org/releases/libunwind/libunwind-${LIBUNWIND}.tar.gz

## No customization should be necessary further down here

PROG="$(basename $0)"

usage () {
cat <<EOF
Usage: $PROG COMPILER-MPILIB

Download and install OpenMPI or MVAPICH software into the
"$HOME/sw/COMPILER-MPILIB" directory.

First argument COMPILER selects the compiler type (gcc443, gcc450,
intel); second argument MPILIB must be one of 'ompi' or 'mvapich'.

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


case "$mpi" in 

    ompi*)
        if [ -z "${OPENMPI}" ]; then
            die 1 "Missing version of OpenMPI to install."
        fi
        (
            _ Installing OpenMPI ${OPENMPI}
            sw="$HOME/sw/${compiler}-${mpi}"
            mkdir -p "$sw"
            cd "$sw"
            mkdir -p build
            cd ${sw}/build/
            PATH=${sw}/bin:$PATH
            LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH
            
            rm -rf openmpi-${OPENMPI}
            #wget -N http://www.open-mpi.org/software/ompi/v1.4/downloads/openmpi-${OPENMPI}.tar.bz2
            mkdir -p openmpi-${OPENMPI}
            cd openmpi-${OPENMPI}
            
            set -x
            # see: /panfs/panfs0.ften.es.hpcn.uzh.ch/share/software/mpi/src/install.info
            ${src_home}/openmpi-${OPENMPI}/configure --prefix=${sw} \
                --disable-ipv6 \
                --disable-mem-debug \
                --disable-mem-profile \
                --disable-mpi-f77 \
                --disable-mpi-f90 \
                --disable-mpi-threads \
                --enable-branch-probabilities \
                --enable-cxx-exceptions \
                --enable-mpirun-prefix-by-default \
                --enable-openib-ibcm \
                --enable-openib-rdmacm \
                --enable-per-user-config-files \
                --enable-peruse \
                --with-sge \
                --without-mpi-param-check \
                CC=${CC} CFLAGS="${cflags} ${std_cflags}" \
                CXX=${CXX} CXXFLAGS="${cxxflags} ${std_cxxflags}";
            $concurrent_make
            make install
            set +x
        )
        ;; # OpenMPI

    mvapich*)
        if [ -z "${MVAPICH2}" ]; then
            die 1 "Missing version of MVAPICH2 to install."
        fi
        (
            _ Installing MVAPICH2 ${MVAPICH2}
            sw="$HOME/sw/${compiler}-${mpi}"
            mkdir -p "$sw"
            cd "$sw"
            mkdir -p build
            cd ${sw}/build/
            PATH=${sw}/bin:$PATH
            LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH
            
            rm -rf mvapich2-${MVAPICH2}
            # svn co https://mvapich.cse.ohio-state.edu/svn/mpi/mvapich2/branches/${MVAPICH2} mvapich2-${MVAPICH2}
            mkdir -p mvapich2-${MVAPICH2}
            cd mvapich2-${MVAPICH2}
            
            set -x
            ${src_home}/mvapich2-${MVAPICH2}/configure --prefix=${sw} \
                --enable-error-checking=runtime \
                --enable-error-messages=all \
                --enable-g=none \
                --enable-fast=ndebug,O3 \
                --enable-cxx \
                --enable-threads=default \
                --with-pm=remshell \
                --with-thread-package=posix \
                --disable-f77 \
                --disable-f90 \
                CC=${CC} CFLAGS="${cflags} ${std_cflags}" \
                CXX=${CXX} CXXFLAGS="${cxxflags} ${std_cxxflags}";
            $concurrent_make
            make install
            set +x
        )
        ;; # MVAPICH2
esac

# libunwind
(   _ Installing libunwind ${LIBUNWIND} ...

    sw="$HOME/sw/${compiler}-${mpi}"
    mkdir -p "$sw"
    cd "$sw"
    mkdir -p build
    cd ${sw}/build/
    PATH=${sw}/bin:$PATH
    LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH
            
    rm -rf libunwind-${LIBUNWIND}
    #wget http://download.savannah.gnu.org/releases/libunwind/libunwind-${LIBUNWIND}.tar.gz
    mkdir -p libunwind-${LIBUNWIND}
    cd libunwind-${LIBUNWIND}
            
    set -x
    env \
        CC=${CC} CFLAGS="${cflags} ${std_cflags}" \
        CXX=${CXX} CXXFLAGS="${cxxflags} ${std_cxxflags}" \
        ${src_home}/libunwind-${LIBUNWIND}/configure \
          --prefix=${sw} --enable-cxx-exceptions
    $concurrent_make
    make install
    set +x
)     

# mpiP
(   _ Installing mpiP ${MPIP}

    sw="$HOME/sw/${compiler}-${mpi}"
    mkdir -p "$sw"
    cd "$sw"
    mkdir -p build
    cd ${sw}/build/
    PATH=${sw}/bin:$PATH
    LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH
            
    rm -rf mpiP-${MPIP}
    #wget http://switch.dl.sourceforge.net/project/mpip/mpiP/mpiP-${MPIP}/mpiP-${MPIP}.tar.gz
    mkdir -p mpiP-${MPIP}
    cd mpiP-${MPIP}
            
    set -x
    # mpiP's ./configure is too old and cannot cope with CFLAGS=... given on the command-line
    env \
        CC=${CC} CFLAGS="${cflags} ${std_cflags}" \
        CXX=${CXX} CXXFLAGS="${cxxflags} ${std_cxxflags}" \
        ${src_home}/mpiP-${MPIP}/configure --prefix=${sw} \
          --disable-mpi-io \
          --enable-demangling=GNU \
          --enable-collective-report-default \
          --with-cc=${CC} --with-cxx=${CXX} 
    $concurrent_make
    make install
    set +x
)

echo === All done. ===
exit $rc

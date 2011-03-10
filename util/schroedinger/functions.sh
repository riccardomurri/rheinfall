#! /bin/sh
#
# This is a collection of functions used to carry out common tasks 
# in Schroedinger helper scripts.  It won't run as a stand-alone script.
#

if test -e $HOME/.bashrc; then
  source $HOME/.bashrc
fi

set_mpi_and_compiler_flavor () {
    flavor="$1"; shift

    if test -d .bzr; then
        current_rev=$( cat .bzr/branch/last-revision | (read revno rest; echo r$revno) );
    elif test -d .svn; then
        current_rev=r$( env LC_ALL=C svn info | grep '^Revision:' | cut -d' ' -f2 );
    elif bzr revno >/dev/null; then
        current_rev=r$( bzr revno );
    else
        current_rev=UNKNOWN
    fi

    rev="$(echo $flavor | cut -d/ -f1)"
    compiler_and_mpilib="$(echo $flavor | cut -d/ -f2)"
    if [ -z "$compiler_and_mpilib" -o "x$compiler_and_mpilib" = "x$rev" ]; then
        # rev can be omitted - try to get it from the BZR/SVN repository
        compiler_and_mpilib="$rev"
        rev="$current_rev"
    fi

    compiler="$(echo $compiler_and_mpilib | (IFS=- read one two; echo $one) )"
    if [ -z "$compiler" ]; then
        compiler=gcc443
    fi
    mpi="$(echo $compiler_and_mpilib |  (IFS=- read one two; echo $two) )"
    if [ -z "$mpi" ]; then
        mpi=ompi
    fi

    # rebuild flavor as a correct filename
    flavor="${rev}-${compiler}-${mpi}"

    ## load modules
    source /panfs/panfs0.ften.es.hpcn.uzh.ch/share/software/Modules/default/init/sh

    supported_compilers='gcc412 gcc434 gcc441 gcc443 gcc450 icc'
    supported_mpilibs='openmpi parastation parastation-mt mvapich intel none'

     # load MPI - must match what the binary was compiled with!
    case "$mpi" in
        ompi|openmpi) # systemwide OpenMPI
            pe=openmpi2
            case "$compiler" in
                gcc450)    module load mpi/openmpi/gcc-4.5.0 ;;
                gcc*)      module load mpi/openmpi/gcc ;;
                icc|intel) module load mpi/openmpi/intel ;;
            esac
            ;;
        ompi*) # My own OpenMPI install in ${sw}/lib etc.
            pe=openmpi2
            ;;
        para*) # systemwide Parastation
            pe=parastation
            export PATH=/opt/parastation/mpi2/bin:/opt/parastation/bin:$PATH
            export LD_LIBRARY_PATH=/opt/parastation/mpi2/lib:$LD_LIBRARY_PATH
            export LD_RUN_PATH=$LD_LIBRARY_PATH
            # apparently required for PSMPI to run,
            # see ChriBo's email to hpcnlist on 2011-02-09
            export PBS_NODEFILE=$TMPDIR/machines
            ;;
        para*mt) # systemwide Parastation w/ threads support
            pe=parastation
            export PATH=/opt/parastation/mpi2-mt/bin:/opt/parastation/bin:$PATH
            export LD_LIBRARY_PATH=/opt/parastation/mpi2-mt/lib:$LD_LIBRARY_PATH
            export LD_RUN_PATH=$LD_LIBRARY_PATH
            # apparently required for PSMPI to run,
            # see ChriBo's email to hpcnlist on 2011-02-09
            export PBS_NODEFILE=$TMPDIR/machines
            ;;
        mvapich) # MVAPICH is not supported by modules, apparently
            pe=mpich_rsh
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
        mvapich*) # My own MVAPICH2 install in ${sw}/lib etc.
            pe=mpich_rsh
            ;;
        impi|intel) # Intel MPI is not supported by modules, apparently
            pe=openmpi2
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
            cflags='-march=nocona'
            cxxflags='-march=nocona'
            std_cflags='-O3'
            std_cxxflags='-O3'
            toolset=gcc
            ;;
        gcc434)
            # GCC 4.3.4 is std on Schroedinger after the 2011-01-29 upgrade
            export CC=gcc-4.3
            export CXX=g++-4.3
            cflags='-march=native'
            cxxflags='-march=native'
            std_cflags='-O3'
            std_cxxflags='-O3'
            toolset=gcc
            ;;
        gcc441)
            module load gcc/4.4.1
            export CC=gcc-4.4.1
            export CXX=g++-4.4.1
            cflags='-march=nocona'
            cxxflags='-march=nocona'
            std_cflags='-O3'
            std_cxxflags='-O3'
            toolset=gcc
            ;;
        gcc443)
            module load gcc/4.4.3
            export CC=/panfs/panfs0.ften.es.hpcn.uzh.ch/share/software/Compilers/gcc-4.4.3-build/bin/gcc
            export CXX=/panfs/panfs0.ften.es.hpcn.uzh.ch/share/software/Compilers/gcc-4.4.3-build/bin/g++
            # however, we need to tell gcc that it has to use the 4.4.3
            # libstdc++, otherwise it will try to use the system default one
            export LD_RUN_PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/share/software/Compilers/gcc-4.4.3-build/lib:$LD_RUN_PATH
            export LD_LIBRARY_PATH=/panfs/panfs0.ften.es.hpcn.uzh.ch/share/software/Compilers/gcc-4.4.3-build/lib:$LD_LIBRARY_PATH
            cflags='-march=native'
            cxxflags='-march=native'
            std_cflags='-O3'
            std_cxxflags='-O3'
            toolset=gcc
            ;;
        gcc450)
            # GCC 4.5 provided by SLES11 package after 2011-02-29 upgrade
            #module load gcc/4.5.0
            export CC=gcc-4.5
            export CXX=g++-4.5
            # however, we need to tell gcc that it has to use the 4.5
            # libstdc++, otherwise it will try to use the 4.3 one
            export LD_RUN_PATH=/usr/lib64/gcc/x86_64-suse-linux/4.5:$LD_RUN_PATH
            export LD_LIBRARY_PATH=/usr/lib64/gcc/x86_64-suse-linux/4.5:$LD_LIBRARY_PATH
            cflags='-march=native'
            cxxflags='-march=native'
            std_cflags='-O3'
            std_cxxflags='-O3'
            toolset=gcc
            ;;
        icc|intel)
            module load intel/comp/11.1.064
            export CC=icc
            export CXX=icpc
            # c(xx)flags are the equivalent of `-fast` w/out the `-static`
            cflags='-xHOST -O3 -ipo -no-prec-div'
            cxxflags='-xHOST -O3 -ipo -no-prec-div'
            std_cflags='-O3'
            std_cxxflags='-O3'
            toolset=intel
            ;;
        *) 
            die 1 "Unknown compiler flavor '${compiler}' - please choose one of: $supported_compilers"
            ;;
    esac
    # MPI wrappers should use the default compiler
    export I_MPI_CC=`which $CC`
    export I_MPI_CXX=`which $CXX`
    export MPICH_CC=`which $CC`
    export MPICH_CXX=`which $CXX`
    export OMPI_CC=`which $CC`
    export OMPI_CXX=`which $CXX`

    # enable local sw
    sw="$HOME/sw/${compiler}-${mpi}"
    PATH=${sw}/bin:$PATH; export PATH
    LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
}


show_build_information () {
    echo === running info ===
    echo flavor: $flavor
    echo revno: $rev
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
}


prepare_openmpi_environment () {
    require_environment_var PE_HOSTFILE
    require_environment_var TMPDIR
    
    cat $PE_HOSTFILE | \
        (while read hostname nslots queue rest; 
            do echo "$hostname slots=$nslots"; done) > $TMPDIR/hostfile

    set -e
    eval `ssh-agent -s`
    ssh-add $HOME/.ssh/id_dsa
    set +e
   
    gmca="--gmca plm_rsh_agent $HOME/bin/qrsh.sh"
    hostfile="--hostfile $TMPDIR/hostfile"
}


## generic functions
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

require_environment_var () {
    if [ -z "$(eval echo '$'$1)" ]; then
        die 1 "Require environment variable '$1' is not defined or has empty value. Aborting."
    fi
}

is_absolute_path () {
    expr match "$1" '/' >/dev/null 2>/dev/null
}

_ () {
    echo
    echo ==== "$@" ...;
}


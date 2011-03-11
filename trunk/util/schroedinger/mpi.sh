#! /bin/bash

flavor="$1"; shift
if expr match "$flavor" 'run/' >/dev/null; then
    flavor=$(echo $flavor | cut -c5-)
fi

bindir="run/$flavor"
if test ! -d "$bindir"; then
    echo 1>&2 "Cannot read directory: '$bindir'"
    exit 1
fi

libdir=$(pwd)/util/schroedinger
source $libdir/functions.sh \
    || { echo 1>&2 "Cannot load 'functions.sh' - aborting."; exit 1; }
alias module=true # work around lack of "modules" on idesl4-g
set_mpi_and_compiler_flavor "$flavor"

flavor="${rev}-${compiler}-${mpi}"

nps='16 32 64' # 128 256 512'
ws='4096 1024 256 64 16' # 4 1'

for np in $nps; do
for w in $ws; do
for prog in \
    rank-int-mpi \
    rank-mod-mpi \
    rank-int-mpi-omp \
    rank-mod-mpi-omp \
    ; 
do
    if test ! -x "${bindir}/${prog}"; then
        echo 1>&2 "Skipping executable '$prog': not found in '$bindir'"
        continue
    fi
    qsub \
        -cwd \
        -S /bin/bash \
        -N "mpi-n${np}-w${w}.${prog}.${flavor}" \
        -l s_rt=$(( 24 * 3600 )) \
        -pe $pe $np \
        -j y \
        run.sh "${bindir}/${prog}" \
        -w $w \
        inputs/M05-D{4,5,6,7,8}.txt \
        inputs/M06-D{5,6,7,8,11,10}.txt \
        ;
done
done
done

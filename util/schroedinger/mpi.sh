#! /bin/bash

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

bindir="run/${flavor}"
if test ! -d "$bindir"; then
    echo 1>&2 "Cannot read directory: '$bindir'"
    exit 1
fi

flavor="${rev}-${compiler}-${mpi}"

nps='16 32 64' # 128 256 512'
ws='4096 1024 256 64 16' # 4 1'

for np in $nps; do
for w in $ws; do
for prog in \
    rank-int-mpi \
    crank-int-mpi \
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
        -l s_cpu=$(( 24 * 3600 )) \
        -pe openmpi2 $np \
        -j y \
        run.sh "${bindir}/${prog}" \
        -w $w \
        inputs/M05-D{4,5,6,7,8}.txt \
        inputs/M06-D{5,6,7,8,11,10}.txt \
        ;
done
done
done

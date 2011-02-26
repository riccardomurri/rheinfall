#! /bin/bash

# no core dumps; especially with MPI, they just fill the home directory
ulimit -c 0

# avoid annoying modulefiles error?
source /panfs/panfs0.ften.es.hpcn.uzh.ch/share/software/Modules/default/init/sh

exec "$@"

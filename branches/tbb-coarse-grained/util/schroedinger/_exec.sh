#! /bin/bash

# no core dumps; especially with MPI, they just fill the home directory
ulimit -c 0

exec "$@"

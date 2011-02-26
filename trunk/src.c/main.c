/**
 * @file   main.c
 *
 * Compute the rank of an integer matrix using the ``Rheinfall''
 * algorithm.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2010 riccardo.murri@gmail.com.  All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 *
 */


#include "config.h"


#ifdef _OPENMP
# include <omp.h>
#endif


#ifdef WITH_MPI
# include <mpi.h>
#endif


#if defined(WITH_MPI) && defined(_OPENMP)
static omp_lock_t mpi_send_lock_;
#endif


#ifdef WITH_MPI
static const char const*
mpi_threading_model_name(const int provided)
{
  if (MPI_THREAD_SINGLE == provided)
    return "MPI_THREAD_SINGLE";
  else if (MPI_THREAD_FUNNELED == provided)
    return "MPI_THREAD_FUNNELED";
  else if (MPI_THREAD_SERIALIZED == provided)
    return "MPI_THREAD_SERIALIZED";
  else if (MPI_THREAD_MULTIPLE == provided)
    return "MPI_THREAD_MULTIPLE";
  else 
    return "an unknown threading model";
}
#endif


int mpi_init(int* argc_p, char*** argv_p)
{
#ifdef WITH_MPI
# ifdef _OPENMP
#  ifdef WITH_MPI_SERIALIZED
  const int required = (1 == omp_get_num_threads() ? MPI_THREAD_SINGLE : MPI_THREAD_SERIALIZED);
#  else
  const int required = (1 == omp_get_num_threads() ? MPI_THREAD_SINGLE : MPI_THREAD_MULTIPLE);
#  endif // defined(WITH_MPI_SERIALIZED)
  int provided;
  const int rc = MPI_Init_thread(argc_p, argv_p, required, &provided);
  if (required > provided) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    fprintf(stderr, 
            "WARNING: MPI rank %d requested %s but the MPI library provided %s\n",
            myrank, mpi_threading_model_name(required), mpi_threading_model_name(provided));
  };
  return rc;
# else
  // no OpenMP, use non-threaded MPI
  return MPI_Init(argc_p, argv_p);
# endif // defined(_OPENMP)
#endif // defined(WITH_MPI)
}



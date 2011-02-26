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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "common.h"
#include "rheinfall.h"
#include "switchboard.h"

#ifdef _OPENMP
# include <omp.h>
#endif

#ifdef WITH_MPI
# include <mpi.h>
#endif

#include <errno.h>
#include <getopt.h> // getopt_long
#include <signal.h> // sigaltstack, sigaction, sig*set
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // memcpy
#include <sys/resource.h> // getrlimit, setrlimit
#include <sys/time.h>
#include <ucontext.h>     // getcontext, setcontext
#include <unistd.h>


#if defined(WITH_MPI_SERIALIZED) || defined(WITH_MPI_MULTIPLE)
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

#ifdef WITH_MPI
int mpi_init(int* argc_p, char*** argv_p)
{
#ifdef _OPENMP
# ifdef WITH_MPI_SERIALIZED
  const int required = (1 == omp_get_num_threads() ? MPI_THREAD_SINGLE : MPI_THREAD_SERIALIZED);
# else
  const int required = (1 == omp_get_num_threads() ? MPI_THREAD_SINGLE : MPI_THREAD_MULTIPLE);
# endif // defined(WITH_MPI_SERIALIZED)
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
#else
  // no OpenMP, use non-threaded MPI
  return MPI_Init(argc_p, argv_p);
#endif // defined(_OPENMP)
}
#endif // WITH_MPI


//
// main
//


int verbose; // FIXME: should be a `rheinfall` parameter


void
usage(FILE* out, int argc, char** argv)
{
  fprintf(out, "Usage: %s matrix_file.sms\n", argv[0]);
}


void 
sigint(int signum)
{
  fprintf(stderr, "*** Terminated upon user request (SIGINT) ***\n");
  exit(signum);
}


static ucontext_t main_loop_ctx;
static bool got_sigfpe = false;

void 
sigfpe(int signum, siginfo_t* siginfo, void* ucp)
{
  fprintf(stderr, "*** Arithmetic error (SIGFPE): ");
  switch(siginfo->si_code) {
  case FPE_INTDIV: fprintf(stderr, "integer divide by zero"); break;
  case FPE_INTOVF: fprintf(stderr, "integer overflow"); break;
  case FPE_FLTDIV: fprintf(stderr, "floating-point divide by zero"); break;
  case FPE_FLTOVF: fprintf(stderr, "floating-point overflow"); break;
  case FPE_FLTUND: fprintf(stderr, "floating-point underflow"); break;
  case FPE_FLTRES: fprintf(stderr, "floating-point inexact result"); break;
  case FPE_FLTINV: fprintf(stderr, "floating-point invalid operation"); break;
  case FPE_FLTSUB: fprintf(stderr, "subscript out of range"); break;
  };
  fprintf(stderr, " ***");

  // try to bounce back into main loop
  got_sigfpe = true;
  setcontext(&main_loop_ctx);

  // if `setcontext` returns, we're in unrecoverable troubles...
  fprintf(stderr, "Cannot continue.\n");
  exit(signum);
}


int
main(int argc, char** argv)
{
  if (1 == argc) {
    usage(stdout, argc, argv);
    return 0;
  }

#ifdef WITH_MPI
  // disable core dumping, as "ulimit -c" is not preserved by "mpirun"
  struct rlimit core;
  core.rlim_cur = 0;
  core.rlim_max = 0;
  setrlimit(RLIMIT_CORE, &core);
#endif

  // set up alternate stack, to gracefully handle out-of-memory errors
  stack_t ss, oss;
  ss.ss_sp = xmalloc(SIGSTKSZ); 
  if (NULL == ss.ss_sp) {
    /* Handle error */
    fprintf(stderr, "Cannot allocate memory for alternate stack. Aborting.\n");
    return 1;
  };
  ss.ss_size = SIGSTKSZ;
  ss.ss_flags = 0;
  if (sigaltstack(&ss, &oss) == -1) {
    fprintf(stderr, "Cannot set alternate stack: sigaltstack() failed."
            " Aborting.\n");
    return 1;
  };


#ifdef WITH_MPI
  int me;
# ifdef _OPENMP
#  ifdef WITH_MPI_SERIALIZED
  const int required = (1 == omp_get_num_threads() ? MPI_THREAD_SINGLE : MPI_THREAD_SERIALIZED);
#  else
  const int required = (1 == omp_get_num_threads() ? MPI_THREAD_SINGLE : MPI_THREAD_MULTIPLE);
#  endif
  int provided;
  MPI_Init_thread(&argc, &argv, required, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  if (required > provided) {
    fprintf(stderr, "WARNING: MPI rank %d requested %s but MPI library provided %s.\n",
            me, mpi_threading_model_name(required), mpi_threading_model_name(provided));
  }
# else
  // no OpenMP, use non-threaded MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
# endif
#else
  const int me = 0;
#endif // WITH_MPI

  // parse command-line options
  static struct option long_options[] = {
    {"help", 0, 0,'h'},
    {"memory", 1, 0, 'm'},
    {"verbose", 0, 0,'v'},
    {0, 0, 0, 0}
  };

  int c;
  while (true) {
    c = getopt_long(argc, argv, "hm:v",
                    long_options, NULL);
    if (-1 == c)
      break;
    else if ('h' == c) {
      usage(stdout, argc, argv);
      return 0;
    }
    else if ('v' == c)
      verbose = 2;
    else if ('m' == c) {
      rlim_t required_memory = atoi(optarg);
      if (required_memory < 1) {
        fprintf(stderr, "Argument to option -m/--memory must be a positive integer,"
                " got '%ld' instead. Aborting.\n", required_memory);
        return 1;
      };
      required_memory *= 1000*1000*1000;
      struct rlimit limit;
      getrlimit(RLIMIT_AS, &limit);
      if (limit.rlim_cur != RLIM_INFINITY && limit.rlim_cur > required_memory) {
        fprintf(stderr, "%s: WARNING: Requested memory limit of %ld bytes"
                " is higher than current system soft limit of %ld bytes."
                " Capping requested memory to system limit.\n",
                argv[0], required_memory, limit.rlim_cur);
        required_memory = limit.rlim_cur;
      };
      if (limit.rlim_max != RLIM_INFINITY && limit.rlim_max > required_memory) {
        fprintf(stderr, "%s: WARNING: Requested memory limit of %ld bytes"
                " is higher than current system hard limit of %ld bytes."
                " Capping requested memory to system limit.\n",
                argv[0], required_memory, limit.rlim_max);
        required_memory = limit.rlim_max;
      };
      limit.rlim_cur = required_memory;
      limit.rlim_max = required_memory;
      setrlimit(RLIMIT_AS, &limit);
    }
    else {
      fprintf(stderr, "Unknown option '%c': type '%s --help' to get usage help.\n",
              c, argv[0]);
      return 1;
    };
  };

  // gracefully handle SIGINT, so that we can interrupt the process
  // and still get profiling information dumped to disk
  struct sigaction sa;
  sa.sa_handler = sigint;
  sa.sa_flags = 0;
  sigemptyset(&sa.sa_mask);
  sigaddset(&sa.sa_mask, SIGINT);
  sigaction (SIGINT, &sa, NULL);

  // start processing
  for (int i = optind; i < argc; ++i)
    {
      FILE* input = fopen(argv[i], "r");
      if (NULL == input) {
        fprintf(stderr, "Cannot open file '%s' for reading: %s."
                " Skipping this one.\n",
                argv[i], strerror(errno));
        continue;
      };

      if (0 == me) {
        printf("%s: file:%s", argv[0], argv[i]);
#ifdef WITH_MPI
        int mpi_comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_comm_size);
        printf(" mpi:%d", mpi_comm_size);
#endif
#ifdef _OPENMP
        printf(" omp:%d", omp_get_max_threads());
#endif
      }

      // read matrix dimensions
      coord_t rows, cols;
      char M;
      int rc = fscanf(input, "%ld %ld %c", &rows, &cols, &M);
      if (rc < 3) {
        fprintf(stderr, "Cannot read SMS file header from file '%s'."
                " Skipping this one.\n", argv[i]);
        continue;
      };
      assert ('M' == M);
      if (0 == me)
        printf(" rows:%ld cols:%ld", rows, cols);

      // possibly transpose matrix so that rows = min(rows, cols)
      // i.e., minimize the number of eliminations to perform
      bool transpose = false;
      if (rows > cols) {
        coord_t tmp = rows;
        rows = cols;
        cols = tmp;
        transpose = true;
      };

#ifdef WITH_MPI
      switchboard_t* rf = switchboard_new(cols, MPI_COMM_WORLD);
#else
      switchboard_t* rf = switchboard_new(cols);
#endif
      coord_t nnz = read_sms_file(rf, input, &rows, &cols, transpose);
      fclose(input);
      if (nnz < 0) {
        fprintf(stderr, "Could not read SMS file '%s': %s.  Skipping it.\n",
                argv[i], (-EILSEQ == nnz ? "malformed SMS header" 
                          : (-EINVAL == nnz ? "columns in header do not match allocated columns"
                             : strerror(-nnz))));
      };
      if (0 == me)
        printf(" nonzero:%ld", nnz);

      // handle SIGFPE: math errors will get back into the main loop and
      // we can proceed with next file
      sa.sa_sigaction = sigfpe;
      sa.sa_flags = SA_SIGINFO|SA_ONSTACK;
      sigemptyset(&sa.sa_mask);
      sigaddset(&sa.sa_mask, SIGFPE);
      sigaction(SIGFPE, &sa, NULL);

      got_sigfpe = false;
      getcontext(&main_loop_ctx);
      // `sigfpe` handler will resume here
      if (got_sigfpe) {
        // clear flag and skip to next command-line arg
        got_sigfpe = false;
        continue;
      };

      // base for computing wall-clock time
      struct timeval wc_t0;
      gettimeofday(&wc_t0, NULL);

      // base for computing CPU time
      struct rusage ru;
      getrusage(RUSAGE_SELF, &ru);
      struct timeval cpu_t0; memcpy(&cpu_t0, &ru.ru_utime, sizeof(struct timeval));

      coord_t rk = rank(rf);
      if (0 == me)
        printf(" rank:%ld", rk);

      // compute CPU time delta
      getrusage(RUSAGE_SELF, &ru);
      struct timeval cpu_t1; memcpy(&cpu_t1, &ru.ru_utime, sizeof(struct timeval));
      struct timeval tdelta; timersub(&cpu_t1, &cpu_t0, &tdelta);
      double consumed = tdelta.tv_sec + (tdelta.tv_usec / 1000000.0);

      // compute wall-clock time delta
      struct timeval wc_t1; gettimeofday(&wc_t1, NULL);
      timersub(&wc_t1, &wc_t0, &tdelta);
      double elapsed = tdelta.tv_sec + (tdelta.tv_usec / 1000000.0);

      if (0 == me)
        printf(" cputime:%.3f wctime:%.3f\n", consumed, elapsed);

      switchboard_free(rf);
    };

#ifdef WITH_MPI
  MPI_Finalize();
#endif

  // free resources
  sigaltstack(&oss, &ss);
  free(ss.ss_sp);

  return 0;
}

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
 * Copyright (c) 2010, 2011 riccardo.murri@gmail.com.  All rights reserved.
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


void
sigsegv(int signum, siginfo_t* siginfo, void* ucp)
{
  fprintf(stderr, "*** Segmentation fault (SIGSEGV): ");
  switch (siginfo->si_code) {
  case SEGV_MAPERR:  fprintf(stderr, "address not mapped"); break;
  case  SEGV_ACCERR: fprintf(stderr, "invalid permissions for mapped object"); break;
  default:           fprintf(stderr, "unknown fault"); break;
  };
  fprintf(stderr, " ***\n");
  exit(signum);
}


#ifdef WITH_GE
void 
sigusr2(int signum)
{
  fprintf(stderr, "*** Terminated upon user request (SIGINT) ***\n");
  exit(signum);
}
#endif // WITH_GE


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
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
#else
  const int me = 0;
#endif // WITH_MPI

  // parse command-line options
  int width = 1;
  bool transpose = false;

  static struct option long_options[] = {
    {"help",      0, 0, 'h'},
    {"verbose",   0, 0, 'v'},
    {"transpose", 0, 0, 't'},
    {"width",     1, 0, 'w'},
    {0, 0, 0, 0}
  };

  int c;
  while (true) {
    c = getopt_long(argc, argv, "hvtw:",
                    long_options, NULL);
    if (-1 == c)
      break;
    else if ('h' == c) {
      usage(stdout, argc, argv);
      return 0;
    }
    else if ('t' == c)
      transpose = true;
    else if ('v' == c)
      verbose = 2;
    else if ('w' == c) {
      width = atoi(optarg);
    }
    else {
      fprintf(stderr, "Unknown option '%c': type '%s --help' to get usage help.\n",
              c, argv[0]);
      return 1;
    };
  };

  // install signal handlers only if not running MPI: OpenMPI does a
  // better job at it than the simple-minded functions here.
  struct sigaction sa;
#ifndef WITH_MPI

  // gracefully handle SIGINT, so that we can interrupt the process
  // and still get profiling information dumped to disk
  sa.sa_handler = sigint;
  sa.sa_flags = 0;
  sigemptyset(&sa.sa_mask);
  sigaddset(&sa.sa_mask, SIGINT);
  sigaction (SIGINT, &sa, NULL);

  // dump backtrace on SIGSEGV
  sa.sa_sigaction = sigsegv;
  sa.sa_flags = SA_SIGINFO|SA_ONSTACK;
  sigemptyset(&sa.sa_mask);
  sigaddset(&sa.sa_mask, SIGSEGV);
  sigaction (SIGSEGV, &sa, NULL);
#endif // WITH_MPI

#ifdef WITH_GE
  // SIGUSR2 is used by GE to notify of a SIGKILL coming shortly,
  // catch it and exit gracefully.
  sa.sa_sigaction = sigsegv;
  sa.sa_flags = SA_SIGINFO|SA_ONSTACK;
  sigemptyset(&sa.sa_mask);
  sigaddset(&sa.sa_mask, SIGSEGV);
  sigaction (SIGSEGV, &sa, NULL);
#endif

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
        printf(" mpi:%d w:%d", mpi_comm_size, width);
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

#ifdef WITH_MPI
      switchboard_t* rf = switchboard_new(cols, width, MPI_COMM_WORLD);
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
        printf(" cputime:%.6f wctime:%.6f\n", consumed, elapsed);

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

/**
 * @file   main.cpp
 *
 * Sample application of the ``rheinfall'' algorithm for computing matrix rank.
 * This file won't compile as-is; you must @c typedef @c val_t and @c coord_t
 * first.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2010 riccardo.murri@gmail.com.  All rights reserved.
 *
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

#if defined(HAVE_GMPXX) and not defined(WITH_GMPXX)
# undef HAVE_GMPXX
#endif


#include <rheinfall.hpp>

#ifdef WITH_MPI
# include <boost/mpi.hpp>
  namespace mpi = boost::mpi;

# include <boost/serialization/export.hpp>
typedef rheinfall::Row<val_t,coord_t> Row_;
typedef rheinfall::SparseRow<val_t,coord_t> SparseRow_;
typedef rheinfall::DenseRow<val_t,coord_t> DenseRow_;

// elminate serialization overhead at the cost of
// never being able to increase the version.
BOOST_CLASS_IMPLEMENTATION(Row_, boost::serialization::object_serializable);
BOOST_CLASS_IMPLEMENTATION(SparseRow_, boost::serialization::object_serializable);
BOOST_CLASS_IMPLEMENTATION(DenseRow_, boost::serialization::object_serializable);

// eliminate object tracking (even if serialized through a pointer)
// at the risk of a programming error creating duplicate objects.
BOOST_CLASS_TRACKING(Row_, boost::serialization::track_never)
BOOST_CLASS_TRACKING(SparseRow_, boost::serialization::track_never)
BOOST_CLASS_TRACKING(DenseRow_, boost::serialization::track_never)
#endif

#ifdef _OPENMP
# include <omp.h>
#endif

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <getopt.h> // getopt_long
#include <signal.h> // sigaltstack, sigaction, sig*set
#include <string.h> // memcpy
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h> // getrlimit, setrlimit
#include <ucontext.h>     // getcontext, setcontext
#include <unistd.h>
#ifdef __GNUC__
// FIXME: `backtrace` and friends are provided by the GNU C library,
// not the GNU compiler, but there is no way to test for the latter
// here...  This should definitely be checked by the `configure`
// script.
# include <execinfo.h> // backtrace, backtrace_symbols_fd
# define HAVE_BACKTRACE_FN 1
#endif

int verbose = 0;

void
usage(std::ostream& out, const int argc, char const* const* argv)
{
  out << "Usage: " << argv[0] << " [options] matrix_file.sms" << std::endl
      << std::endl
      << "Options:" << std::endl
      << "  -m NUM  Limit memory consumption to NUM GBs;" << std::endl
      << "          if more memory is needed, the program will abort." << std::endl
#ifdef WITH_MPI
      << "  -w NUM  Divide matrix in bands of NUM columns each and distribute" << std::endl
      << "          bands to MPI ranks in a round-robin fashion. (Default: NUM=1)" << std::endl
#endif
      << "  -v      Verbosely report about rank computation." << std::endl
      << "  -h      Print this help text." << std::endl;
}


#if defined(WITH_MPI_SERIALIZED) or defined(WITH_MPI_MULTIPLE)
static std::string
mpi_threading_model_name(int provided)
{
  if (MPI_THREAD_SINGLE == provided)
    return "MPI_THREAD_SINGLE";
  else if (MPI_THREAD_FUNNELED == provided)
    return "MPI_THREAD_FUNNELED";
  else if (MPI_THREAD_SERIALIZED == provided)
    return "MPI_THREAD_SERIALIZED";
  else if (MPI_THREAD_MULTIPLE == provided)
    return "MPI_THREAD_MULTIPLE";
  else {
    std::ostringstream o;
    o << "unknown threading model (" << provided <<")";
    return o.str();
  };
}
#endif


#ifdef HAVE_BACKTRACE_FN
static const int  BTSIZE = 64;
static   void    *bt[64];    /* Array to store backtrace symbols */
static   size_t   btsize;    /* To store the exact no of values stored */

static inline void
print_backtrace()
{
  btsize = backtrace(bt, BTSIZE);
  backtrace_symbols_fd(bt, btsize, STDERR_FILENO);
}
#else
static inline void
print_backtrace()
{
  // no way to print backtrace with std C++ library
}
#endif


extern "C"
void 
sigint(int signum)
{
  std::cerr << "*** Terminated upon user request (SIGINT) ***" << std::endl;
  exit(signum);
}


static ucontext_t main_loop_ctx;
static bool got_sigfpe = false;

extern "C"
void 
sigfpe(int signum, siginfo_t* siginfo, void* ucp)
{
  std::cerr << "*** Arithmetic error (SIGFPE): ";
  switch(siginfo->si_code) {
  case FPE_INTDIV: std::cerr << "integer divide by zero"; break;
  case FPE_INTOVF: std::cerr << "integer overflow"; break;
  case FPE_FLTDIV: std::cerr << "floating-point divide by zero"; break;
  case FPE_FLTOVF: std::cerr << "floating-point overflow"; break;
  case FPE_FLTUND: std::cerr << "floating-point underflow"; break;
  case FPE_FLTRES: std::cerr << "floating-point inexact result"; break;
  case FPE_FLTINV: std::cerr << "floating-point invalid operation"; break;
  case FPE_FLTSUB: std::cerr << "subscript out of range"; break;
  };
  std::cerr << " ***" << std::endl;
  print_backtrace();

  // try to bounce back into main loop
  got_sigfpe = true;
  setcontext(&main_loop_ctx);

  // if `setcontext` returns, we're in unrecoverable troubles...
  std::cerr << "Cannot continue." << std::endl;
  exit(signum);
}


extern "C"
void
sigsegv(int signum, siginfo_t* siginfo, void* ucp)
{
  std::cerr << "*** SIGSEGV: "
            << (siginfo->si_code == SEGV_MAPERR?
                "address not mapped"
                : (siginfo->si_code == SEGV_ACCERR? 
                   "invalid permissions for mapped object"
                   : "unknown fault"))
            << " ***" << std::endl;
  print_backtrace();
  exit(signum);
}


int
main(int argc, char** argv)
{
#ifdef WITH_MPI
  // disable core dumping, as "ulimit -c" is not preserved by "mpirun"
  struct rlimit core;
  core.rlim_cur = 0;
  core.rlim_max = 0;
  setrlimit(RLIMIT_CORE, &core);
#endif

  // set up alternate stack, to gracefully handle out-of-memory errors
  stack_t ss, oss;
  ss.ss_sp = malloc(SIGSTKSZ); 
  if (ss.ss_sp == NULL) {
    /* Handle error */
    std::cerr << "Cannot allocate memory for alternate stack."
              << " Aborting." << std::endl;
    exit(1);
  };
  ss.ss_size = SIGSTKSZ;
  ss.ss_flags = 0;
  if (sigaltstack(&ss, &oss) == -1) {
    std::cerr << "Cannot set alternate stack: sigaltstack() failed."
              << " Aborting." << std::endl;
    exit(1);
  };


#ifdef WITH_MPI
# ifdef _OPENMP
#  ifdef WITH_MPI_SERIALIZED
  const int required = (1 == omp_get_num_threads() ? MPI_THREAD_SINGLE : MPI_THREAD_SERIALIZED);
#  else  
  const int required = (1 == omp_get_num_threads() ? MPI_THREAD_SINGLE : MPI_THREAD_MULTIPLE);
#  endif // WITH_MPI_SERIALIZED
  int provided;
  MPI_Init_thread(&argc, &argv, required, &provided);
  if (required > provided) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    std::cerr << "WARNING: MPI rank " << myrank << ":"
              << " requested " << mpi_threading_model_name(required) << ","
              << " but MPI library provided " << mpi_threading_model_name(provided)
              << std::endl;
  }
# else
  // no OpenMP, use non-threaded MPI
  MPI_Init(&argc, &argv);
# endif
  mpi::environment env(argc, argv);
  mpi::communicator world;
  int myid = world.rank();
#else
  const int myid = 0;
#endif // WITH_MPI

  // parse command-line options
  static struct option long_options[] = {
    {"help",    0, 0, 'h'},
    {"memory",  1, 0, 'm'},
    {"verbose", 0, 0, 'v'},
    {"width",   1, 0, 'w'},
    {0, 0, 0, 0}
  };

  coord_t width = 1;

  int c;
  while (true) {
    c = getopt_long(argc, argv, "hm:vw:",
                    long_options, NULL);
    if (-1 == c)
      break;
    else if ('h' == c) {
      usage(std::cout, argc, argv);
      return 0;
    }
    else if ('v' == c)
      verbose = 1;
    else if ('w' == c) {
      std::istringstream arg(optarg);
      arg >> width;
      if (width < 1) {
        std::cerr << "Argument to option -w/--width must be a positive integer."
                  << " Aborting." << std::endl;
        return 1;
      };
    }
    else if ('m' == c) {
      rlim_t required_memory;
      std::istringstream arg(optarg);
      arg >> required_memory;
      if (required_memory < 1) {
        std::cerr << "Argument to option -m/--memory must be a positive integer."
                  << " Aborting." << std::endl;
        return 1;
      };
      required_memory *= 1000*1000*1000;
      struct rlimit limit;
      getrlimit(RLIMIT_AS, &limit);
      if (limit.rlim_cur != RLIM_INFINITY and limit.rlim_cur > required_memory) {
        std::cerr << argv[0] 
                  << " WARNING: Requested memory limit of " 
                  << required_memory 
                  << " bytes is higher than current system soft limit of "
                  << limit.rlim_cur << " bytes."
                  << " Capping requested memory to system limit."
                  << std::endl;
        required_memory = limit.rlim_cur;
      };
      if (limit.rlim_max != RLIM_INFINITY and limit.rlim_max > required_memory) {
        std::cerr << argv[0] 
                  << " WARNING: Requested memory limit of " 
                  << required_memory 
                  << " bytes is higher than current system hard limit of "
                  << limit.rlim_max << " bytes."
                  << " Capping requested memory to system limit."
                  << std::endl;
        required_memory = limit.rlim_max;
      };
      limit.rlim_cur = required_memory;
      limit.rlim_max = required_memory;
      setrlimit(RLIMIT_AS, &limit);
    }
    else {
      std::cerr << "Unknown option; type '" << argv[0] << " --help' to get usage help."
                << std::endl;
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

#ifdef HAVE_BACKTRACE_FN
  // dump backtrace on SIGSEGV
  sa.sa_sigaction = sigsegv;
  sa.sa_flags = SA_SIGINFO|SA_ONSTACK;
  sigemptyset(&sa.sa_mask);
  sigaddset(&sa.sa_mask, SIGSEGV);
  sigaction (SIGSEGV, &sa, NULL);
#endif

  // start processing
  for (int i = optind; i < argc; ++i)
    {
      std::ifstream input(argv[i]);
      if (input.fail()) {
        std::cerr << "Cannot open file '" <<argv[i]<< "' for reading - skipping." 
                  << std::endl;
        continue;
      };

      if (0 == myid) {
        std::cout << argv[0] << " file:"<<argv[i];
#ifdef WITH_MPI
        std::cout << " mpi:"<< world.size();
#endif
#ifdef _OPENMP
        std::cout << " omp:"<< omp_get_max_threads();
#endif
      }

      // read matrix dimensions
      long rows, cols;
      char M;
      input >> std::skipws >> rows >> cols >> M;
      assert ('M' == M);
      if (0 == myid) {
        std::cout << " rows:" << rows;
        std::cout << " cols:" << cols;
      };

      // possibly transpose matrix so that rows = min(rows, cols)
      // i.e., minimize the number of eliminations to perform
      bool transpose = false;
      if (rows > cols) {
        std::swap(rows, cols);
        transpose = true;
      };

#ifdef WITH_MPI
      rheinfall::Rheinfall<val_t, coord_t> w(world, cols, width);
#else
      rheinfall::Rheinfall<val_t, coord_t> w(cols, width);
#endif
      long nnz = w.read(input, rows, cols, true, transpose);
      input.close();
      if (0 == myid)
        std::cout << " nonzero:" << nnz;

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

      int rank = w.rank();
      if (0 == myid)
        std::cout << " rank:" << rank;

      // compute CPU time delta
      getrusage(RUSAGE_SELF, &ru);
      struct timeval cpu_t1; memcpy(&cpu_t1, &ru.ru_utime, sizeof(struct timeval));
      struct timeval tdelta; timersub(&cpu_t1, &cpu_t0, &tdelta);
      double consumed = tdelta.tv_sec + (tdelta.tv_usec / 1000000.0);

      // compute wall-clock time delta
      struct timeval wc_t1; gettimeofday(&wc_t1, NULL);
      timersub(&wc_t1, &wc_t0, &tdelta);
      double elapsed = tdelta.tv_sec + (tdelta.tv_usec / 1000000.0);

      if (0 == myid)
        std::cout << " cputime:" << std::fixed << std::setprecision(3) << consumed;
      if (0 == myid)
        std::cout << " wctime:" << std::fixed << std::setprecision(3) << elapsed << std::endl;
    };

#ifdef WITH_MPI
  MPI_Finalize();
#endif

  // free other malloc'd resources (to please Valgrind)
  free(ss.ss_sp);
  sigaltstack(&oss, NULL);

  return 0;
}

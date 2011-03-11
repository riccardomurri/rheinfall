/**
 * @file   rank.cpp
 *
 * Sample application of the ``Rheinfall'' algorithm for computing matrix rank.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2010, 2011 riccardo.murri@gmail.com.  All rights reserved.
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

#include <config.h>

// matrix dimensions should fit into an `long` integer type
typedef long coord_t;

// Select what type will be used for the matrix coefficients.
#if defined(WITH_MODULAR_VALUES)

# include <types/modular.hpp>
# if defined(HAVE_LONG_LONG_INT)
typedef long long mod_int_t;
# else
typedef long mod_int_t;
# endif
typedef modular::Modular<mod_int_t> val_t;

#elif defined(WITH_MODULAR_INT64_VALUES)

# include <types/modular.hpp>
typedef int64_ mod_int_t;
typedef modular::Modular<mod_int_t> val_t;

#elif defined(WITH_INT_VALUES)

// always use widest standard integer type available
# if defined(HAVE_LONG_LONG_INT)
typedef long long val_t;
# else
typedef long val_t;
# endif

#elif defined(WITH_INT32_VALUES)

typedef int32_t val_t;

#elif defined(WITH_INT64_VALUES)

typedef int64_t val_t;

#elif defined(WITH_DOUBLE_VALUES)

// always select the widest floating-point type available
# ifdef HAVE_LONG_DOUBLE
typedef long double val_t;
# else
typedef double val_t;
# endif 

#elif defined(WITH_MPQ_VALUES)

# ifndef HAVE_GMPXX
#  error This source requires GMP to compile; messed up autoconf settings?
# else
#  define WITH_GMPXX
# endif
# include <gmpxx.h>
# include <types/gmpxx.hpp>
typedef mpq_class val_t;

#elif defined(WITH_MPZ_VALUES)

# ifndef HAVE_GMPXX
#  error This source requires GMP to compile; messed up autoconf settings?
# else
#  define WITH_GMPXX
# endif
# include <gmpxx.h>
# include <types/gmpxx.hpp>
typedef mpz_class val_t;

#elif defined(WITH_XINT_VALUES)

# include <boost/xint/integer.hpp>
# include <types/xint.hpp>
typedef boost::xint::integer val_t;
#else

# error Please define one of: WITH_INT_VALUES, WITH_MODULAR_VALUES, WITH_DOUBLE_VALUES, WITH_MPZ_VALUES, WITH_MPQ_VALUES or WITH_XINT_VALUES

#endif // WITH_..._VALUES 


// use WITH_MPI_SERIALIZED / WITH_MPI_MULTIPLE to select which
// threading model to request to the MPI-2 library
#if (defined(WITH_MPI_SERIALIZED) or defined (WITH_MPI_MULTIPLE)) and not defined(WITH_MPI)
# define WITH_MPI
#endif


// XXX: Undefine `HAVE_GMPXX` if we're not using it.
#if defined(HAVE_GMPXX) and not defined(WITH_GMPXX)
# undef HAVE_GMPXX
#endif


#include <rheinfall.hpp>
#include <types.hpp>

#ifdef WITH_MPI
# include <boost/mpi.hpp>
  namespace mpi = boost::mpi;

// # include <boost/serialization/export.hpp>
// typedef rheinfall::Row<val_t,coord_t> Row_;
// typedef rheinfall::SparseRow<val_t,coord_t> SparseRow_;
// typedef rheinfall::DenseRow<val_t,coord_t> DenseRow_;

// // elminate serialization overhead at the cost of
// // never being able to increase the version.
// BOOST_CLASS_IMPLEMENTATION(Row_, boost::serialization::object_serializable);
// BOOST_CLASS_IMPLEMENTATION(SparseRow_, boost::serialization::object_serializable);
// BOOST_CLASS_IMPLEMENTATION(DenseRow_, boost::serialization::object_serializable);

// // eliminate object tracking (even if serialized through a pointer)
// // at the risk of a programming error creating duplicate objects.
// BOOST_CLASS_TRACKING(Row_, boost::serialization::track_never)
// BOOST_CLASS_TRACKING(SparseRow_, boost::serialization::track_never)
// BOOST_CLASS_TRACKING(DenseRow_, boost::serialization::track_never)
#endif // WITH_MPI

#ifdef _OPENMP
# include <omp.h>
#endif

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
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
#include <execinfo.h> // backtrace, backtrace_symbols_fd

int verbose = 0;

void
usage(std::ostream& out, const int argc, char const* const* argv)
{
  out << "Usage: " << argv[0] << " [options] matrix_file.sms" << std::endl
      << std::endl
      << "Options:" << std::endl
#ifdef WITH_MODULAR_VALUES
      << "  -p NUM  Perform computations modulo NUM (default: 2038076783)." << std::endl
#endif
      << "  -m NUM  Limit memory consumption to NUM GBs;" << std::endl
      << "          if more memory is needed, the program will abort." << std::endl
#if defined(WITH_MPI) or defined(_OPENMP)
      << "  -w NUM  Divide matrix in bands of NUM columns each and distribute" << std::endl
      << "          bands to MPI ranks / OpenMP tasks in a round-robin fashion." << std::endl
      << "          (Default: NUM=1)" << std::endl
#endif
      << "  -v      Verbosely report about rank computation." << std::endl
      << "  -h      Print this help text." << std::endl;
}


#ifdef WITH_MPI
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


#if defined(HAVE_BACKTRACE) && defined(HAVE_BACKTRACE_SYMBOLS_FD)
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


#ifdef WITH_GE
extern "C"
void
sigusr2(int signum, siginfo_t* siginfo, void* ucp)
{
  std::cerr << "*** SIGUSR2: Notification of SIGKILL from GridEngine, exiting now. ***" << std::endl;
  exit(signum);
}
#endif


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
# if defined(_OPENMP) or defined(WITH_MPI_SERIALIZED) or defined(WITH_MPI_MULTIPLE)
#  ifdef WITH_MPI_SERIALIZED
  const int required = (1 == omp_get_num_threads() ? MPI_THREAD_SINGLE : MPI_THREAD_SERIALIZED);
#  else  // WITH_MPI_MULTIPLE
  const int required = (1 == omp_get_num_threads() ? MPI_THREAD_SINGLE : MPI_THREAD_MULTIPLE);
#  endif 
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
  // use non-threaded MPI
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
#ifdef WITH_MODULAR_VALUES
    {"modulus", 1, 0, 'p'},
#endif
    {"tranpose",0, 0, 't'},
    {"verbose", 0, 0, 'v'},
    {"width",   1, 0, 'w'},
    {0, 0, 0, 0}
  };

  coord_t width = 1;
  bool transpose = false;

#ifdef WITH_MODULAR_VALUES
  modular::Modular<mod_int_t>::global_set_modulus(2038076783);
#endif

  int c;
  while (true) {
    c = getopt_long(argc, argv, "hp:tvw:",
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
#ifdef WITH_MODULAR_VALUES
    else if ('p' == c) {
      mod_int_t p;
      std::istringstream arg(optarg);
      arg >> p;
      modular::Modular<mod_int_t>::global_set_modulus(p);
    }
#endif
    else if ('t' == c) {
      transpose = true;
    }
    else {
      std::cerr << "Unknown option; type '" << argv[0] << " --help' to get usage help."
                << std::endl;
      return 1;
    };
  };

  // install signal handlers only if not running MPI: OpenMPI does a
  // better job at it than the simple-minded functions here.
  struct sigaction sa;
#if not defined(WITH_MPI)
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
  sa.sa_sigaction = sigusr2;
  sa.sa_flags = SA_SIGINFO|SA_ONSTACK;
  sigemptyset(&sa.sa_mask);
  sigaddset(&sa.sa_mask, SIGUSR2);
  sigaction (SIGUSR2, &sa, NULL);
#endif // WITH_GE

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
#if defined(WITH_MPI) or defined(_OPENMP)
        std::cout << " w:" << width;
#endif
      }

      // read matrix dimensions
      coord_t rows, cols;
      char M;
      input >> std::skipws >> rows >> cols >> M;
      assert ('M' == M);
      if (0 == myid) {
        std::cout << " rows:" << rows;
        std::cout << " cols:" << cols;
      };
      if (transpose)
        std::swap(rows, cols);

#ifdef WITH_DOUBLE_VALUES
      rheinfall::is_zero_traits<val_t>::tolerance = 
        rows * cols * std::numeric_limits<val_t>::epsilon();
#endif

#ifdef WITH_MPI
      rheinfall::Rheinfall<val_t, coord_t> rf(world, cols, width);
#else
      rheinfall::Rheinfall<val_t, coord_t> rf(cols, width);
#endif

      coord_t nnz = rf.read(input, rows, cols, true, transpose);
      input.close();
      if (0 == myid) {
        std::cout << " nonzero:" << nnz;
      };

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

      int rank = rf.rank();
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
        std::cout << " cputime:" << std::fixed << std::setprecision(6) << consumed;
      if (0 == myid)
        std::cout << " wctime:" << std::fixed << std::setprecision(6) << elapsed << std::endl;
    };

#ifdef WITH_MPI
  MPI_Finalize();
#endif

  // free other malloc'd resources (to please Valgrind)
  free(ss.ss_sp);
  sigaltstack(&oss, NULL);

  return 0;
}
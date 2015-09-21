# Introduction #

Rheinfall consists of a C++ library and example programs.

The C++ library is a header-only library; no compilation is necessary
and the installation procedure just copies the header files into a
directory `rheinfall/` on your compiler's _include_ path.

The example programs compute the rank of an input matrix in
[SMS format][2] The same program is offered in C and C++ versions.
The C++ version demoes the features available through the library; the
C version implements the same algorithm directly, so you will have to
copy and modify the code to adapt it, e.g., to use a different data
type.

At the moment there is no packaged and released code;
the only way to build Rheinfall is to install from the
[Git sources][1].

[1]: https://github.com/riccardomurri/rheinfall
[2]: http://ljk.imag.fr/membres/Jean-Guillaume.Dumas/simc.html


# Quickstart instructions #

The following instructions provide the simplest route to compiling all
the example programs in Rheinfall and trying them out.  See section
["How to build Rheinfall from sources"](#How_to_build_Rheinfall_from_sources.md)
for a detailed description of options.

In addition to the standard Linux commands and a C/C++ compiler, you
will need the commands `make`, `svn`, `autoreconf`, `wget` installed.  (See
["How to build Rheinfall from sources"](#How_to_build_Rheinfall_from_sources.md)
for a table with download URLs for this software.)

If you have the dependencies in place, the following commands will
place the Rheinfall sources in a directory called `rheinfall` and
compile the example programs:
```
  # check out the Rheinfall sources
  svn co http://rheinfall.googlecode.com/svn/trunk rheinfall

  # install the required libraries into ./sw
  cd rheinfall/
  ./util/prereq.sh sw

  # build the example programs
  autoreconf -i
  ./configure --with-boost=`pwd`/sw --with-gmp=`pwd`/sw
  make
```

The test suite will download some matrices of the _Mgn_ group from the online
[Sparse Integer Matrices Collection](http://ljk.imag.fr/membres/Jean-Guillaume.Dumas/simc.html)
and use the compiled programs to compute their ranks (and check if
this is correct).  To run the test suite:
```
  make check
```


# How to build Rheinfall from sources #

At the moment there is no packaged and released code; the only way to
build Rheinfall is to install from the SVN sources.

This section provides complete instructions on how to build Rheinfall
example programs and the various compilation options.

Rheinfall uses the GNU automake/autoconf mechanism for generating
Makefiles and driving compilation. The following sections detail the
important options available to the `./configure` script and their
effect.


## What is actually built ##

Since the C++ library is a header-only library, no compilation is necessary.

C++ example programs are built in the following flavors:

  * `rank-int`, `rank-int32`, `rank-64`: integer-based exact computation, using fixed-width integer types: `rank-int` uses the widest integer type available in the compiler, `rank-int32` uses 32-bit integers and `rank-int64` uses 64-bit integers (if available);
  * `rank-mod`: modular arithmetic (32-bit modulus);
  * `rank-double`: double precision;
  * `rank-mpq`, `rank-mpz`: arbitrary precision rational and integer types (the [GNU GMP library](http://www.gnu.org/software/gmp) is required for building and running these).

The C example programs are built in the following flavors:
  * `crank-int`: integer based exact computation, using the widest integer type available from the compiler;
  * `crank-double`: double precision.

### MPI support ###

There are MPI-based variants of the example programs, available only
for some types (those that have a direct mapping to MPI types):
  * `rank-int-mpi`, `rank-mod-mpi`, `rank-double-mpi`: parallel variants of the corresponding C++ example programs
  * `crank-int-mpi`, `crank-double-mpi`: parallel variants of the C example programs

MPI-based variants of the example programs can be built iff the
following software can be found by the `./configure` script:
  * an MPI library and the corresponding compiler wrappers `mpicc` and `mpicxx`.
  * a patched version of the Boost.MPI library (see below).

Check that the MPI library and compiler wrapper found by `./configure`
is the one used to build Boost.MPI; if there is a mistmatch, errors
will arise at compilation, link, or run time.

You can disable compiling MPI code altogether (and not even check for
MPI presence) by passing option `--with-mpi=no` to `./configure`.

### OpenMP support ###

If your compiler has OpenMP support, programs `rank-*-omp` and
`crank-*-omp` will be built, that implement the parallel version of
Rheinfall using OpenMP threads.  Use of the OpenMP-based versions is
_presently discouraged:_ the OpenMP code does not have any significant
performance advantage over the serial one.  There is no way, at
present, to disable compilation of the OpenMP variants.


## Required Software ##

In addition to the standard Linux commands and a C/C++ compiler,
the following software needs to be installed prior to building Rheinfall:

| **command name** | **home page** | **package name** |
|:-----------------|:--------------|:-----------------|
| `svn`            | SubVersion, http://subversion.tigris.org/         | `subversion` ([Debian](http://packages.debian.org/squeeze/subversion), [Ubuntu](http://packages.ubuntu.com/subversion), [CentOS/Fedora/RHEL and other RPM-based distributions](http://rpmfind.net/linux/rpm2html/search.php?query=subversion)) |
| `make`           | GNU make, http://www.gnu.org/software/make/       | `make` ([Debian](http://packages.debian.org/squeeze/make), [Ubuntu](http://packages.ubuntu.com/make), [CentOS/Fedora/RHEL and other RPM-based distributions](http://rpmfind.net/linux/rpm2html/search.php?query=make)) |
| `wget`           | GNU wget, http://www.gnu.org/software/wget/       | `wget` ([Debian](http://packages.debian.org/squeeze/wget), [Ubuntu](http://packages.ubuntu.com/wget), [CentOS/Fedora/RHEL and other RPM-based distributions](http://rpmfind.net/linux/rpm2html/search.php?query=wget)) |
| `automake`       | GNU automake _(at least version 1.10),_ http://www.gnu.org/software/automake/ | `automake` ([Debian](http://packages.debian.org/squeeze/automake), [Ubuntu](http://packages.ubuntu.com/automake), [CentOS/Fedora/RHEL and other RPM-based distributions](http://rpmfind.net/linux/rpm2html/search.php?query=automake)) |
| `autoreconf`     | GNU autoconf _(at least version 2.65),_ http://www.gnu.org/software/autoconf/ | `automake` ([Debian](http://packages.debian.org/squeeze/automake), [Ubuntu](http://packages.ubuntu.com/automake), [CentOS/Fedora/RHEL and other RPM-based distributions](http://rpmfind.net/linux/rpm2html/search.php?query=automake)) |
| `libtoolize`     | GNU libtool, http://www.gnu.org/software/libtool/ | `libtool` ([Debian](http://packages.debian.org/squeeze/libtool), [Ubuntu](http://packages.ubuntu.com/libtool), [CentOS/Fedora/RHEL and other RPM-based distributions](http://rpmfind.net/linux/rpm2html/search.php?query=libtool)) |
|                  | Boost libraries, http://boost.org/       | _compile from sources, almost no distribution has recent enough packages_ |

The [GNU GMP](http://www.gnu.org/software/gmp) and the
[Boost.MPI](http://www.boost.org/doc/libs/1_46_1/doc/html/mpi.html)
libraries are needed to enable, respectively, arbitrary precision
support and distribued-memory parallel processing support.  The helper script
[util/prereq.sh](http://rheinfall.googlecode.com/svn/trunk/util/prereq.sh)
installs all the needed software into directory `sw`
(located within the directory where you invoke `util/prereq.sh`).
The following sections give more details if you want to compile
GMP and/or Boost.MPI without the helper script.

To change the installation directory for the auxiliary software, pass
the directory path as the sole argument to the `util/prereq.sh`
script, e.g.:
```
  cd /the/rheinfall/sources

  # install dependent software into ./mysw
  ./util/prereq.sh "mysw"
```
Invoking `util/prereq.sh` without any argument is equivalent to the
following:
```
  ./util/prereq.sh ./sw
```

### GNU GMP ###

In order to get the arbitrary precision support, you need to have GMP
installed (with C++ support!).

If the GMP library is not installed in the default compiler search
path, you need to pass the `--with-gmp=PATH` option to Rheinfall's
`./configure`.  The GMP library `PATH` should contain a `lib/`
subdirectory with the GMP and GMPxx libraries in it.

To disable building the arbitrary-precision example programs and not
check for GMP at all, pass option `--with-gmp=no` to Rheinfall's
`./configure` script.

### Boost.MPI ###

Rheinfall can use MPI to distribute the load and memory usage across a
cluster of machines.  It needs the Boost.MPI C++ library, patched to
implement the MPI "synchronous send" feature.
([The patch has been submitted](https://svn.boost.org/trac/boost/ticket/5245)
to the Boost issue tracker, but has not yet been incorporated into
the mainstream distribution.)

Note that the Boost library is always required
by Rheinfall; only the Boost.MPI support is optional!

To build the patched version of the Boost.MPI library, proceed as follows:

0. Download the Boost sources and unpack them into a directory `boost`;
1. Save the [boost.issend.patch](http://rheinfall.googlecode.com/svn/trunk/util/boost_1_45_0.issend.patch) file into the top-level directory of the Boost sources;
2. Apply the patch to the sources (the patch is known to apply cleanly to Boost.MPI version 1.45 and 1.46):
```
  patch -i boost_1_45_0.issend.patch
```
3. Now proceed to build and install Boost as per [instructions in the Boost manual](http://www.boost.org/doc/libs/1_46_1/more/getting_started/unix-variants.html#prepare-to-use-a-boost-library-binary); for example:
```
  bjam threading=multi variant=release --with-libraries=test,serialization,mpi ...
```

If the Boost library is not installed in the default compiler search
path, you need to pass the `--with-boost=PATH` option to Rheinfall's
`./configure`.  The `PATH` should contain a `lib/` subdirectory with
the Boost libraries (static or shared) in it, and an `include/boost`
subdirectory with the header files.

To disable building the MPI example programs and not check for
Boost.MPI at all, pass option `--with-mpi=no` to Rheinfall's
`./configure` script.  Note that the Boost library is always required
by Rheinfall; only the Boost.MPI support is optional!


## Compilation ##

Before compiling, ensure you have installed all the required
software; see section [Required Software](#Required_software.md) above.

Rheinfall uses the GNU automake/autoconf mechanism for generating
Makefiles and driving compilation.
Assuming you have installed the required software into directory `./sw`
(e.g., by using the supplied `util/prereq.sh` script, see above), then
you can proceed as with any automake program.
```
  # create the Makefile
  autoreconf -i
  ./configure --with-boost=./sw --with-gmp=./sw

  # compile
  make

  # download test matrices and run checks
  make check

  # copy library header files to their destination
  make install
```

Besides the standard options to the `./configure` script, the
following Rheinfall-specific options are recognized:

  * `--with-boost=PATH`: set the PATH where the Boost libraries (required) will be searched for; the shared or static libraries are expected to be in `PATH/lib`, and the header files in `PATH/include/boost`.  The Boost libraries are required by Rheinfall, so if they cannot be found compilation will stop.

  * `--with-gmp=PATH`: set the PATH where the GNU GMP libraries are searched for; the shared libraries are expected to be in `PATH/lib`.  If `PATH` is the string `yes`, then GMP is searched for in the default compiler search path; if `PATH` is `no` then GMP support is disabled.  If GMP is not found, then a wanrning message is printed and GMP support will be disabled, but compilation will continue.

  * `--with-mpi` / `--without-mpi`: Enable/disable MPI support.  MPI support requires the Boost.MPI library, see above.

  * `--with-ge` / `--without-ge`: If built with `--with-ge=yes`, the example programs will exit gracefully if they receive a SIGUSR2 (which is used by SGE/OGS to notify jobs of an impending kill).

This is a repository for implementations of the
["Rheinfall" parallel algorithm for Gaussian Elimination][1].

The sample implementation consists of a C++ library and
example programs:

  * The C++ library is a header-only library; it can run the Rheinfall
    algorithm in single-threaded, OpenMP or MPI programs.

  * The example programs compute the rank of an input matrix in
    [SMS format][2]. The same program is offered in C and C++
    versions.  The C++ version demoes the features available through
    the library; the C version is a standalone one and has no
    dependencies (save for the MPI library if you choose to compile
    with MPI support).

The most complete description of the Rheinfall algorithm to-date
can be found in the [arXiv paper 1105.4136][1].

If you want to try out Rheinfall on your own, please proceeed to the
[installation instructions](INSTALL.md).  Documentation is still quite
terse; feel free to [write me][3] with any question or feedback.


[1]: http://arxiv.org/abs/1105.4136
[2]: http://ljk.imag.fr/membres/Jean-Guillaume.Dumas/simc.html
[3]: mailto:riccardo.murri@gmail.com

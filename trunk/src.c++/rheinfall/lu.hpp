/**
 * @file   lu.hpp
 *
 * Interface of the rheinfall class for computing LU decompositions.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2010, 2011 riccardo.murri@gmail.com. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * See http://www.gnu.org/licenses/gpl.html for licence details.
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

#ifndef LU_HPP
#define LU_HPP 1


#include "config.hpp"
#include "stats.hpp"
#include "row.hpp"
#include "sparserow.hpp"
#include "denserow.hpp"
#include "types.hpp"

#ifdef _OPENMP
# include <omp.h>
# include <boost/optional.hpp>
#endif

#include <bits/allocator.h>
#include <iostream>
#include <list>
#include <map>
#include <numeric>
#include <vector>


namespace rheinfall {

  /** Implements the `Rheinfall` algorithm for the computation of a
      matrix LU decomposition. */
  template <typename val_t, typename coord_t = int, 
            template<typename T> class allocator = std::allocator >
  class LU {

public:

    /** Constructor. If compiled with MPI, will use the default
        communicator @c MPI_COMM_WORLD */
    LU(const coord_t ncols, 
       const val_t pivoting_threshold,
       const coord_t width = 1,
       const float dense_threshold = 40.0);
    
#ifdef WITH_MPI
    /** Constructor, taking explicit MPI communicator. */
    LU(mpi::communicator& comm, 
       const coord_t ncols,
       const val_t pivoting_threshold,
       const coord_t width = 1,
       const float dense_threshold = 40.0);
#endif // WITH_MPI

    /** Destructor. When using OpenMP and MPI serialization, destroys
        the @c mpi_send_lock_; otherwise, does nothing. */
    ~LU();
  
    typedef SparseRow<val_t,coord_t,allocator> sparserow_t;
    typedef DenseRow<val_t,coord_t,allocator> denserow_t;

    /** Read a matrix stream into the processors.  Return number of
        nonzero entries read.

        Each line of the input stream must consist of a three numbers
        separated by white space: row index, column index and the
        entry value.  Row- and column indices are 1-based; the entry
        value must be in a format that C++'s @c operator>> can map to
        the @c val_t type.

        No assumption on the order of entries in the input stream is
        made, therefore all entries have to be read into memory before
        @c Row objects can be assembled (and possibly dispatched to
        other MPI ranks, if using MPI and @a local_only is @c false ).

        If arguments @a nrows and @a ncols are null, then the first
        line read from stream @a input is assumed to be a SMS format
        header line, from which the number of rows and columns is
        extracted (and returned in @a nrows and @a ncols).  Otherwise,
        if @a nrows is not zero, no header line is read and @a nrows,
        @a ncols must contain the exact number of rows and columns in
        the matrix.
        
        Reading from the stream is terminated by EOF or after the SMS
        footer line (consisting of three 0's in a row) is read.  In
        any case, the @a input stream is not closed.
        
        If @c local_only is @c true (default), only rows assigned to the
        local MPI rank are retained and other are discarded.  If @c
        local_only is @c false, rows are sent to the destination MPI
        rank as soon as they are read; only one rank should do the I/O.
        
        If @c transpose is @c true, then row and column values read from
        the stream are exchanged, i.e., the transpose of the matrix is
        read into memory. 
    */
    coord_t read_triples(std::istream& input, coord_t& nrows, coord_t& ncols, 
                         const bool local_only=true, const bool transpose=false);
    
    /** Read a matrix stream in SMS format into the processors.
        Return number of nonzero entries read.
        
        The @c input stream should be in SMS format, see
        http://www-ljk.imag.fr/membres/Jean-Guillaume.Dumas/simc.html
        for details.  In particular, this implies that entries
        belonging to one row are not interleaved with entries from
        other rows; i.e., that it can consider reading row @a i
        complete as soon as it finds a row index @a j!=i. 
        
        If arguments @a nrows and @a ncols are null, then the first
        line read from stream @a input is assumed to be the SMS header
        line, from which the number of rows and columns is extracted
        (and written to @a nrows and @a ncols); otherwise, if @a nrows
        is not zero, no header line is read and @a nrows, @a ncols
        must contain the exact number of rows and columns in the
        matrix.

        Reading from the stream is terminated by EOF or after the SMS
        footer line (consisting of three 0's in a row) is read.  In
        any case, the @a input stream is not closed.
        
        If @c local_only is @c true (default), only rows assigned to the
        local MPI rank are retained and other are discarded.  If @c
        local_only is @c false, rows are sent to the destination MPI
        rank as soon as they are read; only one rank should do the I/O.
        
        If @c transpose is @c true, then row and column values read from
        the stream are exchanged, i.e., the transpose of the matrix is
        read into memory. */
    coord_t read_sms(std::istream& input, coord_t& nrows, coord_t& ncols, 
                     const bool local_only=true, const bool transpose=false);

    /** Generic type of a matrix row. */
    typedef Row<val_t,coord_t,allocator> row_t;

    /** Compute LU decomposition and store the result into @a l and @a
     u.  When using MPI, only the locally-computed rows are stored
     into @a l and @a u; it is the caller's responsibility to collect
     them all. */
    void lu(std::vector<row_t*>& l, std::vector<row_t*>& u);
    
    /** A row will switch its storage to dense format when the percent
        of nonzero entries w.r.t. total length exceeds this
        value. Default is to switch to dense storage at 40% fill-in. */
    float dense_threshold;

    /** Coefficient for threshold pivoting. */
    const val_t pivoting_threshold;

#ifdef RF_ENABLE_STATS
    /** Keep counts of various kinds of operations. */
    Stats stats;

    /** Copy the global collected statistics into @a global_stats. */
    void get_global_stats(Stats& global_stats) const;

    /** Copy the statistics of the local processing units into @a local_stats. 
        This is the same as taking a copy of the @c stats member attribute directly. */
    void get_local_stats(Stats& local_stats) const;
#endif


  protected:
    /** A single processing element. */
    class Processor {
      
    public:
      /** Constructor, taking parent instance and index of the columns
          to process. */
      Processor(LU& parent, const coord_t column);
      /** Destructor. Releases OpenMP lock and frees up remaining memory. */
      ~Processor();
      
    protected:     
      /** Parent instance. */
      LU<val_t,coord_t,allocator>& parent_;
      /** Index of matrix column to process. */
      const coord_t column_;
      friend class LU<val_t,coord_t,allocator>; // XXX: see rank()
      
      typedef LU<val_t,coord_t,allocator>::row_t row_t;
      typedef SparseRow<val_t,coord_t,allocator> sparserow_t;
      typedef DenseRow<val_t,coord_t,allocator> denserow_t;

      /** Row used for elimination */
      row_t* u0;
      /** Corresponding row in the L matrix. */
      row_t* l0;
      
      /** Processor state: it is @c running when Gaussian elimination
          operations are being carried out regularly; it turns @c ending
          when a @c TAG_END message is received; it is @c done when the
          final contribution to the rank has been computed and the @c
          step() method should be called no more. */
      enum { running, ending, done } phase;
      
      /** A block of rows. */
      typedef std::list< row_t*, 
                         allocator<row_t*> > row_block;
      /** The block of rows that will be eliminated next time @c step() is called. */
      row_block u_rows;
      /** The corresponding block of rows from the L matrix. */
      row_block l_rows;

      /** List of incoming rows from other processors. */
      row_block inbox_u, inbox_l;
#ifdef _OPENMP
      /** Lock for accessing the @c inbox */
      omp_lock_t inbox_lock_;
#endif

#ifdef WITH_MPI
      /** Type used for storing the list of rows sent to other
          processors and the associated MPI request. */
      // XXX: use boost::tuple ?
      typedef std::list< std::pair< mpi::request, std::pair<row_t*,row_t*> >,
                         allocator< std::pair< mpi::request, std::pair<row_t*,row_t*> > > > outbox_t;
      /** List of rows sent to other processors and the associated MPI request. 
          We need to keep track of these in order to free the resources when we're done. */
      outbox_t outbox;
#endif

#ifdef _OPENMP
    private:
      /** Lock for running the @c step() function. */
      omp_lock_t processing_lock_;
#endif

#ifdef RF_ENABLE_STATS
      Stats *stats_ptr;
#endif
      
    public: 

      /** Initialize a new matrix row with the given content. */
      void recv_row(row_t* new_row);

      /** Add the given pair to @c inbox. */
      void recv_pair(row_t* u, row_t* l);
      
      /** Make a pass over the block of rows and return either 0 or 1
          (depending whether elimination took place or not). */
      void step();
      
      /** Switch processor to @c ending state: after one final round of
          elimination, a @c TAG_END message will be sent to the next
          column VPU and no further VPU activity will follow: each
          invocation to @c step() will return always the same result. */
      void end_phase();
      
      /** Return @a true if the processor is in @c done state. */
      bool is_done() const;
    };


    /** Number of matrix columns. */
    const coord_t ncols_;         
    /** Type for grouping all local processors. */
    typedef std::vector< Processor*, allocator<Processor*> > vpu_array_t;
    /** Local processors. */
    vpu_array_t vpus;

#ifdef WITH_MPI
    mpi::communicator comm_; /**< MPI communicator. */
    const int me_;           /**< MPI rank. */
    const int nprocs_;       /**< Total number of ranks in MPI communicator. */
#else
    // simulate running MPI on 1 rank only
    static const int me_ = 0;
    static const int nprocs_ = 1;
#endif

    //
    // implement distribution of columns to MPI ranks
    //

    /** Number of consecutive columns forming each band assigned to this process. */
    const coord_t w_;

    /** Return @c true if the VPU processing column @a c is in the local @c vpus array. */
    bool       is_local(const coord_t c) const; 
    /** Return index of the @c Processor instance processing column @a c
        within local @a vpus array. */
    coord_t    vpu_index_for_column(const coord_t c) const;
    /** Return pointer to the @c Processor instance processing column @a c. */
    Processor* vpu_for_column(const coord_t c) const;
#ifdef WITH_MPI
    /** Return MPI rank where VPU assigned to column @a c resides. */
    int        owner(const coord_t c) const;
#endif 

    // 
    // communication primitives
    //

    /** MPI tags used to discriminate various types of messages. */
    enum { 
      TAG_END=0, 
      TAG_ROW_SPARSE=1, 
      TAG_ROW_DENSE=2,
      TAG_PAIR_SPARSE_SPARSE=3,
      TAG_PAIR_SPARSE_DENSE=4,
      TAG_PAIR_DENSE_SPARSE=5,
      TAG_PAIR_DENSE_DENSE=6
    };
    
    /** Send the termination signal to the processor owning the block
        starting at @c col. */
    void send_end(const Processor& origin, const coord_t column) const;

    /** Send the row data pointed by @c row to the processor owning the
        block starting at @c col. Frees @c row */
    void send_row(Processor& source, row_t* row);

    /** Send the row data pointed by @c row to the processor owning the
        block starting at @c col. Frees @c row */
    void send_pair(Processor& source, row_t* u, row_t* l);

#ifdef WITH_MPI
    /** Receive messages that have arrived during a computation cycle. */
    void do_receive(); 

# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    /** Lock used to to serialize all outgoing MPI calls
        (receives are already done by the master thread only) */
    // XXX: being an instance variable, this won't serialize
    // MPI calls if there are two concurrent `LU` operating...
    mutable omp_lock_t mpi_send_lock_;
# endif
#endif // WITH_MPI

  private:
    /// only used inside `read` to keep rows in memory until the whole file has been loaded
    typedef std::map< coord_t, SparseRow<val_t,coord_t,allocator>*,
                      std::less<coord_t>,
                      allocator< std::pair< coord_t, 
                                            SparseRow<val_t,coord_t,allocator>* > > > _rowmap_t;

};


//
// implementation
//

// ------- inline methods -------

  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  LU<val_t,coord_t,allocator>::
  LU(const coord_t ncols, 
     const val_t pivoting_threshold,
     const coord_t width,
     const float dense_threshold)
    : ncols_(ncols)
    , vpus()
    , w_(width)
#ifdef WITH_MPI
    , comm_(MPI_COMM_WORLD, mpi::comm_attach)
    , me_(comm_.rank())
    , nprocs_(comm_.size())
#endif
    , dense_threshold(dense_threshold)
    , pivoting_threshold(pivoting_threshold)
#ifdef RF_ENABLE_STATS
    , stats()
#endif
  {  
    // setup the array of processors
#ifdef WITH_MPI
    const int nmemb = width * (1 + ((ncols / width) / nprocs_));
#else 
    const int nmemb = ncols;
#endif
    vpus.reserve(nmemb);
    for (int c = 0; c < ncols_; ++c)
      if (is_local(c))
        vpus.push_back(new LU<val_t,coord_t,allocator>::Processor(*this, c));
    // internal check that the size of the data structure is actually correct
    assert(vpus.size() <= nmemb);
#ifdef RF_ENABLE_STATS
    for (typename vpu_array_t::iterator vpu = vpus.begin(); vpu != vpus.end(); ++vpu)
      (*vpu)->stats_ptr = &stats;
#endif
  };


#ifdef WITH_MPI
  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  LU<val_t,coord_t,allocator>::
  LU(mpi::communicator& comm, 
     const coord_t ncols, 
     const val_t pivoting_threshold,
     const coord_t width,
     const float dense_threshold)
  : ncols_(ncols)
  , vpus()
  , comm_(comm)
  , me_(comm.rank())
  , nprocs_(comm.size())
  , w_(width)
  , dense_threshold(dense_threshold)
  , pivoting_threshold(pivoting_threshold)
#ifdef RF_ENABLE_STATS
  , stats()
#endif
{  
  // setup the array of processors
  vpus.reserve(1 + ncols_ / nprocs_);
    for (int c = 0; c < ncols_; ++c)
      if (is_local(c))
        vpus.push_back(new LU<val_t,coord_t,allocator>::Processor(*this, c));
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_init_lock(& mpi_send_lock_); 
# endif // _OPENMP && WITH_MPI_SERIALIZED
#ifdef RF_ENABLE_STATS
    for (typename vpu_array_t::iterator vpu = vpus.begin(); vpu != vpus.end(); ++vpu)
      (*vpu)->stats_ptr = &stats;
#endif
};
#endif


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  LU<val_t,coord_t,allocator>::~LU()
  {
#if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_destroy_lock(& mpi_send_lock_); 
#endif // _OPENMP && WITH_MPI_SERIALIZED
  };

  
  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  coord_t
  LU<val_t,coord_t,allocator>::
  read_sms(std::istream& input,  
           coord_t& nrows, coord_t& ncols, 
           const bool local_only,
           const bool transpose) 
  {
    // count non-zero items read
    coord_t nnz = 0;
    
    // only read rows whose leading column falls in our domain
    SparseRow<val_t,coord_t,allocator>* row = 
      new SparseRow<val_t,coord_t,allocator>(ncols_ - 1);
    
    if (0 == nrows) {
      // read header
      char M;
      input >> nrows >> ncols >> M;
      if (input.fail() or 'M' != M)
        throw std::domain_error("Cannot read SMS header");
      if (transpose)
        std::swap(nrows, ncols);
    };
    
#ifdef WITH_MPI
    std::list< mpi::request, allocator<mpi::request> > outbox;
#endif
    coord_t last_row_index = -1;
    coord_t i, j;
    val_t value;
    while (not input.eof()) {
      input >> i >> j >> value;
      assert(0 <= i and i <= nrows);
      assert(0 <= j and j <= ncols);
      // ignore zero entries in matrix -- they shouldn't be here in the first place
      if (0 == value and 0 != i and 0 != j) // '0 0 0' is end-of-stream marker
        continue; 
      if (transpose)
        std::swap(i, j);
      // SMS indices are 1-based
      --i;
      --j;
      // if row index changes, then a new row is starting, so commit the old one.
      if (i != last_row_index) {
        // set initial row number
        row->h0 = last_row_index;
        row = row->adjust();
        if (NULL != row) {
          const coord_t starting_column = row->first_nonzero_column();
          assert (starting_column < ncols_);
          if (is_local(starting_column))
            vpu_for_column(starting_column)->recv_row(row);
          else {
#ifdef WITH_MPI
            if (not local_only) {
              mpi::request req;
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
              omp_set_lock(& mpi_send_lock_);
# endif 
              req = comm_.isend(owner(starting_column), TAG_ROW_SPARSE, row);
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
              omp_unset_lock(& mpi_send_lock_);
# endif
              outbox.push_back(req);
            }; // if not local_only
#endif // WITH_MPI
            // discard non-local and null rows
            delete row;
          }; // ! is_local(starting_column)
        };
        row = new SparseRow<val_t,coord_t,allocator>(ncols_ - 1);
        last_row_index = i;
      };
      if (-1 == i and -1 == j and 0 == value)
        break; // end of matrix stream
      ++nnz;
      row->set(j, value);
    }; // while (! eof)
    
#ifdef WITH_MPI
    // wait for all sent rows to arrive
    mpi::wait_all(outbox.begin(), outbox.end());
#endif
    return nnz;
  };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  coord_t
  LU<val_t,coord_t,allocator>::
  read_triples(std::istream& input,  coord_t& nrows, coord_t& ncols, 
               const bool local_only, const bool transpose)
  {
    if (0 == nrows) {
      // read header
      char M;
      input >> nrows >> ncols >> M;
      if (input.fail() or 'M' != M)
        throw std::domain_error("Cannot read SMS header");
      if (transpose)
        std::swap(nrows, ncols);
    };
    
    // count non-zero items read
    coord_t nnz = 0;
    
    // need to keep rows in memory until we reach end of file
    typedef std::pair< const coord_t,val_t > _coord_and_val;
    typedef std::map< coord_t, val_t, std::less<coord_t>, 
                      allocator<_coord_and_val> > _simplerow;
    typedef std::pair< const coord_t,_simplerow > _coord_and_simplerow;
    typedef std::map< coord_t, _simplerow, std::less<coord_t>, 
                      allocator<_coord_and_simplerow> > _simplerows;
    _simplerows m;

    coord_t i, j;
    val_t value;
    while (not input.eof()) {
      input >> i >> j >> value;
      if (0 == i and 0 == j and value == 0)
        break; // end of matrix stream
      if (transpose)
        std::swap(i, j);
      // ignore zero entries in matrix -- they shouldn't be there in the first place
      if (value == 0)
        continue; 
      ++nnz;
      assert(0 <= i and i <= nrows);
      assert(0 <= j and j <= ncols);
      // SMS indices are 1-based
      --i;
      --j;
      m[i][j] = value;
    }; // while(not eof)

#ifdef WITH_MPI
    std::list< mpi::request, allocator<mpi::request> > outbox;
#endif
    for (typename _simplerows::iterator it = m.begin(); 
         it != m.end(); 
         ++it) {
      if (it->second.begin() != it->second.end()) { // non-null matrix row
        SparseRow<val_t,coord_t,allocator>* row = 
          SparseRow<val_t,coord_t,allocator>::new_from_range(it->second.begin(), 
                                                             it->second.end(), 
                                                             ncols-1);
        if (NULL == row)
          continue; // with next `it`
        // set initial row number
        row->h0 = it->first;
        // commit row
        const coord_t starting_column = row->first_nonzero_column();
        if (is_local(starting_column))
          vpu_for_column(starting_column)->recv_row(row);
        else {
#ifdef WITH_MPI
          if (not local_only) {
            mpi::request req;
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
            omp_set_lock(& mpi_send_lock_);
# endif 
            req = comm_.isend(owner(starting_column), TAG_ROW_SPARSE, row);
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
            omp_unset_lock(& mpi_send_lock_);
# endif
            outbox.push_back(req);
          }; // if not local_only
#endif // WITH_MPI
          // discard non-local and null rows
          delete row;
        }; // ! is_local(starting_column)
      }; // it->second.begin() != it->second.end() 
    }; // for (it = m.begin(); ...)
    
#ifdef WITH_MPI
    // wait for all sent rows to arrive
    mpi::wait_all(outbox.begin(), outbox.end());
#endif
    return nnz;
  };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  void 
  LU<val_t,coord_t,allocator>::
  lu(std::vector<row_t*> &l, std::vector<row_t*> &u) 
  {
    // kickstart termination signal
    if (is_local(0)) {
      vpu_for_column(0)->end_phase();
    };

    coord_t n0 = 0;
    while(n0 < vpus.size()) {
      // n0 incremented each time a processor goes into `done` state 
      while (n0 < vpus.size() and vpus[n0]->is_done()) {
        //delete vpus[n0];
        //vpus[n0] = NULL; // faster SIGSEGV
        ++n0;
      };
      if (n0 >= vpus.size())
        break; // exit `while(n0 < vpus.size())` loop
#ifdef _OPENMP
# pragma omp parallel for schedule(static,w_)
#endif // _OPENMP
      for (coord_t n = n0; n < vpus.size(); ++n)
        vpus[n]->step();
#ifdef WITH_MPI
      while (boost::optional<mpi::status> status = comm_.iprobe()) { 
        switch(status->tag()) {
        case TAG_END: {
          // "end" message received; no new rows will be coming.
          // But some other rows could have arrived or could
          // already be in the `inbox`, so we need to make another
          // pass anyway.  All this boils down to: set a flag,
          // make another iteration, and end the loop next time.
          coord_t column = -1;
          comm_.recv(status->source(), status->tag(), column);
          assert(column >= 0 and is_local(column));
          Processor* vpu = vpu_for_column(column);
          assert(NULL != vpu);
          vpu->end_phase();
          break;
        };
        case TAG_ROW_SPARSE: {
          SparseRow<val_t,coord_t,allocator>* new_row = 
            SparseRow<val_t,coord_t,allocator>::new_from_mpi_message(comm_, status->source(), status->tag());
          assert(is_local(new_row->first_nonzero_column()));
          Processor* vpu = vpu_for_column(new_row->first_nonzero_column());
          assert(NULL != vpu);
          vpu->recv_row(new_row);
          break;
        };
        case TAG_ROW_DENSE: {
          DenseRow<val_t,coord_t,allocator>* new_row = 
            DenseRow<val_t,coord_t,allocator>::new_from_mpi_message(comm_, status->source(), status->tag());
          assert(is_local(new_row->first_nonzero_column()));
          Processor* vpu = vpu_for_column(new_row->first_nonzero_column());
          assert(NULL != vpu);
          vpu->recv_row(new_row);
          break;
        };
        case TAG_PAIR_SPARSE_SPARSE: {
          std::pair<sparserow_t*, sparserow_t*> payload;
          payload.first = new  SparseRow<val_t,coord_t,allocator>();
          payload.second = new  SparseRow<val_t,coord_t,allocator>();
          comm_.recv(status->source(), status->tag(), payload);
          assert(is_local(payload.first->first_nonzero_column()));
          Processor* vpu = vpu_for_column(payload.first->first_nonzero_column());
          assert(NULL != vpu);
          vpu->recv_pair(payload.first, payload.second);
          break;
        };
        case TAG_PAIR_SPARSE_DENSE: {
          std::pair<sparserow_t*, denserow_t*> payload;
          payload.first = new  SparseRow<val_t,coord_t,allocator>();
          payload.second = new  DenseRow<val_t,coord_t,allocator>();
          comm_.recv(status->source(), status->tag(), payload);
          assert(is_local(payload.first->first_nonzero_column()));
          Processor* vpu = vpu_for_column(payload.first->first_nonzero_column());
          assert(NULL != vpu);
          vpu->recv_pair(payload.first, payload.second);
          break;
        };
        case TAG_PAIR_DENSE_SPARSE: {
          std::pair<denserow_t*, sparserow_t*> payload;
          payload.first = new  DenseRow<val_t,coord_t,allocator>();
          payload.second = new  SparseRow<val_t,coord_t,allocator>();
          comm_.recv(status->source(), status->tag(), payload);
          assert(is_local(payload.first->first_nonzero_column()));
          Processor* vpu = vpu_for_column(payload.first->first_nonzero_column());
          assert(NULL != vpu);
          vpu->recv_pair(payload.first, payload.second);
          break;
        };
        case TAG_PAIR_DENSE_DENSE: {
          std::pair<denserow_t*, denserow_t*> payload;
          payload.first = new  DenseRow<val_t,coord_t,allocator>();
          payload.second = new  DenseRow<val_t,coord_t,allocator>();
          comm_.recv(status->source(), status->tag(), payload);
          assert(is_local(payload.first->first_nonzero_column()));
          Processor* vpu = vpu_for_column(payload.first->first_nonzero_column());
          assert(NULL != vpu);
          vpu->recv_pair(payload.first, payload.second);
          break;
        };
        }; // switch(status->tag())
      }; // while(iprobe)
#endif // WITH_MPI
    }; // end while(n0 < vpus.size())

#ifdef WITH_MPI
    // wait until all processors have done running
    comm_.barrier();
#endif // WITH_MPI

    // append local results to u, l
    for (coord_t n = 0; n < vpus.size(); ++n)
      u.push_back(vpus[n]->u0);
    for (coord_t n = 0; n < vpus.size(); ++n)
      l.push_back(vpus[n]->l0);

    return;
  };


#ifdef RF_ENABLE_STATS
  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline void
  LU<val_t,coord_t,allocator>::
  get_global_stats(Stats& global_stats) const 
  { 
# ifdef WITH_MPI
    // sum stats count from all processes
    mpi::all_reduce(comm_, stats.ops_count, global_stats.ops_count, std::plus<int>());
    mpi::all_reduce(comm_, stats.sparserow_count, global_stats.sparserow_count, std::plus<int>());
    mpi::all_reduce(comm_, stats.sparserow_elts, global_stats.sparserow_elts, std::plus<int>());
    mpi::all_reduce(comm_, stats.denserow_count, global_stats.denserow_count, std::plus<int>());
    mpi::all_reduce(comm_, stats.denserow_elts, global_stats.denserow_elts, std::plus<int>());
# else
    global_stats = stats;
# endif
  };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline void
  LU<val_t,coord_t,allocator>::
  get_local_stats(Stats& local_stats) const
  { 
    local_stats = stats;
  };
#endif // RF_ENABLE_STATS


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline bool 
  LU<val_t,coord_t,allocator>::
  is_local(const coord_t c) const 
  { 
#ifdef WITH_MPI
    return (owner(c) == me_); 
#else
    return true;
#endif
  };
  

  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline coord_t 
  LU<val_t,coord_t,allocator>::
  vpu_index_for_column(const coord_t c) const
  { 
    assert(c >= 0 and c < ncols_);
    return (w_ * ((c/w_) / nprocs_) + (c % w_));
  };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline typename LU<val_t,coord_t,allocator>::Processor*
  LU<val_t,coord_t,allocator>::
  vpu_for_column(const coord_t c) const
  { 
    return vpus[vpu_index_for_column(c)]; 
  };


#ifdef WITH_MPI
  template <typename val_t, typename coord_t, template<typename T> class allocator>
  inline int 
  LU<val_t,coord_t,allocator>::owner(const coord_t c) const
  { 
    assert(c >= 0 and c < ncols_);
    return ((c/w_) % nprocs_); 
  };
#endif


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline void
  LU<val_t,coord_t,allocator>::
  send_row(Processor& source, Row<val_t,coord_t,allocator>* row) 
{
  coord_t column = row->first_nonzero_column();
  assert(column >= 0 and column < ncols_);
  if (is_local(column)) {
    Processor* vpu = vpu_for_column(column);
    assert(NULL != vpu);
    vpu->recv_row(row);
  }
#ifdef WITH_MPI
  else { // ship to remote process
    mpi::request req;
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_set_lock(& mpi_send_lock_);
# endif
    if (row->kind == Row<val_t,coord_t,allocator>::sparse) 
      req = comm_.issend(owner(column), TAG_ROW_SPARSE, 
                        *(static_cast<SparseRow<val_t,coord_t,allocator>*>(row)));
    else if (row->kind == Row<val_t,coord_t,allocator>::dense)
      req = comm_.issend(owner(column), TAG_ROW_DENSE, 
                        *(static_cast<DenseRow<val_t,coord_t,allocator>*>(row)));
    else 
      // should not happen!
      throw std::logic_error("Unhandled row kind in LU::send_row()");
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_unset_lock(& mpi_send_lock_);
# endif
    source.outbox.push_back(std::make_pair(req, std::make_pair(row,NULL)));
  };
#endif // WITH_MPI
};


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline void
  LU<val_t,coord_t,allocator>::
  send_pair(Processor& source, row_t* u_row, row_t* l_row) 
{
  coord_t column = u_row->first_nonzero_column();
  assert(column >= 0 and column < ncols_);
  if (is_local(column)) {
    Processor* vpu = vpu_for_column(column);
    assert(NULL != vpu);
    vpu->recv_pair(u_row, l_row);
  }
#ifdef WITH_MPI
  else { // ship to remote process
    mpi::request req;
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_set_lock(& mpi_send_lock_);
# endif
    if (u_row->kind == row_t::sparse and l_row->kind == row_t::sparse) 
      req = comm_.issend(owner(column), TAG_PAIR_SPARSE_SPARSE, 
                         std::make_pair(*(static_cast<sparserow_t*>(u_row)),
                                        *(static_cast<sparserow_t*>(l_row))));
    else if (u_row->kind == row_t::sparse and l_row->kind == row_t::dense) 
      req = comm_.issend(owner(column), TAG_PAIR_SPARSE_DENSE, 
                         std::make_pair(*(static_cast<sparserow_t*>(u_row)),
                                        *(static_cast<denserow_t*>(l_row))));
    else if (u_row->kind == row_t::dense and l_row->kind == row_t::sparse) 
      req = comm_.issend(owner(column), TAG_PAIR_DENSE_SPARSE, 
                         std::make_pair(*(static_cast<denserow_t*>(u_row)),
                                        *(static_cast<sparserow_t*>(l_row))));
    else if (u_row->kind == row_t::dense and l_row->kind == row_t::sparse) 
      req = comm_.issend(owner(column), TAG_PAIR_DENSE_DENSE, 
                         std::make_pair(*(static_cast<denserow_t*>(u_row)),
                                        *(static_cast<denserow_t*>(l_row))));
    else 
      // should not happen!
      throw std::logic_error("Unhandled row kind in LU::send_pair()");
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_unset_lock(& mpi_send_lock_);
# endif
    source.outbox.push_back(std::make_pair(req, std::make_pair(u_row,l_row)));
  };
#endif // WITH_MPI
};


template <typename val_t, typename coord_t, 
          template<typename T> class allocator>
inline void 
LU<val_t,coord_t,allocator>::
send_end(Processor const& origin, const coord_t column) const
{
  if (column >= ncols_)
    return;

  if (is_local(column)) {
    Processor* vpu = vpu_for_column(column);
    assert(NULL != vpu);
    vpu->end_phase();
  }
#ifdef WITH_MPI
  else { // ship to remote process
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_set_lock(&mpi_send_lock_);
# endif
    comm_.send(owner(column), TAG_END, column);
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_unset_lock(&mpi_send_lock_);
# endif
  };
#endif
}; // send_end
  


// ------- Processor -------

  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  LU<val_t,coord_t,allocator>::Processor::
  Processor(LU<val_t,coord_t,allocator>& parent, 
            const coord_t column)
      : parent_(parent)
      , column_(column)
      , u0(NULL), l0(NULL)
      , phase(running)
      , u_rows(), l_rows()
      , inbox_u(), inbox_l()
#ifdef WITH_MPI
      , outbox()
#endif
    { 
#ifdef _OPENMP
      omp_init_lock(&inbox_lock_);
      omp_init_lock(&processing_lock_);
#endif
    };


template <typename val_t, typename coord_t, 
          template<typename T> class allocator>
LU<val_t,coord_t,allocator>::Processor::
~Processor()
{
#ifdef _OPENMP
      omp_destroy_lock(&inbox_lock_);
      omp_destroy_lock(&processing_lock_);
#endif
      assert(u_rows.empty());
      assert(l_rows.empty());
      assert(inbox_u.empty());
      assert(inbox_l.empty());
}


template <typename val_t, typename coord_t, 
          template<typename T> class allocator>
inline void 
LU<val_t,coord_t,allocator>::Processor::
recv_row(Row<val_t,coord_t,allocator>* new_u_row) 
  {
    assert(running == phase);
    assert((Row<val_t,coord_t,allocator>::sparse == new_u_row->kind 
            or Row<val_t,coord_t,allocator>::dense == new_u_row->kind));
    assert(column_ == new_u_row->starting_column_);
    // make a corresponding row from the identity matrix
    // note: the starting column in L rows is always 0, since
    // L is lower triangular
    SparseRow<val_t,coord_t,allocator>* new_l_row = 
      new SparseRow<val_t,coord_t,allocator>(0, new_u_row->ending_column_, 0);
    new_l_row->set(new_u_row->h0, 1);
    new_l_row->h0 = new_u_row->h0;
#ifdef RF_ENABLE_STATS
    if (NULL != this->stats_ptr) {
      new_u_row->stats_ptr->ops_count = this->stats_ptr->ops_count;
      new_l_row->stats_ptr->ops_count = this->stats_ptr->ops_count;
    }
#endif
#ifdef _OPENMP
    omp_set_lock(&inbox_lock_);
#endif
    inbox_u.push_back(new_u_row);
    inbox_l.push_back(new_l_row);
#ifdef _OPENMP
    omp_unset_lock(&inbox_lock_);
#endif
  };


template <typename val_t, typename coord_t, 
          template<typename T> class allocator>
inline void 
LU<val_t,coord_t,allocator>::Processor::
recv_pair(Row<val_t,coord_t,allocator>* new_u_row, 
          Row<val_t,coord_t,allocator>* new_l_row) 
  {
    assert(running == phase);
    assert((Row<val_t,coord_t,allocator>::sparse == new_u_row->kind 
            or Row<val_t,coord_t,allocator>::dense == new_u_row->kind));
    assert((Row<val_t,coord_t,allocator>::sparse == new_l_row->kind 
            or Row<val_t,coord_t,allocator>::dense == new_l_row->kind));
    assert(column_ == new_u_row->starting_column_);
    assert(0 == new_l_row->starting_column_);
#ifdef RF_ENABLE_STATS
    if (NULL != this->stats_ptr) {
      new_u_row->stats_ptr->ops_count = this->stats_ptr->ops_count;
      new_l_row->stats_ptr->ops_count = this->stats_ptr->ops_count;
    }
#endif
#ifdef _OPENMP
    omp_set_lock(&inbox_lock_);
#endif
    inbox_u.push_back(new_u_row);
    inbox_l.push_back(new_l_row);
#ifdef _OPENMP
    omp_unset_lock(&inbox_lock_);
#endif
  };


template <typename val_t, typename coord_t, 
          template<typename T> class allocator>
void
LU<val_t,coord_t,allocator>::Processor::
step() 
{
#ifdef _OPENMP
  omp_set_lock(&processing_lock_);
#endif

  if (is_done())
    return;

  // receive new rows
#ifdef _OPENMP
  omp_set_lock(&inbox_lock_);
#endif
  std::swap(inbox_u, u_rows);
  assert(inbox_u.empty());
  std::swap(inbox_l, l_rows);
  assert(inbox_l.empty());
#ifdef _OPENMP
  omp_unset_lock(&inbox_lock_);
#endif

  if (not u_rows.empty()) {
    // ensure there is one row for elimination
    if (NULL == u0) {
      u0 = u_rows.front();
      u_rows.pop_front();
      l0 = l_rows.front();
      l_rows.pop_front();
    }
    assert (NULL != u0);

    assert(u_rows.size() == l_rows.size());
#if 0
    // find the row with "best" pivot
    val_t best = u0->leading_term_;
    for (typename row_block::iterator it = u_rows.begin(); it != u_rows.end(); ++it) {
      if (first_is_better_pivot<val_t> ((*it)->leading_term_, best)) {
        best = (*it)->leading_term_;
      };
    }; // end for (it = rows.begin(); ...)

    // now choose the "more sparse" row among those whose leading term
    // is within a certain factor of the absolute best pivot
    // (threshold pivoting)
    val_t new_pivot = u0->leading_term_;
    double new_pivot_row_weight = u0->weight();
    boost::optional<typename row_block::iterator> new_u0_row_loc;
    boost::optional<typename row_block::iterator> new_l0_row_loc;
    typename row_block::iterator iter_u, iter_l;
    for (iter_u = u_rows.begin(), iter_l = l_rows.begin(); 
         iter_u != u_rows.end(); 
         ++iter_u, ++iter_l) 
      {
        if (not good_enough_pivot(best, parent_.pivoting_threshold, 
                                  (*iter_u)->leading_term_))
          continue;
        // pivot for sparsity and break ties by usual algorithm
        if (((*iter_u)->weight() < new_pivot_row_weight)
            or ((*iter_u)->weight() == new_pivot_row_weight 
                and first_is_better_pivot<val_t>((*iter_u)->leading_term_, new_pivot))) 
          {
            new_pivot = (*iter_u)->leading_term_;
            new_pivot_row_weight = (*iter_u)->weight();
            new_u0_row_loc = iter_u;
            new_l0_row_loc = iter_l;
          };
      }; // end for (it = rows.begin(); ...)
    if (!!new_u0_row_loc) 
      {
        assert(!!new_l0_row_loc);
        std::swap(u0, *(new_u0_row_loc.get()));
        std::swap(l0, *(new_l0_row_loc.get()));
      };
#endif // if 0

    typename row_block::iterator iter_u, iter_l;
    assert(u_rows.size() == l_rows.size());
    for (iter_u = u_rows.begin(), iter_l = l_rows.begin(); 
         iter_u != u_rows.end(); 
         ++iter_u, ++iter_l) 
      {
#if 1
        // pivot for sparsity and break ties by usual algorithm
        if (((*iter_u)->weight() < u0->weight())
            or ((*iter_u)->weight() == u0->weight()
                and first_is_better_pivot<val_t>((*iter_u)->leading_term_, u0->leading_term_)))
          {
            std::swap(u0, *iter_u);
            std::swap(l0, *iter_l);
          };
#endif // if 1

        // perform elimination -- return NULL in case resulting row is full of zeroes
        val_t a, b;
        get_row_multipliers<val_t>(u0->leading_term_, (*iter_u)->leading_term_, 
                                   a, b);
        row_t* new_u_row = 
          u0->linear_combination(*iter_u, a, b, true, parent_.dense_threshold);
        // ship reduced rows to other processors
        if (NULL != new_u_row) 
          {
            assert(new_u_row->starting_column_ > this->column_);
            *iter_l = l0->linear_combination(*iter_l, a, b, false);
            assert(NULL != *iter_l);
            parent_.send_pair(*this, new_u_row, (*iter_l));
          };
      }; // end for (iter_u = u_rows.begin(); ...)
    u_rows.clear();
    l_rows.clear();
  }; // end if not rows.empty()
  assert(u_rows.empty());
  assert(l_rows.empty());

  if (running == phase) {
#ifdef WITH_MPI
    if (outbox.size() > 0) {
      // check if some test messages have arrived and free
      // corresponding resources 
      // XXX: this is basically a rewrite of the Boost.MPI code for
      // mpi::test_some(), but I could find no way of associating a
      // request with the corresponding payload data...
      //outbox.erase(mpi::test_some(outbox.begin(), outbox.end()), outbox.end());
      typename outbox_t::iterator current(outbox.begin());
      typename outbox_t::iterator completed(outbox.end());
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
      omp_set_lock(&(parent_.mpi_send_lock_));
# endif
      while (current != completed) {
        if (!! current->first.test()) {
          --completed;
          std::iter_swap(current, completed);
          continue;
        }
        else
          ++current;
      }; // while current != completed
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
      omp_unset_lock(&(parent_.mpi_send_lock_));
# endif
      // now the range [completed, end) contains the completed requests
      for (current = completed; current != outbox.end(); ++current) {
        // free message payloads
        delete current->second.first;
        delete current->second.second;
      };
      // finally, erase all completed requests
      outbox.erase(completed, outbox.end());
    };
#endif
  }
  else { // `phase == ending`: end message already received
#ifdef WITH_MPI        
    if (outbox.size() > 0) {
      // wait until all sent messages have arrived
      std::vector< mpi::request, allocator<mpi::request> > reqs;
      reqs.reserve(outbox.size());
      for (typename outbox_t::iterator it = outbox.begin();
           it != outbox.end();
           ++it)
        reqs.push_back(it->first);
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
      omp_set_lock(&(parent_.mpi_send_lock_));
# endif
      mpi::wait_all(reqs.begin(), reqs.end());
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
      omp_unset_lock(&(parent_.mpi_send_lock_));
# endif
      // finally, free message payloads and erase all requests
      for (typename outbox_t::iterator it = outbox.begin();
           it != outbox.end();
           ++it) {
        // free message payloads
        delete it->second.first;
        delete it->second.second;
      };
      outbox.clear();
    };
#endif

    // pass end message along
    parent_.send_end(*this, column_ + 1);

    // all done
    phase = done;
  };

#ifdef _OPENMP
  omp_unset_lock(&processing_lock_);
#endif
  return;
}


template <typename val_t, typename coord_t, 
          template<typename T> class allocator>
inline void 
LU<val_t,coord_t,allocator>::Processor::
end_phase() 
{ 
  phase = ending; 
};


template <typename val_t, typename coord_t, 
          template<typename T> class allocator>
inline bool 
LU<val_t,coord_t,allocator>::Processor::
is_done() const 
{ 
  return (done == phase); 
};


}; // namespace rheinfall


#endif // LU_HPP

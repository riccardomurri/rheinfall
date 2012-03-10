/**
 * @file   rank.hpp
 *
 * Interface of the rheinfall class.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2010, 2011, 2012 riccardo.murri@gmail.com. All rights reserved.
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

#ifndef RANK_HPP
#define RANK_HPP 1


#include "config.hpp"
#include "stats.hpp"
#include "types.hpp"
#include "row.hpp"
#include "sparserow.hpp"
#include "denserow.hpp"


#include <boost/optional.hpp>

#ifdef RF_USE_TBB
# include <tbb/concurrent_vector.h>
# include <tbb/mutex.h>
# include <tbb/parallel_do.h>
# include <tbb/task.h>
#endif


#include <bits/allocator.h>
#include <iostream>
#include <list>
#include <map>
#include <numeric>
#include <vector>


namespace rheinfall {

  /** Implements the `Rheinfall` algorithm for the computation of matrix
      rank. */
  template <typename val_t, typename coord_t = int, 
            template<typename T> class allocator = std::allocator >
  class Rank {

public:

    /** Constructor. If compiled with MPI, will use the default
        communicator @c MPI_COMM_WORLD */
    Rank(const coord_t ncols, 
         const val_t pivoting_threshold,
         const coord_t width = 1,
         const float dense_threshold = 40.0);
    
#ifdef WITH_MPI
    /** Constructor, taking explicit MPI communicator. */
    Rank(mpi::communicator& comm, 
         const coord_t ncols,
         const val_t pivoting_threshold,
         const coord_t width = 1,
         const float dense_threshold = 40.0);
#endif // WITH_MPI

    /** Destructor. When using OpenMP and MPI serialization, destroys
        the @c mpi_send_mutex_; otherwise, does nothing. */
    ~Rank();
  

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
        
        If @c local_only is @c true (default), only rows assigned to
        the local MPI rank are retained and other are discarded.  If
        @c local_only is @c false, rows are sent to the destination
        MPI rank as soon as they are read; in this case, only one rank
        should do the I/O.
        
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
        other rows; i.e., the algorithm can consider reading row @a i
        done as soon as it finds a row index @a j!=i. 

        For the rest (including meaning of the arguments), this
        function behaves exacly as @c read_triples (which see).
    */
    coord_t read_sms(std::istream& input, coord_t& nrows, coord_t& ncols, 
                     const bool local_only=true, const bool transpose=false);

    /** Return rank of matrix after in-place destructive computation. */
    coord_t rank();
    
    /** A row will switch its storage to dense format when the percent
        of nonzero entries w.r.t. total length exceeds this
        value. Default is to switch to dense storage at 40% fill-in. */
    const float dense_threshold;

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
#endif // RF_ENABLE_STATS

  protected:
    /** A single processing element. */
    class Processor {
      
    public:
      /** Constructor, taking parent instance and index of the columns
          to process. */
      Processor(Rank& parent, const coord_t column);
      /** Destructor. Releases OpenMP lock and frees up remaining memory. */
      ~Processor();
      
    protected:     
      /** Parent instance. */
      Rank<val_t,coord_t,allocator>& parent_;
      /** Index of matrix column to process. */
      const coord_t column_;
      friend class Rank<val_t,coord_t,allocator>; // XXX: see rank()
      
      /** Row used for elimination */
      Row<val_t,coord_t,allocator>* u;
      
      /** Processor state: it is @c running when Gaussian elimination
          operations are being carried out regularly; it turns @c ending
          when a @c TAG_END message is received; it is @c done when the
          final contribution to the rank has been computed and the @c
          step() method should be called no more. */
      enum { running, ending, done } phase;
      
      /** A block of rows. */
#ifndef RF_USE_TBB
      typedef std::list< 
        Row<val_t,coord_t,allocator>*, 
        allocator<Row<val_t,coord_t,allocator>*> > row_block;
#else
      typedef tbb::concurrent_vector< 
        Row<val_t,coord_t,allocator>*, 
        allocator<Row<val_t,coord_t,allocator>*> > row_block;
#endif
      /** The block of rows that will be eliminated next time @c step() is called. */
      row_block rows;
      /** List of incoming rows from other processors. */
      row_block inbox;

#ifdef RF_USE_TBB
      /** Locked when a task to process the rows in @c inbox has been already spawned. */
      tbb::mutex processing_;
#endif

#ifdef WITH_MPI
      /** Type used for storing the list of rows sent to other
          processors and the associated MPI request. */
      typedef std::list< std::pair< mpi::request, Row<val_t,coord_t,allocator>* >,
                         allocator< std::pair< mpi::request, Row<val_t,coord_t,allocator>* > > > outbox_t;
      /** List of rows sent to other processors and the associated MPI request. 
          We need to keep track of these in order to free the resources when we're done. */
      outbox_t outbox;
#endif

      /** Stores the result of the last invocation of @c step(). */
      coord_t result_;

#ifdef RF_ENABLE_STATS
      Stats *stats_ptr;
#endif
      
    public: 

      /** Append the given row to @c inbox */
      void recv_row(Row<val_t,coord_t,allocator>* new_row);
      
      /** Make a pass over the block of rows and return either 0 or 1
          (depending whether elimination took place or not). */
      coord_t step();
      
      /** Switch processor to @c ending state: after one final round of
          elimination, a @c TAG_END message will be sent to the next
          column VPU and no further VPU activity will follow: each
          invocation to @c step() will return always the same result. */
      void end_phase();
      
      /** Return @a true if the processor is in @c done state. */
      bool is_done() const;
    }; // class Processor


#ifdef RF_USE_TBB
    /** Wrap a @c Processor::step() invocation so that it can be
        executed as a TBB task. */
    class ProcessorStepTask : public tbb::task {
    public:
      /** Constructor, taking pointer to a @c rheinfall::Processor instance. */
      ProcessorStepTask(Processor* const processor)
        : p_(processor) { }; 
      
      /** Copy constructor. */
      ProcessorStepTask(ProcessorStepTask& other)
        : p_(other->p_) { }; 
      
      /** Invoke @c Processor::step() and then release the @c
          processing_ lock on the @c Processor instance @c p_. */
      tbb::task* execute() 
      { p_->step(); lock_.release(); return NULL; };

    private:
      Processor* const p_;
      
      tbb::mutex::scoped_lock lock_;

      friend class Processor;
    }; // class ProcessorStepTask


    /** Wrap a @c Processor::step() invocation so that it can be
        executed in TBB @c parallel_for. */
    class ProcessorStepApply {
    public:
        /** Invoke @c Processor::step() */
        void operator() (Processor* p) const
        { p->step(); };
    }; // class ProcessorStepApply
#endif // RF_USE_TBB


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
    enum { TAG_END=0, TAG_ROW_SPARSE=1, TAG_ROW_DENSE=2 };
    
    /** Send the termination signal to the processor owning the block
        starting at @c col. */
    void send_end(const Processor& origin, const coord_t column) const;

    /** Send the row data pointed by @c row to the processor owning the
        block starting at @c col. Frees up @c row */
    void send_row(Processor& source, Row<val_t,coord_t,allocator>* row);

#ifdef WITH_MPI
    /** Receive messages that have arrived during a computation cycle. */
    void do_mpi_receive(); 

# if defined(RF_USE_TBB) and defined(WITH_MPI_SERIALIZED)
    /** Lock used to to serialize all outgoing MPI calls
        (receives are already done by the master thread only) */
    // XXX: being an instance variable, this won't serialize
    // MPI calls if there are two concurrent `Rank` operating...
    mutable tbb::mutex mpi_send_mutex_;
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
  Rank<val_t,coord_t,allocator>::
  Rank(const coord_t ncols, 
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
        vpus.push_back(new Rank<val_t,coord_t,allocator>::Processor(*this, c));
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
  Rank<val_t,coord_t,allocator>::
  Rank(mpi::communicator& comm, 
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
# ifdef RF_ENABLE_STATS
  , stats()
# endif
# if defined(RF_USE_TBB) and defined(WITH_MPI_SERIALIZED)
  , mpi_send_mutex_()
# endif // RF_USE_TBB && WITH_MPI_SERIALIZED
{  
  // setup the array of processors
  vpus.reserve(1 + ncols_ / nprocs_);
    for (int c = 0; c < ncols_; ++c)
      if (is_local(c))
        vpus.push_back(new Rank<val_t,coord_t,allocator>::Processor(*this, c));
# ifdef RF_ENABLE_STATS
    for (typename vpu_array_t::iterator vpu = vpus.begin(); vpu != vpus.end(); ++vpu)
      (*vpu)->stats_ptr = &stats;
# endif
};
#endif // WITH_MPI


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  Rank<val_t,coord_t,allocator>::~Rank()
  {
    // nothing to do
  };

  
  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  coord_t
  Rank<val_t,coord_t,allocator>::
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
        throw std::runtime_error("Cannot read SMS header");
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
      // SMS indices are 1-based
      --i;
      --j;
      if (transpose)
        std::swap(i, j);
      // check validity of the data we read; since '0 0 0' is the
      // end-of-stream marker, we have to case for it separately
      if (not (-1 == i and -1 == j and 0 == value)) { 
        if (not(0 <= i and i < nrows)) {
          std::ostringstream msg; 
          msg << "Invalid row index '" << (i+1) << "',"
              << " should be >0 and <" << nrows;
          throw std::runtime_error(msg.str());
        };
        if (not (0 <= j and j < ncols)) {
          std::ostringstream msg; 
          msg << "Invalid column index '" << (j+1) << "'"
              << " should be >0 and <" << ncols;
          throw std::runtime_error(msg.str());
        };
        // ignore zero entries in matrix -- they shouldn't be here in the first place
        if (0 == value) 
          continue; 
      };
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
# if defined(RF_USE_TBB) and defined(WITH_MPI_SERIALIZED)
              mpi_send_mutex_.lock();
# endif 
              req = comm_.isend(owner(starting_column), TAG_ROW_SPARSE, row);
# if defined(RF_USE_TBB) and defined(WITH_MPI_SERIALIZED)
              mpi_send_mutex_.unlock();
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
      row->set(j, value);
      ++nnz;
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
  Rank<val_t,coord_t,allocator>::
  read_triples(std::istream& input,  coord_t& nrows, coord_t& ncols, 
               const bool local_only, const bool transpose)
  {
    if (0 == nrows) {
      // read header
      char M;
      input >> nrows >> ncols >> M;
      if (input.fail() or 'M' != M)
        throw std::runtime_error("Cannot read SMS header");
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
      // SMS indices are 1-based
      --i;
      --j;
      // transpose if needed
      if (transpose)
        std::swap(i, j);
      // check validity of the data we read
      if (i < 0 or i >= nrows) {
        std::ostringstream msg; 
        msg << "Invalid row index '" << (i+1) << "',"
            << " should be >0 and <" << nrows;
        throw std::runtime_error(msg.str());
      };
        if (j < 0 or j >= ncols) {
        std::ostringstream msg; 
        msg << "Invalid column index '" << (j+1) << "'"
            << " should be >0 and <" << ncols;
        throw std::runtime_error(msg.str());
      };
      // ignore zero entries in matrix -- they shouldn't be here in the first place
      if (0 == value) 
        continue; 
      m[i][j] = value;
      ++nnz;
    }; // while(not eof)

#ifdef WITH_MPI
    std::list< mpi::request, allocator<mpi::request> > outbox;
#endif
    for (typename _simplerows::iterator it = m.begin(); 
         it != m.end(); 
         ++it) 
      {
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
# if defined(RF_USE_TBB) and defined(WITH_MPI_SERIALIZED)
              mpi_send_mutex_.lock();
# endif 
              req = comm_.isend(owner(starting_column), TAG_ROW_SPARSE, row);
# if defined(RF_USE_TBB) and defined(WITH_MPI_SERIALIZED)
              mpi_send_mutex_.unlock();
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
  coord_t 
  Rank<val_t,coord_t,allocator>::
  rank() 
  {
    // kickstart termination signal
    if (is_local(0)) {
      vpu_for_column(0)->end_phase();
    };

    // where to collect (partial) ranks
    std::vector<coord_t, allocator<coord_t> > r(vpus.size(), 0);

#ifdef RF_USE_TBB
    //tbb::parallel_do(vpus.begin(), vpus.end(), ProcessorStepApply());
#endif
  
    coord_t n0 = 0;
    while(n0 < vpus.size()) {
      for (coord_t n = n0; n < vpus.size(); ++n) {
#ifdef RF_USE_TBB
        // wait for the lock to be released by other processing tasks
        tbb::mutex::scoped_lock processing_lock(vpus[n]->processing_);
#endif
        r[n] = vpus[n]->step();
#ifdef RF_USE_TBB
        processing_lock.release();
#endif
      };
#ifdef WITH_MPI
      do_mpi_receive();
#endif
      // n0 incremented each time a processor goes into `done` state 
      while (n0 < vpus.size() and vpus[n0]->is_done()) {
        delete vpus[n0];
        vpus[n0] = NULL; // faster SIGSEGV
        ++n0;
      };
    } while(n0 < vpus.size());

    // the partial rank is computed as the sum of all ranks computed
    // by local processors
    coord_t local_rank = std::accumulate(r.begin(), r.end(), 0);

#ifdef WITH_MPI
    // wait until all processors have done running
    comm_.barrier();

    // collect the partial ranks for all processes
    coord_t rank = 0;
    mpi::all_reduce(comm_, local_rank, rank, std::plus<int>());
#else
    const coord_t rank = local_rank;
#endif

    return rank;
  };


#ifdef RF_ENABLE_STATS
  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline void
  Rank<val_t,coord_t,allocator>::
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
  Rank<val_t,coord_t,allocator>::
  get_local_stats(Stats& local_stats) const
  { 
    local_stats = stats;
  };
#endif // RF_ENABLE_STATS


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline bool 
  Rank<val_t,coord_t,allocator>::
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
  Rank<val_t,coord_t,allocator>::
  vpu_index_for_column(const coord_t c) const
  { 
    assert(c >= 0 and c < ncols_);
    return (w_ * ((c/w_) / nprocs_) + (c % w_));
  };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline typename Rank<val_t,coord_t,allocator>::Processor*
  Rank<val_t,coord_t,allocator>::
  vpu_for_column(const coord_t c) const
  { 
    return vpus[vpu_index_for_column(c)]; 
  };


#ifdef WITH_MPI
  template <typename val_t, typename coord_t, template<typename T> class allocator>
  inline int 
  Rank<val_t,coord_t,allocator>::owner(const coord_t c) const
  { 
    assert(c >= 0 and c < ncols_);
    return ((c/w_) % nprocs_); 
  };
#endif


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline void
  Rank<val_t,coord_t,allocator>::
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
# if defined(RF_USE_TBB) and defined(WITH_MPI_SERIALIZED)
    mpi_send_mutex_.lock();
# endif
    if (row->kind == Row<val_t,coord_t,allocator>::sparse) 
      req = comm_.issend(owner(column), TAG_ROW_SPARSE, 
                        *(static_cast<SparseRow<val_t,coord_t,allocator>*>(row)));
    else if (row->kind == Row<val_t,coord_t,allocator>::dense)
      req = comm_.issend(owner(column), TAG_ROW_DENSE, 
                        *(static_cast<DenseRow<val_t,coord_t,allocator>*>(row)));
    else 
      // should not happen!
      throw std::logic_error("Unhandled row kind in Rank::send_row()");
# if defined(RF_USE_TBB) and defined(WITH_MPI_SERIALIZED)
    mpi_send_mutex_.unlock();
# endif
    source.outbox.push_back(std::make_pair(req,row));
  };
#endif // WITH_MPI
};


template <typename val_t, typename coord_t, 
          template<typename T> class allocator>
inline void 
Rank<val_t,coord_t,allocator>::
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
# if defined(RF_USE_TBB) and defined(WITH_MPI_SERIALIZED)
    mpi_send_mutex_.lock();
# endif
    comm_.send(owner(column), TAG_END, column);
# if defined(RF_USE_TBB) and defined(WITH_MPI_SERIALIZED)
    mpi_send_mutex_.unlock();
# endif
  };
#endif // WITH_MPI
}; // send_end


#ifdef WITH_MPI
template <typename val_t, typename coord_t, 
          template<typename T> class allocator>
inline void 
Rank<val_t,coord_t,allocator>::
do_mpi_receive()
{
  while (boost::optional<mpi::status> status = comm_.iprobe()) { 
    switch(status->tag()) {
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
    }; // switch(status->tag())
  }; // while(iprobe)
}; // do_mpi_receive()
#endif // WITH_MPI


// ------- Processor -------

  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  Rank<val_t,coord_t,allocator>::Processor::
  Processor(Rank<val_t,coord_t,allocator>& parent, 
            const coord_t column)
      : parent_(parent)
      , column_(column)
      , u(NULL)
      , phase(running)
      , rows()
      , inbox()
#ifdef RF_USE_TBB
      , processing_()
#endif
#ifdef WITH_MPI
      , outbox()
#endif
      , result_(0)
    { };


template <typename val_t, typename coord_t, 
          template<typename T> class allocator>
Rank<val_t,coord_t,allocator>::Processor::
~Processor()
{
      if (NULL != u)
        delete u;
      assert(rows.empty());
      assert(inbox.empty());
}


template <typename val_t, typename coord_t, 
          template<typename T> class allocator>
inline void 
Rank<val_t,coord_t,allocator>::Processor::
recv_row(Row<val_t,coord_t,allocator>* new_row) 
  {
    assert(running == phase);
    assert((Row<val_t,coord_t,allocator>::sparse == new_row->kind 
            or Row<val_t,coord_t,allocator>::dense == new_row->kind));
    assert(column_ == new_row->starting_column_);
#ifdef RF_ENABLE_STATS
    if (NULL != this->stats_ptr)
      new_row->stats_ptr = this->stats_ptr;
#endif
    inbox.push_back(new_row);
#ifdef RF_USE_TBB
    ProcessorStepTask* t = new(tbb::task::allocate_root()) ProcessorStepTask(this);
    // if the lock is held, then a task is already scheduled for
    // stepping the same @c Processor instance.
    if (t->lock_.try_acquire(this->processing_))
      tbb::task::spawn(*t);
    else
      tbb::task::destroy(*t);
#endif
  };


template <typename val_t, typename coord_t, 
          template<typename T> class allocator>
coord_t
Rank<val_t,coord_t,allocator>::Processor::
step() 
{
  if (is_done())
    return result_;

  // receive new rows
  std::swap(inbox, rows);
#ifndef RF_USE_TBB
  // this is only true when executing sequentially
  assert(inbox.empty());
#endif

  bool skip_front = false;
  if (not rows.empty()) {
    // ensure there is one row for elimination
    if (NULL == u) {
      u = rows.front();
      // tbb::concurrent_vector has no erase operation, 
      // so let's just rebase the vector
      skip_front = true;
    }
    assert (NULL != u);

    typename row_block::iterator it;
#if (RF_PIVOT_STRATEGY == RF_PIVOT_THRESHOLD)
    // find the row with "best" pivot
    val_t best = u->leading_term_;
    it = rows.begin();
    if (skip_front)
      ++it;
    for (/* it */; it != rows.end(); ++it) {
      if (first_is_better_pivot<val_t> ((*it)->leading_term_, best)) {
        best = (*it)->leading_term_;
      };
    }; // end for (it = rows.cbegin(); ...)

    // now choose the "more sparse" row among those whose leading term
    // is within a certain factor of the absolute best pivot
    // (threshold pivoting)
    val_t new_pivot = u->leading_term_;
    double new_pivot_row_weight = u->weight();
    boost::optional<typename row_block::iterator> new_pivot_row_loc;
    it = rows.begin();
    if (skip_front)
      ++it;
    for (/* it */; it != rows.end(); ++it) {
      if (not good_enough_pivot(best, parent_.pivoting_threshold, (*it)->leading_term_))
        continue;
      // pivot for sparsity and break ties by usual algorithm
      if (((*it)->weight() < new_pivot_row_weight)
          or ((*it)->weight() == new_pivot_row_weight 
              and first_is_better_pivot<val_t>((*it)->leading_term_, new_pivot))) {
        new_pivot = (*it)->leading_term_;
        new_pivot_row_weight = (*it)->weight();
        new_pivot_row_loc = it;
      };
    }; // end for (it = rows.begin(); ...)
    if (!!new_pivot_row_loc) {
      std::swap(u, *(new_pivot_row_loc.get()));
    };
#endif // if RF_PIVOT_STRATEGY == RF_PIVOT_THRESHOLD

    // perform elimination -- return NULL in case resulting row is full of zeroes
    it = rows.begin();
    if (skip_front)
      ++it;
    for (/* it */; it != rows.end(); ++it) {
#if (RF_PIVOT_STRATEGY == RF_PIVOT_WEIGHT)
      // pivot for sparsity and break ties by usual algorithm
      if (((*it)->weight() < u->weight())
          or ((*it)->weight() == u->weight()
              and first_is_better_pivot<val_t>((*it)->leading_term_, u->leading_term_)))
        std::swap(u, *it);
#elif (RF_PIVOT_STRATEGY == RF_PIVOT_SPARSITY)
      // swap `u` and the new row if the new row is shorter, or has "better" leading term
      if (Row<val_t,coord_t,allocator>::sparse == (*it)->kind) {
        SparseRow<val_t,coord_t,allocator>* s = 
          static_cast<SparseRow<val_t,coord_t,allocator>*>(*it);
        if (Row<val_t,coord_t,allocator>::sparse == u->kind) {
          // if `*it` (new row) is shorter, it becomes the new pivot row
          if (s->size() < u->size())
            std::swap(u, *it);
        }
        else if (Row<val_t,coord_t,allocator>::dense == u->kind)
          std::swap(u, *it);
        else 
          assert(false); // forgot one kind in chained `if`s?
      }
      else if (Row<val_t,coord_t,allocator>::dense == (*it)->kind) {
        if (Row<val_t,coord_t,allocator>::dense == u->kind) {
          // swap `u` and `row` iff `row`'s leading term is "better"
          if (first_is_better_pivot<val_t>
              (static_cast<DenseRow<val_t,coord_t,allocator>*>(*it)->leading_term_,
               static_cast<DenseRow<val_t,coord_t,allocator>*>(u)->leading_term_))
            // swap `u` and `row`
            std::swap(u, *it);
        };
        // else, if `u` is a sparse row, no need to check further
      }
      else
        assert(false); // forgot a row kind in `if` chain?
#elif (RF_PIVOT_STRATEGY == RF_PIVOT_NONE)
      // no row exchanges
# elif (RF_PIVOT_STRATEGY == RF_PIVOT_THRESHOLD)
      // already done above
#else
# error "RF_PIVOT_STRATEGY should be one of: RF_PIVOT_NONE, RF_PIVOT_SPARSITY, RF_PIVOT_THRESHOLD, RF_PIVOT_WEIGHT"
#endif // if RF_PIVOT_STRATEGY == ...
      val_t a, b;
      get_row_multipliers<val_t>((u)->leading_term_, (*it)->leading_term_, 
                                 a, b);
      Row<val_t,coord_t,allocator>* new_row = 
        u->linear_combination(*it, a, b, true, parent_.dense_threshold);
      // ship reduced rows to other processors
      if (NULL != new_row) {
        assert(new_row->starting_column_ > this->column_);
        parent_.send_row(*this, new_row);
      };
    }; // end for (it = rows.begin(); ...)
    rows.clear();
    result_ = (NULL == u? 0 : 1);
  }; // end if not rows.empty()
  assert(rows.empty());

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
# if defined(RF_USE_TBB) and defined(WITH_MPI_SERIALIZED)
      parent_.mpi_send_mutex_.lock();
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
# if defined(RF_USE_TBB) and defined(WITH_MPI_SERIALIZED)
      parent_.mpi_send_mutex_.unlock();
# endif
      // now the range [completed, end) contains the completed requests
      for (current = completed; current != outbox.end(); ++current) {
        // free message payloads
        delete current->second;
      };
      // finally, erase all completed requests
      outbox.erase(completed, outbox.end());
    };
#endif
  }
  else { // `phase == ending`: end message already received
#ifdef WITH_MPI        
    if (outbox.size() > 0) {
      // wait untill all sent messages have arrived
      std::vector< mpi::request, allocator<mpi::request> > reqs;
      reqs.reserve(outbox.size());
      for (typename outbox_t::iterator it = outbox.begin();
           it != outbox.end();
           ++it)
        reqs.push_back(it->first);
# if defined(RF_USE_TBB) and defined(WITH_MPI_SERIALIZED)
      parent_.mpi_send_mutex_.lock();
# endif
      mpi::wait_all(reqs.begin(), reqs.end());
# if defined(RF_USE_TBB) and defined(WITH_MPI_SERIALIZED)
      mpi_send_mutex_.unlock();
# endif
      // finally, free message payloads and erase all requests
      for (typename outbox_t::iterator it = outbox.begin();
           it != outbox.end();
           ++it)
        delete it->second;
      outbox.clear();
    };
#endif

    // pass end message along
    parent_.send_end(*this, column_ + 1);

    // all done
    phase = done;
  };

  // exit with 0 only if we never processed any row
  return result_;
}


template <typename val_t, typename coord_t, 
          template<typename T> class allocator>
inline void 
Rank<val_t,coord_t,allocator>::Processor::
end_phase() 
{ 
  phase = ending; 
};


template <typename val_t, typename coord_t, 
          template<typename T> class allocator>
inline bool 
Rank<val_t,coord_t,allocator>::Processor::
is_done() const 
{ 
  return (done == phase); 
};


}; // namespace rheinfall


#endif // RANK_HPP

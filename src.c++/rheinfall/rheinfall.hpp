/**
 * @file   rheinfall.hpp
 *
 * Interface of the rheinfall class.
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

#ifndef RHEINFALL_HPP
#define RHEINFALL_HPP


#include "config.hpp"
#include "row.hpp"
#include "sparserow.hpp"
#include "denserow.hpp"
#include "types.hpp"

#ifdef _OPENMP
# include <omp.h>
# include <boost/optional.hpp>
#endif

#include <iostream>
#include <list>
#include <map>
#include <numeric>
#include <vector>


namespace rheinfall {

  /** Implements the `Rheinfall` algorithm for the computation of matrix
      rank. */
  template <typename val_t, typename coord_t = int>
  class Rheinfall {

public:

    /** Constructor. If compiled with MPI, will use the default
        communicator @c MPI_COMM_WORLD */
    Rheinfall(const coord_t ncols, 
              const coord_t width = 1,
              const float dense_threshold = 40.0);
    
#ifdef WITH_MPI
    /** Constructor, taking explicit MPI communicator. */
    Rheinfall(mpi::communicator& comm, 
              const coord_t ncols,
              const coord_t width = 1,
              const float dense_threshold = 40.0);
#endif // WITH_MPI

    /** Destructor. When using OpenMP and MPI serialization, destroys
        the @c mpi_send_lock_; otherwise, does nothing. */
    ~Rheinfall();
  

    /** Read a matrix stream into the processors. Does not make any
        assumption on the order of entries in the input stream,
        therefore all entries have to be read into memory. Return
        number of nonzero entries read.

        The @c input stream should be in SMS format, see
        http://www-ljk.imag.fr/membres/Jean-Guillaume.Dumas/simc.html
        for details.  However, this function allows matrix entries to
        appear in any order in the input stream (contrary to the SMS
        specification, which states that entries should be
        lexicographically sorted).
        
        If arguments @a nrows and @a ncols are null, then the first
        line read from stream @a input is assumed to be the header
        line, from which the number of rows and columns is extracted
        (and returned in @a nrows and @a ncols); otherwise, if @a
        nrows is not zero, no header line is read and @a nrows, @a
        ncols must contain the exact number of rows and columns in the
        matrix.
        
        If @c local_only is @c true (default), only rows assigned to the
        local MPI rank are retained and other are discarded.  If @c
        local_only is @c false, rows are sent to the destination MPI
        rank as soon as they are read; only one rank should do the I/O.
        
        If @c transpose is @c true, then row and column values read from
        the stream are exchanged, i.e., the transpose of the matrix is
        read into memory. 

        Finally, if @p ac is not @c NULL, every triple (row, column,
        value) read from the stream is passed to @c ac->process.  This
        could be used to compute, e.g., the matrix norm.
    */
    coord_t read(std::istream& input, coord_t& nrows, coord_t& ncols, 
                 const bool local_only=true, const bool transpose=false);
    
    /** Read a matrix stream into the
        processors. Assumes that entries belonging to one row are not
        interleaved with entries from other rows; i.e., that it can
        consider reading row @a i complete as soon as it finds a row
        index @a j!=i. Return number of nonzero entries read.
        
        The @c input stream should be in SMS format, see
        http://www-ljk.imag.fr/membres/Jean-Guillaume.Dumas/simc.html
        for details.
        
        If arguments @a nrows and @a ncols are null, then the first
        line read from stream @a input is assumed to be the SMS header
        line, from which the number of rows and columns is extracted
        (and returned in @a nrows and @a ncols); otherwise, if @a
        nrows is not zero, no header line is read and @a nrows, @a
        ncols must contain the exact number of rows and columns in the
        matrix.
        
        If @c transpose is @c true, then row and column values read from
        the stream are exchanged, i.e., the transpose of the matrix is
        read into memory. */
    coord_t read_noninterleaved(std::istream& input, coord_t& nrows, coord_t& ncols, 
                                const bool transpose=false);

    /** Return rank of matrix after in-place destructive computation. */
    coord_t rank();
    
    /** A row will switch its storage to dense format when the percent
        of nonzero entries w.r.t. total length exceeds this
        value. Default is to switch to dense storage at 40% fill-in. */
    float dense_threshold;


  protected:
    /** A single processing element. */
    class Processor {
      
    public:
      /** Constructor, taking parent instance and index of the columns
          to process. */
      Processor(Rheinfall& parent,
                const coord_t column);
      /** Destructor. Releases OpenMP lock and frees up remaining memory. */
      ~Processor();
      
    protected:     
      Rheinfall<val_t,coord_t>& parent_; /**< Parent instance. */
      const coord_t column_;             /**< Index of matrix column to process. */
      friend class Rheinfall<val_t,coord_t>; // XXX: see rank()
      
      /** Row used for elimination */
      Row<val_t,coord_t>* u;
      
      /** Processor state: it is @c running when Gaussian elimination
          operations are being carried out regularly; it turns @c ending
          when a @c TAG_END message is received; it is @c done when the
          final contribution to the rank has been computed and the @c
          step() method should be called no more. */
      enum { running, ending, done } phase;
      
      /** A block of rows. */
      typedef std::list< Row<val_t,coord_t>* > row_list;
      /** The block of rows that will be eliminated next time @c step() is called. */
      row_list rows;
      /** List of incoming rows from other processors. */
      row_list inbox;
#ifdef _OPENMP
      /** Lock for accessing the @c inbox */
      omp_lock_t inbox_lock_;
#endif

#ifdef WITH_MPI
      /** Type used for storing the list of rows sent to other
          processors and the associated MPI request. */
      typedef std::list< std::pair< mpi::request,Row<val_t,coord_t>* > > outbox_t;
      /** List of rows sent to other processors and the associated MPI request. 
          We need to keep track of these in order to free the resources when we're done. */
      outbox_t outbox;
#endif

      /** Stores the result of the last invocation of @c step(). */
      coord_t result_;

#ifdef _OPENMP
    private:
      /** Lock for running the @c step() function. */
      omp_lock_t processing_lock_;
#endif
      
    public: 

      /** Append the given row to @c inbox */
      void recv_row(Row<val_t,coord_t>* new_row);
      
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
    };


    const coord_t ncols_;         /**< Number of matrix columns. */
    std::vector<Processor*> vpus; /**< Local processors. */

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
    void send_row(Processor& source, Row<val_t,coord_t>* row);

#ifdef WITH_MPI
    /** Receive messages that have arrived during a computation cycle. */
    void do_receive(); 

# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    /** Lock used to to serialize all outgoing MPI calls
        (receives are already done by the master thread only) */
    // XXX: being an instance variable, this won't serialize
    // MPI calls if there are two concurrent `Rheinfall` operating...
    mutable omp_lock_t mpi_send_lock_;
# endif
#endif // WITH_MPI

  private:
    /// only used inside `read` to keep rows in memory until the whole file has been loaded
    typedef std::map< coord_t, SparseRow<val_t,coord_t>* > _rowmap_t;
};


//
// implementation
//

// ------- inline methods -------

  template <typename val_t, typename coord_t>
  Rheinfall<val_t,coord_t>::Rheinfall(const coord_t ncols, 
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
        vpus.push_back(new Processor(*this, c));
    // internal check that the size of the data structure is actually correct
    assert(vpus.size() <= nmemb);
  };


#ifdef WITH_MPI
  template <typename val_t, typename coord_t>
  Rheinfall<val_t,coord_t>::Rheinfall(mpi::communicator& comm, 
                                      const coord_t ncols, 
                                      const coord_t width,
                                      const float dense_threshold)
  : ncols_(ncols)
  , vpus()
  , comm_(comm)
  , me_(comm.rank())
  , nprocs_(comm.size())
  , w_(width)
  , dense_threshold(dense_threshold)
{  
  // setup the array of processors
  vpus.reserve(1 + ncols_ / nprocs_);
    for (int c = 0; c < ncols_; ++c)
      if (is_local(c))
        vpus.push_back(new Processor(*this, c));
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_init_lock(& mpi_send_lock_); 
# endif // _OPENMP && WITH_MPI_SERIALIZED
};
#endif


  template <typename val_t, typename coord_t>
  Rheinfall<val_t,coord_t>::~Rheinfall()
  {
#if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_destroy_lock(& mpi_send_lock_); 
#endif // _OPENMP && WITH_MPI_SERIALIZED
  };

  
  template <typename val_t, typename coord_t>
  coord_t
  Rheinfall<val_t,coord_t>::read_noninterleaved(std::istream& input,  
                                                coord_t& nrows, coord_t& ncols, 
                                                const bool transpose) 
  {
    // count non-zero items read
    coord_t nnz = 0;
    
    // only read rows whose leading column falls in our domain
    SparseRow<val_t,coord_t>* row = new SparseRow<val_t,coord_t>(ncols_ - 1);
    
    if (0 == nrows) {
      // read header
      char M;
      input >> nrows >> ncols >> M;
      if (input.fail() or 'M' != M)
        throw std::domain_error("Cannot read SMS header");
      if (transpose)
        std::swap(nrows, ncols);
    };
    
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
        row = row->adjust();
        if (NULL != row) {
          const coord_t starting_column = row->first_nonzero_column();
          if (starting_column < ncols_ and is_local(starting_column))
            vpu_for_column(starting_column)->recv_row(row);
          else
            // discard null rows and rows belonging to other MPI ranks
            delete row;
        };
        row = new SparseRow<val_t,coord_t>(ncols_ - 1);
        last_row_index = i;
      };
      if (-1 == i and -1 == j and 0 == value)
        break; // end of matrix stream
      ++nnz;
      row->set(j, value);
    }; // while (! eof)
    
    return nnz;
  };


  template <typename val_t, typename coord_t>
  coord_t
  Rheinfall<val_t,coord_t>::read(std::istream& input,  coord_t& nrows, coord_t& ncols, 
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
    std::map< coord_t, std::map< coord_t,val_t > > m;
    
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
    std::list< mpi::request > outbox;
#endif
    for (typename std::map< coord_t, std::map< coord_t,val_t > >::iterator it = m.begin(); 
         it != m.end(); 
         ++it) {
      if (it->second.begin() != it->second.end()) { // non-null matrix row
        SparseRow<val_t,coord_t>* row = 
          SparseRow<val_t,coord_t>::new_from_range(it->second.begin(), 
                                                   it->second.end(), 
                                                   ncols-1);
        if (NULL == row)
          continue; // with next `it`
        const coord_t starting_column = row->first_nonzero_column();
        // commit row
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


  template <typename val_t, typename coord_t>
  coord_t 
  Rheinfall<val_t,coord_t>::rank() 
  {
#ifdef _OPENMP
    boost::optional <coord_t> ending_next;
#endif

    // kickstart termination signal
    if (is_local(0)) {
      vpu_for_column(0)->end_phase();
#ifdef _OPENMP
      ending_next = 0;
#endif
    };

    // collect (partial) ranks
    std::vector<int> r(vpus.size(), 0);
  
    coord_t n0 = 0;
    while(n0 < vpus.size()) {
      // n0 incremented each time a processor goes into `done` state 
      while (n0 < vpus.size() and vpus[n0]->is_done()) {
        delete vpus[n0];
        vpus[n0] = NULL; // faster SIGSEGV
        ++n0;
      };
      if (n0 >= vpus.size())
        break; // exit `while(n0 < vpus.size())` loop
#ifdef _OPENMP
      size_t total_size, chunk_size, bottom, top;
#     pragma omp parallel private(total_size,chunk_size,bottom,top) default(shared)
      {
        if (!! ending_next) {
          assert(*ending_next == 0);
          if (omp_get_thread_num() == 0) {
            // master thread takes care of the ending ranks
            coord_t col = *ending_next;
            while (col < ncols_ and is_local(col)) {
              const coord_t n = vpu_index_for_column(col);
              r[n] = vpus[n]->step();
              col++;
            };
            if (col < ncols_) 
              n0 = vpu_index_for_column(col);
            else
              n0 = vpus.size();
            ending_next.reset();
          }
          else {
            // other threads loop over all VPUs with a non-empty inbox
            total_size = vpus.size() - n0;
            chunk_size = total_size / std::max(1, omp_get_num_threads()-1);
            bottom = n0 + chunk_size * (omp_get_thread_num() - 1);
            top = std::min(bottom + chunk_size, vpus.size());
            for (coord_t n = top-1; n >= bottom; --n) {
              // XXX: all this locking is definitely slowing us down,
              // but neither OpenMP nor Boost (as of 1.43) provide any
              // way to do concurrent-access data structures...
              omp_set_lock(&(vpus[n]->inbox_lock_));
              // FIXME: could gcc optimize this out of the lock/unlock pair?
              bool data_in_vpu_inbox = not (vpus[n]->inbox.empty());
              omp_unset_lock(&(vpus[n]->inbox_lock_));
              if (data_in_vpu_inbox) {
#               pragma omp task untied firstprivate(n)
                r[n] = vpus[n]->step();
              };
            };
#           pragma omp taskwait
            { /* no-op for taskwait */ };
          };
        } // end if (!! ending_next)
        else {
          // run parallel loop over all VPUs with a non-empty inbox
#         pragma omp for schedule(static)
          for (coord_t n = n0; n < vpus.size(); ++n)
            if (not vpus[n]->inbox.empty()) 
              r[n] = vpus[n]->step();
        }; // end else if (!! ending_next)
      }; // end omp parallel
#else // no OpenMP, loop over all VPUs
      for (coord_t n = n0; n < vpus.size(); ++n)
        r[n] = vpus[n]->step();
#endif // _OPENMP
#ifdef WITH_MPI
      while (boost::optional<mpi::status> status = comm_.iprobe()) { 
        switch(status->tag()) {
        case TAG_ROW_SPARSE: {
          SparseRow<val_t,coord_t>* new_row = 
            SparseRow<val_t,coord_t>::new_from_mpi_message(comm_, status->source(), status->tag());
          assert(is_local(new_row->first_nonzero_column()));
          Processor* vpu = vpu_for_column(new_row->first_nonzero_column());
          assert(NULL != vpu);
          vpu->recv_row(new_row);
          break;
        };
        case TAG_ROW_DENSE: {
          DenseRow<val_t,coord_t>* new_row = 
            DenseRow<val_t,coord_t>::new_from_mpi_message(comm_, status->source(), status->tag());
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
# ifdef _OPENMP
          ending_next = column;
# endif
          break;
        };
        }; // switch(status->tag())
      }; // while(iprobe)
#endif // WITH_MPI
    }; // end while(n0 < vpus.size())

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


  template <typename val_t, typename coord_t>
  inline bool 
  Rheinfall<val_t,coord_t>::is_local(const coord_t c) const 
  { 
#ifdef WITH_MPI
    return (owner(c) == me_); 
#else
    return true;
#endif
  };
  

  template <typename val_t, typename coord_t>
  inline coord_t 
  Rheinfall<val_t,coord_t>::vpu_index_for_column(const coord_t c) const
  { 
    assert(c >= 0 and c < ncols_);
    return (w_ * ((c/w_) / nprocs_) + (c % w_));
  };


  template <typename val_t, typename coord_t>
  inline typename Rheinfall<val_t,coord_t>::Processor*
  Rheinfall<val_t,coord_t>::vpu_for_column(const coord_t c) const
  { 
    return vpus[vpu_index_for_column(c)]; 
  };


#ifdef WITH_MPI
  template <typename val_t, typename coord_t>
  inline int 
  Rheinfall<val_t,coord_t>::owner(const coord_t c) const
  { 
    assert(c >= 0 and c < ncols_);
    return ((c/w_) % nprocs_); 
  };
#endif


  template <typename val_t, typename coord_t>
  inline void
  Rheinfall<val_t,coord_t>::send_row(Processor& source, Row<val_t,coord_t>* row) 
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
    if (row->kind == Row<val_t,coord_t>::sparse) 
      req = comm_.issend(owner(column), TAG_ROW_SPARSE, 
                        *(static_cast<SparseRow<val_t,coord_t>*>(row)));
    else if (row->kind == Row<val_t,coord_t>::dense)
      req = comm_.issend(owner(column), TAG_ROW_DENSE, 
                        *(static_cast<DenseRow<val_t,coord_t>*>(row)));
    else 
      // should not happen!
      throw std::logic_error("Unhandled row kind in Rheinfall::send_row()");
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_unset_lock(& mpi_send_lock_);
# endif
    source.outbox.push_back(std::make_pair(req,row));
  };
#endif // WITH_MPI
};


template <typename val_t, typename coord_t>
inline void 
Rheinfall<val_t,coord_t>::send_end(Processor const& origin, const coord_t column) const
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

  template <typename val_t, typename coord_t>
  Rheinfall<val_t,coord_t>::Processor::Processor(Rheinfall<val_t,coord_t>& parent, 
                                                 const coord_t column)
      : parent_(parent),
        column_(column), 
        u(NULL), 
        phase(running), 
        rows(), 
        inbox(), 
#ifdef WITH_MPI
        outbox(), 
#endif
        result_(0)
    { 
#ifdef _OPENMP
      omp_init_lock(&inbox_lock_);
      omp_init_lock(&processing_lock_);
#endif
    };


template <typename val_t, typename coord_t>
Rheinfall<val_t,coord_t>::Processor::~Processor()
{
#ifdef _OPENMP
      omp_destroy_lock(&inbox_lock_);
      omp_destroy_lock(&processing_lock_);
#endif
      if (NULL != u)
        delete u;
      assert(rows.empty());
      assert(inbox.empty());
}


  template <typename val_t, typename coord_t>
  inline void 
  Rheinfall<val_t,coord_t>::Processor::recv_row(Row<val_t,coord_t>* new_row) 
  {
    assert(running == phase);
    assert((Row<val_t,coord_t>::sparse == new_row->kind 
            or Row<val_t,coord_t>::dense == new_row->kind));
    assert(column_ == new_row->starting_column_);
#ifdef _OPENMP
    omp_set_lock(&inbox_lock_);
#endif
    inbox.push_back(new_row);
#ifdef _OPENMP
    omp_unset_lock(&inbox_lock_);
#endif
  };


template <typename val_t, typename coord_t>
coord_t
Rheinfall<val_t,coord_t>::Processor::step() 
{
#ifdef _OPENMP
  omp_set_lock(&processing_lock_);
#endif

  if (is_done())
    return result_;

  // receive new rows
#ifdef _OPENMP
  omp_set_lock(&inbox_lock_);
#endif
  std::swap(inbox, rows);
  assert(inbox.empty());
#ifdef _OPENMP
  omp_unset_lock(&inbox_lock_);
#endif

  if (not rows.empty()) {
    // ensure there is one row for elimination
    if (NULL == u) {
      u = rows.front();
      rows.pop_front();
    }

    assert (NULL != u);
    for (typename row_list::iterator it = rows.begin(); it != rows.end(); ++it) {
      // swap `u` and the new row if the new row is shorter, or has "better" leading term
      if (Row<val_t,coord_t>::sparse == (*it)->kind) {
        SparseRow<val_t,coord_t>* s = static_cast<SparseRow<val_t,coord_t>*>(*it);
        if (Row<val_t,coord_t>::sparse == u->kind) {
          // if `*it` (new row) is shorter, it becomes the new pivot row
          if (s->size() < u->size())
            std::swap(u, *it);
        }
        else if (Row<val_t,coord_t>::dense == u->kind)
          std::swap(u, *it);
        else 
          assert(false); // forgot one kind in chained `if`s?
      }
      else if (Row<val_t,coord_t>::dense == (*it)->kind) {
        if (Row<val_t,coord_t>::dense == u->kind) {
          // swap `u` and `row` iff `row`'s leading term is "better"
          if (first_is_better_pivot<val_t>
              (static_cast<DenseRow<val_t,coord_t>*>(*it)->leading_term_,
               static_cast<DenseRow<val_t,coord_t>*>(u)->leading_term_))
            // swap `u` and `row`
            std::swap(u, *it);
        };
        // else, if `u` is a sparse row, no need to check further
      }
      else
        assert(false); // forgot a row kind in `if` chain?

      // perform elimination -- return NULL in case resulting row is full of zeroes
      Row<val_t,coord_t>* new_row = u->gaussian_elimination(*it, parent_.dense_threshold);
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
      std::vector< mpi::request > reqs;
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

#ifdef _OPENMP
  omp_unset_lock(&processing_lock_);
#endif
  // exit with 0 only if we never processed any row
  return result_;
}


template <typename val_t, typename coord_t>
inline void 
Rheinfall<val_t,coord_t>::Processor::end_phase() 
{ 
  phase = ending; 
};


template <typename val_t, typename coord_t>
inline bool 
Rheinfall<val_t,coord_t>::Processor::is_done() const 
{ 
  return (done == phase); 
};


}; // namespace rheinfall


#endif // RHEINFALL_HPP

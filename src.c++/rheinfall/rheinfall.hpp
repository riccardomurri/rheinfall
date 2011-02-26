/**
 * @file   rheinfall.hpp
 *
 * Interface of the rheinfall class.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2010 riccardo.murri@gmail.com. All rights reserved.
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <list>
#include <map>
#include <vector>


namespace rheinfall {

  /** Implements the `Rheinfall` algorithm for the computation of matrix
      rank. */
  template <typename val_t, typename coord_t = int>
  class Rheinfall {

public:

  /** Constructor. If compiled with MPI, will use the default
      communicator @c MPI_COMM_WORLD */
    Rheinfall(const coord_t ncols, const double dense_threshold = 40.0);

#ifdef WITH_MPI
  /** Constructor, taking explicit MPI communicator. */
    Rheinfall(const coord_t ncols, mpi::communicator& comm, 
              const double dense_threshold = 40.0);
#endif // WITH_MPI

    /** Destructor. When using OpenMP and MPI serialization, destroys
        the @c mpi_send_lock_; otherwise, does nothing. */
    ~Rheinfall();
  

  /** Read a matrix stream into the processors. Does not make any
      assumption on the order of entries in the input stream,
      therefore all entries have to be read into memory. Return number
      of nonzero entries _read_.  

      If arguments @a nrows and @a ncols are null, then the first line
      read from stream @a input is assumed to be a header line, from
      which the number of rows and columns is extracted (and returned
      in @a nrows and @a ncols); otherwise, if @a nrows is not zero,
      no header line is read and @a nrows, @a ncols must contain the
      exact number of rows and columns in the matrix.

      If @c local_only is @c true (default), only rows assigned to the
      local MPI rank are retained and other are discarded.  If @c
      local_only is @c false, rows are sent to the destination MPI
      rank as soon as they are read; only one rank should do the I/O.

      If @c transpose is @c true, then row and column values read from
      the stream are exchanged, i.e., the transpose of the matrix is
      read into memory. */
    coord_t read(std::istream& input, coord_t& nrows, coord_t& ncols, 
                 const bool local_only=true, const bool transpose=false);

  /** Read a matrix stream into the processors. Assumes that entries
      belonging to one row are not interleaved with entries from other
      rows; i.e., that it can consider reading row @i i complete as
      soon as it finds a row index @i j != i. Return number of nonzero
      entries _read_. 

      If arguments @a nrows and @a ncols are null, then the first line
      read from stream @a input is assumed to be a header line, from
      which the number of rows and columns is extracted (and returned
      in @a nrows and @a ncols); otherwise, if @a nrows is not zero,
      no header line is read and @a nrows, @a ncols must contain the
      exact number of rows and columns in the matrix.

      If @c transpose is @c true, then row and column values read from
      the stream are exchanged, i.e., the transpose of the matrix is
      read into memory. */
    coord_t read_noninterleaved(std::istream& input, coord_t& nrows, coord_t& ncols, 
                                const bool transpose=false);

  /** Return rank of matrix after in-place destructive computation. */
  coord_t rank();

public:
  /** A single processing element. */
  class Processor {

  public:
    Processor(Rheinfall& parent,
              const coord_t column,
              const double dense_threshold = 40.0);
    /** Destructor. Releases OpenMP lock and frees up remaining memory. */
    ~Processor();

  protected:     
    Rheinfall<val_t,coord_t>& parent_;
    const coord_t column_;
    friend class Rheinfall<val_t,coord_t>; // XXX: see rank()

    /** Row used for elimination */
    Row<val_t,coord_t>* u;

    /** Processor state: it is `running` when Gaussian elimination
        operations are being carried out regularly; it turns @c ending
        when a @c TAG_END message is received; it is @c done when the
        final contribution to the rank has been computed and the @c
        step() method should be called no more. */
    enum { running, ending, done } phase;
    
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
    /** List of rows sent to other processors and the associated MPI request. 
        We need to keep track of these in order to free the resources when we're done. */
    typedef std::list< mpi::request > outbox_t;
    outbox_t outbox;
#endif

#ifdef _OPENMP
    /** Used by the adaptive scheduler to collect thread running time. */
    double runtime;
#endif

    /** A row will switch its storage to dense format when the percent of
        nonzero entries w.r.t. total length exceeds this value. */
    const double dense_threshold_;

  public: 

    /** Append the given row to @c inbox */
    void recv_row(Row<val_t,coord_t>* new_row);

    /** Make a pass over the block of rows and return either 0 or 1
        (depending whether elimination took place or not). */
    int step();

    /** Switch processor to @c ending state: after one final round of
        elimination, a @c TAG_END message will be sent to the next
        column VPU and no further VPU activity will follow: each
        invocation to @c step() will return always the same result. */
    void end_phase();

    /** Return @a true if the processor is in @c done state. */
    bool is_done() const;
  };


protected:
  const coord_t ncols_;
  std::vector<Processor*> vpus;   /**< Local processors. */

#ifdef WITH_MPI
    mpi::communicator comm_;
    const int me_;     /**< MPI rank. */
    const int nprocs_; /**< Total number of ranks in MPI communicator. */
#else
    // simulate running MPI on 1 rank only
    static const int me_ = 0;
    static const int nprocs_ = 1;
#endif

  //
  // implement distribution of columns to MPI ranks
  //
  
  /** Return @c true if the VPU processing column @a c is in the local @c vpus array. */
  bool       is_local(const coord_t c) const; 
  /** Return @c Processor instance processing column @a c. */
  Processor& local_owner(const coord_t c) const;
#ifdef WITH_MPI
  /** Return MPI rank where VPU assigned to column @a c resides. */
  int        remote_owner(const coord_t c) const;
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
    /// Lock used to to serialize all outgoing MPI calls
    /// (receives are already done by the master thread only)
    // XXX: being an instance variable, this won't serialize
    // MPI calls if there are two concurrent `Rheinfall` operating...
    mutable omp_lock_t mpi_send_lock_;
# endif
#endif // WITH_MPI

#ifdef _OPENMP
    mutable omp_lock_t stderr_lock_;
#endif

  private:
    /// only used inside `read` to keep rows in memory until the whole file has been loaded
    typedef std::map< coord_t, SparseRow<val_t,coord_t>* > _rowmap_t;
};


//
// implementation
//

// ------- inline methods -------

  template <typename val_t, typename coord_t>
  Rheinfall<val_t,coord_t>::Rheinfall(const coord_t ncols, const double dense_threshold)
    : ncols_(ncols)
    , vpus()
#ifdef WITH_MPI
    , comm_(MPI_COMM_WORLD, mpi::comm_attach)
    , me_(comm_.rank())
    , nprocs_(comm_.size())
#endif
  {  
    // setup the array of processors
    vpus.reserve(1 + ncols_ / nprocs_);
    for (int n = 0; (me_ + n * nprocs_) < ncols_; ++n)
      vpus.push_back(new Processor(*this, me_ + n * nprocs_, ncols_-1));
  };


#ifdef WITH_MPI
  template <typename val_t, typename coord_t>
  Rheinfall<val_t,coord_t>::Rheinfall(const coord_t ncols, mpi::communicator& comm, 
                       const double dense_threshold)
  : ncols_(ncols)
  , vpus()
  , comm_(comm)
  , me_(comm.rank())
  , nprocs_(comm.size())
{  
  // setup the array of processors
  vpus.reserve(1 + ncols_ / nprocs_);
  for (int n = 0; (me_ + n * nprocs_) < ncols_; ++n)
    vpus.push_back(new Processor(*this, me_ + n * nprocs_, ncols_-1));
#if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_init_lock(& mpi_send_lock_); 
#endif // _OPENMP && WITH_MPI_SERIALIZED
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
            local_owner(starting_column).recv_row(row);
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
      if (0 == i and 0 == j and 0 == value)
        break; // end of matrix stream
      if (transpose)
        std::swap(i, j);
      // ignore zero entries in matrix -- they shouldn't be there in the first place
      if (0 == value) 
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
        SparseRow<val_t,coord_t>* row = new SparseRow<val_t,coord_t>(it->second.begin(), it->second.end(), 
                                                                     ncols-1);
        const coord_t starting_column = row->first_nonzero_column();
        // commit row
        if (is_local(starting_column))
          local_owner(starting_column).recv_row(row);
        else {
#ifdef WITH_MPI
          if (not local_only) {
            mpi::request req;
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
            omp_set_lock(& mpi_send_lock_);
# endif 
            req = comm_.isend(remote_owner(starting_column), TAG_ROW_SPARSE, row);
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
    // kickstart termination signal
    if (is_local(0))
      local_owner(0).end_phase();

    // collect (partial) ranks
    std::vector<int> r(vpus.size(), 0);
    size_t n0 = 0;
    
    while(n0 < vpus.size()) {
      // n0 incremented each time a processor goes into `done` state 
      while (n0 < vpus.size() and vpus[n0]->is_done()) {
        delete vpus[n0];
        ++n0;
      };
      if (n0 >= vpus.size())
        break; // out of the while(n0 < vpus.size()) loop
#ifdef _OPENMP
# pragma omp parallel for schedule(guided)        
#endif // _OPENMP
      for (size_t n = n0; n < vpus.size(); ++n)
        r[n] = vpus[n]->step();
#ifdef WITH_MPI
      while (boost::optional<mpi::status> status = comm_.iprobe()) { 
        switch(status->tag()) {
        case TAG_ROW_SPARSE: {
          SparseRow<val_t,coord_t>* new_row = 
            SparseRow<val_t,coord_t>::new_from_mpi_message(comm_, status->source(), status->tag());
          assert(is_local(new_row->first_nonzero_column()));
          local_owner(new_row->first_nonzero_column()).recv_row(new_row);
          break;
        };
        case TAG_ROW_DENSE: {
          DenseRow<val_t,coord_t>* new_row = 
            DenseRow<val_t,coord_t>::new_from_mpi_message(comm_, status->source(), status->tag());
          assert(is_local(new_row->first_nonzero_column()));
          local_owner(new_row->first_nonzero_column()).recv_row(new_row);
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
          local_owner(column).end_phase();
          break;
        };
        }; // switch(status->tag())
      }; // while(iprobe)
#endif // WITH_MPI
    }; // while(n0 < vpus.size())

    // the partial rank is computed as the sum of all ranks computed
    // by local processors
    int local_rank = 0;
    for (size_t n = 0; n < r.size(); ++n)
      local_rank += r[n];

#ifdef WITH_MPI
    // wait until all processors have done running
    comm_.barrier();

    // collect the partial ranks for all processes
    int rank = 0;
    mpi::all_reduce(comm_, local_rank, rank, std::plus<int>());
#else
    int rank = local_rank;
#endif

    return rank;
  };


  template <typename val_t, typename coord_t>
  inline bool 
  Rheinfall<val_t,coord_t>::is_local(const coord_t c) const 
  { 
    assert(c >= 0 and c < ncols_);
    return (me_ == c % nprocs_); 
  };
  

  template <typename val_t, typename coord_t>
  inline typename Rheinfall<val_t,coord_t>::Processor& 
  Rheinfall<val_t,coord_t>::local_owner(const coord_t c) const
  { 
    assert(c >= 0 and c < ncols_);
    return *(vpus[c / nprocs_]); 
  };


#ifdef WITH_MPI
  template <typename val_t, typename coord_t>
  inline int 
  Rheinfall<val_t,coord_t>::remote_owner(const coord_t c) const
  { 
    assert(c >= 0 and c < ncols_);
    return (c % nprocs_); 
  };
#endif


  template <typename val_t, typename coord_t>
  inline void
  Rheinfall<val_t,coord_t>::send_row(Processor& source, Row<val_t,coord_t>* row) 
{
  coord_t column = row->first_nonzero_column();
  if (is_local(column))
    local_owner(column).recv_row(row);
#ifdef WITH_MPI
  else { // ship to remote process
    mpi::request req;
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_set_lock(& mpi_send_lock_);
# endif
    if (row->kind == Row<val_t,coord_t>::sparse)
      req = comm_.isend(remote_owner(column), TAG_ROW_SPARSE, 
                        *(static_cast<SparseRow<val_t,coord_t>*>(row)));
    else if (row->kind == Row<val_t,coord_t>::dense)
      req = comm_.isend(remote_owner(column), TAG_ROW_DENSE, 
                        *(static_cast<DenseRow<val_t,coord_t>*>(row)));
    else
      // should not happen!
      throw std::logic_error("Unhandled row type in Processor::send_row()");
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_unset_lock(& mpi_send_lock_);
# endif
    source.outbox.push_back(req);
    delete row; // no longer needed
  };
#endif // WITH_MPI
};


template <typename val_t, typename coord_t>
inline void 
Rheinfall<val_t,coord_t>::send_end(Processor const& origin, const coord_t column) const
{
  if (column >= ncols_)
    return;

  if (is_local(column))
    local_owner(column).end_phase();
#ifdef WITH_MPI
  else { // ship to remote process
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_set_lock(&mpi_send_lock_);
# endif
    comm_.send(remote_owner(column), TAG_END, column);
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
    omp_unset_lock(&mpi_send_lock_);
# endif
  };
#endif
}; // send_end
  


// ------- Processor -------

  template <typename val_t, typename coord_t>
  Rheinfall<val_t,coord_t>::Processor::Processor(Rheinfall<val_t,coord_t>& parent, 
                                                 const coord_t ending_column,
                                                 const double dense_threshold)
      : parent_(parent),
        column_(ending_column), 
        u(NULL), 
        phase(running), 
        rows(), 
        inbox(), 
#ifdef WITH_MPI
        outbox(), 
#endif
#ifdef _OPENMP
        runtime(0),
#endif
        dense_threshold_(dense_threshold)
    { 
#ifdef _OPENMP
      omp_init_lock(&inbox_lock_);
#endif
    };


template <typename val_t, typename coord_t>
Rheinfall<val_t,coord_t>::Processor::~Processor()
{
#ifdef _OPENMP
      omp_destroy_lock(&inbox_lock_);
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
            or Row<val_t,coord_t>::dense == new_row->kind)
           and (column_ == new_row->starting_column_));
#ifdef _OPENMP
    omp_set_lock(&inbox_lock_);
#endif
    inbox.push_back(new_row);
#ifdef _OPENMP
    omp_unset_lock(&inbox_lock_);
#endif
  };


template <typename val_t, typename coord_t>
int 
Rheinfall<val_t,coord_t>::Processor::step() 
{
  assert(not is_done()); // cannot be called once finished

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
      // swap `u` and the new row if the new row is shorter, or has < leading term
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
          // swap `u` and `row` iff `row`'s leading term is less
          if (std::abs(static_cast<DenseRow<val_t,coord_t>*>(*it)->leading_term_)
              < std::abs(static_cast<DenseRow<val_t,coord_t>*>(u)->leading_term_))
            // swap `u` and `row`
            std::swap(u, *it);
        };
        // else, if `u` is a sparse row, no need to check further
      }
      else
        assert(false); // forgot a row kind in `if` chain?

      // perform elimination -- return NULL in case resulting row is full of zeroes
      Row<val_t,coord_t>* new_row = u->gaussian_elimination(*it, dense_threshold_);
      // ship reduced rows to other processors
       if (NULL != new_row)
        parent_.send_row(*this, new_row);
    };
    rows.clear();
  };
  assert(rows.empty());

  if (running == phase) {
#ifdef WITH_MPI
    if (outbox.size() > 0) {
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
      omp_set_lock(&(parent_.mpi_send_lock_));
# endif
      // check if some test messages have arrived
      // and free corresponding resources
      outbox.erase(mpi::test_some(outbox.begin(), outbox.end()), outbox.end());
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
      omp_unset_lock(&(parent_.mpi_send_lock_));
# endif
    };
#endif
  }
  else { // `phase == ending`: end message already received
    // pass end message along
    parent_.send_end(*this, column_ + 1);

#ifdef WITH_MPI        
    if (outbox.size() > 0) {
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
      omp_set_lock(&(parent_.mpi_send_lock_));
# endif
      // wait untill all sent messages have arrived
      mpi::wait_all(outbox.begin(), outbox.end());
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
      omp_unset_lock(&(parent_.mpi_send_lock_));
# endif
    };
#endif
    // all done
    phase = done;
  };

  // exit with 0 only if we never processed any row
  return (NULL == u? 0 : 1);
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

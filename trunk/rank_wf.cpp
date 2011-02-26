/**
 * @file   rank_wf.cpp
 *
 * Sample implementation of the ``waterfall'' algorithm for computing matrix rank.
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

// TODO:
// === first working version ===
//   - have read_*() return the number of nonzero entries
//   - check Processor phases to optimize communication/computation overlap


#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <list>
#include <stdexcept>
#include <utility>
#include <vector>

#include <sys/time.h>
#include <sys/resource.h>

#include <boost/mpi.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
namespace mpi = boost::mpi;


enum { TAG_END=0, TAG_ROW_SPARSE=1, TAG_ROW_DENSE=2 };


typedef int coord_t;
typedef long long val_t;


/** Implements the `Waterfall` algorithm for distributed computation
    of matrix rank. */
class Waterfall {

public:
  /** Constructor, taking row and column dimension; uses the @c boost::mpi
      default communicator (@c MPI_COMM_WORLD ). */
  Waterfall(const coord_t nrows, const coord_t ncols);
  
  /** Read a matrix stream into the processors. Assumes that entries
      belonging to one row are not interleaved with entries from other
      rows; i.e., that it can consider reading row @i i complete as
      soon as it finds a row index @i j != i. Return number of nonzero
      entries _read_. */
  size_t read_noninterleaved(std::istream& input);

  /** Read a matrix stream into the processors. This should be run
      from one MPI rank only: it will then distribute the data to
      other MPI ranks. Return total number of nonzero entries read. */
  size_t read_and_distribute(std::istream& input) { throw "not implemented"; };

  /** Return rank of matrix after in-place destructive computation. */
  int rank();

  // forward declarations
protected:
  class Processor;
  

private:
  mpi::communicator comm;
  const int mpi_id; /**< MPI rank. */
  const int mpi_nprocs; /**< Total number of ranks in MPI communicator. */
  
  const coord_t nrows;
  const coord_t ncols;
  
  std::vector<Processor*> procs;   /**< Local processors. */

protected:
  // implement cyclic distribution of columns to MPI ranks
  bool is_local(const coord_t c) const; 
  Processor& local_owner(const coord_t c);
  int remote_owner(const coord_t c) const;

  /** A single processing element. */
  class Processor {

  public:
    Processor(Waterfall& parent_, 
              const coord_t starting_column_, const coord_t ending_column_,
              const double dense_threshold_ = 40.0) 
      : parent(parent_),
        starting_column(starting_column_), ending_column(ending_column_), rows(), 
        phase(running), inbox(), outbox(), 
        dense_threshold(dense_threshold_) 
    { };

  public:
    // forward declarations
    class Row;
    class SparseRow;
    class DenseRow;
    typedef Row* row_ptr;
    typedef SparseRow* sparse_row_ptr;
    typedef const SparseRow* const const_sparse_row_ptr;
    typedef DenseRow* dense_row_ptr;
    typedef const DenseRow* const const_dense_row_ptr;

  protected:     
    Waterfall& parent;
    const coord_t starting_column;
    const coord_t ending_column;

    /** The block of rows owned by this processor */
    typedef std::vector<row_ptr> block;
    block rows;

    /** Processor state: it is `running` when Gaussian elimination
        operations are being carried out regularly; it turns @c ending
        when a @c TAG_END message is received; it is @c done when the
        final contribution to the rank has been computed and the @c
        step() method should be called no more. */
    enum { running, ending, done } phase;
    
    /** List of incoming messages from other processors; each message
        consists of a message tag (an int) and (optionally) a rows. */
    typedef std::list< row_ptr > inbox_t;
    inbox_t inbox;

    /** List of rows sent to other processors and the associated MPI request. 
        We need to keep track of these in order to free the resources when we're done. */
    typedef std::list< mpi::request > outbox_t;
    outbox_t outbox;

    /** A row will switch its storage to dense format when the percent of
        nonzero entries w.r.t. total length exceeds this value. */
    const double dense_threshold;

  public:
    /** Abstract base class for matrix rows. */
    class Row 
    {
    public:
      typedef enum { sparse, dense } kind_t;
      const kind_t kind;
      /** Constructor. */
      Row(const kind_t kind_, const coord_t starting_column, const coord_t ending_column)
        : kind(kind_), starting_column_(starting_column), ending_column_(ending_column)
      { };
    public:
      /** Virtual destructor (for actual use in subclasses). */
      virtual ~Row();
      /** Return index of first column that has a nonzero entry. */
      virtual coord_t first_nonzero_column() const;
      /** Return number of allocated (i.e., non-zero) entries. */
      virtual size_t size() const = 0;
      /** Return a reference to the element stored at column @c col */
      virtual val_t& operator[](const coord_t col) = 0;
      /** Return a copy of the element stored at column @c col */
      virtual const val_t operator[](const coord_t col) const = 0;
      /** Perform Gaussian elimination: sum a multiple of this row to
          (a multiple of) row @c r so that the combination has its
          first nonzero entry at a column index strictly larger than
          the one of both rows. Return pointer to the combined row, which
          could possibly be this row if in-place update took place. */
      row_ptr gaussian_elimination(row_ptr other, const double dense_threshold = 40.0) const;
    protected:
      coord_t starting_column_; // would-be `const`: can only be modified by ctor and serialization
      coord_t ending_column_; // would-be `const`: can only be modified by ctor and serialization
      val_t leading_term_; // would-be `const`: can only be modified by ctor and serialization
    }; // class Row

    /** Store a matrix row as a vector of (column, entry) pairs. */
    class SparseRow : public Row
    {
    public:
      SparseRow(const coord_t starting_column, const coord_t ending_column) 
        : Row(sparse, starting_column, ending_column), storage() 
      { };
      /** Constructor initializing an invalid row; should call
          serialize immediately after!  */
      SparseRow() 
        : Row(sparse, -1 , -1), storage() 
      { };
      /** Return fill-in percentage, that is the number of actual
          nonzero entries divided by the number of potential entries.  */
      double fill_in() const { return 100.0 * size() / (ending_column_ - starting_column_ + 1); };
      /** Return number of allocated (i.e., non-zero) entries. */
      virtual size_t size() const;
      /** Return a reference to the element stored at column @c col */
      virtual val_t& operator[](const coord_t col);
      /** Return a copy of the element stored at column @c col */
      virtual const val_t operator[](const coord_t col) const;
      /** Perform Gaussian elimination: sum a multiple of this row to
          (a multiple of) row @c r so that the combination has its
          first nonzero entry at a column index strictly larger than
          the one of both rows. Return pointer to the combined row, which
          could possibly be this row if in-place update took place. */
      sparse_row_ptr gaussian_elimination(sparse_row_ptr other) const;
      /** Perform Gaussian elimination: sum a multiple of this row to
          (a multiple of) row @c r so that the combination has its
          first nonzero entry at a column index strictly larger than
          the one of both rows. Return pointer to the combined row, which
          could possibly be this row if in-place update took place. 
          Elimination on a dense row returns a dense row. */
      dense_row_ptr gaussian_elimination(dense_row_ptr other) const;
    protected:
      typedef std::vector< std::pair<coord_t,val_t> > storage_t;
      storage_t storage;
      friend class DenseRow;
      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive& ar, const unsigned int version);
    }; // class SparseRow

    /** Store a matrix row in a linear vector of consecutive entries. */
    class DenseRow : public Row
    {
    public:
      /** Constructor, copying entries from a @c SparseRow instance. */
      DenseRow(const sparse_row_ptr r);
      /** Construct a null row. */
      DenseRow(const coord_t starting_column_, const size_t ending_column_);
      /** Construct an invalid row; exists only for the purpose of
          calling @c serialize immediately after! */
      DenseRow(); // FIXME: make private?
      ~DenseRow();
      /** Return number of allocated entries. */
      virtual size_t size() const;
      /** Return a reference to the element stored at column @c col */
      virtual val_t& operator[](const coord_t col);
      /** Return a copy of the element stored at column @c col */
      virtual const val_t operator[](const coord_t col) const;
      dense_row_ptr gaussian_elimination(const sparse_row_ptr other) const;
      dense_row_ptr gaussian_elimination(const dense_row_ptr other) const;
    protected:
      /// Row length.  This is of course the number of elements in the @c storage array.
      size_t size_;
      /// Plain array of entries.
      val_t* storage; 
      friend class SparseRow;
      // serialization support (needed for MPI)
      friend class boost::serialization::access;
      template<class Archive> void serialize(Archive& ar, const unsigned int version);
      /** Adjust @c _size and @c starting_column to reflect the actual
          contents of the stored row. */
      size_t _resize();
    }; // class DenseRow

  public: 
    void recv_row(row_ptr new_row);
    /** Send the row data pointed by @c row to the processor owning the
        block starting at @c col. Frees up @c row */
    void send_row(row_ptr row);

    /** Make a pass over the block of rows and return either 0 or 1
        (depending whether elimination took place or not). */
    int step();

    void end_phase();

    bool is_done() const;

  protected:
    /** Send the termination signal to the processor owning the block
        starting at @c col. */
    void send_end(const coord_t column) const;
  };
};


//
// implementation
//

/* Needed arithmetic operations */
static inline
val_t _gcd(long long const n, long long const m) 
{ 
  return m==0? n : _gcd(m, n%m);
};

static inline
val_t gcd(long long const n, long long const m) 
{ 
  return _gcd(std::abs(n), std::abs(m)); 
};



// ------- Waterfall -------

Waterfall::Waterfall(const coord_t nrows, const coord_t ncols)
  : comm(mpi::communicator()), mpi_id(comm.rank()), mpi_nprocs(comm.size()),
    nrows(nrows), ncols(ncols),
    procs()
{  
  // setup the array of processors
  procs.reserve(1 + ncols / mpi_nprocs);
  for (int n = 0; n < 1 + (ncols / mpi_nprocs); ++n)
    procs.push_back(new Processor(*this, mpi_id + n * mpi_nprocs, ncols-1));
};


size_t
Waterfall::read_noninterleaved(std::istream& input) 
{
    // count non-zero items read
    size_t nnz = 0;

    // only read rows whose leading column falls in our domain
    std::vector< std::pair<coord_t,val_t> > row;

    coord_t current_row_index = -1;
    coord_t i, j;
    val_t value;
    while (not input.eof()) {
      input >> i >> j >> value;
      if (0 == i and 0 == j and 0 == value)
        break; // end of matrix stream
      ++nnz;
      // SMS indices are 1-based
      --i;
      --j;
      // new row?
      if (i != current_row_index) { // yes, commit old row
        // find starting column
        coord_t starting_column = ncols;
        for (std::vector< std::pair<coord_t,val_t> >::const_iterator it = row.begin();
             it != row.end();
             ++it)
          if (it->first < starting_column)
            starting_column = it->first;
        // commit row
        if (starting_column < ncols // otherwise, row is empty -- ignore it
            and is_local(starting_column)) {
          Processor::SparseRow* r = new Processor::SparseRow(starting_column, ncols-1);
          for (std::vector< std::pair<coord_t,val_t> >::const_iterator it = row.begin();
               it != row.end();
               ++it)
            {
              (*r)[it->first] = it->second;
            }; 
          local_owner(starting_column).recv_row(r);
        }; // if is_local(starting_column)
        row.clear();
        current_row_index = i;
      }; // if i != current_row_index
      row.push_back(std::make_pair(j, value));
    }; // while (! eof)
    
    return nnz;
};


int 
Waterfall::rank() 
{
  // kickstart termination signal
  if (is_local(0))
    local_owner(0).end_phase();
  // collect (partial) ranks
  std::vector<int> r(procs.size(), 0);
  size_t n0 = 0;
  while(true) {
    // n0 incremented each time a processor goes into `done` state 
    while (procs[n0]->is_done())
      ++n0;
    if (procs.size()-1 == n0)
      break;
    // step all processors
    for (size_t n = n0; n < procs.size(); ++n) {
      r[n] = procs[n]->step();
    };
    // manage arrival of messages
    while (boost::optional<mpi::status> status = comm.iprobe()) { 
      switch(status->tag()) {
      case TAG_ROW_SPARSE: {
        Processor::sparse_row_ptr new_row = new Processor::SparseRow();
        comm.recv(status->source(), status->tag(), new_row);
        assert(is_local(new_row->first_nonzero_column()));
        Processor::row_ptr new_row_ptr(new_row);
        local_owner(new_row->first_nonzero_column()).recv_row(new_row_ptr);
        break;
      };
      case TAG_ROW_DENSE: {
        Processor::dense_row_ptr new_row = new Processor::DenseRow();
        comm.recv(status->source(), status->tag(), 
                  *boost::get<Processor::DenseRow*>(new_row));
        assert(is_local(new_row->first_nonzero_column()));
        Processor::row_ptr new_row_ptr(new_row);
        local_owner(new_row->first_nonzero_column()).recv_row(new_row_ptr);
        break;
      };
      case TAG_END: {
        // "end" message received; no new rows will be coming.
        // But some other rows could have arrived or could
        // already be in the `inbox`, so we need to make another
        // pass anyway.  All this boils down to: set a flag,
        // make another iteration, and end the loop next time.
        coord_t column = -1;
        comm.recv(status->source(), status->tag(), column);
        assert(column >= 0 and is_local(column));
        local_owner(column).end_phase();
      };
      }; // switch(status->tag())
    };
  };

  // the partial rank is computed as the sum of all ranks computed
  // by local processors
  int partial_rank = 0;
  for (size_t n = 0; n < r.size(); ++n)
    partial_rank += r[n];
  
  // wait until all processors have done running
  comm.barrier();
  
  // collect the partial ranks for all processes
  int rank = 0;
  mpi::all_reduce(comm, partial_rank, rank, std::plus<int>());
  return rank;
};


inline bool 
Waterfall::is_local(const coord_t c) const 
{ 
  return (mpi_id == c % mpi_nprocs); 
};


inline Waterfall::Processor& 
Waterfall::local_owner(const coord_t c)
{ 
  return *(procs[c / mpi_nprocs]); 
};


inline int 
Waterfall::remote_owner(const coord_t c) const
{ 
  return (c % mpi_nprocs); 
};


// ------- Processor -------

inline void 
Waterfall::Processor::recv_row(row_ptr new_row) 
{
  // FIXME: requires locking for MT operation
  inbox.push_back(new_row);
  // DEBUG
  std::cerr << "DEBUG: Processor " <<this 
            << " received row"
            << " of size " << new_row->size()
            << " starting at column " << new_row->first_nonzero_column()
            << " (row_ptr " <<new_row<< ")"
            << std::endl;
};


inline void 
Waterfall::Processor::send_row(row_ptr row) 
{
  coord_t column = row->first_nonzero_column();
  if (parent.is_local(column))
    parent.local_owner(column).recv_row(row);
  else { // ship to remote process
    mpi::request req;
    if (row->kind == Row::sparse)
      req = parent.comm.isend(parent.remote_owner(column), TAG_ROW_SPARSE, row);
    else if (row->kind == Row::dense)
      req = parent.comm.isend(parent.remote_owner(column), TAG_ROW_DENSE, row);
    else
      // should not happen!
      throw std::logic_error("Unhandled row type in Processor::send_row()");
    outbox.push_back(req);
    delete row; // no longer needed
  };
};


int 
Waterfall::Processor::step() 
{
  assert(not is_done()); // cannot be called once finished
  std::cerr << "DEBUG: stepping Processor " << this << std::endl;

  // receive new rows
  for (inbox_t::iterator it = inbox.begin(); it != inbox.end(); ++it) {
    row_ptr new_row = *it;
    rows.push_back(new_row);
    if (new_row->size() < rows.front()->size())
      // less nonzero entries: we have a new leading row
      std::swap(rows.back(), rows.front());
  };
  inbox.clear();
          
  if (rows.size() > 0) {
    row_ptr first = rows.front();
    for (block::iterator second = rows.begin()+1; second < rows.end(); ++second) {
      // perform elimination
      row_ptr new_row = first->gaussian_elimination(*second, dense_threshold);
          
      // ship reduced rows to other processors
      if (new_row->size() > 0) {
        send_row(new_row);
      }
    };
    // just keep the "first" row
    rows.clear();
    rows.push_back(first);
  };

  if (running == phase) {
    if (outbox.size() > 0) {
      // check if some test messages have arrived
      // and free corresponding resources
      outbox.erase(mpi::test_some(outbox.begin(), outbox.end()), outbox.end());
    };
  }
  else { // `phase == ending`: end message already received
    // pass end message along
    if (starting_column < ending_column)
      send_end(starting_column + 1);
        
    if (outbox.size() > 0) {
      // wait untill all sent messages have arrived
      mpi::wait_all(outbox.begin(), outbox.end());
    };
        
    // all done
    phase = done;
  };

  // exit with 0 only if we never processed any row
  assert(rows.size() <= 1);
  return rows.size();
};


inline void 
Waterfall::Processor::end_phase() 
{ 
  phase = ending; 
};


inline bool 
Waterfall::Processor::is_done() const 
{ 
  return (done == phase); 
};


inline void 
Waterfall::Processor::send_end(const coord_t column) const
{
  if (parent.is_local(column))
    parent.local_owner(column).end_phase();
  else { // ship to remote process
    parent.comm.send(parent.remote_owner(column), TAG_END, column);
  };
}; // send_end



// ------- Row -------

Waterfall::Processor::Row::~Row()
{
  // override in sub-classes
};


inline coord_t 
Waterfall::Processor::Row::first_nonzero_column() const 
{
  assert(0 != leading_term_);
  return starting_column_;
};


inline Waterfall::Processor::row_ptr
Waterfall::Processor::Row::gaussian_elimination(row_ptr other, const double dense_threshold) const
{
  if (sparse == this->kind and sparse == other->kind) {
    const_sparse_row_ptr s1 = static_cast<const_sparse_row_ptr>(this);
    sparse_row_ptr s2 = static_cast<sparse_row_ptr>(other);
    // assume that `this` has already been checked for fill-in
    // in the main loop, so we won't attempt to convert it here;
    // just test `other` for too much fill-in
    if (s2->fill_in() > dense_threshold) {
      // FIXME: could merge ctor+elimination in one funcall
      dense_row_ptr dense_other(new DenseRow(s2));
      delete other;
      return s1->gaussian_elimination(dense_other);
    }
    else { // `other` kept sparse
    return s1->gaussian_elimination(s2);
    }; // if (fill_in > ...)
  }
  else if (sparse == this->kind and dense == other->kind) {
    const_sparse_row_ptr s = static_cast<const_sparse_row_ptr>(this);
    dense_row_ptr d = static_cast<dense_row_ptr>(other);
    return s->gaussian_elimination(d);
  }
  else if (dense == this->kind and sparse == other->kind) {
    const_dense_row_ptr d = static_cast<const_dense_row_ptr>(this);
    sparse_row_ptr s = static_cast<sparse_row_ptr>(other);
    return d->gaussian_elimination(s);
  }
  else if (dense == this->kind and dense == other->kind) {
    const_dense_row_ptr d1 = static_cast<const_dense_row_ptr>(this);
    dense_row_ptr d2 = static_cast<dense_row_ptr>(other);
    return d1->gaussian_elimination(d2);
  }
  else
    // should not happen!
    throw std::logic_error("Unhandled row type combination in Row::gaussian_elimination()");
};



// ------- SparseRow -------

Waterfall::Processor::sparse_row_ptr 
Waterfall::Processor::SparseRow::gaussian_elimination(sparse_row_ptr other) const
{
  assert(this->starting_column_ == other->starting_column_);

  // compute:
  //   `a`: multiplier for `this` row
  //   `b`: multiplier for `other` row
  const val_t GCD = gcd(leading_term_, other->leading_term_);
  val_t a = - (other->leading_term_ / GCD);
  val_t b = leading_term_ / GCD;

  sparse_row_ptr result = NULL; // XXX: use boost::optional<...> instead?

  storage_t::const_iterator this_i = this->storage.begin(); 
  storage_t::const_iterator other_i = other->storage.begin();
  // loop while one of the two indexes is still valid
  while(this_i != this->storage.end() or other_i != other->storage.end()) {
    coord_t this_col;
    if (this_i != this->storage.end()) {
      this_col = this_i->first;
      assert (this_col > this->starting_column_);
    }
    else // this_i reached end of vector, use out-of-range value
      this_col = this->ending_column_ + 1;

    coord_t other_col;
    if (other_i != other->storage.end()) {
      other_col = other_i->first;
      assert (other_col > other->starting_column_);
    }
    else 
      other_col = other->ending_column_ + 1;

    bool nonzero_found = false;
    coord_t coord;
    val_t entry;
    if (other_col < this_col) {
      entry = b * other_i->second;
      if (0 != entry) {
        coord = other_col;
        nonzero_found = true;
      };
      ++other_i;
    }
    else if (other_col == this_col) {
      entry = a * this_i->second + b * other_i->second;
      if (0 != entry) {
        coord = this_col;
        nonzero_found = true;
      };
      ++this_i;
      ++other_i;
    }
    else if (other_col > this_col) {
      entry = a * this_i->second;
      if (0 != entry) {
        coord = this_col;
        nonzero_found = true;
      };
      ++this_i;
    };
    
    if (nonzero_found) {
      // XXX: the following code makes assumptions about how the
      // storage is organised within `SparseRow`; it duplicates code
      // from SparseRow::operator[] for efficiency
      if (NULL == result) {
        // allocate new SparseRow
        result = new SparseRow(coord, ending_column_);
        result->leading_term_ = entry;
        result->storage.reserve(storage.size() + other->storage.size());
      }
      else {
        // FIXME: change storage order so that this is push_back()!
        result->storage.insert(result->storage.begin(), std::make_pair(coord, entry));
      }; 
    }; 
  }; // while(this_i < storage.rend() ...)
#ifndef NDEBUG
  if (NULL != result) {
    assert(result->starting_column_ > this->starting_column_);
    assert(result->starting_column_ > other->starting_column_);
    assert(result->ending_column_ == this->ending_column_);
    assert(result->size() <= this->size() + other->size());
  };
#endif
  delete other; // release old storage
  if (NULL == result) // row full of zeroes
    return new SparseRow(ending_column_, ending_column_);
  else 
    return result;
}; // sparse_row_ptr gaussian_elimination(sparse_row_ptr r)


Waterfall::Processor::dense_row_ptr 
Waterfall::Processor::SparseRow::gaussian_elimination(dense_row_ptr other) const
{
  // compute:
  //   `a`: multiplier for `this` row
  //   `b`: multiplier for `other` row
  const val_t GCD = gcd(leading_term_, other->leading_term_);
  val_t a = - (other->leading_term_ / GCD);
  val_t b = leading_term_ / GCD;

  for (size_t j = 0; j < other->size(); ++j) 
    other->storage[j] *= b;
  for (storage_t::const_iterator it = storage.begin();
       it != storage.end();
       ++it)
    {
      assert(it->first > starting_column_ and it->first <= ending_column_);
      other->storage[it->first] += a*it->second;
    };

  other->_resize(); // content update done, adjust size and starting column
  return other;
}; // dense_row_ptr gaussian_elimination(dense_row_ptr other)


inline val_t& 
Waterfall::Processor::SparseRow::operator[](const coord_t col) 
{
  assert(col >= starting_column_ and col <= ending_column_);
  if (col == starting_column_) 
    return leading_term_;
  // else, fast-forward to place where element is/would be stored
  size_t jj = 0;
  while (jj < storage.size() and storage[jj].first < col)
    ++jj;
  if (storage.size() == jj) {
    // end of list reached, `col` is larger than any index in this row
    storage.push_back(std::make_pair<size_t,val_t>(col,0));
    return storage.back().second;
  }
  else if (col == storage[jj].first) 
    return storage[jj].second;
  else { // storage[jj].first > j > storage[jj+1].first
    // insert new pair before `jj+1`
    storage.insert(storage.begin()+jj+1, 
                   std::make_pair<size_t,val_t>(col,0));
    return storage[jj+1].second;
  };
}; // SparseRow::operator[](...)


inline const val_t
Waterfall::Processor::SparseRow::operator[](const coord_t col) const
{
  assert(col >= starting_column_ and col <= ending_column_);
  if (col == starting_column_) 
    return leading_term_;
  if (storage.size() == 0)
    return 0;
  // else, fast-forward to place where element is/would be stored
  size_t jj = 0;
  while (jj < storage.size() and storage[jj].first < col)
    ++jj;
  if (storage.size() == jj) {
    // end of list reached, `col` is larger than any index in this row
    return 0;
  }
  else if (col == storage[jj].first) 
    return storage[jj].second;
  else { // storage[jj].first > j > storage[jj+1].first
    return 0;
  };
}; // SparseRow::operator[](...) const 


template<class Archive>
void 
Waterfall::Processor::SparseRow::serialize(Archive& ar, const unsigned int version) 
{
  assert(version == 0);
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  ar & starting_column_ & ending_column_ & leading_term_ & storage;
  assert(starting_column_ >= 0);
  assert(ending_column_ >= 0);
}; // SparseRow::serialize(...)


inline size_t 
Waterfall::Processor::SparseRow::size() const 
{ 
  return storage.size(); 
};



// ------- DenseRow -------

Waterfall::Processor::DenseRow::DenseRow(const sparse_row_ptr r) 
  : Row(dense, r->starting_column_, r->ending_column_),
    size_(ending_column_ - starting_column_ + 1),
    storage(new val_t[size_]) 
{ 
  assert(std::distance(r->storage.begin(), r->storage.end()) <= size_);
  std::fill_n(storage, size_, 0);
  for (SparseRow::storage_t::const_iterator it = r->storage.begin();
       it != r->storage.end();
       ++it) 
    {
      storage[it->first] = it->second;
    };
};


inline
Waterfall::Processor::DenseRow::DenseRow(const coord_t starting_column, 
                                         const size_t ending_column) 
  : Row(dense, starting_column, ending_column),
    size_(ending_column_ - starting_column_ + 1),
    storage(new val_t[size_]) 
{ 
  std::fill_n(storage, size_, 0);
};


inline
Waterfall::Processor::DenseRow::DenseRow() 
  : Row(dense, -1, -1),
    size_(0),
    storage(NULL) 
{ };


Waterfall::Processor::DenseRow::~DenseRow() 
{ 
  delete [] storage; 
};


inline Waterfall::Processor::dense_row_ptr 
Waterfall::Processor::DenseRow::gaussian_elimination(const sparse_row_ptr other) const
{
  // convert `other` to dense storage upfront: adding the
  // non-zero entries from `this` would made it pretty dense
  // anyway
  dense_row_ptr dense_other(new DenseRow(other));
  delete other;
  return this->gaussian_elimination(dense_other);
}; // dense_row_ptr gaussian_elimination(sparse_row_ptr other)


inline Waterfall::Processor::dense_row_ptr 
Waterfall::Processor::DenseRow::gaussian_elimination(const dense_row_ptr other) const
{
  assert(this->size() == other->size());
  // compute:
  //   `a`: multiplier for `this` row
  //   `b`: multiplier for `other` row
  const val_t GCD = gcd(leading_term_, other->leading_term_);
  val_t a = - (other->leading_term_ / GCD);
  val_t b = leading_term_ / GCD;

  for (size_t j = 0; j < size(); ++j) {
    // XXX: is it faster to allocate new storage and fill it with `a*x+b*y`?
    other->storage[j] *= a;
    other->storage[j] += b * storage[j];
  };
  other->_resize(); // update done, adjust size and string column
  return other;
}; // dense_row_ptr gaussian_elimination(dense_row_ptr other)


inline val_t& 
Waterfall::Processor::DenseRow::operator[](const coord_t col) 
{
  assert(col >= starting_column_ and col <= ending_column_);
  if (col == starting_column_)
    return leading_term_;
  else
    return storage[size_ - col - 1];
};


inline const val_t
Waterfall::Processor::DenseRow::operator[](const coord_t col) const
{
  assert(col >= starting_column_ and col <= ending_column_);
  if (col == starting_column_)
    return leading_term_;
  else
    return storage[size_ - col - 1];
};


inline size_t
Waterfall::Processor::DenseRow::_resize()
{
  // compute new starting column
  coord_t new_starting_column = ending_column_;
  assert(size_ > 0);
  for (size_t j = 0; j < size_; ++j)
    if ((*this)[j] != 0)
      new_starting_column = starting_column_ + j;

  size_ -= new_starting_column - starting_column_;
  starting_column_ = new_starting_column;
  return size_;
};


template<class Archive>
void 
Waterfall::Processor::DenseRow::serialize(Archive& ar, const unsigned int version) 
{
  assert(version == 0);
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  ar & starting_column_ & ending_column_ & leading_term_ & size_;
  assert(starting_column_ >= 0);
  assert(ending_column_ >= 0);
  assert(size_ >= 0);
  if (NULL == storage)
    storage = new val_t[size_];
  for (size_t j = 0; j < size_; ++j) {
    ar & storage[j];
  };
};


inline size_t 
Waterfall::Processor::DenseRow::size() const 
{ 
  return size_; 
};


//
// main
//

int
main(int argc, char** argv)
{
  if (1 == argc 
      or std::string("--help") == argv[1] 
      or std::string("-h") == argv[1]) 
    {
      std::cout << "Usage: " << argv[0] << " matrix_file.sms" << std::endl;
      return 1;
    }

  mpi::environment env(argc, argv);
  mpi::communicator world;
  int myid = world.rank();

  for (int i = 1; i < argc; ++i)
    {
      std::ifstream input(argv[i]);
      if (input.fail()) {
        std::cerr << "Cannot open file '" <<argv[i]<< "' for reading - skipping." 
                  << std::endl;
        continue;
      };

      if (0 == myid)
        std::cout << argv[0] << " file:"<<argv[i];

      // read matrix dimensions
      size_t rows, cols;
      char M;
      input >> std::skipws >> rows >> cols >> M;
      assert ('M' == M);
      if (0 == myid) {
        std::cout << " rows:" << rows;
        std::cout << " cols:" << cols;
      };

      Waterfall w(rows, cols);
      size_t nnz = w.read_noninterleaved(input);
      input.close();
      if (0 == myid)
        std::cout << " nonzero:" << nnz;

      struct rusage ru;
      getrusage(RUSAGE_SELF, &ru);
      struct timeval t0; memcpy(&t0, &ru.ru_utime, sizeof(struct timeval));

      if (0 == myid)
        std::cout << " rank:" << w.rank();

      getrusage(RUSAGE_SELF, &ru);
      struct timeval t1; memcpy(&t1, &ru.ru_utime, sizeof(struct timeval));
      struct timeval tdelta; timersub(&t1, &t0, &tdelta);
      double elapsed = tdelta.tv_sec + (tdelta.tv_usec / 1000000.0);
      if (0 == myid)
        std::cout << " time:" << std::fixed << std::setprecision(2) << elapsed << std::endl;
    }

  return 0;
}



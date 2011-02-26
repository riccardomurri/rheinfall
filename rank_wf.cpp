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
#include <stdexcept>
#include <iostream>
#include <list>
#include <utility>
#include <vector>

#include <boost/mpi.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/variant.hpp>
namespace mpi = boost::mpi;


enum { TAG_END=0, TAG_ROW_SPARSE=1, TAG_ROW_DENSE=2 };


typedef int coord_t;
typedef long long val_t;


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


/** Implements the `Waterfall` algorithm for distributed computation
    of matrix rank. */
class Waterfall {

public:
  /** Constructor, taking row and column dimension; uses the @c boost::mpi
      default communicator (@c MPI_COMM_WORLD ). */
  Waterfall(const coord_t nrows, const coord_t ncols)
    : comm(mpi::communicator()), mpi_id(comm.rank()), mpi_nprocs(comm.size()), procs(),
      nrows(nrows), ncols(ncols)
 {  };

  
  /** Read a matrix stream into the processors. Assumes that entries
      belonging to one row are not interleaved with entries from other
      rows; i.e., that it can consider reading row @i i complete as
      soon as it finds a row index @i j != i. Return number of nonzero
      entries _read_. */
  size_t read_noninterleaved(std::istream& input) {
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
        Processor::SparseRow* r = new Processor::SparseRow(starting_column, ncols);
        for (std::vector< std::pair<coord_t,val_t> >::const_iterator it = row.begin();
             it != row.end();
             ++it)
          {
            (*r)[it->first] = it->second;
          }; // for it = row.begin()
        boost::variant<Processor::SparseRow*,Processor::DenseRow*> r_(r);
        local_owner(starting_column).recv_row(r_);
        row.clear();
        current_row_index = i;
      }; // it i != current_row_index
    }; // while (! eof)

    return nnz;
  };

  /** Read a matrix stream into the processors. This should be run
      from one MPI rank only: it will then distribute the data to
      other MPI ranks. Return total number of nonzero entries read. */
  size_t read_and_distribute(std::istream& input) {
    throw "not implemented";
  };


  /** Return rank of matrix after in-place destructive computation. */
  int rank() {
    // setup the array of processors
    for (int n = 0; n < 1 + (nrows / mpi_nprocs); ++n)
      procs.push_back(Processor(*this, mpi_id + n * mpi_nprocs, ncols));
    std::vector<int> r(0);
    for (size_t n0 = 0; n0 < procs.size(); 
         /* n0 incremented each time a processor goes into `done` state */ ) {
      if (procs[n0].is_done())
        ++n0;
      // step all processors
      for (size_t n = n0; n < procs.size(); ++n) {
        r[n] = procs[n].step();
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


  // forward declarations
protected:
  class Processor;
  

private:
  mpi::communicator comm;
  const int mpi_id; /**< MPI rank. */
  const int mpi_nprocs; /**< Total number of ranks in MPI communicator. */
  
  std::vector<Processor> procs;   /**< Local processors. */

  const coord_t nrows;
  const coord_t ncols;
  
protected:
  // implement cyclic distribution of columns to MPI ranks
  bool is_local(const coord_t c) const { return (mpi_id == c % mpi_nprocs); };
  Processor& local_owner(const coord_t c) { return procs[c / mpi_nprocs]; };
  int remote_owner(const coord_t c) { return (c % mpi_nprocs); };

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
    class SparseRow;
    class DenseRow;
    typedef SparseRow* sparse_row_ptr;
    typedef DenseRow* dense_row_ptr;
    typedef boost::variant< sparse_row_ptr,dense_row_ptr > row_ptr;

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
    class SparseRow 
    {
    public:
      SparseRow(const coord_t starting_column_, const coord_t ending_column_) 
        : starting_column(starting_column_), ending_column(ending_column_), storage() 
      { };
      /** Constructor initializing an invalid row; should call
          serialize immediately after!  */
      SparseRow() 
        : starting_column(-1), ending_column(-1), storage() 
      { };
      /** Return index of first column that has a nonzero entry. */
      coord_t first_nonzero_column() const {
        assert(storage.size() > 0);
        return storage.back().first;
      };
      /** Return number of allocated (i.e., non-zero) entries. */
      size_t size() const { return storage.size(); };
      /** Return a reference to the element stored at column @c col */
      val_t& operator[](const coord_t col) {
#ifndef NDEBUG
        if (col < starting_column or col > ending_column)
          throw std::out_of_range("DenseRow::operator[] was passed a column index out of valid range");
#endif
        if (col == starting_column) 
          return leading_term;
        // else, fast-forward to place where element is/would be stored
        int jj = storage.size() - 1;
        while (jj >= 0 and storage[jj].first < col)
          --jj;
        if (0 > jj) {
          // end of list reached, `col` is larger than any index in this row
          storage.insert(storage.begin(), std::make_pair<size_t,val_t>(col,0));
          return storage[0].second;
        }
        else if (col == storage[jj].first) 
          return storage[jj].second;
        else { // storage[jj].first > j > storage[jj+1].first
          // insert new pair before `jj+1`
          storage.insert(storage.begin()+jj+1, 
                         std::make_pair<size_t,val_t>(col,0));
          return storage[jj+1].second;
        };
      }; // operator[](...)
      /** Perform Gaussian elimination: sum a multiple of this row and
          (a multiple of) row @c r so that the combination has its
          first nonzero entry at a column index strictly larger than
          the one of both rows. Return pointer to the combined row, which
          could possibly be this row if in-place update took place. */
      sparse_row_ptr gaussian_elimination(sparse_row_ptr other) {
        assert(this->starting_column == other->starting_column);

        sparse_row_ptr result = new SparseRow(starting_column, ending_column);
        // XXX: is it worth to compute the exact size here?
        result->storage.reserve(storage.size() + other->storage.size());

        // compute:
        //   `a`: multiplier for `this` row
        //   `b`: multiplier for `other` row
        const val_t GCD = gcd(leading_term, other->leading_term);
        val_t a = - (other->leading_term / GCD);
        val_t b = leading_term / GCD;

        size_t ix = 0;
        size_t iy = 0;
        // loop while one of the two indexes is still valid
        while(ix < storage.size() or iy < other->storage.size()) {
          coord_t jx;
          if (ix < storage.size()) {
            jx = storage[ix].first;
            assert (jx > starting_column);
          }
          else // ix reached end of vector
            jx = -1;

          coord_t jy;
          if (iy < other->storage.size()) {
            jy = other->storage[iy].first;
            assert (jy > other->starting_column);
          }
          else 
            jy = -1;

          if (jy > jx) {
            result->storage.push_back(std::make_pair(jy, b*other->storage[iy].second));
            ++iy;
          }
          else if (jy == jx) {
            result->storage.push_back(std::make_pair(jx, a*storage[ix].second 
                                                         + b*other->storage[iy].second));
            ++ix;
            ++iy;
          }
          else if (jy < jx) {
            result->storage.push_back(std::make_pair(jx, a*storage[ix].second));
            ++ix;
          };
        };
        delete other; // release old storage
        return result;
      }; // sparse_row_ptr gaussian_elimination(sparse_row_ptr r)
      /** Perform Gaussian elimination: sum a multiple of this row and
          (a multiple of) row @c r so that the combination has its
          first nonzero entry at a column index strictly larger than
          the one of both rows. Return pointer to the combined row, which
          could possibly be this row if in-place update took place. 
          Elimination on a dense row returns a dense row. */
      dense_row_ptr gaussian_elimination(dense_row_ptr other) {
        // compute:
        //   `a`: multiplier for `this` row
        //   `b`: multiplier for `other` row
        const val_t GCD = gcd(leading_term, other->leading_term);
        val_t a = - (other->leading_term / GCD);
        val_t b = leading_term / GCD;

        for (size_t j = 0; j < other->size(); ++j) 
          other->storage[j] *= b;

        for (storage_t::const_iterator it = storage.begin();
             it != storage.end();
             ++it)
          {
            assert(it->first > starting_column);
            other->storage[it->first] += a*it->second;
          };
        return other;
      }; // dense_row_ptr gaussian_elimination(dense_row_ptr other)
    protected:
      coord_t starting_column; // would-be `const`: can only be modified by ctor and serialization
      coord_t ending_column; // would-be `const`: can only be modified by ctor and serialization
      val_t leading_term; // would-be `const`: can only be modified by ctor and serialization
      typedef std::vector< std::pair<coord_t,val_t> > storage_t;
      storage_t storage;
      friend class DenseRow;
      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive& ar, const unsigned int version) {
        assert(version == 0);
        // When the class Archive corresponds to an output archive, the
        // & operator is defined similar to <<.  Likewise, when the class Archive
        // is a type of input archive the & operator is defined similar to >>.
        ar & starting_column & ending_column & leading_term & storage;
        assert(starting_column >= 0);
        assert(ending_column >= 0);
      }; // serialize(...)
    }; // class SparseRow

    class DenseRow
    {
    public:
      /** Constructor, copying entries from a @c SparseRow instance. */
      DenseRow(const sparse_row_ptr r) 
        : starting_column(r->starting_column), 
          ending_column(r->ending_column),
          _size(ending_column - starting_column + 1),
          storage(new val_t[_size]) { 
        assert(std::distance(r->storage.begin(), r->storage.end()) <= _size);
        std::fill_n(storage, _size, 0);
        for (SparseRow::storage_t::const_iterator it = r->storage.begin();
             it != r->storage.end();
             ++it) {
          storage[it->first] = it->second;
        };
      };
      /** Construct a null row. */
      DenseRow(const coord_t starting_column_, const size_t ending_column_) 
        : starting_column(starting_column_), 
          ending_column(ending_column_),
          _size(ending_column - starting_column + 1),
          storage(new val_t[_size]) { 
        std::fill_n(storage, _size, 0);
      };
      /** Construct an invalid row; exists only for the purpose of
          calling @c serialize immediately after! */
      DenseRow() 
        : starting_column(-1), 
          ending_column(-1),
          _size(0),
          storage(NULL) 
      { };
      ~DenseRow() { 
        delete [] storage;
      };
      /** Return index of first column that has a nonzero entry. */
      coord_t first_nonzero_column() const {
        assert(_size > 0);
        for (size_t j = 0; j < _size; ++j)
          if (storage[j] != 0)
            return starting_column + j;
        assert(false); // `first_nonzero_column` should not be called on zero rows
      };
      /** Return number of allocated entries. */
      size_t size() const { return _size; };
      /** Return a reference to the element stored at column @c col */
      val_t& operator[](const coord_t col) {
#ifndef NDEBUG
        if (col < starting_column or col > ending_column)
          throw std::out_of_range("DenseRow::operator[] was passed a column index out of valid range");
#endif
        if (col == starting_column)
          return leading_term;
        else
          return storage[col];
      };
      dense_row_ptr gaussian_elimination(const sparse_row_ptr other) {
        // convert `other` to dense storage upfront: adding the
        // non-zero entries from `this` would made it pretty dense
        // anyway
        dense_row_ptr dense_other(new DenseRow(other));
        delete other;
        return this->gaussian_elimination(dense_other);
      }; // dense_row_ptr gaussian_elimination(sparse_row_ptr other)
      dense_row_ptr gaussian_elimination(const dense_row_ptr other) {
        assert(this->size() == other->size());
        // compute:
        //   `a`: multiplier for `this` row
        //   `b`: multiplier for `other` row
        const val_t GCD = gcd(leading_term, other->leading_term);
        val_t a = - (other->leading_term / GCD);
        val_t b = leading_term / GCD;
        for (size_t j = 0; j < size(); ++j) {
          // XXX: is it faster to allocate new storage and fill it with `a*x+b*y`?
          other->storage[j] *= a;
          other->storage[j] += b * storage[j];
        };
        return other;
      }; // dense_row_ptr gaussian_elimination(dense_row_ptr other)
    protected:
      coord_t starting_column; 
      coord_t ending_column; 
      val_t leading_term; 
      /// Row length.  This is of course the number of elements in the @c storage array.
      size_t _size;
      /// Plain array of entries.
      val_t* storage; 
      friend class SparseRow;
      // serialization support (needed for MPI)
      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive& ar, const unsigned int version) {
        assert(version == 0);
        // When the class Archive corresponds to an output archive, the
        // & operator is defined similar to <<.  Likewise, when the class Archive
        // is a type of input archive the & operator is defined similar to >>.
        ar & starting_column & ending_column & leading_term & _size;
        assert(starting_column >= 0);
        assert(ending_column >= 0);
        assert(_size >= 0);
        if (NULL == storage)
          storage = new val_t[_size];
        for (size_t j = 0; j < _size; ++j) {
          ar & storage[j];
        };
      };
    }; // class DenseRow

  public: 
    void recv_row(row_ptr& new_row) {
      // FIXME: requires locking for MT operation
      inbox.push_back(new_row);
    };

    /** Send the row data pointed by @c row to the processor owning the
        block starting at @c col. Frees up @c row */
    void send_row(const coord_t column, row_ptr& row) {
      if (parent.is_local(column))
        parent.local_owner(column).recv_row(row);
      else { // ship to remote process
        class do_send : public boost::static_visitor<mpi::request> {
          Waterfall& wf;
          const coord_t column;
        public:
          do_send(Waterfall& wf_, const coord_t column_) : wf(wf_), column(column_) { };
          mpi::request operator() (sparse_row_ptr& r) {
            return wf.comm.isend(wf.remote_owner(column), TAG_ROW_SPARSE, r);
            delete r;
          };
          mpi::request operator() (dense_row_ptr& r) {
            return wf.comm.isend(wf.remote_owner(column), TAG_ROW_DENSE, r);
            delete r;
          };
        };
        mpi::request req = boost::apply_visitor(do_send(parent,column), row);
        outbox.push_back(req);
      };
    };

    /** Make a pass over the block of rows and return either 0 or 1
        (depending whether elimination took place or not). */
    int step() {
      assert(not is_done()); // cannot be called once finished

      // get row size with boost::variant
      class __row_size : public boost::static_visitor<size_t> {
      public:
        size_t operator() (const sparse_row_ptr& r) const {
          return r->size();
        };
        size_t operator() (const dense_row_ptr& r) const {
          return r->size();
        };
      }; // class __row_size
      class _row_size {
      public:
        size_t operator() (row_ptr& r) {
          return boost::apply_visitor(__row_size(), r);
        };
      } const row_size; 

      // get first nonzero column with boost::variant
      class __first_nonzero_column : public boost::static_visitor<coord_t> {
      public:
        coord_t operator() (sparse_row_ptr& r) const {
          return r->first_nonzero_column();
        };
        coord_t operator() (dense_row_ptr& r) const {
          return r->first_nonzero_column();
        };
      }; // class __first_nonzero_column
      class _first_nonzero_column {
      public:
        coord_t operator() (row_ptr& r) {
          return boost::apply_visitor(__first_nonzero_column(), r);
        };
      } first_nonzero_column;

      // receive new rows
      for (inbox_t::iterator it = inbox.begin(); it != inbox.end(); ++it) {
        row_ptr new_row = *it;
        rows.push_back(new_row);
        if (row_size(new_row) < row_size(rows[0]))
          // less nonzero entries: we have a new leading row
          std::swap(rows[0], rows[rows.size()-1]);
      };
      inbox.clear();
          
      if (rows.size() > 0) {
        row_ptr first = rows[0];
        for (block::iterator second = rows.begin()+1; second != rows.end(); ++second) {
          // perform elimination
          row_ptr new_row = gaussian_elimination(first, second);
          
          // ship reduced rows to other processors
          if (row_size(new_row) > 0) {
            send_row(first_nonzero_column(new_row), new_row);
          };
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
        send_end(starting_column + 1);
        
        if (outbox.size() > 0) {
          // wait untill all sent messages have arrived
          mpi::wait_all(outbox.begin(), outbox.end());
        };
        
        // all done, exit with 0 only if we never processed any row
        phase = done;
      };
      
      assert(rows.size() <= 1);
      return rows.size();
    };

    void end_phase() { phase = ending; };

    bool is_done() const { return (done == phase); };

  protected:
    /** Send the termination signal to the processor owning the block
        starting at @c col. */
    void send_end(const coord_t column) {
      if (parent.is_local(column))
        parent.local_owner(column).end_phase();
      else { // ship to remote process
        parent.comm.send(parent.remote_owner(column), TAG_END, column);
      };
    }; // send_end


    /** Perform Gaussian elimination on row `second`, using the entry
        in the leading column of `first` as a pivot. */
    row_ptr gaussian_elimination(row_ptr first, row_ptr second) {
      class _ge_op
        : public boost::static_visitor<row_ptr>
      {
        const coord_t starting_column;
        const coord_t ending_column;
        const double dense_threshold;
      public:
        _ge_op(const coord_t starting_column_, const coord_t ending_column_, 
               const double dense_threshold_)
          : starting_column(starting_column_), ending_column(ending_column_),
            dense_threshold(dense_threshold_)
        { };
        row_ptr operator() (sparse_row_ptr first, sparse_row_ptr second) const {
          // assume that `first` has already been checked for fill-in
          // in the main loop, so we won't attempt to convert it here;
          // just test `second` for too much fill-in
          const double fill_in = 100.0 * second->size()
                                 / (ending_column - starting_column + 1);
          if (fill_in > dense_threshold) {
            // FIXME: could merge ctor+elimination in one funcall
            dense_row_ptr dense_second(new DenseRow(second));
            delete second;
            return first->gaussian_elimination(dense_second);
          }
          else {
            return first->gaussian_elimination(second);
          }; // if (fill_in > dense_threshold)
        };

        row_ptr operator() (dense_row_ptr first, sparse_row_ptr second) const {
          return first->gaussian_elimination(second);
        };

        row_ptr operator() (sparse_row_ptr first, dense_row_ptr second) const {
          return first->gaussian_elimination(second);
        };

        row_ptr operator() (dense_row_ptr first, dense_row_ptr second) const {
          return first->gaussian_elimination(second);
        };
        } ge_op(starting_column, ending_column);
      return boost::apply_visitor(ge_op, first, second);
    }; // gaussian_elimination
  };
};

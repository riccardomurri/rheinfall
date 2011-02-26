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
//   - serialization support for SparseRow / DenseRow
//   - add a way to actually populate SparseRow / DenseRow with data
//   - store leading_{column,term}
//   - have read_*() return the number of nonzero entries
//   - check Processor phases to optimize communication/computation overlap


#include <cassert>
#include <iostream>
#include <list>
#include <utility>
#include <vector>

#include <boost/mpi.hpp>
#include <boost/optional.hpp>
#include <boost/variant.hpp>
namespace mpi = boost::mpi;

enum { TAG_END=0, TAG_ROW_SPARSE=1, TAG_ROW_DENSE=2 };

typedef int coord_t;
typedef long long val_t;


/* Needed arithmetic operations */
static inline
val_t gcd(long long const n, long long const m) 
{ 
  return _gcd(std::abs(n), std::abs(m)); 
}

static inline
val_t _gcd(long long const n, long long const m) 
{ 
  return m==0? n : _gcd(m, n%m)    
};


/** Implements the `Waterfall` algorithm for distributed computation
    of matrix rank. */
class Waterfall {

public:
  /** Constructor, taking row and column dimension, plus a @c
      boost::mpi::communicator instance to use. */
  Waterfall(const coord_t nrows, const coord_t ncols, mpi::communicator& comm_ = mpi::communicator()) 
    : comm(comm_), mpi_id(comm.rank()), mpi_nprocs(comm.size()), procs(),
      nrows(nrows), ncols(ncols)
  {  };

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
  size_t read(std::istream& input) {
    // only read rows whose leading column falls in our domain
    throw "not implemented";
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
    for (size_t n = 0; n < 1 + (nrows / mpi_nprocs); ++n)
      procs.push_back(Processor(mpi_id + n * mpi_nprocs, ncols));

    
    std::vector<int> r(0);
    for (size_t n0 = 0; n0 < procs.size(); 
         /* n0 incremented each time a processor goes into `done` state */ ) {
      if (procs[n0].done())
        ++n0;
      // step all processors
      for (size_t n = n0; n < procs.size(); ++n) {
        r[n] = procs[n].step();
      };
      // manage arrival of messages
      while (boost::optional<mpi::status> status = mpi::iprobe()) { 
        switch(status.tag()) {
        case TAG_ROW_SPARSE:
          row_ptr new_row = new SparseRow();
          mpi::recv(status.source(), status.tag(), *new_row);
          assert(is_local(new_row->first_nonzero_column()));
          local_owner(new_row->first_nonzero_column()).recv_row(new_row);
          break;
        case TAG_ROW_DENSE:
          row_ptr new_row = new DenseRow();
          mpi::recv(status.source(), status.tag(), *new_row);
          assert(is_local(new_row->first_nonzero_column()));
          local_owner(new_row->first_nonzero_column()).recv_row(new_row);
          break;
        case TAG_END:
          // "end" message received; no new rows will be coming.
          // But some other rows could have arrived or could
          // already be in the `inbox`, so we need to make another
          // pass anyway.  All this boils down to: set a flag,
          // make another iteration, and end the loop next time.
          coord_t column = -1;
          mpi::recv(status.source(), status.tag(), column);
          assert(column >= 0 and is_local(column));
          local_owner(column).end_phase();
        }; // switch(status.tag())
      };
    };

    // the partial rank is computed as the sum of all ranks computed
    // by local processors
    int partial_rank = 0;
    for (size_t n = 0; n < r.size(); ++n)
      partial_rank += r[n];

    // wait until all processors have done running
    comm::barrier();

    // collect the partial ranks for all processes
    int rank = 0;
    mpi::all_reduce(comm, partial_rank, rank, std::plus<int,int>);
    return rank;
  };

private:
  const coord_t nrows;
  const coord_t ncols;
  mpi::communicator comm;
  const int mpi_id; /**< MPI rank. */
  const int mpi_nprocs; /**< Total number of ranks in MPI communicator. */
  
  std::vector<Processor> procs;   /**< Local processors. */
  
protected:
  // implement cyclic distribution of columns to MPI ranks
  bool is_local(const coord_t c) { return (mpi_id == c % mpi_nprocs); } const;
  Processor& local_owner(const coord_t c) { return procs[c / mpi_procs]; };
  int remote_owner(const coord_t c) { return (c % mpi_nprocs); };

  /** A single processing element. */
  class Processor {

  public:
    Processor(const coord_t leading_column, const size_t ncols) 
      : leading_column(leading_column), ncols(ncols), rows(), 
        phase(running), inbox(), outbox(), 
        dense_threshold(40.0) 
    { };

  public:
    
    class SparseRow;
    class DenseRow;
    typedef SparseRow* sparse_row_ptr;
    typedef DenseRow* dense_row_ptr;
    typedef boost::variant< sparse_row_ptr,dense_row_ptr > row_ptr;

    class SparseRow 
    {
    public:
      SparseRow(const coord_t leading_column_) 
        : leading_column(leading_column_), storage() { };
      /** Return index of first column that has a nonzero entry. */
      coord_t first_nonzero_column() {
        assert(storage.size() > 0);
        return storage.back().first;
      } const;
      /** Return number of allocated (i.e., non-zero) entries. */
      size_t size() { return storage.size(); } const;
      /** Perform Gaussian elimination: sum a multiple of this row and
          (a multiple of) row @c r so that the combination has its
          first nonzero entry at a column index strictly larger than
          the one of both rows. Return pointer to the combined row, which
          could possibly be this row if in-place update took place. */
      sparse_row_ptr merge(sparse_row_ptr other) {
        assert(this->leading_column == other->leading_column);

        sparse_row_ptr result = new SparseRow(leading_column);
        // XXX: is it worth to compute the exact size here?
        result->storage.reserve(storage.size() + other->storage.size());

        // compute:
        //   `a`: multiplier for `this` row
        //   `b`: multiplier for `other` row
        const val_t GCD = gcd(leading_term, other->leading_term);
        a = - (other->leading_term / GCD);
        b = leading_term / GCD;

        size_t ix = 0;
        size_t iy = 0;
        // loop while one of the two indexes is still valid
        while(ix < storage.size() or iy < other->storage.size()) {
          coord_t jx;
          if (ix < storage.size()) {
            jx = storage[ix].first;
            assert (jx > leading_column);
          }
          else // ix reached end of vector
            jx = -1;

          coord_t jy;
          if (iy < other->storage.size()) {
            jy = other->storage[iy].first;
            assert (jy > other->leading_column);
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
        delete r; // release old storage
        return result;
      }; // sparse_row_ptr merge(sparse_row_ptr r)
      /** Perform Gaussian elimination: sum a multiple of this row and
          (a multiple of) row @c r so that the combination has its
          first nonzero entry at a column index strictly larger than
          the one of both rows. Return pointer to the combined row, which
          could possibly be this row if in-place update took place. 
          Elimination on a dense row returns a dense row. */
      dense_row_ptr merge(dense_row_ptr other) {
        // compute:
        //   `a`: multiplier for `this` row
        //   `b`: multiplier for `other` row
        const val_t GCD = gcd(leading_term, other->leading_term);
        a = - (other->leading_term / GCD);
        b = leading_term / GCD;

        for (size_t j = 0; j < other->size(); ++j) 
          other->storage[j] *= b;

        for (typename storage_t::const_iterator it = storage.begin();
             it != storage.end();
             ++it)
          {
            assert(it->first > leading_column);
            other->storage[it->first] += a*it->second;
          };
        return other;
      }; // dense_row_ptr merge(dense_row_ptr other)
    protected:
      const coord_t leading_column; 
      const val_t leading_term;
      typedef std::vector< std::pair<coord_t,val_t> > storage_t;
      storage_t storage;
      friend class DenseRow;
    }; // class SparseRow

    class DenseRow
    {
    public:
      /** Constructor, copying entries from a @c SparseRow instance. */
      DenseRow(const sparse_row_ptr r) 
        : leading_column(r->leading_column), 
          _size(ncols - r->leading_column),
          storage(new val_t[_size]) { 
        assert(std::distance(r->storage.begin(), r->storage.end()) <= _size);
        std::fill_n(storage, _size, 0);
        for (typename SparseRow::storage_t::const_iterator it = r.storage.begin();
             it != r.storage.end();
             ++it) {
          storage[it->first] = it->second;
        };
      };
      /** Construct a null row. */
      DenseRow(const coord_t leading_columns, const size_t length) 
        : leading_column(r->leading_column), 
          _size(length),
          storage(new val_t[_size]) { 
        std::fill_n(storage, _size, 0);
      };
      ~DenseRow() { 
        delete [] storage;
      };
      /** Return index of first column that has a nonzero entry. */
      coord_t first_nonzero_column() {
        assert(_size > 0);
        for (size_t j = 0; j < _size; ++j)
          if (storage[j] != 0)
            return leading_column + j;
        assert(false); // `first_nonzero_column` should not be called on zero rows
      } const;
      /** Return number of allocated entries. */
      size_t size() { return _size; } const;
      dense_row_ptr merge(const sparse_row_ptr other) {
        // convert `other` to dense storage upfront: adding the
        // non-zero entries from `this` would made it pretty dense
        // anyway
        dense_row_ptr dense_other(new DenseRow(other));
        delete other;
        return this->merge(dense_other);
      }; // dense_row_ptr merge(sparse_row_ptr other)
      dense_row_ptr merge(const dense_row_ptr other) {
        assert(this->size() == other->size());
        for (size_t j = 0; j < size(); ++j) {
          // XXX: is it faster to allocate new storage and fill it with `a*x+b*y`?
          other->storage[j] *= a;
          other->storage[j] += b * storage[j];
        };
        return other;
      }; // dense_row_ptr merge(dense_row_ptr other)
    protected:
      const coord_t leading_column; 
      const val_t leading_term; 
      /// Row length.  This is of course the number of elements in the @c storage array.
      const size_t _size;
      /// Plain array of entries.
      val_t* storage; 
    }; // class DenseRow


    /** The block of rows owned by this processor */
  private:     
    const coord_t leading_column;

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
    std::list< row_ptr > inbox;

    /** List of rows sent to other processors and the associated MPI request. 
        We need to keep track of these in order to free the resources when we're done. */
    std::list< std::pair<mpi::request, row_ptr> > outbox;
    /** Adapt the @c outbox iterator to provide the kind of
        bidirectional iterator over @c mpi::request objects that @c
        boost::mpi::test_some() requires. */
    class request_iterator : public outbox::iterator {
    public:
      request_iterator(outbox::iterator& parent_)
        : outbox::iterator(parent_);
      virtual ~request_iterator();
      mpi::request& operator->() { return static_cast<outbox::iterator>(this)->first; }
    };

    /** A row will switch its storage to dense format when the percent of
        nonzero entries w.r.t. total length exceeds this value. */
    const double dense_threshold;

  public: 
    void recv_row(row_ptr& new_row) {
      // FIXME: requires locking for MT operation
      inbox.push_back(new_row);
    };

    /** Send the row data pointed by @c row to the processor owning the
        block starting at @c col. */
    void send_row(const coord_t column, row_ptr& row) {
      assert(NULL != row);
      if (is_local(column))
        local_owner(column).recv_row(row);
      else { // ship to remote process
        class _do_send : public boost::static_visitor<mpi::request> {
          template <>
          mpi::request operator() (sparse_row_ptr r) {
            return mpi::isend(remote_owner(column), TAG_ROW_SPARSE, r);
          };
          template <>
          mpi::request operator() (dense_row_ptr r) {
            return mpi::isend(remote_owner(column), TAG_ROW_DENSE, r);
          };
        } do_send;
        mpi::request req = boost::apply_visitor(do_send, new_row);
        outbox.push_back(req, new_row);
      };
    };

    /** Make a pass over the block of rows and return either 0 or 1
        (depending whether elimination took place or not). */
    int step() {
      assert(not done()); // cannot be called once finished

      // receive new rows
      for (inbox::iterator it = inbox.begin(); it != inbox.end(); ++it) {
        row_ptr new_row = *it;
        rows.push_back(new_row);
        if (new_row->size() < rows[0]->size())
          // less nonzero entries: we have a new leading row
          std::swap(rows[0], rows[rows.size()-1]);
      };
      inbox.clear();
          
      if (rows.size() > 0) {
        row_ptr first = rows[0];
        for (block::iterator second = rows.begin()+1; second != rows.end(); ++second) {
          // perform elimination
          row_ptr new_row = gaussian_elimination_step(first, second);
          
          // ship reduced rows to other processors
          if (new_row->size() > 0) {
            send(new_row->first_nonzero_column(), new_row);
          };
        };
      };

      // just keep the "first" row
      rows.clear();
      rows.push_back(first);
      
      if (running == phase) {
        if (outbox.size() > 0) {
          // check if some test messages have arrived...
          request_iterator first = outbox.begin();
          request_iterator last = outbox.end();
          request_iterator first_done = mpi::test_some(first, last);
          
          // ...and free corresponding resources
          for (outbox::iterator it = first_done.parent();
               it != outbox.end(); 
               ++it) {
            delete it->second;
            outbox.erase(it);
          };
        };
      }
      else { // `phase == ending`: end message already received
        // pass end message along
        send_end(leading_column + 1);
        
        if (outbox.size() > 0) {
          // wait untill all sent messages have arrived
          request_iterator first = outbox.begin();
          request_iterator last = outbox.end();
          mpi::wait_all(first, last);
          
          // free resources
          for (outbox::iterator it = outbox.begin(); 
               it != outbox.end(); 
               ++it)
            delete it->second;
        };
        
        // all done, exit with 0 only if we never processed any row
        phase = done;
      };
      
      assert(rows.size() <= 1);
      return rows.size();
    };

    void end_phase() { phase = ending; };

    bool done() { return (done == phase); } const;

  protected:
    /** Send the termination signal to the processor owning the block
        starting at @c col. */
    void send_end(const coord_t column) {
      if (is_local(column))
        local_owner(column).recv(TAG_END, NULL);
      else { // ship to remote process
        mpi::send(remote_owner(column), TAG_END, column);
      };
    }; // send_end


    /** Perform Gaussian elimination on row `second`, using the entry
        in the leading column of `first` as a pivot. */
    row_ptr gaussian_elimination(row_ptr first, row_ptr second) {
      class _ge_op
        : public boost::static_visitor<row_ptr>
      {
      public:
        template <>
        row_ptr operator() (sparse_row_ptr first, sparse_row_ptr second) const {
          // assume that `first` has already been checked for fill-in
          // in the main loop, so we won't attempt to convert it here;
          // just test `second` for too much fill-in
          const double fill_in = 100.0 * second.size() / (ncols - leading_column);
          if (fill_in > dense_threshold) {
            // FIXME: could merge ctor+elimination in one funcall
            dense_row_ptr dense_second(new DenseRow(second));
            delete second;
            return first->merge(dense_second);
          }
          else {
            return first->merge(second);
          }; // if (fill_in > dense_threshold)
        };

        template <>
        row_ptr operator() (dense_row_ptr first, sparse_row_ptr second) const {
          return first->merge(second);
        };

        template <>
        row_ptr operator() (sparse_row_ptr first, dense_row_ptr second) const {
          return first->merge(second);
        };

        template <>
        row_ptr operator() (dense_row_ptr first, dense_row_ptr second) const {
          return first->merge(second);
        };
      } ge_op;
      return boost::apply_visitor(ge_op, first, second);
    }; // gaussian_elimination
  };
};

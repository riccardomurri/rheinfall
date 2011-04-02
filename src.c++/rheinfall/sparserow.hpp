/**
 * @file   sparserow.hpp
 *
 * Interface of the sparserow class.
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

#ifndef RF_SPARSEROW_HPP
#define RF_SPARSEROW_HPP 1


#include "config.hpp"
#include "row.hpp"

#ifdef WITH_MPI
# include <boost/mpi.hpp>
  namespace mpi = boost::mpi;
# include <boost/serialization/access.hpp>
# include <boost/serialization/utility.hpp>
# include <boost/serialization/vector.hpp>
# ifdef WITH_GMPXX
#  include <gmpxx_boost_serialization.hpp>
# endif
#endif

#include <cassert>
#include <iostream>
#include <iterator>
#include <memory> // std::allocator
#include <stdexcept>
#include <utility>
#include <vector>


namespace rheinfall {

  // forward declaration of needed classes
  // (otherwise we have circular dependency problems)
  template <typename val_t, typename coord_t, template<typename T> class allocator> class Rheinfall;
  template <typename val_t, typename coord_t, template<typename T> class allocator> class Row;
  template <typename val_t, typename coord_t, template<typename T> class allocator> class SparseRow;
  template <typename val_t, typename coord_t, template<typename T> class allocator> class DenseRow;


  template <typename val_t, typename coord_t,
            template<typename T> class allocator = std::allocator >
  /** Store a matrix row as a vector of (column, entry) pairs. */
  class SparseRow : public Row<val_t,coord_t,allocator>
  {
  public:

    /** Construct a row with the given starting and ending column. The
        entry at @a starting_column is set to @a leading_term, any
        other entry is null. */
    SparseRow(const coord_t starting_column, 
              const coord_t ending_column, 
              const val_t leading_term);

    /** Make a @a SparseRow instance from coordinate/value pairs
        gotten from interval [p0, p1). */
    template <typename ForwardIter>
    static
    SparseRow<val_t,coord_t,allocator>* new_from_range(ForwardIter p0, ForwardIter p1, 
                                                       const coord_t ending_column);

#ifdef WITH_MPI
    /** Make a @a SparseRow instance from an MPI message payload.
        The MPI message is identified by the triple
        communicator/source/tag, that is passed unchanged to @c
        boost::mpi::communicator::recv(). */
    static
    SparseRow<val_t,coord_t,allocator>* new_from_mpi_message(mpi::communicator& comm, 
                                                             const int source, const int tag);
#endif

    /** Return fill-in percentage, that is the number of actual
        nonzero entries divided by the number of potential entries.  */
    float fill_in() const;

    /** Reset starting column and leading term to the first nonzero
        entry in the row.  Return pointer to this row. If this row
        is a null row, then delete it and return @c NULL.  Used for
        two-phase construction. */
    SparseRow* adjust();

    /** Return number of allocated (i.e., non-zero) entries. */
    virtual size_t size() const;

    /** Set the element stored at column @c col */
    virtual void set(const coord_t col, const val_t value);

    /** Return a copy of the element stored at column @c col */
    virtual val_t get(const coord_t col) const;

    /** Perform Gaussian elimination: sum a multiple of this row to
        (a multiple of) row @c r so that the combination has its
        first nonzero entry at a column index strictly larger than
        the one of both rows. Return pointer to the combined row, which
        could possibly be this row if in-place update took place. */
    SparseRow<val_t,coord_t,allocator>* gaussian_elimination(SparseRow<val_t,coord_t,allocator>* other) const;

    /** Perform Gaussian elimination: sum a multiple of this row to
        (a multiple of) row @c r so that the combination has its
        first nonzero entry at a column index strictly larger than
        the one of both rows. Return pointer to the combined row, which
        could possibly be this row if in-place update took place. 
        Elimination on a dense row returns a dense row. */
    DenseRow<val_t,coord_t,allocator>* gaussian_elimination(DenseRow<val_t,coord_t,allocator>* other) const;

  protected:

    /** Nickname for base class; used to shorten templatized expressions. */
    typedef Row<val_t,coord_t,allocator> Row_; 

    /** Type used for storing the (coordinate, value) pairs that make up a sparse row. */
    typedef std::vector< std::pair<coord_t,val_t>,
                         allocator<std::pair<coord_t,val_t> > > storage_t;
    /** Actual storage of the (coordinate, value) pairs that make up a
        sparse row.  The initial nonzero term is stored separately
        (see @ref Row::leading_term_). */
    storage_t storage;

    /** Print a textual representation of the row to stream @c o. */
    virtual void print_on(std::ostream& o) const;

    friend class Rheinfall<val_t,coord_t,allocator>; // read*() needs to call the above ctor
    friend class DenseRow<val_t,coord_t,allocator>;

    /** Default constructor, needed by boost::serialize.
        Initializes the row with invalid values; should call
        serialize immediately after!  */
    SparseRow();

#ifdef WITH_MPI
    friend class boost::serialization::access;
    /// boost::serialization support
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version);
#endif

#ifndef NDEBUG
  private:
    /// Check that the internal data structure are consistent; only
    /// used in assert() calls for debugging.
    bool __ok() const;
#endif
  }; // class SparseRow

}; // namespace rheinfall


  // ------- inline methods -------

#include "types.hpp"
#include "denserow.hpp"

namespace rheinfall {

  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline
  SparseRow<val_t,coord_t,allocator>::
  SparseRow(const coord_t starting_column, 
            const coord_t ending_column, 
            const val_t leading_term) 
    : Row_::Row(Row_::sparse, starting_column, ending_column, leading_term),
      storage() 
  { 
    // init-only ctor
  };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline
  SparseRow<val_t,coord_t,allocator>::SparseRow() 
    : Row_::Row(Row_::sparse, -1, -1, 0),
      storage() 
  { 
    // init-only ctor
  };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  template <typename ForwardIter>
  inline
  SparseRow<val_t,coord_t,allocator>*
  SparseRow<val_t,coord_t,allocator>::
  new_from_range(ForwardIter p0,
                 ForwardIter p1,
                 const coord_t ending_column)
  { 
    SparseRow<val_t,coord_t,allocator>* row = new SparseRow();
    row->ending_column_ = ending_column;
    coord_t starting_column = ending_column;
    val_t leading_term = 0;
    for (; p0 != p1; ++p0) {
      const coord_t coord = p0->first;
      const val_t value = p0->second;
      if(is_zero(value))
        continue;
      if (coord < starting_column) {
        if (not is_zero(leading_term))
          row->set(starting_column, leading_term);
        starting_column = coord;
        leading_term = value;
      }
      else {
        row->set(coord, value);
      };
    };
    row->Row_::starting_column_ = starting_column;
    row->Row_::leading_term_ = leading_term;

    if (row->starting_column_ == row->ending_column_) {
      delete row;
      return NULL;
    } 
    else {
      assert(row->__ok());
      return row;
    };
  };


#ifdef WITH_MPI
  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline
  SparseRow<val_t,coord_t,allocator>*
  SparseRow<val_t,coord_t,allocator>::
  new_from_mpi_message(mpi::communicator& comm, 
                       const int source, const int tag)
  {
    SparseRow<val_t,coord_t,allocator>* row = new  SparseRow<val_t,coord_t,allocator>();
    comm.recv(source, tag, *row);
    return row;
  };
#endif // WITH_MPI


#ifndef NDEBUG
  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline bool
  SparseRow<val_t,coord_t,allocator>::
  __ok() const
  {
    assert(0 <= Row_::starting_column_);
    assert(Row_::starting_column_ <= Row_::ending_column_);
    assert(not is_zero(Row_::leading_term_));
    // entries in `this->storage` are *always* ordered by increasing column index
    coord_t s1 = Row_::starting_column_;
    for (typename storage_t::const_iterator it = storage.begin(); 
         it < storage.end(); 
         ++it) {
      assert(s1 < it->first);
      assert(it->first <= Row_::ending_column_);
      assert(not is_zero(it->second));
      s1 = it->first;
    };
    return true;
  };
#endif


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline SparseRow<val_t,coord_t,allocator>*
  SparseRow<val_t,coord_t,allocator>::
  adjust()
  {
    if (storage.size() > 0) {
      Row_::starting_column_ = storage.front().first;
      assert(0 <= Row_::starting_column_ and Row_::starting_column_ < Row_::ending_column_);
      Row_::leading_term_ = storage.front().second;
      assert(not is_zero(Row_::leading_term_));
      storage.erase(storage.begin());
      assert(this->__ok());
      return this;
    }
    else {
      // nothing in the row
      delete this;
      return NULL;
    };
    assert(false); // should not happen!
  };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  float
  SparseRow<val_t,coord_t,allocator>::fill_in() const 
  { 
    return 100.0 * size() / (Row_::ending_column_ - Row_::starting_column_); 
  };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  SparseRow<val_t,coord_t,allocator>* 
  SparseRow<val_t,coord_t,allocator>::
  gaussian_elimination(SparseRow<val_t,coord_t,allocator>* restrict other) 
    const restrict_this
  {
    assert(this->__ok());
    assert(other->__ok());
    assert(this->starting_column_ == other->starting_column_);

    val_t a, b;
    rheinfall::get_row_multipliers<val_t>(this->Row_::leading_term_,
                                          other->Row_::leading_term_,
                                          a, b);

    SparseRow<val_t,coord_t,allocator>* restrict result = NULL; // XXX: use boost::optional<...> instead?

    typename storage_t::const_iterator this_i = this->storage.begin(); 
    typename storage_t::const_iterator other_i = other->storage.begin();
    // loop while one of the two indexes is still valid
    while(this_i != this->storage.end() or other_i != other->storage.end()) {
      const coord_t this_col =
        (this_i != this->storage.end())?
        // use coord from row element
        this_i->first
        // this_i reached end of vector, use out-of-range value
        : this->ending_column_ + 1;

      const coord_t other_col = 
        (other_i != other->storage.end())?
        // use coord from row element
        other_i->first
        // else, other_i reached end of vector, use out-of-range value
        : other->ending_column_ + 1;

      bool nonzero_found = false;
      coord_t coord;
      val_t entry;
      if (other_col < this_col) {
        entry = b * other_i->second;
        assert(not is_zero(entry));
        coord = other_col;
        nonzero_found = true;
        ++other_i;
      }
      else if (other_col == this_col) {
        entry = a * this_i->second + b * other_i->second;
        if (not is_zero(entry)) {
          coord = this_col;
          nonzero_found = true;
        };
        ++this_i;
        ++other_i;
      }
      else if (other_col > this_col) {
        entry = a * this_i->second;
        assert(not is_zero(entry));
        coord = this_col;
        nonzero_found = true;
        ++this_i;
      }
      else 
        // should not happen!
        throw std::logic_error("Unhandled case in SparseRow::gaussian_elimination(SparseRow*)");
    
      if (nonzero_found) {
        // XXX: the following code makes assumptions about how the
        // storage is organised within `SparseRow`; it duplicates code
        // from SparseRow::operator[] for efficiency
        if (NULL == result) {
          // allocate new SparseRow
          result = new SparseRow(coord, Row_::ending_column_, entry);
          assert(storage.size() + other->storage.size() <= result->storage.max_size());
          result->storage.reserve(storage.size() + other->storage.size());
          assert(0 == result->storage.size());
        }
        else {
          result->storage.push_back(std::make_pair(coord, entry));
        }; 
      }; 
    }; // while(this_i < storage.end() ...)
#ifndef NDEBUG
    if (NULL != result) {
      assert(result->__ok());
      assert(result->starting_column_ > this->starting_column_);
      assert(result->starting_column_ > other->starting_column_);
      assert(result->ending_column_ == this->ending_column_);
      assert(result->size() <= this->size() + other->size());
    };
#endif
    delete other; // release old storage
    return result;
  }; // SparseRow<val_t,coord_t,allocator>* gaussian_elimination(SparseRow<val_t,coord_t,allocator>* r)


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  DenseRow<val_t,coord_t,allocator>* 
  SparseRow<val_t,coord_t,allocator>::
  gaussian_elimination(DenseRow<val_t,coord_t,allocator>* restrict other) 
    const restrict_this
  {
    assert(this->starting_column_ == other->starting_column_);
    assert(this->__ok());
    assert(other->__ok());
    assert(not is_zero(other->leading_term_));

    val_t a, b;
    rheinfall::get_row_multipliers(this->Row_::leading_term_,
                                   other->Row_::leading_term_,
                                   a, b);

    for (size_t j = 0; j < other->size(); ++j) 
      other->storage[j] *= b;
    for (typename storage_t::const_iterator it = storage.begin();
         it != storage.end();
         ++it)
      {
        assert(it->first > Row_::starting_column_ and it->first <= Row_::ending_column_);
        other->storage[other->size()-1 - (it->first - (other->starting_column_ + 1))] += a * it->second;
      };
    other->Row_::leading_term_ = 0; // by construction

    DenseRow<val_t,coord_t,allocator>* result = other->adjust();
    assert(NULL == result 
           or result->Row_::starting_column_ > this->Row_::starting_column_);
    return result; // content update done, adjust size and starting column
  }; // DenseRow<val_t,coord_t,allocator>* gaussian_elimination(DenseRow<val_t,coord_t,allocator>* other)


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline val_t
  SparseRow<val_t,coord_t,allocator>::
  get(const coord_t col) const
  {
    assert(col >= Row_::starting_column_ and col <= Row_::ending_column_);
    if (col == Row_::starting_column_) 
      return Row_::leading_term_;
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
  }; // SparseRow::get(...) const 


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline void
  SparseRow<val_t,coord_t,allocator>::
  print_on(std::ostream& out) const 
  {
    out << "{ "
        << Row_::starting_column_ <<":"<< Row_::leading_term_;
    for (typename storage_t::const_iterator it = storage.begin(); 
         it != storage.end(); 
         ++it) {
      out <<" "<< it->first <<":"<< it->second;
    }
    out << " }";
  };


#ifdef WITH_MPI
  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  template<class Archive>
  inline void 
  SparseRow<val_t,coord_t,allocator>::
  serialize(Archive& ar, const unsigned int version) 
  {
    assert(version == 0);
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    ar 
      & boost::serialization::base_object< Row<val_t,coord_t,allocator> >(*this) 
      & storage;
    assert(this->__ok());
  }; // SparseRow::serialize(...)
#endif // WITH_MPI


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline void
  SparseRow<val_t,coord_t,allocator>::
  set(const coord_t col, const val_t value) 
  {
    assert(col >= Row_::starting_column_ and col <= Row_::ending_column_);
    if (col == Row_::starting_column_) {
      Row_::leading_term_ = value;
      return;
    };
    // else, fast-forward to place where element is/would be stored
    size_t jj = 0;
    while (jj < storage.size() and storage[jj].first < col)
      ++jj;
    // set element accordingly
    if (storage.size() == jj) {
      // end of list reached, `col` is larger than any index in this row
      storage.push_back(std::make_pair<size_t,val_t>(col,value));
    }
    else if (col == storage[jj].first) 
      storage[jj].second = value;
    else { // storage[jj].first > col, insert new pair before `jj`
      storage.insert(storage.begin()+jj, 
                     std::make_pair<size_t,val_t>(col,value));
    };
    //assert(this->__ok());
  }; // SparseRow::set(...)


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline size_t 
  SparseRow<val_t,coord_t,allocator>::size() const 
  { 
    return storage.size(); 
  };


}; // namespace rheinfall


#endif // RF_SPARSEROW_HPP

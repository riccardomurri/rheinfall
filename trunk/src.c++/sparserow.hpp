/**
 * @file   sparserow.hpp
 *
 * Interface of the sparserow class.
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

#ifndef RF_SPARSEROW_HPP
#define RF_SPARSEROW_HPP


#include "row.hpp"
#include "denserow.hpp"
#include "rheinfall.hpp"

#ifdef WITH_MPI
# include <boost/mpi.hpp>
# include <boost/optional.hpp>
# include <boost/serialization/access.hpp>
# include <boost/serialization/utility.hpp>
# include <boost/serialization/vector.hpp>
  namespace mpi = boost::mpi;
#endif

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>


namespace rheinfall {

  // forward declaration of needed classes
  // (otherwise we have circular dependency problems)
  template <typename val_t, typename coord_t> class Rheinfall;
  template <typename val_t, typename coord_t> class Row;
  template <typename val_t, typename coord_t> class SparseRow;
  template <typename val_t, typename coord_t> class DenseRow;


  template <typename val_t, typename coord_t>
  /** Store a matrix row as a vector of (column, entry) pairs. */
  class SparseRow : public Row<val_t,coord_t>
  {
  public:

    SparseRow(const coord_t starting_column, 
              const coord_t ending_column, 
              const val_t leading_term);

    /** Return fill-in percentage, that is the number of actual
        nonzero entries divided by the number of potential entries.  */
    double fill_in() const;

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
    SparseRow<val_t,coord_t>* gaussian_elimination(SparseRow<val_t,coord_t>* other) const;

    /** Perform Gaussian elimination: sum a multiple of this row to
        (a multiple of) row @c r so that the combination has its
        first nonzero entry at a column index strictly larger than
        the one of both rows. Return pointer to the combined row, which
        could possibly be this row if in-place update took place. 
        Elimination on a dense row returns a dense row. */
    DenseRow<val_t,coord_t>* gaussian_elimination(DenseRow<val_t,coord_t>* other) const;

  protected:

    typedef Row<val_t,coord_t> Row_; //< Nickname for base class; used to shorten templatized expressions

    typedef std::vector< std::pair<coord_t,val_t> > storage_t;
    storage_t storage;

    /** Print a textual representation of the row to stream @c o. */
    virtual void print_on(std::ostream& o) const;

    /** Constructor initializing an invalid row; should fill the row
        with values and then call @c adjust(). */
    SparseRow(const coord_t ending_column);

    friend class Rheinfall<val_t,coord_t>; // read*() needs to call the above ctor
    friend class DenseRow<val_t,coord_t>;

#ifdef WITH_MPI
    // serialization support
    /** Default constructor, needed by boost::serialize.
        Initializes the row with invalid values; should call
        serialize immediately after!  */
    SparseRow();
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version);
#endif
#ifndef NDEBUG
  private:
    bool __ok() const;
#endif
  }; // class SparseRow

}; // namespace rheinfall


  // ------- inline methods -------

#include "mathdefs.hpp"

namespace rheinfall {

  template <typename val_t, typename coord_t>
  inline
  SparseRow<val_t,coord_t>::SparseRow(const coord_t starting_column, 
                                      const coord_t ending_column, 
                                      const val_t leading_term) 
    : Row_::Row(Row_::sparse, 
                              starting_column, ending_column, leading_term, true),
      storage() 
  { 
    // init-only ctor
  };

  template <typename val_t, typename coord_t>
  inline
  SparseRow<val_t,coord_t>::SparseRow(const coord_t ending_column) 
    : Row_::Row(Row_::sparse, -1, ending_column, 0, false),
      storage() 
  { 
    // init-only ctor
  };

  // template <typename val_t, typename coord_t>
  // inline
  // SparseRow<val_t,coord_t>::SparseRow() 
  //   : Row_::Row(Row_::sparse, -1, -1, 0, false),
  //     storage() 
  // { 
  //   // init-only ctor
  // };


#ifndef NDEBUG
  template <typename val_t, typename coord_t>
  inline bool
  SparseRow<val_t,coord_t>::__ok() const
  {
    if (Row_::initialized_) {
      assert(0 <= Row_::starting_column_);
      assert(Row_::starting_column_ <= Row_::ending_column_);
      assert(0 != Row_::leading_term_);
    };
    // entries in `this->storage` are *always* ordered by increasing column index
    coord_t s1 = (Row_::initialized_? Row_::starting_column_ : -1);
    for (typename storage_t::const_iterator it = storage.begin(); 
         it < storage.end(); 
         ++it) {
      assert(s1 < it->first);
      assert(it->first <= Row_::ending_column_);
      assert(0 != it->second);
      s1 = it->first;
    };
    return true;
  };
#endif


  template <typename val_t, typename coord_t>
  inline SparseRow<val_t,coord_t>*
  SparseRow<val_t,coord_t>::adjust()
  {
    if (storage.size() > 0) {
      Row_::starting_column_ = storage.front().first;
      assert(0 <= Row_::starting_column_ and Row_::starting_column_ < Row_::ending_column_);
      Row_::leading_term_ = storage.front().second;
      assert(0 != Row_::leading_term_);
      storage.erase(storage.begin());
      Row_::initialized_ = true;
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


  template <typename val_t, typename coord_t>
  double 
  SparseRow<val_t,coord_t>::fill_in() const 
  { 
    return 100.0 * size() / (Row_::ending_column_ - Row_::starting_column_); 
  };


  template <typename val_t, typename coord_t>
  SparseRow<val_t,coord_t>* 
  SparseRow<val_t,coord_t>::gaussian_elimination(SparseRow<val_t,coord_t>* other) const
  {
    assert(this->__ok());
    assert(other->__ok());
    assert(this->starting_column_ == other->starting_column_);

    val_t a, b;
    rheinfall::get_row_multipliers<val_t>(this->Row_::leading_term_,
                                          other->Row_::leading_term_,
                                          a, b);

    SparseRow<val_t,coord_t>* result = NULL; // XXX: use boost::optional<...> instead?

    typename storage_t::const_iterator this_i = this->storage.begin(); 
    typename storage_t::const_iterator other_i = other->storage.begin();
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
        assert(0 != entry);
        coord = other_col;
        nonzero_found = true;
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
        assert(0 != entry);
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
  }; // SparseRow<val_t,coord_t>* gaussian_elimination(SparseRow<val_t,coord_t>* r)


  template <typename val_t, typename coord_t>
  DenseRow<val_t,coord_t>* 
  SparseRow<val_t,coord_t>::gaussian_elimination(DenseRow<val_t,coord_t>* other) const
  {
    assert(this->starting_column_ == other->starting_column_);
    assert(this->__ok());
    assert(0 != other->leading_term_);

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
        other->storage[other->size() - (it->first - other->starting_column_)] += a*it->second;
      };

    return other->adjust(); // content update done, adjust size and starting column
  }; // DenseRow<val_t,coord_t>* gaussian_elimination(DenseRow<val_t,coord_t>* other)


  template <typename val_t, typename coord_t>
  inline val_t
  SparseRow<val_t,coord_t>::get(const coord_t col) const
  {
    assert((not Row_::initialized_) or 
           (col >= Row_::starting_column_ and col <= Row_::ending_column_));
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


  template <typename val_t, typename coord_t>
  inline void
  SparseRow<val_t,coord_t>::print_on(std::ostream& out) const 
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
  template <typename val_t, typename coord_t>
  template<class Archive>
  inline void 
  SparseRow<val_t,coord_t>::serialize(Archive& ar, const unsigned int version) 
  {
    assert(version == 0);
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    ar & boost::serialization::base_object<Row>(*this) & storage;
    initialized_ = true;
    assert(this->__ok());
  }; // SparseRow::serialize(...)
#endif // WITH_MPI


  template <typename val_t, typename coord_t>
  inline void
  SparseRow<val_t,coord_t>::set(const coord_t col, const val_t value) 
  {
    assert((not Row_::initialized_) or 
           (col >= Row_::starting_column_ and col <= Row_::ending_column_));
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


  template <typename val_t, typename coord_t>
  inline size_t 
  SparseRow<val_t,coord_t>::size() const 
  { 
    return storage.size(); 
  };


}; // namespace rheinfall


#endif // RF_SPARSEROW_HPP

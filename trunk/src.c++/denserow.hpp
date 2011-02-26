/**
 * @file   denserow.hpp
 *
 * Interface of the denserow class.
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

#ifndef RF_DENSEROW_HPP
#define RF_DENSEROW_HPP


#include "row.hpp"
#include "sparserow.hpp"

#ifdef WITH_MPI
# include <boost/mpi.hpp>
# include <boost/optional.hpp>
# include <boost/serialization/access.hpp>
# include <boost/serialization/vector.hpp>
  namespace mpi = boost::mpi;
#endif

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>


namespace rheinfall {

  // declaration of Row-based classes
  template <typename val_t, typename coord_t> class Row;
  template <typename val_t, typename coord_t> class SparseRow;
  template <typename val_t, typename coord_t> class DenseRow;


  template <typename val_t, typename coord_t>
    /** Store a matrix row in a linear vector of consecutive entries. */
  class DenseRow : public Row<val_t,coord_t>
    {
    public:
      /** Constructor, copying entries from a @c SparseRow instance. */
      DenseRow(const SparseRow<val_t,coord_t>* r);

      /** Construct a null row. */
      DenseRow(const coord_t starting_column, const size_t ending_column);

      /** Return number of allocated entries. */
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
      DenseRow<val_t,coord_t>* gaussian_elimination(const SparseRow<val_t,coord_t>* other) const;

    /** Perform Gaussian elimination: sum a multiple of this row to
        (a multiple of) row @c r so that the combination has its
        first nonzero entry at a column index strictly larger than
        the one of both rows. Return pointer to the combined row, which
        could possibly be this row if in-place update took place. */
      DenseRow<val_t,coord_t>* gaussian_elimination(const DenseRow<val_t,coord_t>* other) const;

    protected:
      typedef Row<val_t,coord_t> Row_; //< Nickname for base class; used to shorten templatized expressions

      /// array of entries.
      typedef std::vector<val_t> storage_t; 
      storage_t storage;
      friend class SparseRow<val_t,coord_t>;

      /** Print a textual representation of the row to stream @c o. */
      virtual void print_on(std::ostream& o) const;

      /** Adjust @c size_ and @c starting_column to reflect the actual
          contents of the stored row.  Return pointer to adjusted row.
          If this is a null row, then delete it and return NULL. */
      DenseRow<val_t,coord_t>* adjust();

      // serialization support (needed by Boost.MPI)
      /** Default constructor, needed by boost::serialize.
          Initializes the row with invalid values; should call
          serialize immediately after!  */
      DenseRow();

#ifdef WITH_MPI
      friend class boost::serialization::access;
      template<class Archive> void serialize(Archive& ar, const unsigned int version);
#endif
    }; // class DenseRow


  // ------- inline methods -------

  template <typename val_t, typename coord_t>
  inline
  DenseRow<val_t,coord_t>::DenseRow(const SparseRow<val_t,coord_t>* r) 
    : Row_::Row(Row_::dense, r->starting_column_, 
                r->ending_column_, r->leading_term_, true),
      storage(Row_::ending_column_ - Row_::starting_column_ + 1, 0) 
  { 
    assert(std::distance(r->storage.begin(), r->storage.end()) <= storage.size());
    for (typename SparseRow<val_t,coord_t>::storage_t::const_iterator it = r->storage.begin();
         it != r->storage.end();
         ++it) 
      {
        assert(it->first > Row_::starting_column_ and it->first <= Row_::ending_column_);
        storage[size() - (it->first - Row_::starting_column_)] = it->second;
      };
  };


  template <typename val_t, typename coord_t>
  inline
  DenseRow<val_t,coord_t>::DenseRow(const coord_t starting_column, 
                     const size_t ending_column) 
    : Row_::Row(Row_::dense, starting_column, ending_column, 0, true),
      storage(ending_column - starting_column + 1, 0)
  { 
    // nothing to do
  };

  template <typename val_t, typename coord_t>
  inline
  DenseRow<val_t,coord_t>::DenseRow() 
    : Row_::Row(Row_::dense, -1, -1, 0, false), storage() 
  { 
    // nothing to do
  };


  template <typename val_t, typename coord_t>
  inline DenseRow<val_t,coord_t>* 
  DenseRow<val_t,coord_t>::gaussian_elimination(const SparseRow<val_t,coord_t>* other) const
  {
    assert(this->starting_column_ == other->starting_column_);
    assert(0 != this->leading_term_);
    assert(0 != other->leading_term_);

    // convert `other` to dense storage upfront: adding the
    // non-zero entries from `this` would made it pretty dense
    // anyway
    DenseRow<val_t,coord_t>* dense_other(new DenseRow<val_t,coord_t>(other));
    delete other;
    return this->gaussian_elimination(dense_other);
  }; // DenseRow<val_t,coord_t>* gaussian_elimination(SparseRow<val_t,coord_t>* other)


  template <typename val_t, typename coord_t>
  inline DenseRow<val_t,coord_t>*
  DenseRow<val_t,coord_t>::gaussian_elimination(const DenseRow<val_t,coord_t>* other) const
  {
    assert(this->starting_column_ == other->starting_column_);
    assert(0 != this->leading_term_);
    assert(0 != other->leading_term_);
    assert(this->size() == other->size());

    val_t a, b;
    rheinfall::get_row_multipliers(Row_::leading_term_, other->Row_::leading_term_, a, b);

    const coord_t col = 1 + std::max(this->Row_::starting_column_,
                                       other->Row_::starting_column_);
    DenseRow<val_t,coord_t>* result = new DenseRow(col, this->Row_::ending_column_);
    size_t this_j = this->size() - (this->Row_::ending_column_ - col) - 1;
    size_t other_j = other->size() - (other->Row_::ending_column_ - col) - 1;
    for (size_t j = 0; j < result->size(); ++j, ++this_j, ++other_j)
       result->storage[j] = a*other->storage[other_j] + b*this->storage[this_j];
    delete other;
    return result->adjust(); // update done, adjust size and starting column
  }; // dense_row_ptr gaussian_elimination(dense_row_ptr other)


  template <typename val_t, typename coord_t>
  inline val_t
  DenseRow<val_t,coord_t>::get(const coord_t col) const
  {
    assert(col >= Row_::starting_column_ and col <= Row_::ending_column_);
    if (col == Row_::starting_column_)
      return Row_::leading_term_;
    else
      return storage[size() - (Row_::ending_column_ - col) - 1];
  };


  template <typename val_t, typename coord_t>
  inline DenseRow<val_t,coord_t>*
  DenseRow<val_t,coord_t>::adjust()
  {
    // compute new starting column
    for (int j = size()-1; j >= 0; --j)
      if (storage[j] != 0) {
        Row_::leading_term_ = storage[j];
        Row_::starting_column_ += (size() - j);
        storage.erase(storage.begin()+j, storage.end());
        return this;
      };
    // no nonzero element found in storage,
    // this is now a null row
    delete this;
    return NULL;
  };


  template <typename val_t, typename coord_t>
  inline void
  DenseRow<val_t,coord_t>::print_on(std::ostream& out) const 
  {
    out << "["
        << Row_::starting_column_ <<":"<< Row_::leading_term_;
    for (int j = storage.size()-1; j >= 0; --j)
      out <<","<< storage[j];
    out << "]";
  };


#ifdef WITH_MPI
  template <typename val_t, typename coord_t>
  template<class Archive>
  inline void 
  DenseRow<val_t,coord_t>::serialize(Archive& ar, const unsigned int version) 
  {
    assert(version == 0);
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    ar & boost::serialization::base_object<Row>(*this) & storage;
    initialized_ = true;
  };
#endif // WITH_MPI


  template <typename val_t, typename coord_t>
  inline void
  DenseRow<val_t,coord_t>::set(const coord_t col, const val_t value) 
  {
    assert(col >= Row_::starting_column_ and col <= Row_::ending_column_);
    if (col == Row_::starting_column_)
      Row_::leading_term_ = value;
    else
      storage[size() - col - 1] = value;
  };


  template <typename val_t, typename coord_t>
  inline size_t 
  DenseRow<val_t,coord_t>::size() const 
  { 
    return storage.size(); 
  };


}; // namespace rheinfall


#endif // RF_DENSEROW_HPP

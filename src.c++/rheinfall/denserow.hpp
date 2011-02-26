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


#include "config.hpp"
#include "row.hpp"

#ifdef WITH_MPI
# include <boost/mpi.hpp>
# include <boost/serialization/access.hpp>
# include <boost/serialization/vector.hpp>
  namespace mpi = boost::mpi;
# ifdef WITH_GMPXX
#  include <gmpxx_boost_serialization.hpp>
# endif
#endif

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>


namespace rheinfall {

  // forward declarations to avoid circular dependency
  template <typename val_t, typename coord_t> class Rheinfall;
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

#ifdef WITH_MPI
      /** Make a @a DenseRow instance from an MPI message payload.
          The MPI message is identified by the triple
          communicator/source/tag, that is passed unchanged to @c
          boost::mpi::communicator::recv(). */
      static
      DenseRow<val_t,coord_t>* new_from_mpi_message(mpi::communicator& comm, 
                                                    const int source, const int tag);
#endif

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

      friend class Rheinfall<val_t,coord_t>;
      friend class SparseRow<val_t,coord_t>;

      /// array of entries.
      typedef std::vector<val_t> storage_t; 
      storage_t storage;

      /** Print a textual representation of the row to stream @c o. */
      virtual void print_on(std::ostream& o) const;

      /** Adjust @c size_ and @c starting_column to reflect the actual
          contents of the stored row.  Return pointer to adjusted row.
          If this is a null row, then delete it and return NULL. */
      DenseRow<val_t,coord_t>* adjust();

#ifdef WITH_MPI
      /** Default constructor, needed by boost::serialize.
          Initializes the row with invalid values; should call
          serialize immediately after!  */
      DenseRow();

      // serialization support (needed by Boost.MPI)
      friend class boost::serialization::access;
      template<class Archive> void serialize(Archive& ar, const unsigned int version);
#endif
    }; // class DenseRow

}; // namespace rheinfall


  // ------- inline methods -------

#include "mathdefs.hpp"
#include "sparserow.hpp"

namespace rheinfall {

  template <typename val_t, typename coord_t>
  inline
  DenseRow<val_t,coord_t>::DenseRow(const SparseRow<val_t,coord_t>* r) 
    : Row_::Row(Row_::dense, r->starting_column_, 
                r->ending_column_, r->leading_term_, true),
      storage(Row_::ending_column_ - Row_::starting_column_, 0) 
  { 
    assert(std::distance(r->storage.begin(), r->storage.end()) <= storage.size());
    for (typename SparseRow<val_t,coord_t>::storage_t::const_iterator it = r->storage.begin();
         it != r->storage.end();
         ++it) 
      {
        assert(it->first > Row_::starting_column_ and it->first <= Row_::ending_column_);
        storage[size() - (it->first - (Row_::starting_column_ + 1)) - 1] = it->second;
      };
  };

  template <typename val_t, typename coord_t>
  inline
  DenseRow<val_t,coord_t>::DenseRow(const coord_t starting_column, 
                                    const size_t ending_column) 
    : Row_::Row(Row_::dense, starting_column, ending_column, 0, true),
      storage(ending_column - starting_column, 0)
  { 
    // nothing to do
  };

#ifdef WITH_MPI
  template <typename val_t, typename coord_t>
  inline
  DenseRow<val_t,coord_t>* 
  DenseRow<val_t,coord_t>::new_from_mpi_message(mpi::communicator& comm, 
                                                const int source, const int tag)
  {
    DenseRow<val_t,coord_t>* row = new DenseRow<val_t,coord_t>();
    comm.recv(source, tag, row);
    return row;
  };

  template <typename val_t, typename coord_t>
  inline
  DenseRow<val_t,coord_t>::DenseRow() 
    : Row_::Row(Row_::dense, -1, -1, 0, false), storage() 
  { 
    // nothing to do
  };
#endif // WITH_MPI


  template <typename val_t, typename coord_t>
  inline DenseRow<val_t,coord_t>*
  DenseRow<val_t,coord_t>::adjust()
  {
    if (0 == Row_::leading_term_) {
      // compute new starting column
      for (int j = storage.size()-1; j >= 0; --j)
        if (storage[j] != 0) {
          Row_::leading_term_ = storage[j];
          Row_::starting_column_ += (storage.size() - j);
          storage.erase(storage.begin()+j, storage.end());
          return this;
        };
      // no nonzero element found in storage,
      // this is now a null row
      delete this;
      return NULL;
    }
    else
      return this;
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
    assert(this->Row_::starting_column_ == other->Row_::starting_column_);
    assert(this->Row_::ending_column_ == other->Row_::ending_column_);
    assert(0 != this->Row_::leading_term_);
    assert(0 != other->Row_::leading_term_);
    assert(this->size() == other->size());

    val_t a, b;
    rheinfall::get_row_multipliers<val_t>(this->Row_::leading_term_, 
                                          other->Row_::leading_term_, 
                                          a, b);

    DenseRow<val_t,coord_t>* result = new DenseRow(1 + this->Row_::starting_column_, 
                                                   this->Row_::ending_column_);
    assert(result->size() == this->size() - 1);
    size_t this_j = 0; 
    size_t other_j = 0; 
    for (size_t j = 0; j < this->size()-1; ++j, ++this_j, ++other_j)
       result->storage[j] = b*other->storage[other_j] + a*this->storage[this_j];
    result->Row_::leading_term_ = b*other->storage[other_j] + a*this->storage[this_j];
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
      return storage[size() - (col - (Row_::starting_column_ + 1)) - 1];
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
    ar & boost::serialization::base_object<Row_>(*this) & storage;
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
      storage[storage.size() - (col - (Row_::starting_column_ + 1)) - 1] = value;
  };


  template <typename val_t, typename coord_t>
  inline size_t 
  DenseRow<val_t,coord_t>::size() const 
  { 
    return storage.size(); 
  };


}; // namespace rheinfall


#endif // RF_DENSEROW_HPP

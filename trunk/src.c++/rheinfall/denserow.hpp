/**
 * @file   denserow.hpp
 *
 * Interface of the denserow class.
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

#ifndef RF_DENSEROW_HPP
#define RF_DENSEROW_HPP


#include "config.hpp"
#include "row.hpp"

#include <boost/mpl/bool.hpp>
namespace mpl = boost::mpl;

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
#include <memory> // std::allocator
#include <stdexcept>
#include <vector>


namespace rheinfall {

  // forward declarations to avoid circular dependency
  template <typename val_t, typename coord_t, template<typename T> class allocator> class LU;
  template <typename val_t, typename coord_t, template<typename T> class allocator> class Rank;
  template <typename val_t, typename coord_t, template<typename T> class allocator> class Row;
  template <typename val_t, typename coord_t, template<typename T> class allocator> class SparseRow;
  template <typename val_t, typename coord_t, template<typename T> class allocator> class DenseRow;


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator = std::allocator>
    /** Store a matrix row in a linear vector of consecutive entries. */
  class DenseRow : public Row<val_t,coord_t,allocator>
    {
    public:
      /** Constructor, copying entries from a @c SparseRow instance. */
      DenseRow(const SparseRow<val_t,coord_t,allocator>* r);

#ifdef WITH_MPI
      /** Make a @a DenseRow instance from an MPI message payload.
          The MPI message is identified by the triple
          communicator/source/tag, that is passed unchanged to @c
          boost::mpi::communicator::recv(). */
      static
      DenseRow<val_t,coord_t,allocator>* new_from_mpi_message(mpi::communicator& comm, 
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

    /** Return a pointer to a linear combination of this row with row
        @c other.  The returned pointer could possibly be the @p other
        row if in-place update took place.  The two coefficients used
        for the linear combination are given by @a a (multiplier for
        this row) and @a b (multiplier for the @a other row).

        If this row is sparse and its fill-in percentage exceeds @a
        dense_threshold, convert it to dense format storage before
        performing elimination. */
      DenseRow<val_t,coord_t,allocator>* linear_combination(const SparseRow<val_t,coord_t,allocator>* other,
                                                            const val_t& a, const val_t& b, const bool adjust) const;

    /** Return a pointer to a linear combination of this row with row
        @c other.  The returned pointer could possibly be the @p other
        row if in-place update took place.  The two coefficients used
        for the linear combination are given by @a a (multiplier for
        this row) and @a b (multiplier for the @a other row).

        If this row is sparse and its fill-in percentage exceeds @a
        dense_threshold, convert it to dense format storage before
        performing elimination. */
      DenseRow<val_t,coord_t,allocator>* linear_combination(DenseRow<val_t,coord_t,allocator>* other, 
                                                            const val_t& a, const val_t& b, const bool adjust) const;

    protected:
      typedef Row<val_t,coord_t,allocator> Row_; //< Nickname for base class; used to shorten templatized expressions

      friend class LU<val_t,coord_t,allocator>;
      friend class Rank<val_t,coord_t,allocator>;
      friend class SparseRow<val_t,coord_t,allocator>;

      /// array of entries.
      typedef std::vector<val_t, allocator<val_t> > storage_t; 
      storage_t storage;

      /** Print a textual representation of the row to stream @c o. */
      virtual void print_on(std::ostream& o) const;

      /** Adjust @c size_ and @c starting_column to reflect the actual
          contents of the stored row.  Return pointer to adjusted row.
          If this is a null row, then delete it and return NULL. */
      DenseRow<val_t,coord_t,allocator>* adjust();

      /** Perform Gaussian Elimination, adding a suitable multiple of
          @a this to @p other, in-place.  Return pointer to @p other. */
      DenseRow<val_t,coord_t,allocator>* linear_combination_impl(mpl::true_ inplace_update, 
                                                                 DenseRow<val_t,coord_t,allocator>* restrict other,
                                                                 const val_t& a, const val_t& b,
                                                                 const bool adjust) const; 
      /** Perform Gaussian Elimination, storing a suitable linear
          combination of @a this and @p other into a newly-allocated
          @c DenseRow.  Return pointer the new row and delete @p other. */
      DenseRow<val_t,coord_t,allocator>* linear_combination_impl(mpl::false_ inplace_update, 
                                                                 const DenseRow<val_t,coord_t,allocator>* restrict other,
                                                                 const val_t& a, const val_t& b,
                                                                 const bool adjust) const;

#ifdef WITH_MPI
      /** Default constructor, needed by boost::serialize.
          Initializes the row with invalid values; should call
          serialize immediately after!  */
      DenseRow();

      // serialization support (needed by Boost.MPI)
      friend class boost::serialization::access;
      template<class Archive> void serialize(Archive& ar, const unsigned int version);
#endif
      
#ifndef NDEBUG
    private:
      /** Correctness checks on data structure. */
      bool __ok() const;
#endif 
  }; // class DenseRow

}; // namespace rheinfall


  // ------- inline methods -------

#include "types.hpp"
#include "sparserow.hpp"

namespace rheinfall {

  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline
  DenseRow<val_t,coord_t,allocator>::
  DenseRow(const SparseRow<val_t,coord_t,allocator>* restrict r) 
    : Row_::Row(Row_::dense, r->starting_column_, r->ending_column_, r->leading_term_),
      storage(this->ending_column_ - this->starting_column_, 0) 
  { 
    assert(std::distance(r->storage.begin(), r->storage.end()) <= storage.size());
    for (typename SparseRow<val_t,coord_t,allocator>::storage_t::const_iterator it = r->storage.begin();
         it != r->storage.end();
         ++it) 
      {
        assert(it->first > this->starting_column_ and it->first <= this->ending_column_);
        storage[size()-1 - (it->first - (this->starting_column_ + 1))] = it->second;
      };
    assert(this->__ok());
  };

  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline
  DenseRow<val_t,coord_t,allocator>::
  DenseRow(const coord_t starting_column, 
           const size_t ending_column) 
    : Row_::Row(Row_::dense, starting_column, ending_column, 0),
      storage(ending_column - starting_column, 0)
  { 
    // nothing to do
  };

#ifdef WITH_MPI
  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline
  DenseRow<val_t,coord_t,allocator>* 
  DenseRow<val_t,coord_t,allocator>::
  new_from_mpi_message(mpi::communicator& comm, 
                       const int source, const int tag)
  {
    DenseRow<val_t,coord_t,allocator>* row = new DenseRow<val_t,coord_t,allocator>();
    comm.recv(source, tag, *row);
    assert(row->__ok());
    return row;
  };

  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline
  DenseRow<val_t,coord_t,allocator>::
  DenseRow() 
    : Row_::Row(Row_::dense, -1, -1, 0), storage() 
  { 
    // nothing to do
  };
#endif // WITH_MPI


#ifndef NDEBUG
  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline bool
  DenseRow<val_t,coord_t,allocator>::
  __ok() const
  {
    assert(0 <= this->starting_column_);
    assert(this->starting_column_ <= this->ending_column_);
    //assert(not is_zero(this->leading_term_));
    assert(storage.size() == this->ending_column_ - this->starting_column_);
    return true;
  };
#endif


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline DenseRow<val_t,coord_t,allocator>*
  DenseRow<val_t,coord_t,allocator>::
  adjust()
  {
    if (is_zero(this->leading_term_)) {
      // compute new starting column
      for (int j = storage.size()-1; j >= 0; --j)
        if (not is_zero(storage[j])) {
          this->leading_term_ = storage[j];
          this->starting_column_ += (storage.size() - j);
          storage.erase(storage.begin()+j, storage.end());
          assert(this->__ok());
          return this;
        };
      // no nonzero element found in storage,
      // this is now a null row
      delete this;
      return NULL;
    }
    else {
      assert(this->__ok());
      return this;
    };
  };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline DenseRow<val_t,coord_t,allocator>* 
  DenseRow<val_t,coord_t,allocator>::
  linear_combination(const SparseRow<val_t,coord_t,allocator>* restrict other,
                     const val_t& a, const val_t& b, const bool adjust) 
    const restrict_this
  {
    assert(this->starting_column_ == other->starting_column_);
    assert(this->ending_column_ == other->ending_column_);
    assert(not is_zero(this->leading_term_));
    assert(not is_zero(other->leading_term_));

    // convert `other` to dense storage upfront: adding the
    // non-zero entries from `this` would made it pretty dense
    // anyway
    DenseRow<val_t,coord_t,allocator>* dense_other(new DenseRow<val_t,coord_t,allocator>(other));
    delete other;
    return this->linear_combination(dense_other, a, b, adjust);
  }; // DenseRow<val_t,coord_t,allocator>* linear_combination(SparseRow<val_t,coord_t,allocator>* other)


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline DenseRow<val_t,coord_t,allocator>*
  DenseRow<val_t,coord_t,allocator>::
  linear_combination(DenseRow<val_t,coord_t,allocator>* restrict other,
                     const val_t& a, const val_t& b, const bool adjust) 
    const restrict_this
  {
    assert(this->starting_column_ == other->starting_column_);
    assert(this->ending_column_ == other->ending_column_);
    assert(not is_zero(this->leading_term_));
    assert(not is_zero(other->leading_term_));
    assert(this->size() == other->size());

    return this->linear_combination_impl(use_inplace_update<val_t>(), other, a, b, adjust);
  }; // dense_row_ptr linear_combination(dense_row_ptr other)


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline DenseRow<val_t,coord_t,allocator>*
  DenseRow<val_t,coord_t,allocator>::
  linear_combination_impl(mpl::true_ inplace_update, 
                            DenseRow<val_t,coord_t,allocator>* restrict other,
                          const val_t& a, const val_t& b, const bool adjust) 
    const restrict_this
  {
    for (size_t j = 0; j <= this->size()-1; ++j) {
      other->storage[j] *= b;
      other->storage[j] += a*this->storage[j];
    };
#ifdef RF_COUNT_OPS
    if (this->ops_count_ptr != NULL)
      *(this->ops_count_ptr) += 3 * this->size();
#endif
    other->leading_term_ = 0;

    DenseRow<val_t,coord_t,allocator>* const result = 
      (adjust?
       other->adjust() // update done, adjust size and starting column
       : other);
    assert(NULL == result 
           or (adjust and result->starting_column_ > this->starting_column_));
    return result;
  }; // dense_row_ptr linear_combination_impl(mpl::true_, dense_row_ptr other)


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline DenseRow<val_t,coord_t,allocator>*
  DenseRow<val_t,coord_t,allocator>::
  linear_combination_impl(mpl::false_ inplace_update, 
                            const DenseRow<val_t,coord_t,allocator>* restrict other,
                          const val_t& a, const val_t& b, const bool adjust) 
    const restrict_this
  {
    DenseRow<val_t,coord_t,allocator>* restrict result = 
      new DenseRow(1 + this->starting_column_, this->ending_column_);
    assert(result->size() == this->size() - 1);
    size_t j = 0;
    for (; j < this->size()-1; ++j) {
       result->storage[j] = b*other->storage[j] + a*this->storage[j];
#ifdef RF_COUNT_OPS
       if ((this->ops_count_ptr) != NULL)
         *(this->ops_count_ptr) += 3;
#endif
    };
    result->leading_term_ = b*other->storage[j] + a*this->storage[j];
#ifdef RF_COUNT_OPS
    if ((this->ops_count_ptr) != NULL)
      *(this->ops_count_ptr) += 3;
#endif
    delete other;

    result = (adjust?
              result->adjust() // update done, adjust size and starting column
              : result);
    assert(NULL == result 
           or (adjust and result->starting_column_ > this->starting_column_));
    return result;
  }; // dense_row_ptr linear_combination_impl(mpl::false_, dense_row_ptr other)


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline val_t
  DenseRow<val_t,coord_t,allocator>::
  get(const coord_t col) const
  {
    assert(col >= this->starting_column_ and col <= this->ending_column_);
    if (col == this->starting_column_)
      return this->leading_term_;
    else
      return storage[size() - (col - (this->starting_column_ + 1)) - 1];
  };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline void
  DenseRow<val_t,coord_t,allocator>::
  print_on(std::ostream& out) const 
  {
    out << "["
        << this->starting_column_ <<":"<< this->leading_term_;
    for (int j = storage.size()-1; j >= 0; --j)
      out <<","<< storage[j];
    out << "]";
  };


#ifdef WITH_MPI
  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  template<class Archive>
  inline void 
  DenseRow<val_t,coord_t,allocator>::
  serialize(Archive& ar, const unsigned int version) 
  {
    assert(version == 0);
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    ar & boost::serialization::base_object<Row_>(*this) & storage;
  };
#endif // WITH_MPI


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline void
  DenseRow<val_t,coord_t,allocator>::
  set(const coord_t col, const val_t value) 
  {
    assert(col >= this->starting_column_ and col <= this->ending_column_);
    if (col == this->starting_column_)
      this->leading_term_ = value;
    else
      storage[storage.size()-1 - (col - (this->starting_column_ + 1))] = value;
  };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline size_t 
  DenseRow<val_t,coord_t,allocator>::
  size() const 
  { 
    return storage.size(); 
  };


}; // namespace rheinfall


#endif // RF_DENSEROW_HPP

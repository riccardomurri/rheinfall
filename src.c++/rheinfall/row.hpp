/**
 * @file   row.hpp
 *
 * The rheinfall::row class.
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

#ifndef RF_ROW_HPP
#define RF_ROW_HPP


#include "config.hpp"

#ifdef WITH_MPI
// class needs to be serializable for Boost.MPI to handle it
# include <boost/serialization/access.hpp>
# ifdef WITH_GMPXX
#  include <gmpxx_boost_serialization.hpp>
# endif
#endif 

#include <cassert>
#include <iostream>
#include <memory> // std::allocator
#include <stdexcept>


namespace rheinfall {

  // declaration of Row-based classes
  template <typename val_t, typename coord_t, template<typename T> class allocator> class Row;
  template <typename val_t, typename coord_t, template<typename T> class allocator> class SparseRow;
  template <typename val_t, typename coord_t, template<typename T> class allocator> class DenseRow;


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator = std::allocator>
  /** Abstract base class for matrix rows. */
  class Row 
  {
  public:
    /** The two Row subclasses uses in the Rheinfall algorithm. */
    typedef enum { sparse, dense } kind_t;
    /** Discriminating field to tell whether an instance of a Row
        subclass is a @ref SparseRow or a @ref DenseRow */
    const kind_t kind;

    /** Original row number. (As read from the file or the original matrix.) */
    coord_t h0;

    /** Row number after the initial reordering step. */
    coord_t h1;

    /** Final row number. */
    coord_t h2;

    /** Constructor. */
    Row(const kind_t kind_, 
        const coord_t starting_column, const coord_t ending_column, 
        const val_t leading_term);

    /** Virtual destructor (for actual use in subclasses). */
    virtual ~Row();

    /** Add @a value to the entry in column @a col. */
    virtual void add_to_entry(const coord_t col, const val_t value) = 0;

    /** Return a copy of the element stored at column @c col */
    virtual val_t get(const coord_t col) const = 0;

    /** Return index of first column that has a nonzero entry. */
    virtual coord_t first_nonzero_column() const;

    /** Set the element stored at column @c col */
    virtual void set(const coord_t col, const val_t value) = 0;

    /** Return number of allocated (i.e., non-zero) entries. */
    virtual size_t size() const = 0;

    /** Perform Gaussian elimination: sum a multiple of this row to
        (a multiple of) row @c r so that the combination has its
        first nonzero entry at a column index strictly larger than
        the one of both rows. Return pointer to the combined row, which
        could possibly be this row if in-place update took place. */
    Row* gaussian_elimination(Row* other, 
                              const float dense_threshold) const;

    /** Print a textual representation of the row to stream @c out.
        See http://www.parashift.com/c++-faq-lite/input-output.html#faq-15.11 */
    friend std::ostream& operator<< (std::ostream& out, 
                                     Row<val_t,coord_t,allocator> const& row) {
      row.print_on(out);
      return out;
    };

    // XXX: the following is public only to allow the test programs to
    // check the internal state of objects.  Adjust the tests so we
    // can have proper protection of the fields.
  public:
    /** Index of the first column where a nonzero element is stored. */
    coord_t starting_column_; // would-be `const`: can only be modified by ctor and serialization
    /** Index of the last column in the matrix. */
    coord_t ending_column_; // would-be `const`: can only be modified by ctor and serialization
    /** Value of the first nonzero element. (Stored separately from
        the rest for speed reasons.) */
    val_t leading_term_; // would-be `const`: can only be modified by ctor and serialization

    /** Print a textual representation of the row to stream @c out.
        Provided to actually implement @c operator<< in derived classes; see
        http://www.parashift.com/c++-faq-lite/input-output.html#faq-15.11 */
    virtual void print_on(std::ostream& out) const = 0;

#ifdef WITH_MPI
    friend class boost::serialization::access;
    /// boost::serialization support
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version);
#endif

  public:
#ifdef RF_COUNT_OPS
    unsigned long *ops_count_ptr;
#endif 

  }; // class Row

}; // namespace rheinfall


  // ------- inline methods -------
  
#include "sparserow.hpp"
#include "denserow.hpp"

namespace rheinfall {

  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline
  Row<val_t,coord_t,allocator>::
  Row(const kind_t kind_, 
      const coord_t starting_column, 
      const coord_t ending_column, 
      const val_t leading_term)
    : kind(kind_)
    , h0(-1)  // XXX: use boost::optional instead
    , h1(-1)
    , h2(-1)
    , starting_column_(starting_column)
    , ending_column_(ending_column)
    , leading_term_(leading_term)
#ifdef RF_COUNT_OPS
    , ops_count_ptr(NULL)
#endif
    { 
      // init-only ctor
    };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline
  Row<val_t,coord_t,allocator>::
  ~Row()
  {
    // override in sub-classes
  };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline coord_t 
  Row<val_t,coord_t,allocator>::
  first_nonzero_column() const 
  {
    assert(not is_zero(leading_term_));
    return starting_column_;
  };


  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  inline Row<val_t,coord_t,allocator>*
  Row<val_t,coord_t,allocator>::
  gaussian_elimination(Row<val_t,coord_t,allocator>* restrict other, 
                       const float dense_threshold) 
    const restrict_this
  {
    if ((sparse == this->kind) and (sparse == other->kind)) {
      SparseRow<val_t,coord_t,allocator> const* s1 = static_cast<SparseRow<val_t,coord_t,allocator> const*>(this);
      SparseRow<val_t,coord_t,allocator>* s2 = static_cast<SparseRow<val_t,coord_t,allocator>*>(other);
      if (s2->fill_in() > dense_threshold) {
        // FIXME: could merge ctor+elimination in one funcall
        DenseRow<val_t,coord_t,allocator>* dense_other(new DenseRow<val_t,coord_t,allocator>(s2));
        delete other;
        return s1->gaussian_elimination(dense_other);
      }
      else { // `other` kept sparse
        return s1->gaussian_elimination(s2);
      }; // if (fill_in > ...)
    }
    else if ((sparse == this->kind) and (dense == other->kind)) {
      SparseRow<val_t,coord_t,allocator> const* s = static_cast<SparseRow<val_t,coord_t,allocator> const*>(this);
      DenseRow<val_t,coord_t,allocator>* d = static_cast<DenseRow<val_t,coord_t,allocator>*>(other);
      return s->gaussian_elimination(d);
    }
    else if ((dense == this->kind) and (sparse == other->kind)) {
      DenseRow<val_t,coord_t,allocator> const* d = static_cast<DenseRow<val_t,coord_t,allocator> const*>(this);
      SparseRow<val_t,coord_t,allocator>* s = static_cast<SparseRow<val_t,coord_t,allocator>*>(other);
      return d->gaussian_elimination(s);
    }
    else if ((dense == this->kind) and (dense == other->kind)) {
      DenseRow<val_t,coord_t,allocator> const* d1 = static_cast<DenseRow<val_t,coord_t,allocator> const*>(this);
      DenseRow<val_t,coord_t,allocator>* d2 = static_cast<DenseRow<val_t,coord_t,allocator>*>(other);
      return d1->gaussian_elimination(d2);
    }
    else
      // should not happen!
      throw std::logic_error("Unhandled row type combination in Row::gaussian_elimination()");
  };


#ifdef WITH_MPI
  template <typename val_t, typename coord_t, 
            template<typename T> class allocator>
  template<class Archive>
  inline void 
  Row<val_t,coord_t,allocator>::
  serialize(Archive& ar, const unsigned int version) 
  {
    assert(version == 0);
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    ar 
      & h0
      & h1
      & h2
      & starting_column_
      & ending_column_ 
      & leading_term_;
    assert(starting_column_ >= 0);
    assert(ending_column_ >= 0);
    assert (starting_column_ <= ending_column_);
    assert(not is_zero(leading_term_));
  }; // Row::serialize(...)
#endif


}; // namespace rheinfall


#endif // RF_ROW_HPP

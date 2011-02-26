/**
 * @file   ge.c
 *
 * Implementation of Gaussian Elimination functions on sparse and dense rows.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2010 riccardo.murri@gmail.com.  All rights reserved.
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

#include "config.h"
#include "row.h"

#include <stdbool.h>


/* Needed arithmetic operations */
static inline
val_t _gcd(intmax_t const n, intmax_t const m) 
{ 
  return m==0? n : _gcd(m, n%m);
};

static inline
val_t gcd(intmax_t const n, intmax_t const m) 
{ 
  return _gcd(val_abs(n), val_abs(m)); 
};



sparse_row_t* 
gaussian_elimination_sparse_with_sparse_pivot(sparse_row_t* const pivot_row, 
                                              sparse_row_t* other)
{
  assert(sparse_row_ok(pivot_row));
  assert(sparse_row_ok(other));
  assert(pivot_row->starting_column_ == other->starting_column_);

  // compute:
  //   `a`: multiplier for `this` row
  //   `b`: multiplier for `other` row
  const val_t GCD = gcd(pivot_row->leading_term_, other->leading_term_);
  val_t a = - (other->leading_term_ / GCD);
  val_t b = pivot_row->leading_term_ / GCD;
  assert (0 != a and 0 != b);

  sparse_row_t result = NULL;

  size_t pivot_i = 0;
  size_t other_i = 0;
  const size_t pivot_row_size = sparse_row_size(pivot_row);
  const size_t other_row_size = sparse_row_size(other);
  // loop while one of the two indexes is still valid
  while(pivot_i < pivot_row_size || other_i < other_row_size) {
    coord_t pivot_col;
    if (pivot_i != pivot_row_size) {
      pivot_col = pivot_row->storage[pivot_i].column;
      assert (pivot_col > pivot_row->starting_column_);
    }
    else // this_i reached end of vector, use out-of-range value
      pivot_col = pivot_row->ending_column_ + 1;

    coord_t other_col;
    if (other_i != other_row_size) {
      other_col = other->storage[other_i].column;
      assert (other_col > other->starting_column_);
    }
    else 
      other_col = other->ending_column_ + 1;

    bool nonzero_found = false;
    coord_t coord;
    val_t value;
    if (other_col < this_col) {
      value = b * other->storage[other_i].value;
      assert(0 != value);
      coord = other_col;
      nonzero_found = true;
      ++other_i;
    }
    else if (other_col == this_col) {
      value = a * pivot_row->storage[pivot_i].value + b * other->storage[other_i].value;
      if (0 != value) {
        coord = this_col;
        nonzero_found = true;
      };
      ++pivot_i;
      ++other_i;
    }
    else if (other_col > this_col) {
      value = a * pivot_row->storage[pivot_i].value;
      assert(0 != value);
      coord = this_col;
      nonzero_found = true;
      ++pivot_i;
    }
    else 
      assert(0); // should not happen!
    
    if (nonzero_found) {
      if (NULL == result) {
        // allocate new SparseRow
        result = sparse_row_new(pivot_row_size + other_row_size);
        result->kind = ROW_SPARSE;
        result->starting_column_ = coord;
        result->ending_column_ = ending_column_;
        result->leading_term_ = value;
      }
      else {
        entry_t* p = sparse_row_extend1(result);
        p->column = coord;
        p->value = value;
      }; 
    }; 
  }; // while(pivot_i < pivot_row_size && ...)
#ifndef NDEBUG
  if (NULL != result) {
    assert(sparse_row_ok(result));
    assert(result->starting_column_ > pivot_row->starting_column_);
    assert(result->starting_column_ > other->starting_column_);
    assert(result->ending_column_ == pivot_row->ending_column_);
    assert(sparse_row_size(result) <= pivot_row_size + other_row_size);
  };
#endif
  sparse_row_free(other); // release old storage
  return result;
} // gaussian_elimination_sparse_with_sparse_pivot


dense_row_t* 
gaussian_elimination_dense_with_sparse_pivot(sparse_row_t* const pivot_row,
                                             dense_row_t* other)
{
  assert(pivot_row->starting_column_ == other->starting_column_);
  assert(sparse_row_ok(pivot_row));
  assert(0 != other->leading_term_);

  // compute:
  //   `a`: multiplier for `this` row
  //   `b`: multiplier for `other` row
  const val_t GCD = gcd(pivot_row->leading_term_, other->leading_term_);
  val_t a = - (other->leading_term_ / GCD);
  val_t b = pivot_row->leading_term_ / GCD;

  const size_t other_row_size = dense_row_size(other);
  for (size_t j = 0; j < other_row_size; ++j) 
    other->storage[j] *= b;
  const size_t pivot_row_size = sparse_row_size(pivot_row);
  for (size_t j = 0; j < pivot_row_size; ++j) {
    assert(pivot_row->storage[j].column > starting_column_ 
           && pivot_row->storage[j].column <= ending_column_);
    other->storage[other_row_size - (pivot_row->storage[j].column - other->starting_column_)] += a*it->second;
  };

  return dense_row_adjust(other); // content update done, adjust size and starting column
} // gaussian_elimination_dense_with_sparse_pivot


dense_row_t* 
gaussian_elimination_sparse_with_dense_pivot(dense_row_t* const pivot_row,
                                             sparse_row_t* other)
{
  assert(pivot_row->starting_column_ == other->starting_column_);
  assert(0 != pivot_row->leading_term_);
  assert(0 != other->leading_term_);

  // convert `other` to dense storage upfront: adding the non-zero
  // entries from `pivot_row` would made it pretty dense anyway
  dense_row_t* dense_other = dense_row_new_from_sparse_row(other);
  sparse_row_free(other);
  return gaussian_elimination_dense_with_dense_pivot(pivot_row, dense_other);
}


dense_row_t* 
gaussian_elimination_dense_with_dense_pivot(dense_row_t* const pivot_row,
                                            dense_row_t* other)
{
  assert(pivot_row->starting_column_ == other->starting_column_);
  assert(0 != pivot_row->leading_term_);
  assert(0 != other->leading_term_);
  assert(dense_row_size(pivot_row) == dense_row_size(other));

  // compute:
  //   `a`: multiplier for `this` row
  //   `b`: multiplier for `other` row
  const val_t GCD = gcd(pivot_row->leading_term_, other->leading_term_);
  val_t a = - (other->leading_term_ / GCD);
  val_t b = pivot_row->leading_term_ / GCD;

  const size_t pivot_row_size = dense_row_size(pivot_row);
  for (size_t j = 0; j < pivot_row_size; ++j) {
    // XXX: is it faster to allocate new storage and fill it with `a*x+b*y`?
    other->storage[j] *= a;
    other->storage[j] += b * storage[j];
  };
  return dense_row_adjust(other); // update done, adjust size and starting column
}

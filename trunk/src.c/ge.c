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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "ge.h"

#include "common.h"
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



/** Perform Gaussian elimination, combining a sparse row @a other with
    a sparse @a pivot_row. Return pointer to newly-allocated dense row,
    and frees pointer to @a other.  @see gaussian_elimination */
static inline sparse_row_t* 
gaussian_elimination_sparse_with_sparse_pivot(sparse_row_t* const restrict pivot_row, 
                                              sparse_row_t* restrict other_row)
{
  assert(sparse_row_ok(pivot_row));
  assert(sparse_row_ok(other_row));
  assert(pivot_row->starting_column_ == other_row->starting_column_);
  assert(pivot_row->ending_column_ == other_row->ending_column_);

  // compute:
  //   `a`: multiplier for `this` row
  //   `b`: multiplier for `other_row` row
  const val_t GCD = gcd(other_row->leading_term_, pivot_row->leading_term_);
  val_t a = - (other_row->leading_term_ / GCD);
  val_t b = pivot_row->leading_term_ / GCD;
  assert (0 != a && 0 != b);

  sparse_row_t* result = NULL;

  entry_t* pivot_row_cur = sparse_row_at(pivot_row, 0);
  entry_t* other_row_cur = sparse_row_at(other_row, 0);
  const entry_t* pivot_row_ub = sparse_row_ub(pivot_row);
  const entry_t* other_row_ub = sparse_row_ub(other_row);
  // loop while one of the two indexes is still valid
  while(pivot_row_cur < pivot_row_ub || other_row_cur < other_row_ub) {
    coord_t pivot_col;
    if (pivot_row_cur < pivot_row_ub) {
      pivot_col = pivot_row_cur->column;
      assert (pivot_col > pivot_row->starting_column_);
    }
    else // this_i reached end of vector, use out-of-range value
      pivot_col = pivot_row->ending_column_ + 1;

    coord_t other_col;
    if (other_row_cur < other_row_ub) {
      other_col = other_row_cur->column;
      assert (other_col > other_row->starting_column_);
    }
    else 
      other_col = other_row->ending_column_ + 1;

    bool nonzero_found = false;
    coord_t coord;
    val_t value;
    if (other_col < pivot_col) {
      value = b * other_row_cur->value;
      assert(0 != value);
      coord = other_col;
      nonzero_found = true;
      ++other_row_cur;
    }
    else if (other_col == pivot_col) {
      value = a * pivot_row_cur->value + b * other_row_cur->value;
      if (0 != value) {
        coord = pivot_col;
        nonzero_found = true;
      };
      ++pivot_row_cur;
      ++other_row_cur;
    }
    else if (other_col > pivot_col) {
      value = a * pivot_row_cur->value;
      assert(0 != value);
      coord = pivot_col;
      nonzero_found = true;
      ++pivot_row_cur;
    }
    else 
      assert(0); // should not happen!
    
    if (nonzero_found) {
      if (NULL == result) {
        // allocate new SparseRow
        result = sparse_row_alloc(sparse_row_size(pivot_row) + sparse_row_size(other_row));
        result->starting_column_ = coord;
        result->ending_column_ = other_row->ending_column_;
        result->leading_term_ = value;
      }
      else {
        entry_t* p = sparse_row_extend1(&result);
        p->column = coord;
        p->value = value;
      }; 
    }; 
  }; // while(pivot_i < pivot_row_size && ...)
#ifndef NDEBUG
  if (NULL != result) {
    assert(sparse_row_ok(result));
    assert(result->starting_column_ > pivot_row->starting_column_);
    assert(result->starting_column_ > other_row->starting_column_);
    assert(result->ending_column_ == pivot_row->ending_column_);
    assert(sparse_row_size(result) <= sparse_row_size(pivot_row) + sparse_row_size(other_row));
  };
#endif
  sparse_row_free(other_row); // release old storage
  return result;
} // gaussian_elimination_sparse_with_sparse_pivot


/** Perform Gaussian elimination, combining a dense row @a other with
    a sparse @a pivot_row. Return pointer to newly-allocated dense row,
    and frees pointer to @a other.  @see gaussian_elimination */
static inline dense_row_t* 
gaussian_elimination_dense_with_sparse_pivot(sparse_row_t* const restrict pivot_row,
                                             dense_row_t* restrict other_row)
{
  assert(pivot_row->starting_column_ >= other_row->starting_column_);
  assert(pivot_row->ending_column_ == other_row->ending_column_);
  assert(sparse_row_size(pivot_row) <= dense_row_size(other_row));
  assert(sparse_row_ok(pivot_row));
  assert(0 != other_row->leading_term_);

  // compute:
  //   `a`: multiplier for `this` row
  //   `b`: multiplier for `other_row` row
  const val_t GCD = gcd(other_row->leading_term_, pivot_row->leading_term_);
  val_t a = - (other_row->leading_term_ / GCD);
  val_t b = pivot_row->leading_term_ / GCD;

  const size_t other_row_size = dense_row_size(other_row);
  for (size_t j = 0; j < other_row_size; ++j) 
    other_row->storage[j] *= b;
  const size_t pivot_row_size = sparse_row_size(pivot_row);
  for (size_t j = 0; j < pivot_row_size; ++j) {
    other_row->storage[other_row_size - (pivot_row->storage[j].column - other_row->starting_column_)] += a * pivot_row->storage[j].value;
  };

  return dense_row_adjust(other_row); // content update done, adjust size and starting column
} // gaussian_elimination_dense_with_sparse_pivot


/** Perform Gaussian elimination, updating a dense row @a other
    in-place with a dense @a pivot_row.  Return pointer to the updated
    row @a other.  @see gaussian_elimination */
static inline dense_row_t* 
gaussian_elimination_dense_with_dense_pivot(dense_row_t* const restrict pivot_row,
                                            dense_row_t* restrict other)
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
    other->storage[j] += b * pivot_row->storage[j];
  };
  return dense_row_adjust(other); // update done, adjust size and starting column
}


/** Perform Gaussian elimination, combining a sparse row @a other with
    a dense @a pivot_row into a newly-allocated dense row.  Return
    pointer to the combined row and frees pointer to @a other.  @see
    gaussian_elimination */
static inline dense_row_t* 
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


row_t*
gaussian_elimination(row_t* const restrict pivot_row, row_t* restrict other_row, 
                     const double dense_threshold)
{
  if ((ROW_SPARSE == pivot_row->kind) && (ROW_SPARSE == other_row->kind)) {
    sparse_row_t* s1 = (sparse_row_t*)(pivot_row->data);
    sparse_row_t* s2 = (sparse_row_t*)(other_row->data);
    if (sparse_row_fill_in(s2) > dense_threshold) {
      // XXX: could merge ctor+elimination in one funcall only
      dense_row_t* dense_other = dense_row_new_from_sparse_row(s2);
      sparse_row_free(s2);
      other_row->kind = ROW_DENSE;
      other_row->data = gaussian_elimination_dense_with_sparse_pivot(s1, dense_other);
    }
    else { 
      // `other_row` kept sparse
      other_row->data = gaussian_elimination_sparse_with_sparse_pivot(s1, s2);
    }; // if (fill_in > ...)
  }
  else if ((ROW_SPARSE == pivot_row->kind) && (ROW_DENSE == other_row->kind)) {
    sparse_row_t* s = (sparse_row_t*)(pivot_row->data);
    dense_row_t* d = (dense_row_t*)(other_row->data);
    other_row->data = gaussian_elimination_dense_with_sparse_pivot(s, d);
  }
  else if ((ROW_DENSE == pivot_row->kind) && (ROW_SPARSE == other_row->kind)) {
    dense_row_t* d = (dense_row_t*)(pivot_row->data);
    sparse_row_t* s = (sparse_row_t*)(other_row->data);
    other_row->kind = ROW_DENSE;
    other_row->data = gaussian_elimination_sparse_with_dense_pivot(d, s);
  }
  else if ((ROW_DENSE == pivot_row->kind) && (ROW_DENSE == other_row->kind)) {
    dense_row_t* d1 = (dense_row_t*)(pivot_row->data);
    dense_row_t* d2 = (dense_row_t*)(other_row->data);
    other_row->data = gaussian_elimination_dense_with_dense_pivot(d1, d2);
  }
  else
    assert(0); // should not happen!
  return other_row;
};

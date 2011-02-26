/**
 * @file   row.h
 *
 * Interface and definitions for sparse and dense matrix rows
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

#ifndef ROW_H
#define ROW_H


#include "config.h"
#include "xarray.h"

#include <stdbool.h>


typedef enum {
  ROW_SPARSE = 0,
  ROW_DENSE = 1
} row_kind_t;


/** A row consists of a pointer to a @c sparse_row_t or a @c
    dense_row_t, together with a type tag for casting the pointee to
    the correct type. */
typedef struct {
  void* data;
  row_kind_t kind;
} row_t;

/** An extensible set of rows. */
XARRAY_DECLARE(rows_list, row_t*, /* no extra data */);


/********************** sparse row definition **************************/

typedef struct {
  coord_t column;
  val_t   value;
} entry_t;
XARRAY_DECLARE(sparse_row, entry_t,      \
               /* extra fields, the same for sparse and dense matrices */ \
               row_kind_t kind_;         \
               coord_t starting_column_; \
               coord_t ending_column_;   \
               val_t leading_term_;      \
               );

/** Return fill-in percentage, that is the number of actual
    nonzero entries divided by the number of potential entries.  */
inline double sparse_row_fill_in(const sparse_row_t const* row) {
  return 100.0 * sparse_row_size(row) / (row->ending_column_ - row->starting_column_ + 1);
};
/** Reset starting column and leading term to the first nonzero
    entry in the row.  Return pointer to this row, or @c NULL if
    the row consists entirely of null values.  Used for
    two-phase construction. */
extern sparse_row_t* sparse_row_adjust(sparse_row_t const* row);
/** Return a copy of the element stored at column @a col */
extern val_t sparse_row_get(const sparse_row_t const* row, const coord_t col);
/** Set the element stored at column @a col */
extern void sparse_row_set(sparse_row_t const* row, const coord_t col, const val_t value);

#ifndef NDEBUG
extern bool sparse_row_ok(sparse_row_t *row);
#endif



/********************** dense row definition ***************************/

XARRAY_DECLARE(dense_row, val_t,         \
               /* extra fields, the same for sparse and dense matrices */ \
               row_kind_t kind_;         \
               coord_t starting_column_; \
               coord_t ending_column_;   \
               val_t leading_term_;      \
               );

/** Allocate and initialize a new @c dense_row_t instance, copying
    entries from a @c sparse_row_t instance. */
extern dense_row_t* dense_row_new_from_sparse_row(sparse_row_t* const row);
/** Adjust size and @c starting_column to reflect the actual
    contents of the stored row.  Return pointer to adjusted row.
    If @a row is a null row, then free @a row and return NULL. */
extern dense_row_t* dense_row_adjust(dense_row_t* const row);
/** Return a copy of the element stored at column @a col */
inline val_t dense_row_get(dense_row_t* const row, const coord_t col) {
  assert(col >= row->starting_column_ && col <= row->ending_column_);
  if (col == row->starting_column_)
    return row->leading_term_;
  else
    return row->storage[dense_row_size(row) - col - 1];
};
/** Set the element stored at column @a col */
extern void dense_row_set(dense_row_t* const row, const coord_t col, const val_t value);



#endif // ROW_H

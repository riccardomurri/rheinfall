/**
 * @file   row.c
 *
 * Implementation of sparse and dense matrix rows.
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


#include <row.h>

#include <assert.h>


/********************** sparse row implementation **************************/

sparse_row_t* sparse_row_adjust(sparse_row_t const* row)
{
  if (sparse_row_size(row) > 0) {
    row->starting_column_ = row->storage[0].column;
    assert(0 <= row->starting_column_ 
           && row->starting_column_ < row->ending_column_);
    row->leading_term_ = row->storage[0].value;
    assert(0 != row->leading_term_);
    sparse_row_erase(row, 1);
    //assert(sparse_row_ok(row));
    return row;
  }
  else {
    // nothing in the row
    sparse_row_free(row);
    return NULL;
  };
}


#ifndef NDEBUG
bool
sparse_row_ok(sparse_row_t *row)
{
  assert(0 <= row->starting_column_);
  assert(row->starting_column_ <= row->ending_column_);
  assert(0 != row->leading_term_);
  // entries in `this->storage` are *always* ordered by increasing column index
  coord_t s1 = row->starting_column_;
  for (size_t n = 0; n < sparse_row_size(row); ++n) {
      assert(s1 < row->storage[n].column);
      assert(row->storage[n].column <= row->ending_column_);
      assert(0 != row->storage[n].value);
      s1 = row->storage[n].column;
    };
  return true;
};
#endif


val_t sparse_row_get(const sparse_row_t const* row, const coord_t col)
{
  assert((col >= row->starting_column_ && col <= row->ending_column_));
  if (col == row->starting_column_) 
    return row->leading_term_;
  if (sparse_row_size(row) == 0)
    return 0;
  // else, fast-forward to place where element is/would be stored
  size_t jj = 0;
  while (jj < sparse_row_size(row) && row->storage[jj].column < col)
    ++jj;
  if (sparse_row_size(row) == jj) {
    // end of list reached, `col` is larger than any index in this row
    return 0;
  }
  else if (col == row->storage[jj].column) 
    return row->storage[jj].value;
  else { // storage[jj].column > j > storage[jj+1].column
    return 0;
  };
}; // sparse_row_get(...)


void sparse_row_set(sparse_row_t const* row, const coord_t col, const val_t value) 
{
  assert((col >= row->starting_column_ && col <= row->ending_column_));
  if (col == row->starting_column_) {
    row->leading_term_ = value;
    return;
  };
  // else, fast-forward to place where element is/would be stored
  size_t jj = 0;
  while (jj < sparse_row_size(row) && row->storage[jj].column < col)
    ++jj;
  // set element accordingly
  if (sparse_row_size(row) == jj) {
    // end of list reached, `col` is larger than any index in this row
    entry_t *new_entry = sparse_row_extend1(row);
    new_entry->column = col;
    new_entry->value = value;
  }
  else if (col == storage[jj].column) 
    storage[jj].value = value;
  else { // storage[jj].first > col, insert new pair before `jj`
    entry_t *new_entry = sparse_row_insert(row, jj+1);
    new_entry->column = col;
    new_entry->value = value;
  };
  assert(sparse_row_ok(row));
}; // SparseRow::set(...)



/********************** dense row implementation **************************/

dense_row_t*
dense_row_new_from_sparse_row(sparse_row_t* const row)
{
  const size_t result_size = row->ending_column_ - row->starting_column_ + 1;
  dense_row_t* result = dense_row_new(result_size);
  result->starting_column_ = row->starting_column_;
  result->ending_column_ = row->ending_column_;
  result->leading_term_ = row->leading_term_;

  val_t *values = dense_row_extend(result, result_size);
  for (size_t n = 0; n < result_size; ++n)
    values[n] = 0;

  const size_t src_row_size = sparse_row_size(row);
  for (size_t n = 0; n < src_row_size; ++n) {
      assert(row->storage[n].column > result->starting_column_ 
             && row->storage[n].column <= result->ending_column_);
      values[result_size - (row->storage[n].column - result->starting_column_)] = row->storage[n].value;
    };
}


dense_row_t* 
dense_row_adjust(dense_row_t* const row)
{
  // compute new starting column
  for (int j = dense_row_size(row)-1; j >= 0; --j)
    if (row->storage[j] != 0) {
      row->leading_term_ = row->storage[j];
      row->starting_column_ += (dense_row_size(row) - j);
      dense_row_shorten(row, (dense_row_size(row) - j)); // XXX
      return row;
    };
  // no nonzero element found in storage,
  // this is now a null row
  dense_row_free(row);
  return NULL;
};


void
dense_row_set(dense_row_t* const row, const coord_t col, const val_t value) 
{
  assert(col >= row->starting_column_ and col <= row->ending_column_);
  if (col == row->starting_column_)
    row->leading_term_ = value;
  else
    row->storage[dense_row_size(row) - col - 1] = value;
};

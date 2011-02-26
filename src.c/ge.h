/**
 * @file   ge.h
 *
 * Interface file for Gaussian Elimination functions on sparse and dense rows.
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

#ifndef GE_H
#define GE_H


/** Perform Gaussian elimination. @see gaussian_elimination */
extern sparse_row_t* gaussian_elimination_sparse_with_sparse_pivot(sparse_row_t* const pivot_row, 
                                                                   sparse_row_t* other);
/** Perform Gaussian elimination. @see gaussian_elimination */
extern dense_row_t* gaussian_elimination_dense_with_sparse_pivot(sparse_row_t* const pivot_row,
                                                                 dense_row_t* other);
/** Perform Gaussian elimination. @see gaussian_elimination */
extern dense_row_t* gaussian_elimination_sparse_with_dense_pivot(dense_row_t* const pivot_row,
                                                                 sparse_row_t* other);
/** Perform Gaussian elimination. @see gaussian_elimination */
extern dense_row_t* gaussian_elimination_dense_with_dense_pivot(dense_row_t* const pivot_row,
                                                                dense_row_t* other);


/** Perform Gaussian elimination: sum a multiple of the pivot row to
    (a multiple of) row @a other so that the combination has its first
    nonzero entry at a column index strictly larger than the one of
    both rows. Return pointer to the combined row. Frees the pointer
    to row @a other, unless it was updated in-place. */
inline void*
gaussian_elimination(void* const pivot_row, void* other, 
                     row_kind_t pivot_row_kind, row_kind_t other_row_kind,
                     const double dense_threshold)
{
  
  if ((ROW_SPARSE == pivot_row_kind) && (ROW_SPARSE == other_row_kind)) {
    sparse_row_t* s1 = (sparse_row_t*)pivot_row;
    sparse_row_t* s2 = (sparse_row_t*)other_row;
    if (sparse_row_filln_in(s2) > dense_threshold) {
      // FIXME: could merge ctor+elimination in one funcall
      dense_row_t* dense_other = dense_row_new_from_sparse_row(s2);
      dense_row_free(other);
      return gaussian_elimination_dense_with_sparse_pivot(s1, dense_other);
    }
    else { // `other` kept sparse
      return gaussian_elimination_sparse_with_sparse_pivot(s1, s2);
    }; // if (fill_in > ...)
  }
  else if ((ROW_SPARSE == pivot_row_kind) && (ROW_DENSE == other_row_kind)) {
    sparse_row_t* s = (sparse_row_t*)pivot_row;
    dense_row_t* d = (dense_row_t*)other;
    return gaussian_elimination_dense_with_sparse_pivot(s, d);
  }
  else if ((ROW_DENSE == pivot_row_kind) && (ROW_SPARSE == other_row_kind)) {
    dense_row_t* d = (dense_row_t*)pivot_row;
    sparse_row_t* s = (sparse_row_t*)other;
    return gaussian_elimination_sparse_with_dense_pivot(d, s);
  }
  else if ((ROW_DENSE == pivot_row_kind) && (ROW_DENSE == other_row_kind)) {
    dense_row_t* d1 = (dense_row_t*)pivot_row;
    dense_row_t* d2 = (dense_row_t*)other;
    return gaussian_elimination_dense_with_dense_pivot(d1, d2);
  }
  else
    assert(0); // should not happen!
};


#endif // GE_H

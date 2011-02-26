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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "row.h"


/** Perform Gaussian elimination: sum a multiple of the pivot row to
    (a multiple of) row @a other so that the combination has its first
    nonzero entry at a column index strictly larger than the one of
    both rows. Update @a other_row in place and return pointer to it.
    In the process, the kind of @a other_row may change from @c
    ROW_SPARSE to @c ROW_DENSE.  In any event, the @c data pointer in
    @a other_row is a pointer to newly-allocated row data; the old one
    is freed using @c sparse_row_free or @c dense_row_free. */
extern row_t* gaussian_elimination(row_t* const pivot_row, row_t* other_row,
                                   const double dense_threshold);

#endif // GE_H

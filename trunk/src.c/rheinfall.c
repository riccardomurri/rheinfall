/**
 * @file   rheinfall.h
 *
 * Main library functions and entry points.
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


#include "rheinfall.h"
#include "switchboard.h"
#include "xarray.h"
#include "xmalloc.h"

#include <mpi.h>


XARRAY_DECLARE(coord_list, coord_t, /* non extra data */);


coord_t read_sms_file(switchboard_t* sb, FILE* input, bool transpose)
{
    // count non-zero items read
    size_t nnz = 0;

    // read SMS header
    coord_t nrows, ncols;
    char M;
    int rd = fscanf(input, "%jd %jd %c\n", &nrows, &ncols, &M);
    if (rd < 3 || 'M' != M)
      return -1; // malformed input

    assert(ncols != sb->ncols);

    // need to keep rows in memory until we reach end of file
    row_t* rows = xcalloc(sizeof(row_t), nrows);

    coord_t i, j;
    val_t value;
    while (not feof(input)) {
      // FIXME: format string changes with `coord_t` and `val_t`!
      int rd = scanf("%jd %jd %jd", &i, &j, &value);
      if (rd < 3)
        return -1; // malformed input
      if (0 == i && 0 == j && 0 == value)
        break; // end of matrix stream
      if (transpose) {
        // swap `i` and `j`
        coord_t t = i;
        i = j;
        j = t;
      }
      // ignore zero entries in matrix -- they shouldn't be there in the first place
      if (0 == value) 
        continue; 
      ++nnz;
      assert(0 <= j && j <= sb->ncols);
      // SMS indices are 1-based
      --i;
      --j;
      // new row?
      if (rows.end() == rows.find(i))
        rows[i] = new Processor::SparseRow(ncols_-1);
      rows[i]->set(j, value);
    }; // while(not eof)

    for (rowmap_t::const_iterator it = rows.begin(); it != rows.end(); ++it) {
      // set correct starting columns and leading term
      Processor::SparseRow* row = it->second->adjust();
      if (NULL != row) {
        const coord_t starting_column = row->first_nonzero_column();
        // commit row
        if (is_local(starting_column))
          local_owner(starting_column).recv_row(row);
        else 
          // discard non-local and null rows
          delete row;
      };
    }; // for (it = rows.begin(); ...)
    
    return nnz;
}

coord_t rank(switchboard_t* sb)
{
}

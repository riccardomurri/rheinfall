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


#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "rheinfall.h"

#include "comm.h"
#include "row.h"
#include "switchboard.h"
#include "vpu.h"
#include "xarray.h"

#include <xalloc.h>

#ifdef WITH_MPI
# include <mpi.h>
#endif

#include <errno.h>
#include <stdio.h>


XARRAY_DECLARE(coord_list, coord_t, /* no extra data */);


coord_t read_sms_file(switchboard_t* sb, FILE* input, coord_t* nrows_p, coord_t* ncols_p, bool transpose)
{
  assert (NULL != sb);
  assert (NULL != input);
  assert (NULL != nrows_p);
  assert (NULL != ncols_p);

  // count non-zero items read
  size_t nnz = 0;

  if (0 == *nrows_p) {
    // read SMS header
    char M;
    int rd = fscanf(input, fmtspec_coord_t " " fmtspec_coord_t " %c\n", nrows_p, ncols_p, &M);
    if (rd < 0)
      return -errno;
    else if (rd < 3 || 'M' != M)
      return -EILSEQ; // malformed input
  };
  if (*ncols_p != sb->ncols)
    return -EINVAL; 

  // need to keep rows in memory until we reach end of file
  row_t** rows = xcalloc(sizeof(row_t*), *nrows_p);

  coord_t i, j;
  val_t value;
  while (! feof(input)) {
    // FIXME: format string changes with `coord_t` and `val_t`!
    int rd = fscanf(input, fmtspec_coord_t " " fmtspec_coord_t " " fmtspec_val_t, &i, &j, &value);
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
    if (NULL == rows[i]) {
      rows[i] = XMALLOC(row_t);
      rows[i]->data = sparse_row_new(-1, *ncols_p - 1, 0);
      rows[i]->kind = ROW_SPARSE;
    }
    sparse_row_set((sparse_row_t**)&(rows[i]->data), j, value);
  }; // while(not feof...)

  for (coord_t n = 0; n < *nrows_p; ++n) {
    // skip empty rows
    if (NULL == rows[n])
      continue;
    // set correct starting columns and leading term
    sparse_row_t* r = sparse_row_adjust((sparse_row_t*) rows[n]->data);
    const coord_t starting_column = r->starting_column_;
    if (is_local(sb, starting_column))  
      // commit row
      vpu_recv_row(local_owner(sb, starting_column), r, ROW_SPARSE);
    else  
      // discard non-local and null rows
      sparse_row_free(rows[n]->data);
    free(rows[n]);
  }; // for (n = 0; n < ncols; ...)

  free(rows);

  return nnz;
}


coord_t rank(switchboard_t* sb) 
{
  // kickstart termination signal
  if (is_local(sb, 0))
    vpu_end_phase(local_owner(sb, 0));

  // collect (partial) ranks
  coord_t* r = xcalloc(sizeof(coord_t), sb->nvpus);

  size_t n0 = 0;
  while(true) {
    // n0 incremented each time a processor goes into `done` state 
    while (n0 < sb->nvpus && vpu_is_done(sb->vpu[n0])) {
      vpu_free(sb->vpu[n0]);
      sb->vpu[n0] = NULL;
      ++n0;
    };
    if (sb->nvpus == n0)
      break; // exit `while(true)` loop
    for (size_t n = n0; n < sb->nvpus; ++n) {
      r[n] = vpu_step(sb->vpu[n], sb);
    };
    comm_receive(sb);
  };

  // the partial rank is computed as the sum of all ranks computed
  // by local processors
  coord_t local_rank = 0;
  for (size_t n = 0; n < sb->nvpus; ++n)
    local_rank += r[n];

#ifdef WITH_MPI
  // wait until all processors have done running
  MPI_Barrier(sb->comm);

  // collect the partial ranks for all processes
  coord_t rank = 0;
  // FIXME: MPI_Datatype depends on `val_t`!
  MPI_Allreduce(&local_rank, &rank, 1, MPI_LONG, MPI_SUM, sb->comm);
#else
  coord_t rank = local_rank;
#endif

  // free resources
  free(r);

  return rank;
};

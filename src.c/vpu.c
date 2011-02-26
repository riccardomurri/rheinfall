/**
 * @file   vpu.h
 *
 * Interface for VPUs
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

#include "vpu.h"


vpu_t* 
vpu_new(const coord_t column)
{
  vpu_t* new_vpu = xmalloc(sizeof(vpu_t));
  new_vpu->column = column;
  new_vpu->phase = VPU_RUNNING;
  new_vpu->workset = rows_list_new();
  new_vpu->inbox = rows_list_new();
  new_vpu->u.row = NULL;
  return new_vpu;
}


void
vpu_free(vpu_t* vpu)
{
  rows_list_free(vpu->workset);
  rows_list_free(vpu->inbox);
  if (NULL != vpu->u.row)
    free(vpu->u.row);
}


long 
vpu_step(vpu_t* self)
{
  if (VPU_DONE == self->phase)
    return (NULL != u);

  /* swap inboxes to avoid race conditions */
#ifdef _OPENMP
  omp_set_lock(& self->inbox_lock);
#endif
  // swap `workset` and `inbox`, so to perform elimination
  // on the rows that have arrived since last call
  rows_list_t* swp = inbox;
  inbox = workset;
  workset = swp;
#ifdef _OPENMP
  omp_unset_lock(& self->inbox_lock);
#endif

  const size = rows_list_size(workset);
  if (size > 0) {
    int n = 0;
    if (NULL = u.row) {
      u.row = workset->storage[0].row;
      u.kind = workset->storage[0].kind;
      ++n;
    };
    for (n; n < size; ++n) {
      // perform elimination -- return NULL in case resulting row is full of zeroes
      row_t* new_row = gaussian_elimination(u.row, workset->storage[n].row, 
                                            u.kind, workset->storage[n].kind);
      // ship reduced rows to other processors
      if (NULL != new_row)
        comm_send_row(new_row);
    };
    // just keep the "first" row
    rows_list_clear(workset);
  };

  if (VPU_RUNNING == self->phase) {
#ifdef WITH_MPI
    comm_remove_completed(& self->outbox);
#endif
  }
  else if (VPU_ENDING == self->phase) {
    // pass end message along
    comm_send_end(sb, &self->outbox, self->starting_column + 1);

#ifdef WITH_MPI        
    if (outbox.size() > 0) {
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
      omp_set_lock(&mpi_send_lock_);
# endif
      // wait untill all sent messages have arrived
      comm_wait_all_and_then_free(& self->outbox);
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
      omp_unset_lock(&mpi_send_lock_);
# endif
    };
#endif
    // all done
    self->phase = VPU_DONE;
  };

  // exit with 0 only if we never processed any row
  return (NULL == self->u.data) ? 0 : 1;
}


void 
vpu_recv_row(struct vpu_s* self, void* row, row_kind_t kind)
{
#ifdef _OPENMP
  omp_set_lock(&(self->inbox_lock));
#endif
  row_t* new_row = rows_list_extend1(&(self->inbox));
  new_row->data = row;
  new_row->kind = kind;
#ifdef _OPENMP
  omp_unset_lock(&(self->inbox_lock));
#endif
}

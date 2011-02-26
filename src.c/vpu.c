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

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "vpu.h"

#include "comm.h"
#include "ge.h"
#include "row.h"
#include "switchboard.h"

#include <xalloc.h>


struct vpu_s* 
vpu_new(const coord_t column)
{
  vpu_t* new_vpu = xmalloc(sizeof(vpu_t));
  new_vpu->column = column;
  new_vpu->phase = VPU_RUNNING;
  new_vpu->workset = rows_list_alloc(0);
  new_vpu->inbox = rows_list_alloc(0);
  new_vpu->u.data = NULL;
  return new_vpu;
}


void
vpu_free(vpu_t* vpu)
{
  rows_list_free(vpu->workset);
  rows_list_free(vpu->inbox);
  if (NULL != vpu->u.data)
    free(vpu->u.data);
  free(vpu);
}


long 
vpu_step(vpu_t* self, switchboard_t* sb)
{
  if (VPU_DONE == self->phase)
    return (NULL != self->u.data);

  /* swap inboxes to avoid race conditions */
#ifdef _OPENMP
  omp_set_lock(& self->inbox_lock);
#endif
  // swap `workset` and `inbox`, so to perform elimination
  // on the rows that have arrived since last call
  rows_list_t* swp = self->inbox;
  self->inbox = self->workset;
  self->workset = swp;
#ifdef _OPENMP
  omp_unset_lock(& self->inbox_lock);
#endif

  const size_t size = rows_list_size(self->workset);
  if (size > 0) {
    int n = 0;
    if (NULL == self->u.data) {
      self->u.data = self->workset->storage[0].data;
      self->u.kind = self->workset->storage[0].kind;
      ++n;
    };
    for (; n < size; ++n) {
      // perform elimination -- return NULL in case resulting row is full of zeroes
      row_t* new_row = gaussian_elimination(&(self->u), &(self->workset->storage[n]), 
                                            DENSE_THRESHOLD);
      // ship reduced rows to other processors
      if (NULL != new_row->data)
        comm_send_row(sb, &(self->outbox), new_row);
    };
    // just keep the "first" row
    rows_list_clear(&(self->workset));
  };

  if (VPU_RUNNING == self->phase) {
#ifdef WITH_MPI
    comm_remove_completed(&(self->outbox));
#endif
  }
  else if (VPU_ENDING == self->phase) {
    // pass end message along
    comm_send_end(sb, self->column + 1);

#ifdef WITH_MPI        
    if (outbox_size(&(self->outbox)) > 0) {
# if defined(_OPENMP) and defined(WITH_MPI_SERIALIZED)
      omp_set_lock(&mpi_send_lock_);
# endif
      // wait untill all sent messages have arrived
      comm_wait_all_and_then_free(&(self->outbox));
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
  omp_set_lock(&inbox_lock);
#endif
  row_t* new_row = rows_list_extend1(&(self->inbox));
  new_row->data = row;
  new_row->kind = kind;
#ifdef _OPENMP
  omp_unset_lock(&inbox_lock);
#endif
}

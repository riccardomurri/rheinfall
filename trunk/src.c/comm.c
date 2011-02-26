/**
 * @file   comm.c
 *
 * Implementation of the communications functions.
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

#include "common.h"
#include "row.h"
#include "switchboard.h"
#include "vpu.h"

#include <error.h>
#include <xalloc.h>

#include <assert.h>
#include <stdio.h>


#ifdef _OPENMP
# include <omp.h>
#endif


#ifdef WITH_MPI
# include <mpi.h>
enum mpi_tags {
  TAG_END        = 0,
  TAG_ROW_SPARSE = 1,
  TAG_ROW_DENSE  = 2,
};
#endif


#if defined(WITH_MPI) && defined(_OPENMP)
static omp_lock_t mpi_send_lock_;
#endif

int comm_send_end(const switchboard_t* sb, coord_t dest)
{
  if (sb->ncols == dest)
    return 0;
  if (is_local(sb, dest)) {
    vpu_end_phase(vpu_for_column(sb, dest));
    return 0;
  }
#ifdef WITH_MPI
  else { // ship to remote process
# if defined(_OPENMP) && defined(WITH_MPI_SERIALIZED)
    omp_set_lock(&mpi_send_lock_);
# endif
    const int rc = MPI_Send((void*)&dest, 1, mpitype_coord_t, 
                            owner(sb, dest), TAG_END, sb->comm);
# if defined(_OPENMP) && defined(WITH_MPI_SERIALIZED)
    omp_unset_lock(&mpi_send_lock_);
# endif
    return rc;
  };
#endif
  assert(false);
  return -1; // should not happen!
}


#ifdef WITH_MPI
static void
mpi_fatal(const switchboard_t* sb, const MPI_Status* statusp, const char* msg)
{
  char mpi_error_string[MPI_MAX_ERROR_STRING];
  int len;
  MPI_Error_string(statusp->MPI_ERROR, mpi_error_string, &len);
  error(0, 0, "%s: %s", msg, mpi_error_string);
  MPI_Abort(sb->comm, statusp->MPI_ERROR);
}
#endif


int comm_send_sparse_row(const switchboard_t* sb, outbox_t* outbox, sparse_row_t* row)
{
  coord_t column = row->starting_column_;
  if (is_local(sb, column)) {
    vpu_recv_row(vpu_for_column(sb, column), row, ROW_SPARSE);
    return 0;
  }
#ifdef WITH_MPI
  else { // ship to remote process
    MPI_Request* req_p = requests_list_extend1(&(outbox->requests));
    row_t* row_p = rows_list_extend1(&(outbox->rows));
    row_p->data = row; /* will free row when MPI req is complete */
    row_p->kind = ROW_SPARSE; 
# if defined(_OPENMP) && defined(WITH_MPI_SERIALIZED)
    omp_set_lock(&mpi_send_lock_);
# endif
    /* XXX: send row as a byte string - assume homogeneous arch here */
    const int rc = MPI_Isend(row, sparse_row_ub(row) - sparse_row_lb(row), MPI_BYTE,
                             owner(sb, column), TAG_ROW_SPARSE, MPI_COMM_WORLD,
                             req_p);
# if defined(_OPENMP) && defined(WITH_MPI_SERIALIZED)
    omp_unset_lock(&mpi_send_lock_);
# endif
    return rc;
  };
#endif // WITH_MPI
  assert(false);
  return -1; // should not happen!
}


int comm_send_dense_row(const switchboard_t* sb, outbox_t* outbox, dense_row_t* row)
{
  coord_t column = row->starting_column_;
  if (is_local(sb, column)) {
    vpu_recv_row(vpu_for_column(sb, column), row, ROW_DENSE);
    return 0;
  }
#ifdef WITH_MPI
  else { // ship to remote process
    MPI_Request* req_p = requests_list_extend1(&(outbox->requests));
    row_t* row_p = rows_list_extend1(&(outbox->rows));
    row_p->data = row; /* will free row when MPI req is complete */
    row_p->kind = ROW_DENSE;
# if defined(_OPENMP) && defined(WITH_MPI_SERIALIZED)
    omp_set_lock(&mpi_send_lock_);
# endif
    /* XXX: send row as a byte string - assume homogeneoys arch here */
    const int rc = MPI_Isend(row, dense_row_ub(row) - dense_row_lb(row), MPI_BYTE,
                             owner(sb, column), TAG_ROW_DENSE, MPI_COMM_WORLD,
                             req_p);
# if defined(_OPENMP) && defined(WITH_MPI_SERIALIZED)
    omp_unset_lock(&mpi_send_lock_);
# endif
    return rc;
  };
#endif // WITH_MPI
  assert(false);
  return -1; // should not happen!
}


int comm_send_row(const switchboard_t* sb, outbox_t* outbox, const row_t* cargo)
{
  if (ROW_SPARSE == cargo->kind)
    return comm_send_sparse_row(sb, outbox, (sparse_row_t*)cargo->data);
  else if (ROW_DENSE == cargo->kind)
    return comm_send_dense_row(sb, outbox, (dense_row_t*)cargo->data);
  // control should never get here!
  assert (false); 
  return -1; 
}


outbox_t* comm_remove_completed(outbox_t* outbox)
{
  assert(NULL != outbox);
#ifdef WITH_MPI
  const int size = requests_list_size(outbox->requests);
  if (size > 0) {
    int count = 0;
    int* completed = xcalloc(sizeof(int), size);
    MPI_Status* statuses = xmalloc(sizeof(MPI_Status)*size);
# if defined(_OPENMP) && defined(WITH_MPI_SERIALIZED)
    omp_set_lock(&mpi_send_lock_);
# endif
    // check if some test messages have arrived...
    MPI_Testsome(size, outbox->requests->storage,
                 &count, completed, statuses);
# if defined(_OPENMP) && defined(WITH_MPI_SERIALIZED)
    omp_unset_lock(&mpi_send_lock_);
# endif
    // ...and copy incomplete requests into a new outbox
    outbox_t* new_outbox = outbox_new(size - count);
    for (int n = 0; n < size; ++n) {
      bool copy = true;
      // do not copy if index is in array of indices returned by MPI_Testsome()
      for (int m = 0; m < count; ++m) {
        if (completed[m] == n) {
          copy = false;
          break;
        };
      }; // for (int m = ...)
      if (copy) {
        MPI_Request* req_p = requests_list_extend1(&(new_outbox->requests));
        *req_p = outbox->requests->storage[n];
        row_t* row = rows_list_extend1(&(new_outbox->rows));
        row->data = outbox->rows->storage[n].data;
        row->kind = outbox->rows->storage[n].kind;
      }
      else { // free rows corresponding to completed requests
        if (ROW_SPARSE == outbox->rows->storage[n].kind)
          sparse_row_free((sparse_row_t*) outbox->rows->storage[n].data);
        else
          dense_row_free((dense_row_t*) outbox->rows->storage[n].data);
      };
    }; // for (int n = ...)
    outbox_free(outbox);
    // finally, free local resources and return
    free(completed);
    free(statuses);
    return new_outbox;
  };
#endif // WITH_MPI
  return outbox;
}


outbox_t* comm_wait_all_and_then_free(outbox_t* outbox)
{
#ifdef WITH_MPI
  const int size = requests_list_size(outbox->requests);
  if (size > 0) {
    MPI_Status* statuses = calloc(sizeof(MPI_Status), size);
# if defined(_OPENMP) && defined(WITH_MPI_SERIALIZED)
    omp_set_lock(&mpi_send_lock_);
# endif
    // check if some test messages have arrived...
    MPI_Waitall(size, outbox->requests->storage, statuses);
# if defined(_OPENMP) && defined(WITH_MPI_SERIALIZED)
    omp_unset_lock(&mpi_send_lock_);
# endif
    // ...free corresponding rows...
    for (int n = 0; n < size; ++n) {
      if (ROW_SPARSE == outbox->rows->storage[n].kind)
        sparse_row_free((sparse_row_t*) outbox->rows->storage[n].data);
      else
        dense_row_free((dense_row_t*) outbox->rows->storage[n].data);
    };
    outbox_free(outbox);
    // finally, free local resources and return
    free(statuses);
  };
#endif // WITH_MPI
  return NULL;
}


void
comm_receive(const switchboard_t* sb)
{
#ifdef WITH_MPI
  MPI_Status status;
  int flag = 0;
  MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
  if (status.MPI_ERROR != MPI_SUCCESS)
    mpi_fatal(sb, &status, "Got error while receiving row");
  while(flag) { /* only execute the body if `flag` is true */
    int size;
    MPI_Get_count(&status, MPI_BYTE, &size);
    void* xa = xmalloc(size);
    switch(status.MPI_TAG) {
    case TAG_ROW_SPARSE: {
      sparse_row_t* new_row = sparse_row_alloc_placed(xa, size);
      assert(NULL != new_row);
      MPI_Recv(new_row, size, MPI_BYTE,
               status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD,
               &status);
      if (status.MPI_ERROR != MPI_SUCCESS)
        mpi_fatal(sb, &status, "Got error while receiving sparse row");
# ifndef NDEBUG
      int recv_size;
      MPI_Get_count(&status, MPI_BYTE, &recv_size);
      assert(size == recv_size);
# endif
      assert(is_local(sb, new_row->starting_column_));
      vpu_recv_row(vpu_for_column(sb, new_row->starting_column_), new_row, ROW_SPARSE);
      break;
    };
    case TAG_ROW_DENSE: {
      dense_row_t* new_row = dense_row_alloc_placed(xa, size);
      assert(NULL != new_row);
      MPI_Recv(new_row, size, MPI_BYTE,
               status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD,
               &status);
      if (status.MPI_ERROR != MPI_SUCCESS)
        mpi_fatal(sb, &status, "Got error while receiving dense row");
# ifndef NDEBUG
      int recv_size;
      MPI_Get_count(&status, MPI_BYTE, &recv_size);
      assert(size == recv_size);
# endif
      assert(is_local(sb, new_row->starting_column_));
      vpu_recv_row(vpu_for_column(sb, new_row->starting_column_), new_row, ROW_DENSE);
      break;
    };
    case TAG_END: {
      // "end" message received; no new rows will be coming.
      // But some other rows could have arrived or could
      // already be in the `inbox`, so we need to make another
      // pass anyway.  All this boils down to: set a flag,
      // make another iteration, and end the loop next time.
      coord_t column;
      MPI_Recv(&column, 1, MPI_LONG,
               status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD,
               &status);
      if (status.MPI_ERROR != MPI_SUCCESS)
        mpi_fatal(sb, &status, "Got error while receiving sparse END tag");
# ifndef NDEBUG
      int recv_size;
      MPI_Get_count(&status, MPI_BYTE, &recv_size);
      assert(size == recv_size);
# endif
      assert(column >= 0 && is_local(sb, column));
      vpu_end_phase(vpu_for_column(sb, column));
      break;
    };
    }; // switch(status.MPI_TAG)
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    if (status.MPI_ERROR != MPI_SUCCESS)
      mpi_fatal(sb, &status, "Got error while receiving sparse row");
  }; // while(MPI_Iprobe)
#endif
}

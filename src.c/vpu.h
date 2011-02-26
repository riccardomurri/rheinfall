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

#ifndef VPU_H
#define VPU_H


#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "common.h"
#include "row.h"
#include "xarray.h"

// we can't `#include "switchboard.h"` because of circular dependencies,
// so just declare a pointer type
typedef struct switchboard_s* switchboard_ptr;


#ifdef _OPENMP
# include <omp.h>
#endif

#ifdef WITH_MPI
# include <mpi.h>
#endif

#include <stdbool.h>


/** Phases of the VPU algorithm.  When in @c VPU_RUNNING state, rows
    are received and eliminated in the @c step() method.  The @c
    VPU_ENDING state indicates that an "end" message has been
    received: no more rows will be coming, all pending communications
    will be waited for and then the state will transition to @c
    VPU_DONE.  In @c VPU_DONE state, it is an error to call any VPU
    method. */
typedef enum {
  VPU_RUNNING = 0,
  VPU_ENDING = 1,
  VPU_DONE = 2
} vpu_phase_t;


#ifdef WITH_MPI
/** An `outbox` is a list of pending MPI requests and associated row
    pointers (the payload of such requests). These lists are kept in
    separate xarrays in order to be able to use the xarray's storage
    directly to MPI_{Test,Wait} calls. */
XARRAY_DECLARE(requests_list, MPI_Request, /* no extra data */);
typedef struct {
  requests_list_t* requests;
  rows_list_t*     rows;
} outbox_t;

/** Allocate a new @c outbox_t structure and return a pointer to it. */
_inline outbox_t* outbox_new(size_t nmemb) {
  outbox_t* ob = xmalloc(sizeof(outbox_t));
  ob->requests = requests_list_alloc(nmemb);
  ob->rows = rows_list_alloc(nmemb);
  return ob;
};

/** Free an @c outbox_t object previously allocated with `outbox_new`. */
_inline void outbox_free(outbox_t *ob) {
  requests_list_free(ob->requests);
  rows_list_free(ob->rows);
  free(ob);
};

_inline size_t outbox_size(outbox_t* ob) {
  assert(requests_list_size(ob->requests) == rows_list_size(ob->rows));
  return rows_list_size(ob->rows);
};
#else
typedef void outbox_t;
#endif // WITH_MPI


/** VPUs perform elimination on a set of rows beginning at the same column. */
struct vpu_s {
  coord_t column;
  vpu_phase_t phase;
  /** Chosen row for performing elimination. */
  row_t u;
#ifdef _OPENMP
  omp_lock_t inbox_lock;
#endif
  /** Block of rows on which elimination is performed. */
  rows_list_t* workset; 
  /** Block of rows arriving from other VPUs. */ 
  rows_list_t* inbox;
#ifdef WITH_MPI
  /** List of sent rows, of which the receiver has not yet acknowledged receipt. */
  outbox_t* outbox;
#else
  void* outbox; // placeholder
#endif
};
//typedef struct vpu_s vpu_t;
#define vpu_t struct vpu_s


struct vpu_s* vpu_new(const coord_t column);
void vpu_free(struct vpu_s* vpu);

long vpu_step(struct vpu_s* self, switchboard_ptr sb);
_inline void vpu_end_phase(struct vpu_s* self);
_inline bool vpu_is_done(struct vpu_s* self);
void vpu_recv_row(struct vpu_s* self, void* row, row_kind_t kind);


// --- inline defs ---

_inline void vpu_end_phase(struct vpu_s* self) {
  self->phase = VPU_ENDING;
};

_inline bool vpu_is_done(struct vpu_s* self) {
  return VPU_DONE == self->phase;
};


#endif // VPU_H

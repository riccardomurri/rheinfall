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

#include "config.h"
#include "comm.h"
#include "row.h"
#include "switchboard.h"
#include "xarray.h"

#ifdef _OPENMP
# include <omp.h>
#endif


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


/** VPUs perform elimination on a set of rows beginning at the same column. */
typedef struct {
  coord_t column;
  vpu_phase_t phase;
  switchboard_t* sb;
#ifdef _OPENMP
  omp_lock_t inbox_lock;
#endif
  /** Block of rows on which elimination is performed. */
  rows_list_t* workset; 
  /** Block of rows arriving from other VPUs. */ 
  rows_list_t* inbox;
  /** Chosen row for performing elimination. */
  row_t u;
} vpu_t;

vpu_t* vpu_new(const coord_t column, switchboard_t* sb);
void vpu_free(vpu_t* vpu);

long vpu_step(vpu_t* self);
inline void vpu_end_phase(vpu_t* self);
inline bool vpu_is_done(vpu_t* self);
void vpu_recv_row(vpu_t* self, void* row, row_kind_t kind);


// --- inline defs ---

inline void vpu_end_phase(vpu_t* self) {
  self->phase = VPU_ENDING;
};

inline bool vpu_is_done(vpu_t* self) {
  return VPU_DONE == self->phase;
};


#endif // VPU_H

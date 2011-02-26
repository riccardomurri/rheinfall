/**
 * @file   switchboard.h
 *
 * Interface file for the switchboard system (connecting VPUs).
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

#ifndef SWITCHBOARD_H
#define SWITCHBOARD_H

#include "config.h"
#include "vpu.h"

#ifdef WITH_MPI
# include <mpi.h>
#endif


typedef struct {
#ifdef WITH_MPI
  MPI_Comm comm;    /**< MPI communicator to use */
  size_t comm_size; /**< MPI communicator size */
  int me;           /**< MPI rank of this process */
#endif
  int ncols;        /**< Total number of columns */
  size_t nvpus;     /**< Number of local VPUs */
  vpu_t* vpu;       /**< Array of local VPUs */
} switchboard_t;


#ifdef WITH_MPI
inline switchboard_t* switchboard_new(const int ncols, MPI_Comm comm_) {
  switchboard_t* sb = xmalloc(sizeof(switchboard_t));
  sb->comm = comm_;
  MPI_Comm_rank(sb->comm, &(sb->me));
  MPI_Comm_size(sb->comm, &(sb->comm_size));
  sb->nvpus = ncols;
  for (int n = 0; n < ncols; ++n) {
    sb->vpu[n] = vpu_new(n, sb);
  };
};
#else
inline switchboard_t* switchboard_new(const int ncols) {
  switchboard_t* sb = xmalloc(sizeof(switchboard_t));
  sb->nvpus = ncols;
  for (int n = 0; n < ncols; ++n) {
    sb->vpu[n] = vpu_new(n, sb);
  };
};
#endif


inline bool is_local(const switchboard_t* sb, const coord_t column) {
#ifdef WITH_MPI
  return (column % sb->comm_size) == sb->me;
#else
  return true;
#endif
};


inline vpu_t* local_owner(const switchboard_t* sb, const coord_t column) {
#ifdef WITH_MPI
  return vpu[column / sb->comm_size]; 
#else
  return vpu[column];
#endif
};


#ifdef WITH_MPI
inline int remote_owner(const switchboard_t* sb, const coord_t column) {
  return (column % sb->comm_size); 
};
#endif


#endif // SWITCHBOARD_H

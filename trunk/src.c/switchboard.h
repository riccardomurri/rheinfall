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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "common.h"
#include "vpu.h"

#include <xalloc.h> // from GNUlib

#include <stdbool.h>

#ifdef WITH_MPI
# include <mpi.h>
#endif


struct switchboard_s {
#ifdef WITH_MPI
  MPI_Comm comm;     /**< MPI communicator to use */
#endif
  int comm_size;     /**< MPI communicator size */
  int me;            /**< MPI rank of this process */
  int ncols;         /**< Total number of columns */
  size_t nvpus;      /**< Number of local VPUs */
  struct vpu_s** vpu;/**< Array of local VPUs */
};
//typedef struct switchboard_s switchboard_t;
#define switchboard_t struct switchboard_s


#ifdef WITH_MPI
switchboard_t* switchboard_new(const int ncols, MPI_Comm comm);
#else
switchboard_t* switchboard_new(const int ncols);
#endif
void switchboard_free(switchboard_t* sb);

_inline bool is_local(const switchboard_t* sb, const coord_t column) {
#ifdef WITH_MPI
  return (column % sb->comm_size) == sb->me;
#else
  return true;
#endif
};


_inline struct vpu_s* local_owner(const switchboard_t* sb, const coord_t column) {
#ifdef WITH_MPI
  return sb->vpu[column / sb->comm_size]; 
#else
  return sb->vpu[column];
#endif
};


#ifdef WITH_MPI
_inline int remote_owner(const switchboard_t* sb, const coord_t column) {
  return (column % sb->comm_size); 
};
#endif


#endif // SWITCHBOARD_H

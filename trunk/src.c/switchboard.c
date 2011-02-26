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


#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "switchboard.h"

#include <xalloc.h>


switchboard_t* 
switchboard_new(const int ncols
#ifdef WITH_MPI
                , MPI_Comm comm
#endif
                ) {
  switchboard_t* sb = xmalloc(sizeof(switchboard_t));
#ifdef WITH_MPI
  sb->comm = comm;
  MPI_Comm_rank(sb->comm, &(sb->me));
  MPI_Comm_size(sb->comm, &(sb->comm_size));
#else
  // simulate running with just 1 MPI rank
  sb->me = 0;
  sb->comm_size = 1;
#endif
  sb->ncols = ncols;
  sb->nvpus = 0;
  sb->vpu = xmalloc(sizeof(vpu_t*) * (1 + ncols / sb->comm_size));
  for (int n = 0; (sb->me + n*sb->comm_size) < ncols; ++n) {
    sb->vpu[n] = vpu_new(sb->me + n*sb->comm_size);
    sb->nvpus += 1;
  };
  return sb;
};


void
switchboard_free(switchboard_t* sb)
{
  for (int n = 0; n < sb->nvpus; ++n)
    if (NULL != sb->vpu[n])
      vpu_free(sb->vpu[n]);
  free(sb->vpu);
  free(sb);
}

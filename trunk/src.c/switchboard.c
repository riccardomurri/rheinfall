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
                , const int width
                , MPI_Comm comm
#endif
                ) {
  switchboard_t* sb = xmalloc(sizeof(switchboard_t));
#ifdef WITH_MPI
  sb->comm = comm;
  MPI_Comm_rank(sb->comm, &(sb->me));
  MPI_Comm_size(sb->comm, &(sb->comm_size));
  sb->w = width;
#else
  // simulate running with just 1 MPI rank
  sb->me = 0;
  sb->comm_size = 1;
  sb->w = 1;
#endif
  sb->ncols = ncols;
  // initialize `vpu` array and its count `nvpus`
#ifdef WITH_MPI
  const int nmemb = sb->w * (1 + ((ncols / sb->w) / sb->comm_size));
#else
  const int nmemb = ncols;
#endif
  sb->vpu = xcalloc(nmemb, sizeof(vpu_t*));
  int n = 0;
  for (int c = 0; c < ncols; ++c)
    if (is_local(sb, c))
      sb->vpu[n++] = vpu_new(c);
  sb->nvpus = n;
  // check that the VPUs actually fit in the array allocated for them
  assert(sb->nvpus <= nmemb);
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

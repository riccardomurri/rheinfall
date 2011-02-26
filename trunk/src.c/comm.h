/**
 * @file   comm.h
 *
 * Interface for communications functions common between VPUs and the
 * master.
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

#ifndef COMM_H
#define COMM_H

#include "config.h"
#include "row.h"
#include "xarray.h"

#ifdef WITH_MPI
# include <mpi.h>
#endif


/** An `outbox` is a list of pending MPI requests and associated row
    pointers (the payload of such requests). These lists are kept in
    separate xarrays in order to be able to use the xarray's storage
    directly to MPI_{Test,Wait} calls. */
XARRAY_DECLARE(requests_list, MPI_Request, /* no extra data */);
typedef {
  requests_list_t* requests;
  rows_list_t*     rows;
} outbox_t;

/** Allocate a new @c outbox_t structure and return a pointer to it. */
inline outbox_t* outbox_new(size_t nmemb) {
  outbox_t* ob = xmalloc(sizeof(outbox_t));
  ob->requests = requests_list_new(nmemb);
  ob->rows = rows_list_new(nmemb);
  return ob;
};

/** Free an @c outbox_t object previously allocated with `outbox_new`. */
inline void outbox_free(outbox_t *ob) {
  requests_list_free(ob->requests);
  rows_list_free(ob->rows);
  free(ob);
};

int comm_send_row(outbox_t* outbox, row_kind_t kind, void* row);
int comm_send_end(outbox_t* outbox, coord_t dest);


#endif // COMM_H

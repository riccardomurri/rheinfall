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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "common.h"
#include "row.h"
#include "switchboard.h"
#include "xarray.h"

#ifdef WITH_MPI
# include <mpi.h>
#endif


int comm_send_end(const switchboard_t* sb, coord_t dest);
int comm_send_row(const switchboard_t* sb, outbox_t* const outbox, const row_t* cargo);
void comm_receive(const switchboard_t* sb);
outbox_t* comm_remove_completed(outbox_t* outbox);
outbox_t* comm_wait_all_and_then_free(outbox_t* outbox);


#endif // COMM_H

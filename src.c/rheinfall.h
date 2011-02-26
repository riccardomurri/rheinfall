/**
 * @file   rheinfall.h
 *
 * Main library functions and entry points.
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

#ifndef RHEINFALL_H
#define RHEINFALL_H

#include "switchboard.h"

#include <mpi.h>

/** Read matrix entries from the given @a input stream, and send them
    to VPUs.  The @a input stream must be in J.-G. Dumas' "Sparse
    Matrix Stream" format. Return number of nonzero entries read. */
size_t read_sms_file(switchboard_t *sb, FILE* input);

/** Compute and return matrix rank. */
coord_t rank(switchboard_t* sb);


#endif // RHEINFALL_H

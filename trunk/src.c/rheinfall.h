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

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "common.h"
#include "switchboard.h"

#ifdef WITH_MPI
# include <mpi.h>
#endif

#include <stdbool.h>
#include <stdio.h>


/** Read matrix entries from the given @a input stream, and send them
    to VPUs.  The @a input stream must be in J.-G. Dumas' "Sparse
    Matrix Stream" format. Return number (>=0) of nonzero entries
    read, or negative if an error occurred.  Number of rows and
    columns in the matrix must be given in the locations pointed to by
    @a nrows_p and @a ncols_p; if those are 0, then the matrix
    dimensions are read from the SMS header and *stored* back into @c
    *nrows_p and @c *ncols_p.  If @c *nrows_p is non-zero, then the @c
    input stream *must not* contain a header line.
    
    @param sb      Pointer to an (already initialized) @c switchboard_t structure
    @param input   C stream from which the matrix data is read
    @param nrows_p Ptr to number of rows.  If @c *nrows_p is 0, then number
                   of rows and columns is read from SMS header and stored back
                   into @c *nrows_p.
    @param ncols_p Ptr to number of columns.  If @c *nrows_p is 0, then number
                   of columns is read from SMS header and stored back
                   into @c *ncols_p.
    @param transpose  If @c true, then the matrix is read into memory transposed.

    @return Number of nonzero entries read from @a input, or a negative number in case of errors:
            -EILSEQ  could not read SMS header from stream
            -EINVAL  number of columns in SMS header does not match number of columns in @a sb
            any other negative value: @c -errno (system error code)
*/
coord_t read_sms_file(switchboard_t* sb, FILE* input, coord_t* nrows_p, coord_t* ncols_p, bool transpose);

/** Compute and return matrix rank. */
coord_t rank(switchboard_t* sb);


#endif // RHEINFALL_H

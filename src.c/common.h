/**
 * @file   config.h
 *
 * Library-wide configuration parameters.
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

#ifndef CONFIG_H
#define CONFIG_H


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <inttypes.h>


typedef long coord_t;
#define mpitype_coord_t MPI_LONG
#define fmtspec_coord_t "%ld"

typedef long long val_t;
#define mpitype_val_t MPI_LONG_LONG
#define fmtspec_val_t "%lld"
#define val_abs(x) imaxabs(x)

#define DENSE_THRESHOLD 40.0

#if (defined(WITH_MPI_SERIALIZED) || defined(WITH_MPI_MULTIPLE)) && ! defined(WITH_MPI)
# define WITH_MPI
#endif

#define _inline static inline

#endif // CONFIG_H

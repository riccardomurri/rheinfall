/**
 * @file   irank.cpp
 *
 * Sample application of the ``rheinfall'' algorithm for computing matrix rank,
 * built with @c{double} floating-point data type.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2010 riccardo.murri@gmail.com.  All rights reserved.
 *
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

#include <config.h>

#ifndef HAVE_GMPXX
# error This source requires GMP to compile; messed up autoconf settings?
#else
#define WITH_GMPXX
#endif


#include <gmpxx.h>
typedef mpq_class val_t;
typedef long coord_t;


#include "main.cpp"

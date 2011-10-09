/**
 * @file   config.hpp
 *
 * Interface of the config class.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2010, 2011 riccardo.murri@gmail.com. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * See http://www.gnu.org/licenses/gpl.html for licence details.
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

#ifndef CONFIG_HPP
#define CONFIG_HPP


// uncomment/define to enable computing the arithmetic operations count
#define RF_ENABLE_STATS


// use WITH_MPI_SERIALIZED / WITH_MPI_MULTIPLE to select which
// threading model to request to the MPI-2 library
#if (defined(WITH_MPI_SERIALIZED) or defined (WITH_MPI_MULTIPLE)) and not defined(WITH_MPI)
# define WITH_MPI
#endif


#ifdef WITH_GMPXX
# include <gmpxx.h>
#endif // WITH_GMPXX


#endif // CONFIG_HPP

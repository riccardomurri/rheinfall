/**
 * @file   stats.hpp
 *
 * Interface of the stats class.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2011 riccardo.murri@gmail.com. All rights reserved.
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

#ifndef RF_STATS_HPP
#define RF_STATS_HPP 1


namespace rheinfall {

  /** Counts of various types of operations performed by the algorithm. */
  struct Stats {
    /// total number of arithmetic operations
    unsigned long ops_count;

    /// total number of @c SparseRow instances processed
    unsigned long sparserow_count;
    /// total number of matrix entries processed by @c SparseRow
    unsigned long sparserow_elts;

    /// total number of @c DenseRow instances processed
    unsigned long denserow_count;
    /// total number of matrix entries processed by @c DenseRow
    unsigned long denserow_elts;

    /// Default ctor: initialize all fields to 0
    Stats () :
      ops_count(0), 
      sparserow_count(0), sparserow_elts(0),
      denserow_count(0), denserow_elts(0)
    { };

  };

};


#endif // RF_STATS_HPP

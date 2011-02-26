/**
 * @file   xarray.h
 *
 * Definitions for module @c xarray
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

#ifndef XARRAY_H
#define XARRAY_H

#include <xalloc.h>

#include <assert.h>
#include <string.h> // memmove
#include <stdlib.h> // malloc, realloc, free

#define _inline static inline

/** @def XARRAY_CREATE(aname, elt_t, extra)
 *
 * Declare an xarray-like type @c aname_t and functions to operate on
 * it.  An xarray is a (malloc'd) contiguous block of memory holding
 * instances of elements of type @a elt_t; functions are provided to
 * allocate a new xarray (@c aname_alloc), append elements at the end of
 * the array, automatically relocating it if needed (@c aname_extend),
 * access/insert/remove elements at specific positions (@c aname_at,
 * @c aname_insert, @c aname_erase), wipe off all elements from an
 * array (@c aname_clear) or dispose of it and return the storage to
 * the system.
 *
 */
#define XARRAY_DECLARE(aname, elt_t, extra)                            \
                                                                       \
  typedef struct {                                                     \
    /**  user defined extra info that is fitted into the array */      \
    extra;                                                             \
    /** number of elements that can be stored in the xarray before     \
     *  realloc() must take place. */                                  \
    size_t allocated;                                                  \
    /** number of elements stored in the xarray */                     \
    size_t count;                                                      \
    /** memory space where xarray elements are stored */               \
    elt_t  storage[];                                                  \
  } aname##_t;                                                         \
                                                                       \
                                                                       \
/** Return a newly-allocated xarray, or @c NULL if the allocation      \
    failed. The xarray is initially sized to contain @c nmemb elements \
    without the need to resize. */                                     \
_inline aname##_t* aname##_alloc(const size_t nmemb) {                 \
  aname##_t* xa = (aname##_t*)xmalloc(sizeof(aname##_t)                \
                                      + sizeof(elt_t)*nmemb);          \
  if (NULL == xa)                                                      \
    return NULL;                                                       \
                                                                       \
  xa->allocated = nmemb;                                               \
  xa->count = 0;                                                       \
                                                                       \
  assert(NULL != xa);                                                  \
  assert(xa->count <= xa->allocated);                                  \
  return xa;                                                           \
};                                                                     \
                                                                       \
                                                                       \
/** Initialize a new xarray placed at address @a xa and fitting within \
    @a size bytes. Return pointer to the start of the xarray, or       \
    @c NULL if there was some error (e.g., @a size is insufficient). */\
_inline aname##_t* aname##_alloc_placed(void* xa, const size_t size) { \
  assert(NULL != xa);                                                  \
  assert(size > sizeof(aname##_t));                                    \
  if(size < sizeof(aname##_t))                                         \
    return NULL;                                                       \
                                                                       \
  aname##_t* xa_ = (aname##_t*)xa;                                     \
  xa_->allocated = (size - sizeof(aname##_t)) / sizeof(elt_t);         \
  xa_->count = 0;                                                      \
                                                                       \
  assert(xa_->count <= xa_->allocated);                                \
  return xa_;                                                          \
};                                                                     \
                                                                       \
                                                                       \
/** Return a pointer to the element in the xarray at position @c pos.  \
 Xarray positions follow the usual C convention, ranging from 0 to @c  \
 xa->count-1. */                                                       \
_inline size_t aname##_size(const aname##_t* const xa) {               \
  assert(NULL != xa);                                                  \
  return xa->count;                                                    \
};                                                                     \
                                                                       \
                                                                       \
 /** Return a pointer to the first byte in the xarray. */              \
_inline void* aname##_lb(const aname##_t* const xa) {                  \
  assert(NULL != xa);                                                  \
  return (void*)xa;                                                    \
};                                                                     \
                                                                       \
                                                                       \
 /** Return a pointer past the last used byte in the xarray. */        \
_inline void* aname##_ub(const aname##_t* const xa) {                  \
  assert(NULL != xa);                                                  \
  return (void*) &(xa->storage[xa->count]);                            \
};                                                                     \
                                                                       \
                                                                       \
/** Return a pointer to the element in the xarray at position @c pos.  \
 Xarray positions follow the usual C convention, ranging from 0 to @c  \
 xa->count-1. */                                                       \
_inline elt_t* aname##_at(aname##_t* const xa, size_t pos) {           \
  assert(NULL != xa);                                                  \
  assert(pos < xa->count);                                             \
  return &(xa->storage[pos]);                                          \
};                                                                     \
                                                                       \
                                                                       \
/** Ensure that xarray @c xa can be extended by appending @c nmemb     \
    elements at the end of the array without incurring in any          \
    relocation; return a pointer to the first newly-added place, or @c \
    NULL on failure. Note that the newly-added positions are not       \
    initialized. */                                                    \
_inline void aname##_reserve(aname##_t** xa, const size_t nmemb) {     \
  assert(NULL != xa);                                                  \
  assert(NULL != *xa);                                                 \
  assert((*xa)->count <= (*xa)->allocated);                            \
                                                                       \
  const long more = nmemb - ((*xa)->allocated - (*xa)->count);         \
  if (more > 0) {                                                      \
    /* resize xarray->storage */                                       \
    aname##_t* new_xa = xrealloc((*xa), sizeof(aname##_t) +            \
                                 sizeof(elt_t) * (more + (*xa)->allocated)); \
    if (NULL == new_xa)                                                \
      return;                                                          \
    (*xa) = new_xa;                                                    \
    (*xa)->allocated += more;                                          \
  }                                                                    \
  assert((*xa)->count <= (*xa)->allocated);                            \
};                                                                     \
                                                                       \
                                                                       \
/** Extend xarray @c xa by adding @c nmemb elements at the end of the  \
    array; return a pointer to the first newly-added place, or @c      \
    NULL on failure. Note that the newly-added positions are not       \
    initialized. */                                                    \
_inline elt_t*  aname##_extend(aname##_t** xa, size_t nmemb) {         \
  assert(NULL != xa);                                                  \
  assert(NULL != *xa);                                                 \
  assert(0 < nmemb);                                                   \
                                                                       \
  aname##_reserve(xa, nmemb);                                          \
  (*xa)->count += nmemb;                                               \
                                                                       \
  assert((*xa)->count <= (*xa)->allocated);                            \
  return aname##_at(*xa, (*xa)->count - nmemb);                        \
};                                                                     \
                                                                       \
                                                                       \
/** Extend xarray @c xa by adding 1 element at the end of the array,   \
    and return a pointer to the newly-added place, or @c NULL on       \
    failure. */                                                        \
_inline elt_t*  aname##_extend1(aname##_t** xa) {                      \
  return aname##_extend(xa, 1);                                        \
};                                                                     \
                                                                       \
                                                                       \
/** Shorten @c xa by removing @a nmemb elements from the end of the array. */\
_inline void aname##_shorten(aname##_t** const xa, const size_t nmemb){\
  assert(NULL != xa);                                                  \
  assert(NULL != *xa);                                                 \
                                                                       \
  (*xa)->count -= nmemb;                                               \
  if((*xa)->count < 0)                                                 \
    (*xa)->count = 0;                                                  \
};                                                                     \
                                                                       \
                                                                       \
/** Shorten @c xa by removing one element from the end of the array. */\
_inline void aname##_shorten1(aname##_t** const xa) {                  \
  assert(NULL != xa);                                                  \
  assert(NULL != *xa);                                                 \
                                                                       \
  if((*xa)->count > 0)                                                 \
    (*xa)->count--;                                                    \
};                                                                     \
                                                                       \
                                                                       \
/** Insert one element in the array at position @c pos and return      \
 pointer to the newly-added place, or @c NULL if some error occurred.  \
 Array positions follow the usual C convention, ranging from 0 to @c   \
 xa->count-1. Note that the new location is not initialized. */        \
_inline elt_t* aname##_insert(aname##_t** const xa, const size_t pos) {\
  assert(NULL != xa);                                                  \
  assert(NULL != *xa);                                                 \
  assert(pos < (*xa)->count);                                          \
                                                                       \
  elt_t* p = aname##_extend1(xa);                                      \
  if(NULL == p)                                                        \
    return NULL;                                                       \
                                                                       \
  elt_t* src = aname##_at(*xa, pos);                                   \
  elt_t* dst = aname##_at(*xa, pos+1);                                 \
  memmove(dst, src, sizeof(elt_t) * ((*xa)->count - pos));             \
  return src;                                                          \
};                                                                     \
                                                                       \
                                                                       \
/** Remove the element at position @c pos from the xarray. Array       \
 positions follow the usual C convention, ranging from 0 to @c         \
 xa->count-1. */                                                       \
_inline void aname##_erase(aname##_t** const xa, const size_t pos) {   \
  assert(NULL != xa);                                                  \
  assert(NULL != *xa);                                                 \
  assert(pos < (*xa)->count);                                          \
                                                                       \
  if (pos == (*xa)->count - 1) {                                       \
    /* remove last element */                                          \
    aname##_shorten1(xa);                                              \
  }                                                                    \
  else {                                                               \
    void* src = aname##_at(*xa, pos+1);                                \
    void* dst = aname##_at(*xa, pos);                                  \
    memmove(dst, src, sizeof(elt_t) * ((*xa)->count - pos - 1));       \
    (*xa)->count -= 1;                                                 \
  }                                                                    \
};                                                                     \
                                                                       \
                                                                       \
/** Forget all the contents and return @c xa to 0 size. The memory     \
    allocated to the array is *not* freed. */                          \
_inline void aname##_clear(aname##_t** const xa) {                     \
  assert(NULL != xa);                                                  \
  assert(NULL != *xa);                                                 \
  (*xa)->count = 0;                                                    \
};                                                                     \
                                                                       \
                                                                       \
/** Free array @c xa and return allocated memory to the system. */     \
_inline void aname##_free(aname##_t *xa) {                             \
  free(xa);                                                            \
};                                                                     \


#endif // XARRAY_H

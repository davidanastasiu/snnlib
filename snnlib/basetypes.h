/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#pragma once

#include <cstdint>
#include <sys/types.h>

#include "BitVector.h"

namespace snnlib {

/* Index type should be a signed int type */
#ifdef IDXTYPE_LONG
#define PRNT_IDXTYPE "%ld"
#define IDX_MAX INT64_MAX
#define IDX_MIN INT64_MIN
typedef int64_t idx_t;
#elif defined IDXTYPE_SHORT
#define IDX_MAX SHRT_MAX
#define IDX_MIN SHRT_MIN
#define PRNT_IDXTYPE "%hi"
typedef short idx_t;
#else
#define IDXTYPE_INT
#define IDX_MAX INT32_MAX
#define IDX_MIN INT32_MIN
#define PRNT_IDXTYPE "%d"
typedef int32_t idx_t;
#endif

/* Pointer type should be a signed int type */
#ifdef PTRTYPE_INT
#define PRNT_PTRTYPE "%d"
#define PTR_MAX INT32_MAX
#define PTR_MIN INT32_MIN
typedef int32_t ptr_t;
#elif defined PTRTYPE_LONG
#define PRNT_PTRTYPE "%ld"
#define PTR_MAX INT64_MAX
#define PTR_MIN INT64_MIN
typedef int64_t ptr_t;
#else
#define PTRTYPE_SSZIE
#define PTR_MAX INT64_MAX
#define PTR_MIN INT64_MIN
#define PRNT_PTRTYPE "%zu"
typedef ssize_t ptr_t;
#endif


/* Value type should be a float or double */
#ifdef VALTYPE_DOUBLE
typedef double val_t;
#define VAL_MAX DBL_MAX
#define VAL_MIN DBL_MIN
#define PRNT_VALTYPE "%g"
#else
#define VALTYPE_FLOAT
#define VAL_MAX FLT_MAX
#define VAL_MIN FLT_MIN
#define PRNT_VALTYPE "%f"
typedef float val_t;
#endif

/* Accumulator type should be a float or double */
#ifdef ACMTYPE_FLOAT
#define PRNT_ACMTYPE "%f"
typedef float acm_t;
#else
#define ACMTYPE_DOUBLE
#define PRNT_ACMTYPE "%g"
typedef double acm_t;
#endif

#define MAXIDX	(1<<8*sizeof(idx_t)-2)

/* Key-Value Pairs */
struct iakv_t {
    idx_t key;
    acm_t val;

    iakv_t& operator=(const iakv_t& a)
    {
        key=a.key;
        val=a.val;
        return *this;
    }

    bool operator==(const iakv_t& rhs) const
    {
        return val == rhs.val;
    }

    bool operator<(const iakv_t& rhs) const
    {
        return val < rhs.val;
    }

    bool operator>(const iakv_t& rhs) const
    {
        return val > rhs.val;
    }

    bool operator<=(const iakv_t& rhs) const
    {
        return val <= rhs.val;
    }

    bool operator>=(const iakv_t& rhs) const
    {
        return val >= rhs.val;
    }
};

struct ivkv_t {
    idx_t key;
    val_t val;

    ivkv_t& operator=(const ivkv_t& a)
    {
        key=a.key;
        val=a.val;
        return *this;
    }

    bool operator==(const ivkv_t& rhs) const
    {
        return val == rhs.val;
    }

    bool operator<(const ivkv_t& rhs) const
    {
        return val < rhs.val;
    }

    bool operator>(const ivkv_t& rhs) const
    {
        return val > rhs.val;
    }

    bool operator<=(const ivkv_t& rhs) const
    {
        return val <= rhs.val;
    }

    bool operator>=(const ivkv_t& rhs) const
    {
        return val >= rhs.val;
    }
};

struct pikv_t {
    ptr_t key;
    idx_t val;

    // assignment operator modifies object, therefore non-const
    pikv_t& operator=(const pikv_t& a)
    {
        key=a.key;
        val=a.val;
        return *this;
    }

    bool operator==(const pikv_t& rhs) const
    {
        return val == rhs.val;
    }

    bool operator<(const pikv_t& rhs) const
    {
        return val < rhs.val;
    }

    bool operator>(const pikv_t& rhs) const
    {
        return val > rhs.val;
    }

    bool operator<=(const pikv_t& rhs) const
    {
        return val <= rhs.val;
    }

    bool operator>=(const pikv_t& rhs) const
    {
        return val >= rhs.val;
    }
};


/**
 * Read-only queue of pre-defined list of values
*/
template<typename T>
struct rqueue {
    T* data;
    size_t f;     // front of list
    ssize_t size; // list size
    
    rqueue(T* data, const size_t size):
        data(data), f(0), size(size){}

    ~rqueue(){}

    void reset(){
        f = 0;
    }

    void reset(T* data, const size_t size){
        this->data = data;
        this->size = size;
        this->f = 0;
    }

    void push(const T item){
        throw std::runtime_error("Read-only queue. Cannot push.");
    }

    T pop(){
        if(f == size){
            throw std::out_of_range("No other items left to pop.");
        }
        return data[f++];
    }
};

} /** end namespace */
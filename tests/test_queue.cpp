/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#include <iostream>     // std::cout
#include <chrono>
#include <gtest/gtest.h>
#include <math.h>

//header file
#include "../snnlib/basetypes.h"

using namespace snnlib;

TEST(rqueue, consume){
    size_t n=10;
    auto data = new idx_t[n];
    for(idx_t i=0; i < n; ++i){
        data[i] = i;
    }
    auto q = rqueue<idx_t>(data, n);
    EXPECT_ANY_THROW(q.push(1));
    for(idx_t i=0; i < n; ++i){
        EXPECT_EQ(q.pop(), i);
    }
    EXPECT_ANY_THROW(q.pop());
    q.reset();
    for(idx_t i=0; i < n; ++i){
        EXPECT_EQ(q.pop(), i);
    }
    delete data;
}

TEST(rqueue, reset){
    size_t n=10;
    auto data = new idx_t[n];
    auto data2 = new idx_t[n];
    for(idx_t i=0; i < n; ++i){
        data[i] = i;
        data2[i] = i * 2;
    }
    auto q = rqueue<idx_t>(data, n);
    EXPECT_ANY_THROW(q.push(1));
    for(idx_t i=0; i < n; ++i){
        EXPECT_EQ(q.pop(), i);
    }
    EXPECT_ANY_THROW(q.pop());
    q.reset(data2, n);
    for(idx_t i=0; i < n; ++i){
        EXPECT_EQ(q.pop(), i * 2);
    }
    delete data;
    delete data2;
}
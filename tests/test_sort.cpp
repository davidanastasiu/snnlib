/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#include <iostream>     // std::cout
#include <limits>
#include <gtest/gtest.h>

//header file
#include "../snnlib/basetypes.h"
#include "../snnlib/sort.h"

using namespace snnlib;


/**
 * Test sorting
 */
TEST(sort, sort){
    auto keys = new idx_t[100];
    auto vals = new val_t[100];
    for(int i = 0; i < 100; ++i){
        keys[i] = std::rand();
        vals[i] = (double) std::rand()/RAND_MAX;
    }
    sorti(100, keys);
    sorti(100, vals);
    for(int i = 1; i < 100; ++i){
        EXPECT_GE(keys[i], keys[i-1]);
        EXPECT_GE(vals[i], vals[i-1]);
    }
    sortd(100, keys);
    sortd(100, vals);
    for(int i = 1; i < 100; ++i){
        EXPECT_LE(keys[i], keys[i-1]);
        EXPECT_LE(vals[i], vals[i-1]);
    }
    delete keys;
    delete vals;
}

TEST(sort, kvsort){
    auto keys = new idx_t[100];
    auto vals = new val_t[100];
    for(int i = 0; i < 100; ++i){
        keys[i] = std::rand();
        vals[i] = (double) std::rand()/RAND_MAX;
    }
    kvsorti(100, keys, vals);
    for(int i = 1; i < 100; ++i){
        EXPECT_GE(vals[i], vals[i-1]);
    }
    kvsortd(100, keys, vals);
    for(int i = 1; i < 100; ++i){
        EXPECT_LE(vals[i], vals[i-1]);
    }
    delete keys;
    delete vals;
}


/**
 * Test select
 */
TEST(sort, select){
    auto vals = new val_t[100];
    val_t v;
    for(int i = 0; i < 100; ++i){
        vals[i] = (double) std::rand()/RAND_MAX;
    }
    selectd(100, vals, 10);
    v = -100000000.0;
    for(int i = 0; i < 10; ++i){
        if(vals[i] > v){
            v = vals[i];
        }
    }
    for(int i = 10; i < 100; ++i){
        EXPECT_LE(vals[i], v);
    }
    selecti(100, vals, 10);
    for(int i = 0; i < 10; ++i){
        if(vals[i] < v){
            v = vals[i];
        }
    }
    for(int i = 10; i < 100; ++i){
        EXPECT_GE(vals[i], v);
    }
    delete vals;
}


/**
 * Test select
 */
TEST(sort, selectsort){
    auto vals = new val_t[100];
    val_t v;
    for(int i = 0; i < 100; ++i){
        vals[i] = (double) std::rand()/RAND_MAX;
    }
    selectsortd(100, vals, 10);
    v = -100000000.0;
    for(int i = 0; i < 10; ++i){
        if(vals[i] > v){
            v = vals[i];
        }
        if(i > 0){
            EXPECT_LE(vals[i], vals[i-1]);
        }
    }
    for(int i = 10; i < 100; ++i){
        EXPECT_LE(vals[i], v);
    }
    selectsorti(100, vals, 10);
    for(int i = 0; i < 10; ++i){
        if(vals[i] < v){
            v = vals[i];
        }
        if(i > 0){
            EXPECT_GE(vals[i], vals[i-1]);
        }
    }
    for(int i = 10; i < 100; ++i){
        EXPECT_GE(vals[i], v);
    }
    delete vals;
}

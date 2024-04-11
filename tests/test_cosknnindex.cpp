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

//header files
#include "../snnlib/CosKnnIndex.h"

using namespace snnlib;

TEST(cosknnindex, knnga){
    auto mat = csr_t::random(100, 837, 0.25);
    auto index = new CosKnnIndex(mat);
    auto graph = index->knng(10, true);
    EXPECT_EQ(mat->nrows, graph->nrows);
    EXPECT_EQ(mat->nrows, graph->ncols);
    for(idx_t i=0; i < graph->nrows; ++i){
        EXPECT_LE(graph->rptr[i+1]-graph->rptr[i], 10);
    }
    delete mat;
    delete index;
    delete graph;
}

TEST(cosknnindex, knngbf){
    auto mat = csr_t::random(100, 837, 0.25);
    auto index = new CosKnnIndex(mat);
    auto graph = index->knngbf(10);
    EXPECT_EQ(mat->nrows, graph->nrows);
    EXPECT_EQ(mat->nrows, graph->ncols);
    for(idx_t i=0; i < graph->nrows; ++i){
        EXPECT_LE(graph->rptr[i+1]-graph->rptr[i], 10);
    }
    delete mat;
    delete index;
    delete graph;
}

TEST(cosknnindex, knnga_pruning){
    auto mat = csr_t::random(10000, 1237, 0.05);
    auto index = new CosKnnIndex(mat);
    auto g1 = index->knng(10, true, 2.0, 1, 0.0, true);
    auto g2 = index->knng(10, true, 20.0, 1, 0.5, true);
    // TODO: Should test recall here, but not yet implemented in C++
    EXPECT_EQ(g1->rptr[g1->nrows], g2->rptr[g2->nrows]);
    delete mat;
    delete index;
    delete g1;
    delete g2;
}
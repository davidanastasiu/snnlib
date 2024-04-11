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
#include "../snnlib/pqueue.h"

using namespace snnlib;


/**
 * Test creation 
 */
TEST(Pqueue, test_create_internal){
    auto pq = new pqueue(10);
    EXPECT_EQ(pq->nnodes, 0);
    EXPECT_EQ(pq->maxnodes, 10);
    for(int i = 0; i < 10; i++){
        EXPECT_EQ(pq->heap[i].key, -1);
        EXPECT_EQ(pq->heap[i].val, 0.0);
    }
    delete pq;
}

/**
 * Test insertion 
 */
TEST(Pqueue, test_insert){
    auto pq = new pqueue(10);
    // insert till full
    int i;
    for(i = 9; i >= 0; i--){
        pq->insert(i, i * 1.1);
    }

    EXPECT_EQ(pq->see_top_val(), 0.0);
    EXPECT_EQ(pq->nnodes, 10);
    EXPECT_TRUE(pq->check());


    auto pq1 = new pqueue(10);
    // insert till full
    for(i = 0; i < 10; i++){
        pq1->insert(i, i * 1.1);
    }

    EXPECT_EQ(pq1->see_top_key(), 0);
    EXPECT_EQ(pq1->see_top_val(), 0.0);
    EXPECT_EQ(pq1->nnodes, 10);
    EXPECT_TRUE(pq1->check());

    // should not insert -- only insert values > min value if queue is full
    pq1->insert(11, 0); // equal with min, so no insert
    EXPECT_EQ(pq1->see_top_key(), 0);
    EXPECT_EQ(pq1->see_top_val(), 0.0);
    EXPECT_TRUE(pq1->check());
    EXPECT_EQ(pq1->nnodes, 10);
    for(i = 11; i < 100; i++){
        pq1->insert(i, i * -1.1);
        EXPECT_EQ(pq1->nnodes, 10);
        EXPECT_EQ(pq1->see_top_key(), 0);
        EXPECT_EQ(pq1->see_top_val(), 0.0);
    }
    EXPECT_TRUE(pq1->check());

    // insert when full -- larger values than min should insert
    pq1->insert(10, 10.5);
    EXPECT_TRUE(pq1->check());
    EXPECT_EQ(pq1->see_top_key(), 1);
    EXPECT_FLOAT_EQ(pq1->see_top_val(), 1.1);
    EXPECT_EQ(pq1->nnodes, 10);
    for(i = 11; i < 100; i++){
        pq1->insert(i, i * 1.1);
        EXPECT_EQ(pq1->nnodes, 10);
    }
    EXPECT_EQ(pq1->see_top_key(), 90);
    EXPECT_FLOAT_EQ(pq1->see_top_val(), 99.0);
    EXPECT_TRUE(pq1->check());

    delete pq;
    delete pq1;
}


/**
 * Test remove
 */
TEST(Pqueue, test_remove_heap){
    auto pq1 = new pqueue(10);
    int i;
    // insert untill full
    for(i = 0; i < 10; i++){
        pq1->insert(i, i * 1.1);
    }

    EXPECT_EQ(pq1->see_top_val(), 0.0);
    EXPECT_EQ(pq1->nnodes, 10);
    EXPECT_TRUE(pq1->check());

    for(int i = 0; i < 10; i++){
        auto top = pq1->see_top_val();
        auto topk = pq1->see_top_key();
        ivkv_t pair = pq1->get_top();
        EXPECT_EQ(top, pair.val);
        EXPECT_EQ(topk, pair.key);
        EXPECT_EQ(pq1->nnodes, 9 - i);
        EXPECT_TRUE(pq1->check());
    }
    delete pq1;
}

/**
 * Test external memory queue
 */
TEST(Pqueue, test_external){
    ivkv_t * heap = (ivkv_t *) malloc(10 * sizeof(ivkv_t));
    auto pq = new pqueue(10, 0, heap);
    EXPECT_EQ(pq->nnodes, 0);
    EXPECT_EQ(pq->maxnodes, 10);
    int i;
    for(i = 0; i < 10; i++){
        pq->insert(i, i * 1.1);
    }

    EXPECT_EQ(pq->see_top_val(), 0.0);
    EXPECT_EQ(pq->nnodes, 10);
    EXPECT_TRUE(pq->check());

    for(i = 19; i >= 10; i--){
        pq->insert(i, i * 1.1);
    }

    EXPECT_EQ(pq->see_top_val(), 11.0);
    EXPECT_EQ(pq->nnodes, 10);
    EXPECT_TRUE(pq->check());

    for(int i = 0; i < 10; i++){
        auto top = pq->see_top_val();
        auto topk = pq->see_top_key();
        ivkv_t pair = pq->get_top();
        EXPECT_EQ(top, pair.val);
        EXPECT_EQ(topk, pair.key);
        EXPECT_EQ(pq->nnodes, 9 - i);
        EXPECT_TRUE(pq->check());
    }

    delete pq;
    free(heap);
}



/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#include <iostream>     // std::cout
#define MAX64 18446744073709551615ULL
#include <limits.h>
#include <chrono>
#include <gtest/gtest.h>

//header file
#include "../snnlib/BitVector.h"

namespace snnlib {

    /**
     * Test if vector is allocated with all 0s
     * @param  BitVector Vector
    */
    TEST(BitVector, zero_vector){
        uint64_t size = 10000;
        BitVector bv = BitVector(size);
        EXPECT_EQ(bv.allocation_size(size), ((size + 63) / 64));
        EXPECT_EQ(bv.isset(0), 0);
        EXPECT_EQ(bv.isset(size-1), 0);
        EXPECT_EQ(bv[0], 0);
    }

    /**
     * Test various setting functions of vector
     * @param  BitVector Vector
     * @return true if expected values are set to 1
     * @return true if expected values are set to 0
    */
    TEST(BitVector, set_vector){
        uint64_t size = 10000;
        BitVector bv = BitVector(size);
        //set to 1
        for(int i = 1; i < size; i <<= 2){
            bv.set(i);
            EXPECT_EQ(bv.isset(i), 1);
        }
        //set some bits to 0
        for(int i = 1; i < size; i <<=4){
            bv.zero(i);
            EXPECT_EQ(bv.isset(i), 0);
        }
        //set all bits to 0
        bv.reset();
        EXPECT_EQ(bv.isset(0), 0);
        EXPECT_EQ(bv.isset(size-1), 0);
        EXPECT_EQ(bv.isset(1024), 0);

    }

    /**
     * Test vector index operator overload
     * @param  BitVector Vector
     * @return true if index operator correctly reads and writes values
    */
    TEST(BitVector, index_operator){
        uint64_t size = 128; //sz = 2
        BitVector bv = BitVector(size);
        bv.set(63);
        bv.set(62);
        bv.set(0);
        //Expect index operator to return the value of first index
        EXPECT_EQ(bv[0], 3+((uint64_t)1<<63));

        //change value of first index
        bv[0] = 0;
        EXPECT_EQ(bv[0], 0);
        EXPECT_EQ(bv.isset(1), 0);

    }

    /**
     * Test vector functionality with out of range inout
     * @param  BitVector Vector
    */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshift-count-overflow"
    TEST(BitVector, out_of_range_input){
        BitVector bv = BitVector(128);
        bv[0] = (uint64_t)1<<65;
        //truncated
        EXPECT_EQ(bv[0], 0);
    }
#pragma GCC diagnostic pop

    /**
     * Test vector functionality with negative value. Also tests edge with max value
     * @param  BitVector Vector
     * @return true if values are properly changed from negative int to uint64
    */
    TEST(BitVector, negative_input){
        BitVector bv = BitVector(128);
        //setting to -1 should cause bv[0] to be max int of 64 bit int
        bv[0] = -1;
        EXPECT_EQ(bv.isset(0), 1);
        EXPECT_EQ(bv.isset(63), 1);
        EXPECT_EQ(bv.isset(1), 1);


        //MAX64 == 18446744073709551615
        //max 64bit uint
        EXPECT_EQ(bv[0], ULONG_MAX);

        bv[1] = -2;
        EXPECT_EQ(bv[1], MAX64-1);
        EXPECT_EQ(bv[0], MAX64-1);

        bv[64] = -10;
        EXPECT_EQ(bv[64], MAX64-9);
        EXPECT_EQ(bv[0], MAX64-1);
    }

    /**
     * Test constructors
     */
    TEST(BitVector, true_constructor){
        BitVector * bv = new BitVector(128, true);
        for(int i = 0; i < 2; i++){
            EXPECT_EQ((*bv)[i], ULONG_MAX);
            for(int j = 0; j < 64; j++){
                EXPECT_TRUE(bv->isset(i * 64 + j));
            }
        }
        delete bv;
    }

    TEST(BitVector, false_constructor){
        BitVector* bv = new BitVector(128, false);
        for(int i = 0; i < 2; i++){
            EXPECT_EQ((*bv)[i], 0);
            for(int j = 0; j < 64; j++)
                EXPECT_FALSE(bv->isset(i * 64 + j));
        }
        delete bv;
    }


    /**
     * Test toggle
     */
    TEST(BitVector, toggle_false){
        BitVector* bv = new BitVector(128, true);
        for(int i = 0; i < 2; i++){
            EXPECT_EQ(bv->operator[](i*64), ULONG_MAX);
            for(int j = 0; j < 64; j++){
                EXPECT_TRUE(bv->isset(i * 64 + j));
                bv->toggle(i * 64 + j);
            }
        }
        for(int i = 0; i < 2; i++){
            EXPECT_EQ(bv->operator[](i*64), 0);
            for(int j = 0; j < 64; j++){
                EXPECT_FALSE(bv->isset(i * 64 + j));
            }
        }
        delete bv;
    }

    TEST(BitVector, toggle_true){
        BitVector* bv = new BitVector(128, false);
        for(int i = 0; i < 2; i++){
            EXPECT_EQ(bv->operator[](i*64), 0);
            for(int j = 0; j < 64; j++){
                EXPECT_FALSE(bv->isset(i * 64 + j));
                bv->toggle(i * 64 + j);
            }
        }
        for(int i = 0; i < 2; i++){
            EXPECT_EQ(bv->operator[](i*64), ULONG_MAX);
            for(int j = 0; j < 64; j++){
                EXPECT_TRUE(bv->isset(i * 64 + j));
            }
        }
        delete bv;
    }

    TEST(BitVector, toggle_half){
        BitVector* bv = new BitVector(128, true);
        EXPECT_EQ(bv->count(), 128);
        for(int j = 0; j < 64; j++){
            EXPECT_TRUE(bv->isset(j));
            bv->toggle(j);
            EXPECT_FALSE(bv->isset(j));
        }
        EXPECT_EQ(bv->count(), 64);
        EXPECT_EQ(bv->operator[](0), 0);
        EXPECT_EQ(bv->operator[](64), ULONG_MAX);
        delete bv;
    }

} /* end namespace */

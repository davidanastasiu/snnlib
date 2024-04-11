/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#include <iostream>     // std::cout
#include <gtest/gtest.h>
#include <chrono>
#include <omp.h>

#include "../snnlib/KVVectorView.h"

using namespace snnlib;

int testvalarr[15] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
int testkeyarr[15] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
TEST(KVV_View, creation_with_keyval_arrays){
    KVVectorView<int, int> ivkv = KVVectorView<int, int>(15,testkeyarr,testvalarr);
    EXPECT_EQ(ivkv.size(), 15);
    for(int i = 0; i < 15; i++){
        EXPECT_EQ(ivkv.key(i), i+1);
        EXPECT_EQ(ivkv.val(i), i+1);
    }
}

TEST(KVV_View, creation_with_val_array){
    KVVectorView<int, int> ivkv = KVVectorView<int, int>(15,testvalarr);
    EXPECT_EQ(ivkv.size(), 15);
    for(int i = 0; i < 15; i++){
        EXPECT_EQ(ivkv.key(i), i);
        EXPECT_EQ(ivkv.val(i), i+1);
    }
}

/**
 * Should throw for all functions that increase the size/ reallocate
 * insert
 * pushback
 * assign
 * shrink
 * reserve
 * emplace
 */
TEST(KVV_View, throw_tests){
    KVVectorView<int, int> ivkv = KVVectorView<int, int>(15,testkeyarr,testvalarr);
    KVVectorView<int, int> ivkv2 = KVVectorView<int, int>(15,testkeyarr,testvalarr);
    int * begin = ivkv.k_begin();
    int * end = ivkv.k_end();
    auto list = {std::make_pair(5, 10), std::make_pair(1, 12), std::make_pair(3, 7)};
    auto vals = {1,23,5,3,2,1};
    //insertion
    EXPECT_ANY_THROW(ivkv.insert(ivkv.k_begin(), 10));
    EXPECT_ANY_THROW(ivkv.insert(ivkv.k_end(),std::make_pair(1,2)));
    EXPECT_ANY_THROW(ivkv.insert(ivkv.k_end(), 5, 10));
    EXPECT_ANY_THROW(ivkv.insert(ivkv.k_end(),5,std::make_pair(1,2)));
    EXPECT_ANY_THROW(ivkv.insert(begin, list));
    EXPECT_ANY_THROW(ivkv.insert(end, vals));
    EXPECT_ANY_THROW(ivkv.insert(begin + 5, ivkv2.k_begin(), ivkv2.k_end(), ivkv2));

    //push back/emplace
    EXPECT_ANY_THROW(ivkv.push_back(std::make_pair(5,10)));
    EXPECT_ANY_THROW(ivkv.push_back(5));
    EXPECT_ANY_THROW(ivkv.push_back(5,3));
    EXPECT_ANY_THROW(ivkv.emplace(begin, 5));
    EXPECT_ANY_THROW(ivkv.emplace_back(10,2));

    //resizing
    EXPECT_ANY_THROW(ivkv.shrink_to_fit());
    EXPECT_ANY_THROW(ivkv.reserve(0));
    EXPECT_ANY_THROW(ivkv.resize(1));
    EXPECT_ANY_THROW(ivkv.resize(120, (int)5));
    EXPECT_ANY_THROW(ivkv.resize(120, (int)2, (int)4));
    EXPECT_ANY_THROW(ivkv.resize(10, std::make_pair(12,2)));
    EXPECT_ANY_THROW(ivkv.swap(ivkv2));

    //assign
    EXPECT_ANY_THROW(ivkv.assign(list));
    EXPECT_ANY_THROW(ivkv.assign(vals));
    EXPECT_ANY_THROW(ivkv.assign(5,1,2));
    EXPECT_ANY_THROW(ivkv.assign(5,std::make_pair(1,2)));
    EXPECT_ANY_THROW(ivkv.assign(5,testvalarr));
}
/**
 * Tests for functions that don't throw errors
 */
TEST(KVVectorView, element_access){
    auto ivkv = KVVectorView<int, int>(15,testkeyarr, testvalarr);
    KVVector<int, int>::kiterator itb = ivkv.k_begin();
    KVVector<int, int>::kiterator ite = ivkv.k_end();
    KVVector<int, int>::reverse_kiterator itrb = ivkv.k_rbegin();
    KVVector<int, int>::reverse_kiterator itre = ivkv.k_rend();

    EXPECT_EQ(&*itb, &*(itre-1));
    EXPECT_EQ(&*(ite-1), &*itrb);

    EXPECT_EQ(&*itb, &ivkv.k_data()[0]);
    EXPECT_EQ(&*(ite-1), &ivkv.k_data()[ivkv.size()-1]);

    KVVector<int, int>::viterator itbv = ivkv.v_begin();
    KVVector<int, int>::viterator itev = ivkv.v_end();
    KVVector<int, int>::reverse_viterator itrbv = ivkv.v_rbegin();
    KVVector<int, int>::reverse_viterator itrev = ivkv.v_rend();

    EXPECT_EQ(&*itbv, &*(itrev-1));
    EXPECT_EQ(&*(itev-1), &*itrbv);

    EXPECT_EQ(&*itbv, &ivkv.v_data()[0]);
    EXPECT_EQ(&*(itev-1), &ivkv.v_data()[ivkv.size()-1]);
}

/**
 * Tests erase function by checking
 * Size should be set to 0
 */
TEST(KVVectorView, erase){
    auto ivkv = KVVectorView<int, int>(15, testkeyarr, testvalarr);
    EXPECT_EQ(ivkv.size(), 15);
    ivkv.erase(ivkv.k_begin());
    EXPECT_EQ(ivkv.size(), 14);
    EXPECT_EQ(ivkv.key(0), 2);
    ivkv.erase(ivkv.k_begin(), ivkv.k_begin()+5);
    EXPECT_EQ(ivkv.key(0), 7);
    EXPECT_EQ(ivkv.size(), 9);
    ivkv.erase(ivkv.k_begin(), ivkv.k_end());
    EXPECT_EQ(ivkv.size(), 0);
}

TEST(KVVectorView, scan){
    int arrk[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    int arrv[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    auto ivkv = KVVectorView<int, int>(15, arrk, arrv);
    EXPECT_NO_THROW(ivkv.scan('>'));
    for(std::size_t i = 0; i < ivkv.size(); i++) {
        EXPECT_EQ(ivkv.key(i), i+1);
        EXPECT_EQ(ivkv.val(i), i+1);
    }
    //should change base array
    EXPECT_NO_THROW(ivkv.scan('+'));
    EXPECT_EQ(arrv[1], 3);
    EXPECT_EQ(arrv[2], 6);
}

TEST(KVVectorView, merge){
    int arrk1[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    int arrv1[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    auto test = KVVector<int, int>(15,arrk1,arrv1);
    test.merge_keys('+');
    for(int i = 0; i < 15; i++){
        EXPECT_EQ(test.key(i), i + 1);
    }
    EXPECT_EQ(test.size(),15);

    int arrk[] = {1,1,1,1,1,1,1,1};
    int arrv[] = {1,1,1,1,1,1,1,1};
    auto ivkv2 = KVVectorView<int,int>(8,arrk,arrv);
    ivkv2.merge_keys('+');
    EXPECT_EQ(ivkv2.key(0),1);
    EXPECT_EQ(ivkv2.val(0),8);
}

TEST(KVVectorView, sort_partial) {
    int arrk1[100];
    double arrv1[100];
    for (int i = 0; i < 100; i++) {
    	arrk1[i] = 100 - i;
    	arrv1[i] = (100 - i) * 1.0;
    }
    KVVectorView<int,double> * ivkv = new KVVectorView<int, double>(100, arrk1, arrv1);
    omp_set_num_threads(5);
    #pragma omp parallel num_threads(5)
    {
        #pragma omp for schedule(static)
        for(int i = 0; i < 5; i++) {
        	ivkv->sort_kvi(20, i * 20); // sort in blocks of 20
        }
    }
    for(int i = 0; i < 5; i++) {
    	for(int j = i * 20; j < (i + 1) * 20 - 1; j++){
    	   EXPECT_TRUE(ivkv->key(j+1) >= ivkv->key(j));
        }
    }
    delete ivkv;
}

/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#include <iostream>     // std::cout
#include <functional>   // std::greater
#include <algorithm>    // std::sort
#include <chrono>
#include <gtest/gtest.h>
#include <memory>


#include "../snnlib/KVVector.h"


using namespace snnlib;

template<typename K, typename V>
void printdata(KVVector<K, V> vec) {
    for (int i = 0; i < vec.size(); i++) {
        std::cout << "(" << vec.key(i) << ", " << vec.at(i) << ") ";
    }
    std::cout << "\n";
}

/**
 * Tests to do:
 * Constructor
 * Comparison Operators
 * Sorting
 * Data
 * Modifiers
 * Element Access
 * Capacity
 * Iterators
 */

/**
 * Constructor tests
 */

/**
 * Test KVVector creation with no input
 * @param void
 * @return true if size and pushed values are correct
*/
TEST(KVVector, creation) {
    // number of key-value pairs to add to the vector
    unsigned int n = 1024;

    // make a KV vector of <int, float>
    KVVector<int, float> ivkv = KVVector<int, float>();
    for(unsigned int i=0; i < n; ++i){
        ivkv.push_back(i, (float)i*2);
    }
    EXPECT_EQ(ivkv.size(), n);
    EXPECT_EQ(ivkv.key(0), 0);
    EXPECT_EQ(ivkv.at(0), 0);
    EXPECT_EQ(ivkv.key(1023), 1023);
    EXPECT_EQ(ivkv.at(1023), 2046);

}

/**
 * Test KVVector creation with size input
 * @param size of Vector
 * @return true if size is correct
 * @return true if values are init as 0
*/
TEST(KVVector, creation_with_size){
    //size
    std::size_t n = 1024;
    //KV vector of <int, int>
    KVVector<int, int> ivkv =  KVVector<int, int>(n);
    EXPECT_EQ(ivkv.capacity(), n);
    EXPECT_EQ(ivkv.size(), 0);

}

/**
 * Test if Vector constructor with size and pairs
 * @param n        size of Vector
 * @param p        Key Value pair
 * @return true    if size is correct
 * @return true    if pairs are same as parameter
*/
TEST(KVVector, creation_with_size_pair){
    //size
    std::size_t n = 1024;
    std::pair<int, float> p(5, 10.0);
    KVVector<int, float> ivkv =  KVVector<int, float>(n, p);

    EXPECT_EQ(ivkv.key(0), 5);
    EXPECT_EQ(ivkv.key(1023), 5);
    EXPECT_NEAR(ivkv.at(0), 10.0, 1e-7);
    EXPECT_NEAR(ivkv.at(1023), 10.0, 1e-7);

}

/**
 * Test Vector construction with size and value
 * @param n        size of Vector
 * @param value    KVV pair value
 * @return true    if size is correct
 * @return true    if values are same as parameter
*/
TEST(KVVector, creation_with_size_value){
    //size
    std::size_t n = 1024;
    float value = 66.65f;
    //KV vector of <int, float>
    KVVector<int, float> ivkv = KVVector<int, float>(n, value);
    EXPECT_EQ(ivkv.size(), n);
    EXPECT_EQ(ivkv.key(0), 0);
    EXPECT_EQ(ivkv.key(1023), 1023);
    EXPECT_NEAR(ivkv.at(0), value, 1e-7);
    EXPECT_NEAR(ivkv.at(1023), value, 1e-7);
    EXPECT_NEAR(ivkv.at(512), value, 1e-7);
}

/**
 * Test Vector constructor with size, key, and value parameters
 * @param n        size of Vector
 * @param key      KVV pair key
 * @param value    KVV pair value
 * @return true    if size is correct
 * @return true    if pairs of values and keys are same as parameter
*/
TEST(KVVector, creation_with_size_key_value){
    std::size_t n = 1024;
    float value = 66.49f;
    int key = 78;
    //KV vector of <int, float>
    KVVector<int, float> ivkv = KVVector<int, float>(n, key, value);
    EXPECT_EQ(ivkv.size(), n);
    EXPECT_EQ(ivkv.key(0), key);
    EXPECT_EQ(ivkv.key(1023), key);
    EXPECT_EQ(ivkv.key(512), key);
    EXPECT_NEAR(ivkv.at(0), value, 1e-7);
    EXPECT_NEAR(ivkv.at(1023), value, 1e-7);
    EXPECT_NEAR(ivkv.at(512), value, 1e-7);
}

/**
 * Test if Vector is constructed with array of values and keys being read in
 * @param n        size of Vector
 * @param keys     array of KVV pair keys
 * @param values   array of KVV pair values
 * @return true    if size is correct
 * @return true    if pairs of values and keys  are same as input arrays
*/
TEST(KVVector, allocation_from_arrays) {
    int n = 5;
    int keys[] = {1, 2, 3, 4, 5};
    double vals[] = {5.0, 4.0, 3.0, 2.0, 1.0};

    // make a KV vector of <int, double>
    auto idkv = KVVector<int, double>(n, keys, vals);

    EXPECT_EQ(idkv.size(), n);
    EXPECT_EQ(idkv.key(0), 1);
    EXPECT_EQ(idkv.at(0), 5.0);
    EXPECT_EQ(idkv.key(n-1), 5);
    EXPECT_EQ(idkv.val(n-1), 1.0);
}

//check if array of values is read into Vector
/**
 * Test if Vector is constructed with array of values and keys that increment from 0
 * @param n        size of Vector
 * @param values   array of KVV pair values
 * @return true    if size is correct
 * @return true    if pairs of values are same as input arrays
 * @return true    if keys incremented correctly
*/
TEST(KVVector, allocation_from_value_array){
    int n = 5;
    double values[] = {1.1, 2.2, 3.3, 4.4, 5.5};
    KVVector<int, double> ivkv = KVVector<int, double>(n, values);
    EXPECT_EQ(ivkv.size(), n);
    EXPECT_EQ(ivkv.key(0), 0);
    EXPECT_EQ(ivkv.at(0), 1.1);
    EXPECT_EQ(ivkv.key(n-1), n-1);
    EXPECT_EQ(ivkv.val(n-1), 5.5);
}

/**
 * Test if vector1 pairs are read into vector2 using iterators
 * @param first    iterator pointing to beginning of Vector 1
 * @param last     iterator pointing to end of Vector 1
 * @return true    if pairs of values and keys are same in Vector1 and Vector2
*/
TEST(KVVector, creation_with_iterators){
    int n = 1024;
    float value= 3.0;
    KVVector<int, float> ivkv1 = KVVector<int, float>(n, value);

    KVVector<int, float>::kiterator first = ivkv1.k_begin();
    KVVector<int, float>::kiterator last = ivkv1.k_end();
    KVVector<int, float> ivkv2 = KVVector<int, float>(first, last, ivkv1);

    EXPECT_EQ(ivkv1.key(0), ivkv2.key(0));
    EXPECT_EQ(ivkv1.key(1023), ivkv2.key(1023));
    EXPECT_NEAR(ivkv1.at(0), ivkv2.at(0) , 1e-7);
    EXPECT_NEAR(ivkv1.val(1023), ivkv2.at(1023), 1e-7);
}

/**
 * Test if KVV constructor reads from list of pairs
 * @param il       init list<int, double>
 * @return true    if pairs of values and keys are same in Vector and init list
*/
TEST(KVVector, creation_with_initializerPair_list){
    auto il = {std::make_pair(5, 10.0), std::make_pair(1, 12.0), std::make_pair(3, 7.0)};
    auto ivkv = KVVector<int, double>(il);
    EXPECT_EQ(ivkv.size(), il.size());
    EXPECT_EQ(ivkv.key(0), 5);
    EXPECT_DOUBLE_EQ(ivkv.at(0), 10.0);
    EXPECT_EQ(ivkv.key(2) , 3);
    EXPECT_DOUBLE_EQ(ivkv.val(2), 7.0);

}

/**
 * Test if KVV constructor reads from list of values and if keys are incremented by 1
 * @param il       init list<int>
 * @return true    if pairs of values are same in Vector and init list of values.
 * @return true    if keys increment by 1 from 0
*/
TEST(KVVector, creation_with_initializerValue_list){
    auto il = {1, 2, 3, 4, 5, 6, 7};
    auto ivkv = KVVector<int, int>(il);
    EXPECT_EQ(ivkv.size(), il.size());
    for(unsigned int i = 0; i< il.size(); i++){
        EXPECT_EQ(ivkv.val(i), i+1);
        EXPECT_EQ(ivkv.key(i), i);
    }
    std::initializer_list<int> il2 = {};
    auto ivkv2 = KVVector<int, int>(il2);
    EXPECT_EQ(il2.size(), ivkv2.size());

}

/**
 * Test if KVV constructor reads from list of values and if keys are incremented by 1
 * @param size     size of vector 1
 * @param ivkv     First vector
 * @return true    if pairs of values and keys and size are same in Vector1 and Vector2
*/
TEST(KVVector, creation_with_other_KVV){
    std::size_t size = 10;
    float value = 5.0;
    auto ivkv = KVVector<int, float>(size, value);
    auto ivkv2 = KVVector<int, float>(ivkv);

    EXPECT_EQ(ivkv.size(), ivkv2.size());
    for(unsigned int i = 0; i < size; i++){
        EXPECT_EQ(ivkv2.key(i), ivkv.key(i));
        EXPECT_FLOAT_EQ(ivkv[i], ivkv2[i]);
    }

}

/**
 * TEST if KVV equals operator overload works with init lists of pairs and values
 * @param il        initializer list of pairs
 * @param ilv       initializer list of values
 */
TEST(KVVector, equals_operator_overload){
    auto il = {std::make_pair(3, 3), std::make_pair(1, 1), std::make_pair(0, 0)};
    KVVector<int, int> ivkv = il;
    EXPECT_EQ(ivkv.key(0), 3);
    EXPECT_EQ(ivkv.at(0), 3);
    EXPECT_EQ(ivkv.size(), 3);

    auto ilv = {5, 4, 3, 2, 1};
    KVVector<int, int> ivkv2 = ilv;
    EXPECT_EQ(ivkv2.key(0), 0);
    EXPECT_EQ(ivkv2[0], 5);
    EXPECT_EQ(ivkv2.size(), 5);
}

/**
 * Testing KVVector assignment with different parameters
 * @param       Key and value
 * @param       Pair
 * @param       value array
 * @param       Key and value array
 * @param       init list of vals and pairs
 */
TEST(KVVector, assign_test){
    std::size_t size = 49;
    float value = 5.0;
    auto ivkv = KVVector<int, float>(size, value);
    EXPECT_EQ(ivkv.at(0), 5.0);
    EXPECT_EQ(ivkv.key(0), 0);
    EXPECT_EQ(ivkv.key(48), 48);
    EXPECT_EQ(ivkv[48], 5.0);

    //assign with key and value
    ivkv.assign(100, 0, 10.0);
    for(int i = 0; i < 100; i++){
        EXPECT_EQ(ivkv.key(i), 0);
        EXPECT_EQ(ivkv[i], 10.0);
    }
    EXPECT_EQ(ivkv.size(), 100);
    EXPECT_EQ(ivkv.capacity(), 100);

    //assign with pair
    ivkv.assign(1, std::make_pair((int)5, (float)6.0));
    EXPECT_EQ(ivkv.size(), 1);
    EXPECT_EQ(ivkv.key(0), 5);
    EXPECT_EQ(ivkv.at(0), 6.0);
    EXPECT_EQ(ivkv.capacity(), 100);

    //assign with value array
    float farr[] = {3.0, 1.0, 4.0, 1.0, 5.0};
    ivkv.assign(5, farr);
    EXPECT_EQ(ivkv.size(), 5);
    for(int i = 0; i < 5; i++){
        EXPECT_EQ(ivkv.key(i), i);
        EXPECT_EQ(ivkv[i], farr[i]);
    }

    //assign with key and value array
    int karr[] = {10, 5, 2, 3, 4};
    ivkv.assign(5, karr, farr);
    for(int i = 0; i < 5; i++){
        EXPECT_EQ(ivkv.key(i), karr[i]);
        EXPECT_EQ(ivkv[i], farr[i]);
    }

    //iterator assign
    KVVector<int, float>::kiterator first = ivkv.k_begin();
    KVVector<int, float>::kiterator last = ivkv.k_end();
    auto ivkv2 = KVVector<int, float>(ivkv);
    ivkv2.sort_vkd();
    EXPECT_FALSE(ivkv2 == ivkv);
    ivkv2.assign(first, last, ivkv);
    EXPECT_TRUE(ivkv == ivkv2);

    //initializer list of pairs
    auto il = {std::make_pair(5, 10.0f), std::make_pair(1, 12.0f), std::make_pair(3, 7.0f)};
    ivkv.assign(il);
    EXPECT_EQ(ivkv.size(), 3);
    EXPECT_EQ(ivkv.key(0), 5);
    EXPECT_EQ(ivkv.at(0), 10.0);
    EXPECT_EQ(ivkv.key(2), 3);
    EXPECT_EQ(ivkv[2], 7.0);

    //init list of values
    auto ilv = {1.0f, 2.0f, 3.0f};
    ivkv.assign(ilv);
    EXPECT_EQ(ivkv.size(), 3);
    for(int i = 0; i < 3; i++){
        EXPECT_EQ(ivkv.key(i), i);
        EXPECT_EQ(ivkv[i], (float)i+1);
    }

}


/** Comparison Operator Tests */

/**
 * Test == operator
 * Test equals equals operator checks size, values and keys
 * of different base types
 * Return false if any are not equal
 */
TEST(KVVector, equals_equals_operator){
    //check with ints
    auto ivkv1 = KVVector<int, int>(1);
    auto ivkv2 = KVVector<int, int>(10, 5);
    EXPECT_FALSE(ivkv1==ivkv2);
    auto ivkv3 = KVVector<int, int>(ivkv2);
    EXPECT_TRUE(ivkv2==ivkv3);

    ivkv3.sort_vkd();
    EXPECT_FALSE(ivkv3==ivkv2);

    ivkv3.sort_vki();
    EXPECT_TRUE(ivkv3==ivkv2);

    ivkv3.pop_back();
    EXPECT_FALSE(ivkv3==ivkv2);

    ivkv3.push_back(5);
    EXPECT_TRUE(ivkv3==ivkv2);

    ivkv3.clear();
    EXPECT_FALSE(ivkv3 == ivkv2);

    //check with floats
    auto ivkvf1 = KVVector<int, float>(10, 5.000002f);
    auto ivkvf2 = KVVector<int, float>(10, 5.000001f);
    EXPECT_FALSE(ivkvf1==ivkvf2);

    //check with doubles
    auto ivkvd1 = KVVector<int, double>(10, 4.999999999999);
    auto ivkvd2 = KVVector<int, double>(10, 5.000000000000);
    EXPECT_FALSE(ivkvd1==ivkvd2);
}

/**
 * Test not equals operator
 * checks size, keys and values of different base types
 * Return false if all are equal
 */
TEST(KVVector, not_equals_operator) {
    auto ivkv1 = KVVector<int, int>(1);
    auto ivkv2 = KVVector<int, int>(10, 5);
    EXPECT_TRUE(ivkv1 != ivkv2);
    auto ivkv3 = KVVector<int, int>(ivkv2);
    EXPECT_FALSE(ivkv2 != ivkv3);

    ivkv3.sort_vkd();
    EXPECT_TRUE(ivkv3 != ivkv2);

    ivkv3.sort_vki();
    EXPECT_FALSE(ivkv3 != ivkv2);

    ivkv3.pop_back();
    EXPECT_TRUE(ivkv3 != ivkv2);

    ivkv3.push_back(5);
    EXPECT_FALSE(ivkv3 != ivkv2);

    ivkv3.clear();
    EXPECT_TRUE(ivkv3 != ivkv2);

    //check with floats
    auto ivkvf1 = KVVector<int, float>(10, 5.000002f);
    auto ivkvf2 = KVVector<int, float>(10, 5.000001f);
    EXPECT_TRUE(ivkvf1 != ivkvf2);

    //check with doubles
    auto ivkvd1 = KVVector<int, double>(10, 4.999999999999);
    auto ivkvd2 = KVVector<int, double>(5, 5.000000000000);
    EXPECT_TRUE(ivkvd1 != ivkvd2);

}

/**
 * Tests less than and greater than operator
 * Checks if values are less/greater than and if size is same or not
 */
TEST(KVVector, less_greater_than_operator){
    size_t n = 5;
    int keys[] = {1, 2, 3, 4, 5};
    int values[] = {6, 7, 8, 9, 10};
    auto ivkv1 = KVVector<int, int>(n, keys, values);
    auto ivkv2 = KVVector<int, int>(n, keys, values);
    //same size so its equal
    EXPECT_FALSE(ivkv1<ivkv2);
    EXPECT_TRUE(ivkv1<=ivkv2);
    EXPECT_FALSE(ivkv1>ivkv2);
    EXPECT_TRUE(ivkv1>=ivkv2);

    //size changes so it is less than
    ivkv2.push_back(11);
    EXPECT_TRUE(ivkv1<ivkv2);
    EXPECT_TRUE(ivkv1<=ivkv2);
    EXPECT_TRUE(ivkv2>ivkv1);
    EXPECT_TRUE(ivkv2>=ivkv1);

    //both 0
    ivkv1.clear();
    ivkv2.clear();
    EXPECT_FALSE(ivkv1 < ivkv1);
    EXPECT_TRUE(ivkv1>=ivkv2);
    EXPECT_FALSE(ivkv1>ivkv2);

    //Testing with sorting
    for(int i = 0; i < 3; i++){
        ivkv1.push_back(i);
        ivkv2.push_back(i*2);
    }
    EXPECT_TRUE(ivkv1<= ivkv2);
    EXPECT_FALSE(ivkv2<ivkv1);
    EXPECT_TRUE(ivkv2>=ivkv1);
    EXPECT_TRUE(ivkv2>ivkv1);
    ivkv1.sort_vkd();
    EXPECT_FALSE(ivkv1<ivkv2);
    ivkv2.sort_vkd();
    EXPECT_TRUE(ivkv1<ivkv2);

    //Test when vector1 is larger but has one smaller number than vector 2
    int keya[] = {0, 5, 7};
    int vala[] = {3, 4, 5};

    int keyb[] = {1, 4, 6, 8};
    int valb[] = {2, 5, 7, 8};
    auto ivkva = KVVector<int, int>(3, keya, vala);
    auto ivkvb = KVVector<int, int>(4, keyb, valb);


    EXPECT_TRUE(ivkvb > ivkva );
    EXPECT_FALSE(ivkvb <= ivkva );

}


/** Sort Tests */

/**
     * Test if all sort functions are accurate
     * @param n        size of vector
     * @param ivkv     Vector
     * @return true    if sort functions sort in correct order
*/
TEST(KVVector, sort) {
    // number of key-value pairs to add to the vector
    int n = 1024;

    // make a KV vector of <int, float>
    KVVector<int, float> ivkv = KVVector<int, float>();
    for(int i=0; i < n; ++i){
        ivkv.push_back(i, (float)i*2);
    }

    ivkv.sort_vkd();
    EXPECT_EQ(ivkv.key(0), 1023);
    EXPECT_NEAR(ivkv.at(0), (float) 2*1023, 1e-7);
    for(int i = 1; i < n-1; i++){
        EXPECT_GE(ivkv.key(i), ivkv.key(i+1));
        EXPECT_GT(ivkv[i], ivkv[i+1]);
    }

    ivkv.sort_vki();
    EXPECT_EQ(ivkv.key(0), 0);
    EXPECT_NEAR(ivkv.at(0), (float) 0, 1e-7);
    for(int i = 1; i < n-1; i++){
        EXPECT_LE(ivkv.key(i), ivkv.key(i+1));
        EXPECT_LT(ivkv[i], ivkv[i+1]);
    }

    KVVector<int, float> ivkv2 = KVVector<int, float>();
    //keys increase, vals decrease
    for(int i=0; i < n; ++i){
        ivkv2.push_back(i, (float)(n-i));
    }
    ivkv2.sort_kvd();
    for(int i = 1; i < n-1; i++){
        EXPECT_GT(ivkv2.key(i), ivkv2.key(i+1));
        EXPECT_LT(ivkv2[i], ivkv2[i+1]);
    }
    ivkv2.sort_kvi();
    for(int i = 1; i < n-1; i++){
        EXPECT_LT(ivkv2.key(i), ivkv2.key(i+1));
        EXPECT_GT(ivkv2[i], ivkv2[i+1]);
    }
    ivkv2.sort_vkd();
    for(int i = 1; i < n-1; i++){
        EXPECT_LT(ivkv2.key(i), ivkv2.key(i+1));
        EXPECT_GT(ivkv2[i], ivkv2[i+1]);
    }

}

/**
 * Test if sort and clear functions work
 * @param n        size of vector
 * @param ivkv     Vector
 * @return true    if sort functions sort in correct order after clearing
*/
TEST(KVVector, sort_clear_sort) {
    // number of key-value pairs to add to the vector
    int n = 1024;

    // make a KV vector of <int, float>
    KVVector<int, float> ivkv = KVVector<int, float>();
    for(int i=0; i < n; ++i){
        ivkv.push_back(i, (float)i*2);
    }

    ivkv.sort_vkd();
    EXPECT_EQ(ivkv.key(0), 1023);

    ivkv.clear();
    for(int i=0; i < n/2; ++i){
        ivkv.push_back(i+1, (float)i+1);
    }

    ivkv.sort_vki();
    EXPECT_EQ(ivkv.key(0), 1);
    EXPECT_EQ((ivkv.k_back()), n/2);

    ivkv.sort_vkd();
    EXPECT_EQ(ivkv.at(0), n/2);
    EXPECT_EQ((ivkv.v_back()), 1);
}

/**
 * Test sort functions that sort with both key and value.
 * Testing whether sort works when key1 > key2 but val1 < val2
 * @param n        size of vector
 * @param keys[]   Array of keys
 * @param vals[]   Array of values
 * @return true    if sort functions that use keys and values sort properly
*/
TEST(KVVector, kv_sort){
    int n = 6;
    int keys[] = {1, 2, 3, 4, 5, 6};
    double vals[] = {5.0, 4.0, 3.0, 2.0, 1.0, 6.0};
    auto ivkv1 = KVVector<int, double>(n, keys, vals);
    auto ivkv2 = KVVector<int, double>(n, keys, vals);

    ivkv1.sort_vki();
    ivkv2.sort_vkd();
    //should be different
    for(int i = 0; i < n; i++) {
        EXPECT_NE(ivkv1.key(i), ivkv2.key(i));
        EXPECT_NE(ivkv1[i], ivkv2[i]);
    }

}

/**Capacity Function Tests*/

/**
 * Test capacity functions: resize, capacity, and size by
 * pushing objects and resizing and
 * expecting correct resize values, size, and capacity
 */
TEST(KVVector, capacity_functions){
    KVVector<int, int> ivkv = KVVector<int, int>(0);
    EXPECT_TRUE(ivkv.empty());
    EXPECT_EQ(ivkv.size(), 0);
    EXPECT_EQ(ivkv.capacity(), 0);
    for(int i = 0; i < 4; i++){
        ivkv.push_back(i);
    }
    EXPECT_EQ(ivkv.size(), ivkv.capacity());
    ivkv.push_back(4);
    EXPECT_EQ(ivkv.capacity(), 8);
    EXPECT_EQ(ivkv.size(), 5);

    ivkv.resize(9);
    EXPECT_EQ(ivkv.capacity(), 9);
    EXPECT_EQ(ivkv.size(), 9);

    //check if resize changed old values
    for(int i = 0; i < 5; i++){
        EXPECT_EQ(ivkv.key(i), i);
        EXPECT_EQ(ivkv[i], i);
    }

    ivkv.resize(1);
    EXPECT_EQ(ivkv.size(), 1);
    EXPECT_EQ(ivkv.capacity(), 1);
    EXPECT_EQ(ivkv.key(0), 0);
    EXPECT_EQ(ivkv.at(0), 0);

    ivkv.resize(5, std::make_pair(0, 1));
    EXPECT_EQ(ivkv.size(), 5);
    for(int i = 1; i < 5; i++){
        EXPECT_EQ(ivkv.key(i), 0);
        EXPECT_EQ(ivkv[i], 1);
    }

    ivkv.resize(7, 10);
    EXPECT_EQ(ivkv.size(), 7);
    for(int i = 5; i < 7; i++){
        EXPECT_EQ(ivkv.key(i), i);
        EXPECT_EQ(ivkv[i], 10);
    }

    ivkv.resize(9, 12, 17);
    EXPECT_EQ(ivkv.size(), 9);
    for (int i = 7; i < 9; i++){
        EXPECT_EQ(ivkv.key(i), 12);
        EXPECT_EQ(ivkv[i], 17);
    }

}

/**Data Access tests */

/**
* Test if data access functions
*/
TEST(KVVector, data_access){
    auto ivkv = KVVector<int, int>(5, 5);
    EXPECT_NO_FATAL_FAILURE(ivkv.k_data());
    EXPECT_NO_FATAL_FAILURE(ivkv.v_data());
}

/**Element Access tests*/

/**
 * Test element access operators
 * Should return element at index
 * Tests index out of bounds throw
 */
TEST(KVVector, array_operator){
    auto ivkv = KVVector<int, int>(5, 10);
    auto ivkv2 = KVVector<int, int>();

    std::pair<int, int> p1(0, 10);
    std::pair<int, int> p2(0, 0);
    EXPECT_EQ(ivkv.v_data()[0], p1.second);

    EXPECT_EQ(ivkv.k_front(), p1.first);
    EXPECT_EQ(ivkv.k_back(), 4);
    EXPECT_EQ(ivkv.k_data()[1], ivkv.key(1));
    EXPECT_EQ(ivkv.at(2), ivkv.val(2));

    EXPECT_ANY_THROW(ivkv[11]);
    EXPECT_ANY_THROW(ivkv.at(11));

}


/** Modifier Tests */

/**
 * Tests emplace_back by checking values and size
 */
TEST(KVVector, emplace_back){
    auto ivkv = KVVector<int, int>(1);
    ivkv.emplace_back(2);
    EXPECT_EQ(ivkv.key(0), 0);
    EXPECT_EQ(ivkv.at(0), 2);
    ivkv.emplace_back(54);
    ivkv.emplace_back( 3);
    ivkv.emplace_back(4);
    for(int i = 0; i < 100; i++){
        ivkv.emplace_back(i);
        EXPECT_EQ(ivkv[ivkv.size()-1], i);
    }
}

    /**
     * Test emplace by checking values and size
     * Large emplace checks if rellacoation during emplace works
     */
TEST(KVVector, emplace){
    auto ivkv = KVVector<int, int>(1);
    for(int i = 0; i < 200; i++){
        ivkv.emplace(ivkv.k_end()-i, i);
        ivkv.emplace(ivkv.k_begin()+i, i);
    }
    EXPECT_EQ(ivkv.size(), 400);
    EXPECT_EQ(ivkv[399], 0);
    EXPECT_EQ(ivkv.at(0), 0);
    EXPECT_EQ(ivkv.key(399), 0);
    EXPECT_EQ(ivkv.key(0), 1);

    for(int i = 0; i < 200; i++){
        EXPECT_EQ(ivkv.key(i), 1 + 2*i);
    }
}

/**
 * Test insertion of pair
 * Check if all values are correct and if large insertion works
 * Check if throw properly throws out of bounds exceptions
 */
TEST(KVVector, insert_pair){
    auto ivkv = KVVector<int, int>(3, 5);
    ivkv.insert(ivkv.k_end(), std::make_pair(0, 1));

    ivkv.insert(ivkv.k_begin(), std::make_pair(11, 7));
    ivkv.insert(ivkv.k_end(), std::make_pair(12, 9));
    ivkv.insert(ivkv.k_begin(), std::make_pair(3, 11));
    ivkv.insert(ivkv.k_end(), std::make_pair(15, 23));
    ivkv.insert(ivkv.k_end(), std::make_pair(99, 97));

    ivkv.insert(ivkv.k_begin(), std::make_pair(14, 13));
    KVVector<int,int>::kiterator it = ivkv.insert(ivkv.k_begin(), std::make_pair(15, 17));
    EXPECT_EQ(*it, 15);
    EXPECT_EQ(ivkv[it-ivkv.k_data()], 17);

    EXPECT_ANY_THROW(ivkv.insert(ivkv.k_begin()-10, std::make_pair(12, 12)));

    EXPECT_EQ(ivkv.size(), 11);

    EXPECT_EQ(ivkv.key(0), 15);
    EXPECT_EQ(ivkv.at(0), 17);

    EXPECT_EQ(ivkv.key(ivkv.size()-1), 99);
    EXPECT_EQ(ivkv[ivkv.size()-1], 97);

    for(int i = 0; i < 100; i++){
        ivkv.insert(ivkv.k_begin()+i, std::make_pair(i, i));
    }
}

/**
 * Checks insertion of single value
 * Checks if key is equal to the key of element inserted behind
 * Check memory reallocation with large insertion
 */
TEST(KVVector, insert_value){
    auto ivkv = KVVector<int, int>(2, 5);
    ivkv.insert(ivkv.k_begin(), 7);
    EXPECT_EQ(ivkv.key(0), ivkv.key(1) -1);
    EXPECT_EQ(ivkv.at(0), 7);

    ivkv.insert(ivkv.k_begin(), 11);

    EXPECT_EQ(ivkv.key(0), ivkv.key(1)-1);
    EXPECT_EQ(ivkv.at(0), 11);
    EXPECT_EQ(ivkv.size(), 4);
    for(int i = 0; i < 10; i++){
        ivkv.insert(ivkv.k_end(), i);
    }

    ivkv.insert(ivkv.k_begin()+4,1);
    for(unsigned int i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv.key(i),i);
    }
    for(int i = 0; i < 15; i++){
        ivkv.push_back(1,3);
    }

    ivkv.insert(ivkv.k_end()-15,1);
    for(int i = 0; i < 16; i++){
        EXPECT_EQ(ivkv.key(i+15), i+1);
    }
}

/**
 * Checks insertion with multiple pairs
 * Large insertion checks proper memory reallocation
 */
TEST(KVVector, insert_mult_pairs){
    auto ivkv = KVVector<int, int>(2, 5);
    KVVector<int, int>::kiterator it = ivkv.insert(ivkv.k_begin(), 100, std::make_pair(5, 12));
    EXPECT_EQ(it, ivkv.k_begin());
    ivkv.insert(ivkv.k_end(), 100, std::make_pair(10, 23));
    EXPECT_EQ(ivkv.size(), 202);
    EXPECT_EQ(ivkv.key(201), 10);
    EXPECT_EQ(ivkv[201], 23);
    EXPECT_EQ(ivkv.key(0), 5);
    EXPECT_EQ(ivkv.at(0), 12);
}

/**
 * Checks insertion with multiple values
 * Large insertion checks memory reallocation
 * Checks if keys iterate properly
 * Check throw if out of bounds iterator
 */
TEST(KVVector, insert_mult_vals){
    auto ivkv = KVVector<int, int>(2, 5);
    std::size_t size = 100;
    int value1 = 13;
    int value2 = 7;
    ivkv.insert(ivkv.k_end(), size, value1);
    ivkv.insert(ivkv.k_begin(), size, value2);
    EXPECT_EQ(ivkv.size(), 202);

    for(unsigned int i = 0; i < size; i++){
        EXPECT_EQ(ivkv.key(i), i);
        EXPECT_EQ(ivkv[i], value2);
    }

    EXPECT_EQ(ivkv.key(size), 100);
    EXPECT_EQ(ivkv[size], 5);
    EXPECT_EQ(ivkv.key(size+1), size+1);
    EXPECT_EQ(ivkv[size+1], 5);


    for(std::size_t i = size+2; i < 2*size+2; i++){
        EXPECT_EQ(ivkv.key(i), i);
        EXPECT_EQ(ivkv[i], value1);
    }

    ivkv.clear();
    for(int i = 0; i < 100; i++){
        ivkv.push_back(3,i);
    }
    ivkv.insert(ivkv.k_begin()+50, (std::size_t)10, 10);
    EXPECT_EQ(ivkv.key(50), 3);
    EXPECT_EQ(ivkv[50], 10);
    for(unsigned int i = 0; i < (ivkv.size()-50); i++){
        EXPECT_EQ(ivkv.key(i+50), 3+i);
    }
    ivkv.insert(ivkv.k_begin(), (std::size_t)1, 5);
    EXPECT_EQ(ivkv.key(0), 3);
    EXPECT_EQ(ivkv.at(0), 5);
    for(unsigned int i = 1; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv.key(i), 3+i);
    }

    EXPECT_ANY_THROW((ivkv.insert(ivkv.k_begin()-2, size, 10)));

    KVVector<int,int>::kiterator it = ivkv.insert(ivkv.k_begin() + 1, size, 3);
    EXPECT_EQ(ivkv[it-ivkv.k_data()], 3);
    EXPECT_EQ(*it, 4);

}

/**
 * Checks insertion with iterators
 * Uses large insertion to check proper memory reallocation
 * Checks throw statement for invalid inputs
 */
TEST(KVVector, insert_with_iterators){
    auto ivkv = KVVector<int, int>(100, 5);
    auto first = ivkv.k_begin();
    auto second = ivkv.k_end();
    auto ivkv2 = KVVector<int, int>(2, 1);
    KVVector<int,int>::kiterator it = ivkv2.insert(ivkv2.k_begin(), first, second, ivkv);
    EXPECT_EQ(it, ivkv2.k_begin());
    ivkv2.insert(ivkv2.k_end(), first, second, ivkv);

    EXPECT_ANY_THROW(ivkv2.insert(ivkv2.k_end(), second, first, ivkv));
    EXPECT_EQ(ivkv2.size(), 202);
    for(int i = 0; i < 100; i++){
        EXPECT_EQ(ivkv2.key(i), ivkv2.key(102+i));
        EXPECT_EQ(ivkv2[i], ivkv2[102+i]);
    }
}

/**
 * Tests insertion with lists of values and pairs
 * Checks if size and values are the same
 */
TEST(KVVector, insert_with_list){
    auto ivkv = KVVector<int, int>(2, 5);
    auto ilv = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    auto ilp = {std::make_pair(1, 2), std::make_pair(2, 3), std::make_pair(3, 4), std::make_pair(4, 5), std::make_pair(5, 6), };
    KVVector<int,int>::kiterator it = ivkv.insert(ivkv.k_begin(), ilp);
    EXPECT_EQ(it, ivkv.k_begin());
    for(unsigned int i = 0; i < ilp.size(); i++){
        EXPECT_EQ(ivkv.key(i), (ilp.begin()+i)->first);
        EXPECT_EQ(ivkv.val(i), (ilp.begin()+i)->second);
    }
    ivkv.insert(ivkv.k_end(), ilv);
    ivkv.insert(ivkv.k_begin(), ilv);
    EXPECT_EQ(ivkv.size(), 37);
    for(unsigned int i = 0; i < ilv.size(); i++){
        EXPECT_EQ(ivkv[i], ivkv[i+22]);
        EXPECT_EQ(ivkv[i], *(ilv.begin()+i));
    }
}

/**
 * Test merge_keys
 */
TEST(KVVector, merge_keys_additive){
    int testvalarr[15] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    int testkeyarr[15] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    auto test = KVVector<int, int>(15,testkeyarr,testvalarr);
    test.merge_keys('+');
    for(int i = 0; i < 15; i++){
        EXPECT_EQ(test.key(i), i + 1);
    }
    EXPECT_EQ(test.size(),15);

    auto ivkv = KVVector<int, int>(10, 1,5);
    ivkv.merge_keys('+');
    EXPECT_EQ(ivkv.size(), 1);
    EXPECT_EQ(ivkv.at(0), 50);
    EXPECT_EQ(ivkv.key(0),1);

    //should not merge anything bc all keys are different
    for(int i = 0; i < 10; i++){
        ivkv.push_back(i%2, i);
    }

    ivkv.merge_keys('+');

    EXPECT_EQ(ivkv.size(), 11);
    EXPECT_EQ(ivkv.at(0), 50);
    for(int i = 1; i < 11; i++){
        EXPECT_EQ(ivkv[i], i-1);
    }

    //should merge after index 10
    for(int i = 0; i < 10; i++){
        ivkv.push_back(10,i);
    }
    ivkv.merge_keys('+');
    EXPECT_EQ(ivkv.size(), 12);
    EXPECT_EQ(ivkv[11], 9 * (9+1) / 2);


    //merge with negatives
    ivkv.clear();
    for(int i = 0; i < 10; i++){
        ivkv.push_back(1, (i % 2 ? i : -1 * i)); //odds are negative, keys are same
    }
    ivkv.merge_keys('+');
    EXPECT_EQ(ivkv.size(),1);
    EXPECT_EQ(ivkv.at(0), (9+1)/2);

    //merge pairs
    ivkv.clear();
    for(int i = 0; i < 10; i++){
        ivkv.push_back(i % 2, i);
        ivkv.push_back(i % 2, i);
    }
    ivkv.merge_keys('+');
    EXPECT_EQ(ivkv.size(), 10);
    for(unsigned int i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv[i], i * 2);
    }

    //merge triplets of huge input
    ivkv.clear();
    for(int i = 0; i < 10000; i++){
        ivkv.push_back(i % 3, i);
        ivkv.push_back(i % 3, i);
        ivkv.push_back(i % 3, i);
    }
    ivkv.merge_keys('+');
    EXPECT_EQ(ivkv.size(), 10000);
    for(unsigned int i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv[i], i * 3);
    }
}

TEST(KVVector, merge_keys_less_than){
    auto ivkv = KVVector<int, int>(10, 1,5);
    ivkv.merge_keys('<');
    EXPECT_EQ(ivkv.size(), 1);
    EXPECT_EQ(ivkv.at(0), 5);

    //should not merge anything bc all keys are different
    for(int i = 0; i < 10; i++){
        ivkv.push_back(i%2, i);
    }
    ivkv.merge_keys('<');
    EXPECT_EQ(ivkv.size(), 11);
    EXPECT_EQ(ivkv.at(0), 5);
    for(int i = 1; i < 11; i++){
        EXPECT_EQ(ivkv[i], i-1);
    }

    //should merge after index 10
    for(int i = 0; i < 10; i++){
        ivkv.push_back(10,i);
    }
    ivkv.merge_keys('<');
    EXPECT_EQ(ivkv.size(), 12);
    EXPECT_EQ(ivkv[11], 0);


    //merge with negatives
    ivkv.clear();
    for(int i = 0; i < 10; i++){
        ivkv.push_back(1, (i % 2 ? i : -1 * i)); //odds are negative, keys are same
    }
    ivkv.merge_keys('<');
    EXPECT_EQ(ivkv.size(),1);
    EXPECT_EQ(ivkv.at(0), -8);

    //merge pairs
    ivkv.clear();
    for(int i = 0; i < 10; i++){
        ivkv.push_back(i % 2, i);
        ivkv.push_back(i % 2, (i * -1));
    }
    ivkv.merge_keys('<');
    EXPECT_EQ(ivkv.size(), 10);
    for(unsigned int i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv[i], (i * -1));
    }

    //merge triplets of huge input
    ivkv.clear();
    for(int i = 0; i < 10000; i++){
        ivkv.push_back(i % 3, i+10);
        ivkv.push_back(i % 3, i+100);
        ivkv.push_back(i % 3, i+1000);
    }
    ivkv.merge_keys('<');
    EXPECT_EQ(ivkv.size(), 10000);
    for(unsigned int i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv[i], (i + 10));
    }
}

TEST(KVVector, merge_keys_greater_than){
    auto ivkv = KVVector<int, int>(10, 1,5);
    ivkv.merge_keys('>');
    EXPECT_EQ(ivkv.size(), 1);
    EXPECT_EQ(ivkv.at(0), 5);

    //should not merge anything bc all keys are different
    for(int i = 0; i < 10; i++){
        ivkv.push_back(i%2, i);
    }
    ivkv.merge_keys('>');
    EXPECT_EQ(ivkv.size(), 11);
    EXPECT_EQ(ivkv.at(0), 5);
    for(int i = 1; i < 11; i++){
        EXPECT_EQ(ivkv[i], i-1);
    }

    //should merge after index 10
    for(int i = 0; i < 10; i++){
        ivkv.push_back(10,i);
    }
    ivkv.merge_keys('>');
    EXPECT_EQ(ivkv.size(), 12);
    EXPECT_EQ(ivkv[11], 9);


    //merge with negatives
    ivkv.clear();
    for(int i = 0; i < 10; i++){
        ivkv.push_back(1, (i % 2 ? i : -1 * i)); //odds are negative, keys are same
    }
    ivkv.merge_keys('>');
    EXPECT_EQ(ivkv.size(),1);
    EXPECT_EQ(ivkv.at(0), 9);

    //merge pairs
    ivkv.clear();
    for(int i = 0; i < 10; i++){
        ivkv.push_back(i % 2, i);
        ivkv.push_back(i % 2, (i * -1));
    }
    ivkv.merge_keys('>');
    EXPECT_EQ(ivkv.size(), 10);
    for(std::size_t i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv[i], i);
    }

    //merge triplets of huge input
    ivkv.clear();
    for(int i = 0; i < 10000; i++){
        ivkv.push_back(i % 3, i+10);
        ivkv.push_back(i % 3, i+100);
        ivkv.push_back(i % 3, i+1000);
    }
    ivkv.merge_keys('>');
    EXPECT_EQ(ivkv.size(), 10000);
    EXPECT_EQ(ivkv[9999], 10999);
    for(std::size_t i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv[i], (i + 1000));
    }
}

TEST(KVVector, merge_keys_absval){
    auto ivkv = KVVector<int, int>(10, 1,5);
    ivkv.merge_keys('^');
    EXPECT_EQ(ivkv.size(), 1);
    EXPECT_EQ(ivkv.at(0), 5);

    //should not merge anything bc all keys are different
    for(int i = 0; i < 10; i++){
        ivkv.push_back(i%2, i);
    }
    ivkv.merge_keys('^');
    EXPECT_EQ(ivkv.size(), 11);
    EXPECT_EQ(ivkv.at(0), 5);
    for(int i = 1; i < 11; i++){
        EXPECT_EQ(ivkv[i], i-1);
    }

    //should merge after index 10
    for(int i = 0; i < 10; i++){
        ivkv.push_back(10,i);
    }
    ivkv.merge_keys('^');
    EXPECT_EQ(ivkv.size(), 12);
    EXPECT_EQ(ivkv[11], 9);

    //merge with negatives
    ivkv.clear();
    for(int i = 0; i < 9; i++){
        ivkv.push_back(1, (i % 2 ? i : -1 * i)); //odds are negative, keys are same
    }
    ivkv.merge_keys('^');
    EXPECT_EQ(ivkv.size(),1);
    EXPECT_EQ(ivkv.at(0), -8);

    //merge pairs
    ivkv.clear();
    for(int i = 0; i < 10; i++){
        ivkv.push_back(i % 2, i);
        ivkv.push_back(i % 2, ((i+1) * -1));
    }
    ivkv.merge_keys('^');
    EXPECT_EQ(ivkv.size(), 10);
    for(std::size_t i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv[i], (i+1)*-1);
    }

    //merge triplets of huge input
    ivkv.clear();
    for(int i = 0; i < 10000; i++){
        ivkv.push_back(i % 3, i+10);
        ivkv.push_back(i % 3, i+100);
        ivkv.push_back(i % 3, i+1000);
    }
    ivkv.merge_keys('^');
    EXPECT_EQ(ivkv.size(), 10000);
    EXPECT_EQ(ivkv[9999], 10999);
    for(unsigned int i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv[i], (i + 1000));
    }
}

TEST(KVVector, merge_keys_mean){
    auto ivkv = KVVector<int, int>(10, 1,5);
    ivkv.merge_keys('m');
    EXPECT_EQ(ivkv.size(), 1);
    EXPECT_EQ(ivkv.at(0), 5);
    //should not merge anything bc all keys are different
    for(int i = 0; i < 10; i++){
        ivkv.push_back(i%2, i);
    }

    ivkv.merge_keys('m');

    EXPECT_EQ(ivkv.size(), 11);
    EXPECT_EQ(ivkv.at(0), 5);
    for(int i = 1; i < 11; i++){
        EXPECT_EQ(ivkv[i], i-1);
    }

    //should merge after index 10
    for(int i = 0; i < 10; i++){
        ivkv.push_back(10,i);
    }
    ivkv.merge_keys('m');
    EXPECT_EQ(ivkv.size(), 12);
    EXPECT_EQ(ivkv[11], 4);

    //merge with negatives
    ivkv.clear();
    for(int i = 0; i < 10; i++){
        ivkv.push_back(1, (i % 2 ? i : -1 * i)); //odds are negative, keys are same
    }
    ivkv.merge_keys('m');
    EXPECT_EQ(ivkv.size(),1);
    EXPECT_EQ(ivkv.at(0), (9+1)/18);

    //merge pairs
    ivkv.clear();
    for(int i = 0; i < 10; i++){
        ivkv.push_back(i % 2, i);
        ivkv.push_back(i % 2, i);
    }
    ivkv.merge_keys('m');
    EXPECT_EQ(ivkv.size(), 10);
    for(unsigned int i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv[i], i);
    }

    //merge triplets of huge input
    ivkv.clear();
    for(int i = 0; i < 10000; i++){
        ivkv.push_back(i % 3, i);
        ivkv.push_back(i % 3, i);
        ivkv.push_back(i % 3, i);
    }
    ivkv.merge_keys('m');
    EXPECT_EQ(ivkv.size(), 10000);
    for(unsigned int i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv[i], i);
    }
}

/**
 * Test scan
 */
TEST(KVVector, scan_additive){
    //adds all the values
    auto ivkv = KVVector<int, int>(10, 1, 3);
    ivkv.scan('+');
    for(unsigned int i = 0; i < ivkv.size(); i++) {
        EXPECT_EQ(ivkv.key(i), 1);
        EXPECT_EQ(ivkv[i], (i+1) * 3);
    }

    //no change
    auto ivkv_dif_key = KVVector<int, int>(100, 55);
    ivkv.scan('+');
    for(unsigned int i = 0; i < ivkv_dif_key.size(); i++){
        EXPECT_EQ(ivkv_dif_key.key(i), i);
        EXPECT_EQ(ivkv_dif_key[i], 55);
    }

    auto ivkv_zero =  KVVector<float, double>(0);
    ivkv.scan('+');

    auto ivkv_neg = KVVector<float, double>();
    for(int i = 0; i < 50; i++){
        ivkv_neg.push_back(std::make_pair((float)1.0, (i % 2 == 0 ? i : i * -1 )));
    }

    ivkv_neg.scan('+');

    for(int i = 0; i < 50; i++){
        EXPECT_FLOAT_EQ(ivkv_neg.key(i),1.0);
        //every two pairs after 0 have same value but different signs
        EXPECT_DOUBLE_EQ(ivkv_neg[i], (i % 2 == 0 ? (i+1)/2 : -(i+1)/2 ) );
    }

    EXPECT_ANY_THROW(ivkv.scan('a'));
}

TEST(KVVector, scan_less_than){
    auto ivkv = KVVector<int, int> (10);
    for(int i = 0; i < 10; i++){
        ivkv.push_back(std::make_pair(i,i));
    }
    ivkv.sort_vkd();
    //should not do anything
    ivkv.scan('<');
    for(int i = 0; i < 10; i++){
        EXPECT_EQ(ivkv[i], 9-i);
    }

    ivkv.sort_vki();
    ivkv.scan('<');
    for(int i = 0; i < 10; i++){
        EXPECT_EQ(ivkv[i], 0);
    }

    ivkv.push_back(-100);
    ivkv.push_back(std::make_pair(1,-5));
    ivkv.scan('<');
    EXPECT_EQ(ivkv[ivkv.size()-1], -100);
    EXPECT_EQ(ivkv.at(0), 0);

    ivkv.sort_vki();
    ivkv.scan('<');
    for(unsigned int i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv[i], -100);
    }

}

TEST(KVVector, scan_greater_than){
    auto ivkv = KVVector<int, int>(10,7);
    for(int i = 0; i < 20; i++){
        ivkv.push_back(i);
    }
    ivkv.scan('>');
    for(int i = 0; i < 17; i++){
        EXPECT_EQ(ivkv[i],7);
    }
    for(int i = 0; i < 4; i++){
        EXPECT_EQ(ivkv[17+i], 7 + i);
    }

    ivkv.sort_vkd();
    ivkv.scan('>');
    for(unsigned int i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv[i], 19);
    }

    auto ivkv_neg = KVVector<int,int>(10);
    for(int i = 1; i < 21; i++){
        ivkv_neg.push_back(i * -1);
    }

    ivkv_neg.scan('>');

    for(unsigned int i = 0; i < ivkv_neg.size(); i++){
        EXPECT_EQ(ivkv_neg[i], -1);
    }
    ivkv_neg.push_back(20);
    ivkv_neg.scan('>');
    EXPECT_EQ(ivkv_neg[20], 20);
}

TEST(KVVector, scan_abs_max){
    auto ivkv = KVVector<int,int>();
    for(int i = 20; i >= 0; i--){
        ivkv.push_back(i);
    }
    ivkv.scan('^');

    for(unsigned int i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv[i], 20);
    }

    auto ivkv_neg = KVVector<int,int>();
    for(int i = 20; i >= -30; i--){
        ivkv_neg.push_back(i);
    }

    ivkv_neg.scan('^');
    for(int i = 0; i < 41; i++){
        EXPECT_EQ(ivkv_neg[i], 20);
    }
    for(int i = 1; i < 11; i++){
        EXPECT_EQ(ivkv_neg[i + 40], -1 * (20+i));
    }

    ivkv_neg.sort_vki();
    EXPECT_EQ(ivkv_neg[0], -30);
    EXPECT_EQ(ivkv_neg[50],20);

    ivkv_neg.scan('^');

    for(unsigned int i = 0; i < ivkv_neg.size(); i++){
        EXPECT_EQ(ivkv_neg[i], -30);
    }
}

TEST(KVVector, scan_mean){
    //scan with all same
    auto ivkv = KVVector<int, int>(10, 5);
    ivkv.scan('m');
    for(std::size_t i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv[i], 5);
    }

    //scan with incrementing numbers
    ivkv.clear();
    for(int i = 0; i < 10; i++){
        ivkv.push_back(i);
    }
    ivkv.scan('m');
    int sum = 0;
    for(int i = 0; i < 10; i++){
        sum += i;
        EXPECT_EQ(ivkv[i], sum/(i+1));
    }

    //scan with negative and positive numbers
    ivkv.clear();
    ivkv.push_back(5);
    ivkv.push_back(-7);
    ivkv.push_back(10);
    ivkv.push_back(-20);
    ivkv.push_back(12);
    ivkv.scan('m');
    EXPECT_EQ(ivkv.at(0), 5);  // 5/1
    EXPECT_EQ(ivkv[1], -1); // -2/2
    EXPECT_EQ(ivkv[2], 2);  // 8/3 truncated
    EXPECT_EQ(ivkv[3], -3);  // -12/4
    EXPECT_EQ(ivkv[4], 0); // 0/5

    //scan with doubles
    auto ivkvd = KVVector<int, double>();
    for(int i = 0; i < 100; i++){
        ivkvd.push_back(i + i/10);
    }
    ivkvd.scan('m');
    double sumd = 0;
    for(int i = 0; i < 100; i++){
        sumd += (i + i/10);
        EXPECT_DOUBLE_EQ(ivkvd[i], sumd/(i+1));
    }

    ivkvd.clear();
    ivkvd.push_back(2);
    ivkvd.push_back(1);
    ivkvd.scan('m');
    EXPECT_NE(ivkvd[1], 1.50001);
}

/**
 * Tests erase function by checking
 * Size should be set to 0
 */
TEST(KVVector, erase){
    auto ivkv = KVVector<int, int>(10, 5);
    EXPECT_EQ(ivkv.size(), 10);
    ivkv.erase(ivkv.k_begin());
    EXPECT_EQ(ivkv.size(), 9);
    EXPECT_EQ(ivkv.key(0), 1);
    ivkv.erase(ivkv.k_begin(), ivkv.k_begin()+5);
    EXPECT_EQ(ivkv.key(0), 6);
    EXPECT_EQ(ivkv.size(), 4);
    ivkv.erase(ivkv.k_begin(), ivkv.k_end());
    EXPECT_EQ(ivkv.size(), 0);
}

/**
 * Tests if swap correctly changes size, value and keys
 */
TEST(KVVector, swap){
    auto ivkv = KVVector<int, int>(10, 5);
    auto ivkv2 = KVVector<int, int>(10, 3, 9);
    for(unsigned int i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv.key(i), i);
        EXPECT_EQ(ivkv[i], 5);
        EXPECT_EQ(ivkv2.key(i), 3);
        EXPECT_EQ(ivkv2[i], 9);
    }

    ivkv.swap(ivkv2);
    for(unsigned int i = 0; i < ivkv.size(); i++){
        EXPECT_EQ(ivkv2.key(i), i);
        EXPECT_EQ(ivkv2[i], 5);
        EXPECT_EQ(ivkv.key(i), 3);
        EXPECT_EQ(ivkv[i], 9);
    }

}


/**Iterator Tests*/
/**
 * Check if iterators point to correct spot
 */
TEST(KVVector, iterators){
    auto ivkv = KVVector<int, int>(5, 10);
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
 * Test select
 */
TEST(KVVector, select){
    auto ivkv = new KVVector<int, double>(1000);
    for(int i = 0; i < 1000; i++){
        ivkv->set_key_unsafe(i, i);
        ivkv->set_val_unsafe(i, i % 3 ? 1045 + i * 1.1 : 1045 - i);
    }
    ivkv->set_size(1000);
    // move 50 smallest values from the whole vector to the start of vector
    ivkv->selecti(50);
    //find minimum value
    //should be 46
    //max should be 193
    double min = ivkv->val_unsafe(0);
    double max = ivkv->val_unsafe(0);
    for(int i = 1; i < 50; i++){
        double val = ivkv->val_unsafe(i);
        min = min < val ? min : val;
        max = max > val ? max : val;
    }
    ASSERT_EQ(min, 46);
    ASSERT_EQ(max, 193);

    //check rest of array
    for(int i = 50; i < 1000; i++){
        ASSERT_TRUE(max < ivkv->val_unsafe(i));
    }

    //set in decreasing order
    for(int i = 0; i < 1000; i++){
        ivkv->set_key(i, 1000 - i);
        ivkv->set_val(i, i % 3 ? 1045 + i * 1.1 : 1045 - i);
    }

    ivkv->selecti(50);
    //find minimum value
    //should be 46
    //max should be 193
    min = ivkv->val_unsafe(0);
    max = ivkv->val_unsafe(0);
    for(int i = 1; i < 50; i++){
        double val = ivkv->val_unsafe(i);
        min = min < val ? min : val;
        max = max > val ? max : val;
    }
    ASSERT_EQ(min, 46);
    ASSERT_EQ(max, 193);

    //check rest of array
    for(int i = 50; i < 1000; i++){
        ASSERT_TRUE(max < ivkv->val_unsafe(i));
    }

    //Try selecting from part of the array
    for(int i = 0; i < 1000; i++){
        ivkv->set_key(i, i);
        ivkv->set_val(i, i * 1.0);
    }
    ivkv->selectd(0, 100, 50);
    //find minimum value
    //should be 50
    //max should be 99
    min = ivkv->val_unsafe(0);
    max = ivkv->val_unsafe(0);
    for(int i = 1; i < 50; i++){
        double val = ivkv->val_unsafe(i);
        min = min < val ? min : val;
        max = max > val ? max : val;
    }
    EXPECT_NEAR(min, 50, 1e-5);
    EXPECT_NEAR(max, 99, 1e-5);

    //check rest of array
    for(int i = 50; i < 100; i++){
        ASSERT_TRUE(min > ivkv->val_unsafe(i));
    }

    ivkv->selectd(300, 400, 50);
    //find minimum value
    //should be 350
    //max should be 399
    min = ivkv->val_unsafe(300);
    max = ivkv->val_unsafe(300);
    for(int i = 301; i < 350; i++){
        double val = ivkv->val_unsafe(i);
        min = min < val ? min : val;
        max = max > val ? max : val;
    }
    EXPECT_NEAR(min, 350, 1e-5);
    EXPECT_NEAR(max, 399, 1e-5);

    //check rest of array
    for(int i = 350; i < 400; i++){
        ASSERT_TRUE(min > ivkv->val_unsafe(i));
    }

    delete ivkv;
}

/**
 * Test median
 */
TEST(KVVector, median){
    auto ivkv = new KVVector<int, double>(1000);
    // [0...1000)
    for(int i=0; i < 1000; ++i){
        ivkv->set_val(i, (double) i);
    }
    EXPECT_NEAR(ivkv->median(), 499.5, 1e-7);
    ivkv->pop_back();
    EXPECT_NEAR(ivkv->median(), 499.0, 1e-7);
    // (1000...0]
    for(int i=0; i < 999; ++i){
        ivkv->set_val(i, (double) (999-i));
    }
    ivkv->push_back(0);
    EXPECT_NEAR(ivkv->median(), 499.5, 1e-7);
    ivkv->pop_back();
    EXPECT_NEAR(ivkv->median(), 499.0, 1e-7);
    delete ivkv;

    auto iikv = new KVVector<int, int>(1000);
    // [0...1000)
    for(int i=0; i < 1000; ++i){
        iikv->set_val(i, i);
    }
    EXPECT_EQ(iikv->median(), 499);
    iikv->pop_back();
    EXPECT_EQ(iikv->median(), 499);
    // (1000...0]
    for(int i=0; i < 999; ++i){
        iikv->set_val(i, (999-i));
    }
    iikv->push_back(0);
    EXPECT_EQ(iikv->median(), 499);
    iikv->pop_back();
    EXPECT_EQ(iikv->median(), 499);
    delete iikv;
}

/**
 * Test mean
 */
TEST(KVVector, mean){
    auto ivkv = new KVVector<int, double>(1000);
    // [0...1000)
    for(int i=0; i < 1000; ++i){
        ivkv->set_val(i, (double) i);
    }
    EXPECT_NEAR(ivkv->mean(), 499.5, 1e-7);
    ivkv->pop_back();
    EXPECT_NEAR(ivkv->mean(), 499.0, 1e-7);
    // (1000...0]
    for(int i=0; i < 999; ++i){
        ivkv->set_val(i, (double) (999-i));
    }
    ivkv->push_back(0);
    EXPECT_NEAR(ivkv->mean(), 499.5, 1e-7);
    ivkv->pop_back();
    EXPECT_NEAR(ivkv->mean(), 500.0, 1e-7);
    delete ivkv;

    auto iikv = new KVVector<int, int>(1000);
    // [0...1000)
    for(int i=0; i < 1000; ++i){
        iikv->set_val(i, i);
    }
    EXPECT_NEAR(iikv->mean(), 499.5, 1e-7);
    iikv->pop_back();
    EXPECT_NEAR(iikv->mean(), 499.0, 1e-7);
    // (1000...0]
    for(int i=0; i < 999; ++i){
        iikv->set_val(i, (999-i));
    }
    iikv->push_back(0);
    EXPECT_NEAR(iikv->mean(), 499.5, 1e-7);
    iikv->pop_back();
    EXPECT_NEAR(iikv->mean(), 500.0, 1e-7);
    delete iikv;
}

/**
 * Test stdev
 */
TEST(KVVector, stdev){
    auto ivkv = new KVVector<int, double>(1000);
    // [0...1000)
    for(int i=0; i < 1000; ++i){
        ivkv->set_val(i, (double) i);
    }
    EXPECT_NEAR(ivkv->stdev(), 288.6749902572095, 1e-7);
    EXPECT_NEAR(ivkv->stdev(1), 288.8194360957494, 1e-7);
    ivkv->pop_back();
    EXPECT_NEAR(ivkv->stdev(), 288.38631497813253, 1e-7);
    EXPECT_NEAR(ivkv->stdev(1), 288.5307609250702, 1e-7);
    // (1000...0]
    for(int i=0; i < 999; ++i){
        ivkv->set_val(i, (double) (999-i));
    }
    ivkv->push_back(0);
    EXPECT_NEAR(ivkv->stdev(), 288.6749902572095, 1e-7);
    EXPECT_NEAR(ivkv->stdev(1), 288.8194360957494, 1e-7);
    ivkv->pop_back();
    EXPECT_NEAR(ivkv->stdev(), 288.38631497813253, 1e-7);
    EXPECT_NEAR(ivkv->stdev(1), 288.5307609250702, 1e-7);
    delete ivkv;

    auto iikv = new KVVector<int, int>(1000);
    // [0...1000)
    for(int i=0; i < 1000; ++i){
        iikv->set_val(i, i);
    }
    EXPECT_NEAR(iikv->stdev(), 288.6749902572095, 1e-7);
    EXPECT_NEAR(iikv->stdev(1), 288.8194360957494, 1e-7);
    iikv->pop_back();
    EXPECT_NEAR(iikv->stdev(), 288.38631497813253, 1e-7);
    EXPECT_NEAR(iikv->stdev(1), 288.5307609250702, 1e-7);
    // (1000...0]
    for(int i=0; i < 999; ++i){
        iikv->set_val(i, (999-i));
    }
    iikv->push_back(0);
    EXPECT_NEAR(iikv->stdev(), 288.6749902572095, 1e-7);
    EXPECT_NEAR(iikv->stdev(1), 288.8194360957494, 1e-7);
    iikv->pop_back();
    EXPECT_NEAR(iikv->stdev(), 288.38631497813253, 1e-7);
    EXPECT_NEAR(iikv->stdev(1), 288.5307609250702, 1e-7);
    delete iikv;
}

/**
 * Test min, max, sum
 */
TEST(KVVector, min_max_sum){
    auto ivkv = new KVVector<int, double>(1000);
    // [0...1000)
    for(int i=0; i < 1000; ++i){
        ivkv->set_val(i, (double) i);
    }
    EXPECT_NEAR(ivkv->min(), 0, 1e-7);
    EXPECT_NEAR(ivkv->max(), 999.0, 1e-7);
    EXPECT_NEAR(ivkv->sum(), 499500.0, 1e-7);
    ivkv->pop_back();
    EXPECT_NEAR(ivkv->min(), 0, 1e-7);
    EXPECT_NEAR(ivkv->max(), 998.0, 1e-7);
    EXPECT_NEAR(ivkv->sum(), 498501.0, 1e-7);
    // (1000...0]
    for(int i=0; i < 999; ++i){
        ivkv->set_val(i, (double) (999-i));
    }
    ivkv->push_back(0);
    EXPECT_NEAR(ivkv->min(), 0, 1e-7);
    EXPECT_NEAR(ivkv->max(), 999.0, 1e-7);
    EXPECT_NEAR(ivkv->sum(), 499500.0, 1e-7);
    ivkv->pop_back();
    EXPECT_NEAR(ivkv->min(), 1, 1e-7);
    EXPECT_NEAR(ivkv->max(), 999.0, 1e-7);
    EXPECT_NEAR(ivkv->sum(), 499500.0, 1e-7);
    delete ivkv;

    auto iikv = new KVVector<int, int>(1000);
    // [0...1000)
    for(int i=0; i < 1000; ++i){
        iikv->set_val(i, i);
    }
    EXPECT_EQ(iikv->min(), 0);
    EXPECT_EQ(iikv->max(), 999);
    EXPECT_EQ(iikv->sum(), 499500);
    iikv->pop_back();
    EXPECT_EQ(iikv->min(), 0);
    EXPECT_EQ(iikv->max(), 998);
    EXPECT_EQ(iikv->sum(), 498501);
    // (1000...0]
    for(int i=0; i < 999; ++i){
        iikv->set_val(i, (999-i));
    }
    iikv->push_back(0);
    EXPECT_EQ(iikv->min(), 0);
    EXPECT_EQ(iikv->max(), 999);
    EXPECT_EQ(iikv->sum(), 499500);
    iikv->pop_back();
    EXPECT_EQ(iikv->min(), 1);
    EXPECT_EQ(iikv->max(), 999);
    EXPECT_EQ(iikv->sum(), 499500);
    delete iikv;
}

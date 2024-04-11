/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#pragma once

#include <chrono>
#include <cstdint>
#include <sstream>

namespace snnlib {

    /**
     * Counter class
     */
    typedef struct Counter {
        uint_fast64_t count = 0;
        bool enabled = true;

        inline void reset(){
            count = 0;
        }

        Counter& operator+= (const Counter& rhs){
            count += rhs.count;
            return *this;
        }

        Counter& operator+= (const uint_fast64_t rhs){
            count += rhs;
            return *this;
        }

        Counter& operator++(){
            count++;
            return *this;
        }

        Counter& operator++(int){
            auto temp = *this;
            ++*this;
            return temp;
        }

        Counter& operator=(const uint64_t rhs){
            count = rhs;
            return *this;
        }

        friend Counter operator+ (Counter lhs, const Counter& rhs){
            lhs.count += rhs.count;
            return lhs;
        }

    } Counter;

    std::ostream& operator<< (std::ostream &os, const Counter &c) {
        return (os << c.count);
    }

    /**
     * Add up the counters from source to target.
     * @param source Source timers array.
     * @param target Target timers array.
     * @param n      Size of the source and target arrays.
     */
    void add_counters(
        const Counter* const source,
        Counter* const target,
        size_t const n
    ){
        for(uint8_t i=0; i < n; ++i){
            target[i].count += source[i].count;
        }
    }

}  /* end namespace */

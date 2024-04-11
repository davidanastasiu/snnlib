/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#pragma once

#include <bits/stdc++.h>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <cstring>

namespace snnlib {

    class BitVector {

    private:
        uint64_t size;
        uint64_t sz;
        bool complement = false;
        uint64_t * vec = nullptr;

    public:

        /**
         * Construct a bit vector with `size` elements and all values
         * set to 1 if `complement` or 0 otherwise.
         * @param size Number of elements in the bit vector
         * @param complement Start with all 1s instead of all 0s
         */
        BitVector(const uint64_t size, bool complement = false){
            this->size = size;
            this->sz = (size + 63) / 64;
            this->vec = (uint64_t *) malloc(sizeof(uint64_t) * this->sz);
            if (complement){
                memset(this->vec, ~0, sizeof(uint64_t) * this->sz);
                this->complement = true;
            } else {
                memset(this->vec, 0, sizeof(uint64_t) * this->sz);
            }
        }

        /**
         * Destructor
         */
        ~BitVector(){
            free(this->vec);
        }

        /**
         * Size of allocation for a bit vector that can hold numbers 0 -- size-1
         * @param size Max # elements held in the bit vector
         * @return Size of allocation
         */
        static uint64_t allocation_size(const uint64_t size){
            return (size + 63) / 64;
        }

        /**
         * Reset the bit vector to 0.
         */
        void reset(){
            if (this->complement){
                memset(this->vec, ~0, sizeof(uint64_t) * this->sz);
            } else {
                memset(this->vec, 0, sizeof(uint64_t) * this->sz);
            }
        }

        /**
         * Reset the vector only for the buckets that the elements in `items` are in.
         * @param n     Size of `items` list
         * @param items IDs that should be reset
         */
        void reset_all(const uint64_t n, const uint64_t * items)
        {
            if(n >= this->size / 4){ // faster to just reset the whole buffer
                this->reset();
                return;
            }
            if (this->complement){
                for(uint64_t i=0; i < n; ++ i){
                    this->set_all(items[i]);
                }
            } else {
                for(uint64_t i=0; i < n; ++ i){
                    this->zero_all(items[i]);
                }
            }
        }

        /**
         * Reset the vector only for the buckets that the elements in `items` are in.
         * @param n     Size of `items` list
         * @param items IDs that should be reset
         */
        void reset_all(const uint64_t n, const uint32_t * items)
        {
            if(n >= this->size / 4){ // faster to just reset the whole buffer
                this->reset();
                return;
            }
            if (this->complement){
                for(uint64_t i=0; i < n; ++ i){
                    this->set_all(items[i]);
                }
            } else {
                for(uint64_t i=0; i < n; ++ i){
                    this->zero_all(items[i]);
                }
            }
        }

        /**
         * Reset the vector only for the buckets that the elements in `items` are in.
         * @param n     Size of `items` list
         * @param items IDs that should be reset
         */
        void reset_all(const uint64_t n, const uint16_t * items)
        {
            if(n >= this->size / 4){ // faster to just reset the whole buffer
                this->reset();
                return;
            }
            if (this->complement){
                for(uint64_t i=0; i < n; ++ i){
                    this->set_all(items[i]);
                }
            } else {
                for(uint64_t i=0; i < n; ++ i){
                    this->zero_all(items[i]);
                }
            }
        }

        /**
         * Reset the vector only for the buckets that the elements in `items` are in.
         * @param n     Size of `items` list
         * @param items IDs that should be reset
         */
        void reset_all(const uint64_t n, const uint8_t * items)
        {
            if(n >= this->size / 4){ // faster to just reset the whole buffer
                this->reset();
                return;
            }
            if (this->complement){
                for(uint64_t i=0; i < n; ++ i){
                    this->set_all(items[i]);
                }
            } else {
                for(uint64_t i=0; i < n; ++ i){
                    this->zero_all(items[i]);
                }
            }
        }


        /** Count the number of set bits. */
        uint64_t count()
        {
            uint64_t cnt = 0;
            for(uint64_t i=0; i < sz; ++i){
                cnt += __builtin_popcountl(this->vec[i]);
            }
            return cnt;
        }

        /**
         * Check value of bit vector for given id.
         * For efficiency, no check is made that id is valid...
         * @param  vec    Vector buffer
         * @param  id     ID checking for
         * @return Whether bit at index id is set or not
         */
        inline bool isset(const uint64_t id){
            return (this->vec[id / 64] & ((uint64_t)1 << (63 - (id % 64)))) ? true: false;
        }

        /**
         * Set bit value at index id to 1.
         * For efficiency, no check is made that id is valid...
         * @param  vec    Vector buffer
         * @param  id     index of element
         */
        inline void set(const uint64_t id){
            this->vec[id / 64] |= ((uint64_t)1 << (63 - (id % 64)));
        }

        /**
         * Set bit value at index id, along with all neighbors in the same vector pod, to 1.
         * For efficiency, no check is made that id is valid...
         * @param  vec    Vector buffer
         * @param  id     index of element
         */
        inline void set_all(const uint64_t id){
            this->vec[id / 64] = ~((uint64_t) 0);
        }

        /**
         * Sets the bit to the opposite value,
         * i.e, if it was 1, it will be set to 0
         * @param id index of element
         */
        inline void toggle(const uint64_t id){
            this->vec[id / 64] ^= ((uint64_t)1 << (63 - (id % 64)));
        }

        /**
         * Set bit value at index id to 0.
         * For efficiency, no check is made that id is valid...
         * @param  vec    Vector buffer
         * @param  id     ID checking for
         */
        inline void zero(const uint64_t id){
            this->vec[id / 64] &= ~((uint64_t)1 << (63 - (id % 64)));
        }

        /**
         * Set bit value at index id, along with all neighbors in the same vector pod, to 0.
         * For efficiency, no check is made that id is valid...
         * @param  vec    Vector buffer
         * @param  id     ID checking for
         */
        inline void zero_all(const uint64_t id){
            this->vec[id / 64] = 0;
        }

        /**
         * Access value at given index in the underlying bitvector buffer.
         * For efficiency, no check is made that id is valid...
         * @return        Value of buffer at that index
         */
        uint64_t & operator[](const uint64_t id)
        {
            return this->vec[id/64];
        }

        /**
         * Access value at given index in the underlying bitvector buffer.
         * For efficiency, no check is made that id is valid...
         * @param index Index in the vector buffer
         * @return Value of buffer at that index
         */

        uint64_t operator[](const uint64_t id) const
        {
            return this->vec[id/64];
        }

        /**
         * Copy constructor
         */
        BitVector operator=(const BitVector& other) {
            this->size = other.size;
            this->sz = other.sz;
            this->vec = (uint64_t *) malloc(sizeof(uint64_t) * sz);
            for(uint64_t i = 0; i < sz; i++){
                this->vec[i] = other.vec[i];
            }
            return *this;
        }

    };

} /* end namespace */


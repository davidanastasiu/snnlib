/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#pragma once

#include <cmath>
#include <cstddef>
#include <cstring>
#include <utility>
#include <iterator>
#include <stdexcept>
#include <limits>
#include <assert.h>
#include "sort.h"
#include <iostream>

/**
 * @namespace snnlib
 *
 * Key-Value Pair vector class with built-in sort routines and ability to build from
 * key & value arrays. Elements are stored as KEY and VALUE, so keys may be
 * accessed via key[i] and values via val[i].
 *
 * Note that, for efficiency, clear and resize do not clear elements in the vector.
 * One should use reserve + push_back for filling the vector sequentially, or
 * resize(n) + vec[i] for i < n. If possible, always pre-allocate the necessary
 * amount of space needed. Growth factor is 2x, starting with max_sz = 4.
 */


namespace snnlib {


    /**
     * Construct a KVVector
     * @tparam K    Key Type
     * @tparam V    Value Type
     */
    template<typename K, typename V>
    class KVVector {
    public:
        // types:
        typedef std::pair<K, V> T;
        typedef V & vreference;
        typedef K & kreference;
        typedef const V & const_vreference;
        typedef const K & const_kreference;
        typedef K * kiterator;
        typedef V * viterator;
        typedef const K * const_kiterator;
        typedef const V * const_viterator;
        typedef std::reverse_iterator<kiterator> reverse_kiterator;
        typedef std::reverse_iterator<viterator> reverse_viterator;
        typedef std::reverse_iterator<const_kiterator> const_reverse_kiterator;
        typedef std::reverse_iterator<const_viterator> const_reverse_viterator;
        typedef ptrdiff_t difference_type;
        typedef unsigned int size_type;

        // constructors/destructor:
        KVVector();

        explicit KVVector(const size_t n);

        KVVector(const size_t n, const T &keyval);

        KVVector(const size_t n, const V &val);

        KVVector(const size_t n, const K &key, const V &val);

        KVVector(const size_t n, const K *keys, const V *vals);

        KVVector(const size_t n, const V *vals);

        KVVector(
            typename KVVector<K, V>::kiterator& first,
            typename KVVector<K, V>::kiterator& last,
            const KVVector<K, V> &other
        );

        KVVector(std::initializer_list<T>);

        KVVector(std::initializer_list<V>);

        KVVector(const KVVector<K, V> &);

        KVVector(KVVector<K, V> &&);

        ~KVVector();

        KVVector<K, V> &operator=(const KVVector<K, V> &);

        KVVector<K, V> &operator=(KVVector<K, V> &&);

        KVVector<K, V> &operator=(std::initializer_list<T>);

        KVVector<K, V> &operator=(std::initializer_list<V>);

        void copy_unsafe(const KVVector<K,V> &other);

        /**
         * Assign
         * Rewrites original values and resizes if necessary
         */
        void assign(const size_t, const T &keyval);

        void assign(const size_t, const K &key, const V &val);

        void assign(const size_t, const K *keys, const V *vals);

        void assign(const size_t, const V *vals);

        void assign(const size_t, const V &val);

        void assign(
            typename KVVector<K, V>::kiterator,
            typename KVVector<K, V>::kiterator,
            const KVVector<K, V> &other
        );

        void assign(std::initializer_list<T>);

        void assign(std::initializer_list<V>);

        void set_key(const size_t idx, const K key);

        void set_val(const size_t idx, const V val);

        inline void set_key_unsafe(const size_t idx, const K key);

        inline void set_val_unsafe(const size_t idx, const V val);

        inline void set_size(const size_t size);

        /**
         * Iterators:
         * Create iterators pointing to array of data
         */
        kiterator k_begin();

        const_kiterator const_k_begin() const;

        viterator v_begin();

        const_viterator const_v_begin() const;

        kiterator k_end();

        const_kiterator const_k_end() const;

        viterator v_end();

        const_viterator const_v_end() const;

        reverse_kiterator k_rbegin();

        const_reverse_kiterator const_k_rbegin() const;

        reverse_viterator v_rbegin();

        const_reverse_viterator const_v_rbegin() const;

        reverse_kiterator k_rend();

        const_reverse_kiterator const_k_rend() const;

        reverse_viterator v_rend();

        const_reverse_viterator const_v_rend() const;

        // capacity:
        inline bool empty() const;

        inline size_t size() const;

        inline size_t max_size() const;

        inline size_t capacity() const;

        /**
         * Resize
         * Changes the size to an inputted value
         * Can decrement or increment
         * Will fill empty slots with inputted parameter
         */
        void resize(const size_t);

        void resize(const size_t, const T &);

        void resize(const size_t, const V &);

        void resize(const size_t, const K &, const V &);

        void resize_unsafe(const size_t sz);

        void reserve(const size_t);

        void shrink_to_fit();

        // element access
        inline vreference operator[](const size_t);

        inline const_vreference operator[](const size_t) const;

        inline vreference at(const size_t);

        inline const_vreference at(const size_t) const;

        inline kreference k_front();

        inline const_kreference k_front() const;

        inline vreference v_front();

        inline const_vreference v_front() const;

        inline kreference k_back();

        inline const_kreference k_back() const;

        inline vreference v_back();

        inline const_vreference v_back() const;

        inline K key(const size_t);

        inline K key(const size_t) const;

        inline K key_unsafe(const size_t);

		inline K key_unsafe(const size_t) const;

        inline V val(const size_t);

        inline V val(const size_t) const;

        inline V val_unsafe(const size_t);

		inline V val_unsafe(const size_t) const;

        // data access:
        K *k_data();

        V *v_data();

        const K *k_data() const;

        const V *v_data() const;

        // modifiers:
        template<class ... Args>
        void emplace_back(Args &&... args);

        /**
         * Push_back
         * Inserts a specified element at the back
         */
        void push_back(const T &);

        void push_back(const K &, const V &);

        void push_back(const V &);

        void push_back(T &&);

        void pop_back();

        inline void key(const size_t, const K &);

        inline void val(const size_t, const V &);

        inline void key_unsafe(const size_t, const K &);

        inline void val_unsafe(const size_t, const V &);

        template<class ... Args>
        kiterator emplace(const_kiterator, Args &&... args);

        /**
         * Insert
         * Inserts a specified element before a specified position
         * Return an iterator pointing to inserted element
         */
        kiterator insert(const_kiterator, const T &);

        kiterator insert(const_kiterator, const V &);

        kiterator insert(const_kiterator, T &&);

        kiterator insert(const_kiterator, V &&);

        kiterator insert(const_kiterator, const size_t, const T &);

        kiterator insert(const_kiterator, const size_t, const V &);

        template<class InputIt>
        kiterator insert(
            const_kiterator,
            InputIt,
            InputIt,
            const KVVector<K,V> &other
        );

        kiterator insert(const_kiterator, std::initializer_list<T>);

        kiterator insert(const_kiterator, std::initializer_list<V>);

        kiterator erase(const_kiterator);

        kiterator erase(const_kiterator, const_kiterator);

        unsigned int merge_keys(const char how);

        void scan(const char how);

        void swap(KVVector<K, V> &);

        void clear();

        // get a copy of the keys array
        K * keys()
        {
            K * keys = (K *) malloc( sizeof(K) * vec_sz);
            memcpy(keys, karr, sizeof(K) * vec_sz);
            return keys;
        }

        // statistics
        V min() const;
        V max() const;
        V sum() const;
        double mean() const;
        double stdev(const uint64_t ddof = 0) const;
        V median();

        // comparison
        bool operator==(const KVVector<K, V> &) const;

        bool operator!=(const KVVector<K, V> &) const;

        bool operator<(const KVVector<K, V> &) const;

        bool operator<=(const KVVector<K, V> &) const;

        bool operator>(const KVVector<K, V> &) const;

        bool operator>=(const KVVector<K, V> &) const;


        /**
         * Sort
         * Sorting in a specified order
         */

        void sort_vki();

        void sort_vki(const size_t sort_sz, const size_t start_index);

        void sort_vkd();

        void sort_vkd(const size_t sort_sz, const size_t start_index);

        void sort_kvi();

        void sort_kvi(const size_t sort_sz, const size_t start_index);

        void sort_kvd();

        void sort_kvd(const size_t sort_sz, const size_t start_index);

        void selectd (const size_t k);

        void selectd (const size_t start, const size_t end, const size_t k);

        void selecti (const size_t k);

        void selecti (const size_t start, const size_t end, const size_t k);

        /**
         * I/O
         */
        void printvec();

        friend std::ostream& operator<< (std::ostream& out, const KVVector<K,V>& v)
        {
            out << "KVVector(" << v.vec_sz << ")" << std::endl;
            for(int i = 0; i < v.vec_sz; i++){
                out << "("<< (K) v.karr[i] << ", " << (V) v.varr[i] << ") ";
            }
            out << std::endl;
            return out;
        }

    protected:
        size_type max_sz = 4;    /* max amount of space allocated */
        size_type vec_sz = 0;    /* current size of the vector */
        K *karr = nullptr;
        V *varr = nullptr;

        inline void reallocate();

    };


    /**
     * Default Constructor
     * No initialization
     * @tparam K    Key Type
     * @tparam V    Value Type
     */
    template<typename K, typename V>
    KVVector<K, V>::KVVector() {
        vec_sz = 0;
        max_sz = 0;
    }

    /**
     * Constructor with size
     * Sets the max size to n, does not initialize pairs
     * @tparam K    Key Type
     * @tparam V    Value Type
     * @param n     size of KVVector
     */
    template<typename K, typename V>
    KVVector<K, V>::KVVector(const size_t n) {
        vec_sz = 0;
        max_sz = n;
        if(n == 0){
            return;
        }
        reallocate();
    }

    /**
     * Constructor with size and pair
     * Initializes vector with KeyValue pair
     * @param n         size of KVVector
     * @tparam K        Key Type
     * @tparam V        Value Type
     * @param keyval    pair of key and value that must have same type as KVV
     */
    template<typename K, typename V>
    KVVector<K, V>::KVVector(const size_t n, const T &keyval) {
        vec_sz = max_sz = n;
        if(n == 0){
            return;
        }
        reallocate();
        for (size_t i = 0; i < n; ++i) {
            //assigns pair to each index
            karr[i] = keyval.first;
            varr[i] = keyval.second;
        }
    }

    /**
     * Constructor with size and value
     * Initializes vector with KeyValue pair
     * @param n         size of KVVector
     * @param Value     Value of pair that must have same type as KVV
     */
    template<typename K, typename V>
    KVVector<K, V>::KVVector(const size_t n, const V &value) {
        vec_sz = max_sz = n;
        if(n == 0){
            return;
        }
        reallocate();
        for (size_t i = 0; i < n; ++i) {
            karr[i] = (K) i;     // set keys to increment (assuming conversion is possible)
            varr[i] = value;
        }
    }

    /**
     * Constructor with size, key, and pair
     * Initializes vector with KeyValue pair made from key and value
     * @param Key         Key of pair with same type as KVV
     * @param Value     Value of pair
     */
    template<typename K, typename V>
    KVVector<K, V>::KVVector(const size_t n, const K &key, const V &value) {
        vec_sz = max_sz = n;
        if(n == 0){
            return;
        }
        reallocate();
        for (size_t i = 0; i < n; ++i) {
            //makes pair with key and value
            karr[i] = key;
            varr[i] = value;
        }
    }

    /**
     * Constructor with size, key array and value array
     * Initializes vector by reading parallel indexes of arrays and creating pairs.
     * @param n         size of KVVector and size of Arrays
     * @param vals        Array of Values
     * @param keys        Array of Keys
     */
    template<typename K, typename V>
    KVVector<K, V>::KVVector(const size_t n, const K *keys, const V *vals) {
        vec_sz = max_sz = n;
        if(n == 0){
            return;
        }
        reallocate();
        for (size_t i = 0; i < n; ++i) {
            karr[i] = keys[i];
            varr[i] = vals[i];
        }
    }

    /**
     * Constructor with size and value array
     * Initializes vector by reading parallel indexes of an array of values
     * And creating a key from the index of vector.
     * @param n         size of KVVector and size of Array
     * @param vals        Array of Values
     */
    template<typename K, typename V>
    KVVector<K, V>::KVVector(const size_t n, const V *vals) {
        vec_sz = max_sz = n;
        if(n == 0){
            return;
        }
        reallocate();
        for (size_t i = 0; i < n; ++i) {
            karr[i] = (K) i;
            varr[i] = vals[i];
        }
    }

    /**
     * Constructor with 2 iterators
     * Initializes vector by reading from a beginning iterator and end iterator
     * @param n         size of KVVector
     * @param first     iterator pointing to begin
     * @param last        iterator pointing to end
     */
    template<typename K, typename V>
    KVVector<K, V>::KVVector(KVVector<K, V>::kiterator& first,
                             KVVector<K, V>::kiterator& last,
                             const KVVector<K, V> &other) {
        size_t i, count = last - first, start = first-other.karr;
        max_sz = vec_sz = count;
        if(count == 0){
            return;
        }
        reallocate();

        viterator parallel_vit = &other.varr[start];
        for (i = 0; i < count; ++i, ++first, ++parallel_vit) {
            karr[i] = *first;
            varr[i] = *parallel_vit;
        }
    }

    /**
     * Constructor with initializer_list of pairs
     * Initializes vector by reading pairs from a list
     * @param lst<T>    List of pairs
     */
    template<typename K, typename V>
    KVVector<K, V>::KVVector(std::initializer_list<T> lst) {
        max_sz = vec_sz = lst.size();
        if(lst.size() == 0){
            return;
        }
        reallocate();
        vec_sz = 0;
        for (auto &item: lst) {
            karr[vec_sz] = item.first;
            varr[vec_sz++] = item.second;
        }
    }

    /**
     * Constructor with initializer_list of Values
     * Initialzies vector by reading values from a list
     * And creating keys from the index of Vector
     * @param lst<V>    List of Values of same type as KVV
     */
    template<typename K, typename V>
    KVVector<K, V>::KVVector(std::initializer_list<V> lst) {
        max_sz = vec_sz = lst.size();
        if(lst.size() == 0){
            return;
        }
        reallocate();
        vec_sz = 0;
        K key = (K) 0;
        for (auto &item: lst) {
            varr[vec_sz] = item;
            karr[vec_sz++] = key++;
        }
    }

    /**
     * Constructor with other KVVector
     * Initialzies vector by copying pairs from other KVV
     * Size will match that of the other vector
     * @param other<K,V>    KVV with same types
     */
    template<typename K, typename V>
    KVVector<K, V>::KVVector(const KVVector<K, V> &other) {
        vec_sz = other.vec_sz;
        if(vec_sz == 0){
            return;
        }
        max_sz = other.max_sz;
        assert(max_sz >= vec_sz);
        reallocate();
        for (size_t i = 0; i < other.vec_sz; ++i) {
            karr[i] = other.karr[i];
            varr[i] = other.varr[i];
        }
    }

    /**
     * Constructor with other KVVector
     * Initialzies vector by copying pairs from other KVV
     * @param &other<K,V>     KVV rvalue with same types
     */
    template<typename K, typename V>
    KVVector<K, V>::KVVector(KVVector<K, V> &&other) {
        vec_sz = other.vec_sz;
        if(vec_sz == 0){
            return;
        }
        max_sz = other.max_sz;
        assert(max_sz >= vec_sz);
        reallocate();
        for (size_t i = 0; i < other.vec_sz; ++i) {
            karr[i] = std::move(other.karr[i]);
            varr[i] = std::move(other.varr[i]);
        }
    }

    /**
     * Deconstructor
     */
    template<typename K, typename V>
    KVVector<K, V>::~KVVector() {
        if(karr){
            free(karr);
        }
        if(varr){
            free(varr);
        }
    }

    /**
     * @param other     KVVector of same type
     * = operator overload. Sets vector fields to other vector
     * By copying size (reallocating if necessary) and copying array
     */
    template<typename K, typename V>
    KVVector<K, V> &KVVector<K, V>::operator=(const KVVector<K, V> &other) {
        if (max_sz < other.vec_sz) {
            max_sz = other.vec_sz;
            reallocate();
        }
        for (size_t i = 0; i < other.vec_sz; ++i) {
            karr[i] = other.karr[i];
            varr[i] = other.varr[i];
        }
        vec_sz = other.vec_sz;
    }

    template<typename K, typename V>
    KVVector<K, V> &KVVector<K, V>::operator=(KVVector<K, V> &&other) {
        if (max_sz < other.vec_sz) {
            max_sz = other.vec_sz;
            reallocate();
        }
        for (size_t i = 0; i < other.vec_sz; ++i) {
            karr[i] = std::move(other.karr[i]);
            varr[i] = std::move(other.varr[i]);
        }
        vec_sz = other.vec_sz;
    }

    /**
     * @param lst     list of pairs
     * = operator overload
     * Copies size (reallocating if necessary) and copying pairs from list
     */
    template<typename K, typename V>
    KVVector<K, V> &KVVector<K, V>::operator=(std::initializer_list<T> lst) {
        if (max_sz < lst.size()) {
            max_sz = lst.size();
            reallocate();
        }
        vec_sz = 0;
        for (auto &item: lst) {
            karr[vec_sz] = item.first;
            varr[vec_sz++] = item.second;
        }
    }

    /**
     * @param lst     list of values
     * = operator overload
     * Copies size (reallocating if necessary) and copying values from list
     * Keys generated from index
     */
    template<typename K, typename V>
    KVVector<K, V> &KVVector<K, V>::operator=(std::initializer_list<V> lst) {
        if (max_sz < lst.size()) {
            max_sz = lst.size();
            reallocate();
        }
        vec_sz = 0;
        for (auto &item: lst) {
            karr[vec_sz] = vec_sz;
            varr[vec_sz++] = item;
        }
    }

    /**
     * Copies other to current
     */
    template<typename K, typename V>
    void KVVector<K,V>::copy_unsafe(const KVVector<K, V> &other) {
        for (size_t i = 0; i < other.vec_sz; ++i) {
            karr[i] = other.karr[i];
            varr[i] = other.varr[i];
        }
        vec_sz = other.vec_sz;
    }

    /**
     * @param value     KeyValue pair
     * @param count     number of pairs to assign
     * Assigns new content to the vector, overwriting previous content
     */
    template<typename K, typename V>
    void KVVector<K, V>::assign(const size_t count, const T &value) {
        if (count > max_sz) {
            max_sz = count;
            reallocate();
        }
        for (size_t i = 0; i < count; ++i) {
            karr[i] = value.first;
            varr[i] = value.second;
        }
        vec_sz = count;
    }

    /**
     * @param key         Key of a keyvalue pair
     * @param value     Value of a keyvalue pair
     * @param count     number of pairs to assign
     * Assigns new content to the vector, overwriting previous content
     */
    template<typename K, typename V>
    void KVVector<K, V>::assign(const size_t count, const K &key, const V &value) {
        if (count > max_sz) {
            max_sz = count;
            reallocate();
        }
        for (size_t i = 0; i < count; ++i) {
            karr[i] = key;
            varr[i] = value;
        }
        vec_sz = count;
    }

    /**
     * Assign with key and value arrays, copying count elements
     * Assigns new content to the vector, overwriting previous content
     */
    template<typename K, typename V>
    void KVVector<K, V>::assign(const size_t count, const K *keys, const V *vals) {
        if (count > max_sz) {
            max_sz = count;
            reallocate();
        }
        for (size_t i = 0; i < count; ++i) {
            karr[i] = keys[i];
            varr[i] = vals[i];
        }
        vec_sz = count;
    }

    /**
     * Assign with value array, copying count elements
     * Keys are the index of array
     * Assigns new content to the vector, overwriting previous content
     */
    template<typename K, typename V>
    void KVVector<K, V>::assign(const size_t count, const V *vals) {
        if (count > max_sz) {
            max_sz = count;
            reallocate();
        }
        for (size_t i = 0; i < count; ++i) {
            karr[i] = (K) i;
            varr[i] = vals[i];
        }
        vec_sz = count;
    }

    /**
     * Assign with value , copying count elements
     * The key is copied from the previous last element, unless there are no elements
     * Assigns new content to the vector, overwriting previous content
     */
    template<typename K, typename V>
    void KVVector<K, V>::assign(const size_t count, const V &value) {
        if (count > max_sz) {
            max_sz = count;
            reallocate();
        }
        K key = vec_sz > 0 ? (K) karr[vec_sz - 1] : 0;
        for (size_t i = 0; i < count; ++i) {
            karr[i] = key++;
            varr[i] = value;
        }
        vec_sz = count;
    }

    /**
     * Assign with two iterators
     * Assigns new content to the vector, overwriting previous content
     */
    template<typename K, typename V>
    void KVVector<K, V>::assign(KVVector<K, V>::kiterator first,
                                KVVector<K, V>::kiterator last, const KVVector<K, V> &other) {
        size_t count = last - first, start = first - other.karr;
        if (count > max_sz) {
            max_sz = count;
            reallocate();
        }
        viterator parallel_vit = &other.varr[start];
        for (size_t i = 0; i < count; ++i, ++first, ++parallel_vit) {
            karr[i] = *first;
            varr[i] = *parallel_vit;
        }
        vec_sz = count;
    }

    /**
     * Assign using initializer list of pairs
     * Assigns new content to the vector, overwriting previous content
     */
    template<typename K, typename V>
    void KVVector<K, V>::assign(std::initializer_list<T> lst) {
        size_t i, count = lst.size();
        if (count > max_sz) {
            max_sz = count;
            reallocate();
        }
        i = 0;
        for (auto &item: lst) {
            karr[i] = item.first;
            varr[i++] = item.second;
        }
        vec_sz = count;
    }

    /**
     * Assign with initializer list of values
     * Keys are the index
     * Assigns new content to the vector, overwriting previous content
     */
    template<typename K, typename V>
    void KVVector<K, V>::assign(std::initializer_list<V> lst) {
        size_t i, count = lst.size();
        if (count > max_sz) {
            max_sz = count;
            reallocate();
        }
        i = 0;
        for (auto &item: lst) {
            karr[i] = i;
            varr[i++] = item;
        }
        vec_sz = count;
    }

    template<typename K, typename V>
    void KVVector<K, V>::set_key(const size_t idx, const K key) {
        if(idx >= max_sz){
            max_sz = max_sz > 1 ? max_sz << 1 : 4;
            if(max_sz <= idx){
                max_sz = idx+1;
            }
            reallocate();
        }
        karr[idx] = key;
        if(idx >= vec_sz){
            vec_sz = idx+1;
        }
    }

    template<typename K, typename V>
    void KVVector<K, V>::set_val(const size_t idx, const V val) {
        if(idx >= max_sz){
            max_sz = max_sz > 1 ? max_sz << 1 : 4;
            if(max_sz <= idx){
                max_sz = idx+1;
            }
            reallocate();
        }
        varr[idx] = val;
        if(idx >= vec_sz){
            vec_sz = idx+1;
        }
    }

    template<typename K, typename V>
    inline void KVVector<K, V>::set_key_unsafe(const size_t idx, const K key) {
        karr[idx] = key;
    }

    template<typename K, typename V>
    inline void KVVector<K, V>::set_val_unsafe(const size_t idx, const V val) {
        varr[idx] = val;
    }

    template<typename K, typename V>
    inline void KVVector<K, V>::set_size(const size_t size) {
        vec_sz = size;
    }
    /**
     * Return the key at a certain index
     * @param idx     index of pair in vector
     * @return K    Key of pair
     */
    template<typename K, typename V>
    inline K KVVector<K, V>::key(const size_t idx) {
        if (idx < vec_sz){
            return karr[idx];
        }
        throw std::out_of_range("accessed position is out of range");
    }

    template<typename K, typename V>
    inline K KVVector<K, V>::key(const size_t idx) const {
        if (idx < vec_sz){
            return karr[idx];
        }
        throw std::out_of_range("accessed position is out of range");
    }

    /**
     * Returns key at idx without checking bounds
     */
	template<typename K, typename V>
	inline K KVVector<K, V>::key_unsafe(const size_t idx) {
		return karr[idx];
	}

	template<typename K, typename V>
    inline K KVVector<K, V>::key_unsafe(const size_t idx) const{
        return karr[idx];
    }

    /**
     * Return the Value at a certain index
     * @param idx     index of pair in vector
     * @return V    Value of pair
     */
    template<typename K, typename V>
    inline V KVVector<K, V>::val(const size_t idx) {
        if (idx < vec_sz){
            return varr[idx];
        }
        throw std::out_of_range("Accessed position is out of range");
    }

    template<typename K, typename V>
    inline V KVVector<K, V>::val(const size_t idx) const {
        if (idx < vec_sz){
            return varr[idx];
        }
        throw std::out_of_range("Accessed position is out of range");
    }

    /**
     * Returns val without checking if idx is less than size
     */
	template<typename K, typename V>
	inline V KVVector<K, V>::val_unsafe(const size_t idx) {
		return varr[idx];
	}

	template<typename K, typename V>
    inline V KVVector<K, V>::val_unsafe(const size_t idx) const{
        return varr[idx];
    }


    /**
     * Sets the key at an index to an inputted value
     * @param idx     index of pair in vector
     * @param key     Key must match original type
     */
    template<typename K, typename V>
    void KVVector<K, V>::key(const size_t idx, const K &key) {
        if (idx < vec_sz){
            karr[idx] = key;
        }
        throw std::out_of_range("Accessed position is out of range");
    }

    /**
     * Sets the value at an index to an inputted value
     * @param idx     index of pair in vector
     * @param val     value must match original type
     */
    template<typename K, typename V>
    void KVVector<K, V>::val(const size_t idx, const V &val) {
        if (idx < vec_sz){
            varr[idx] = val;
        }
        throw std::out_of_range("Accessed position is out of range");
    }

    /**
     * Return iterator pointing to beginning of key vector
     * @return kiterator
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::kiterator KVVector<K, V>::k_begin() {
        return karr;
    }

    template<typename K, typename V>
    inline typename KVVector<K, V>::const_kiterator KVVector<K, V>::const_k_begin() const {
        return karr;
    }

    /**
     * Return iterator pointing to beginning of value vector
     * @return viterator
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::viterator KVVector<K, V>::v_begin() {
        return varr;
    }

    template<typename K, typename V>
    inline typename KVVector<K, V>::const_viterator KVVector<K, V>::const_v_begin() const {
        return varr;
    }

    /**
     *    Return iterator pointing to the end of key vector
     * @return kiterator
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::kiterator KVVector<K, V>::k_end() {
        return karr + vec_sz;
    }

    template<typename K, typename V>
    inline typename KVVector<K, V>::const_kiterator KVVector<K, V>::const_k_end() const {
        return karr + vec_sz;
    }

    /**
     *    Return iterator pointing to the end of val vector
     * @return viterator
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::viterator KVVector<K, V>::v_end() {
        return varr + vec_sz;
    }

    template<typename K, typename V>
    inline typename KVVector<K, V>::const_viterator KVVector<K, V>::const_v_end() const {
        return varr + vec_sz;
    }

    /**
     *    Return a reverse iterator pointing to the end of key vector
     * @return reverse_kiterator
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::reverse_kiterator KVVector<K, V>::k_rbegin() {
        return reverse_kiterator(karr + vec_sz);
    }

    template<typename K, typename V>
    inline typename KVVector<K, V>::const_reverse_kiterator KVVector<K, V>::const_k_rbegin() const {
        return reverse_kiterator(karr + vec_sz);
    }

    /**
     *    Return a reverse iterator pointing to the end of value vector
     * @return reverse_viterator
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::reverse_viterator KVVector<K, V>::v_rbegin() {
        return reverse_viterator(varr + vec_sz);
    }

    template<typename K, typename V>
    inline typename KVVector<K, V>::const_reverse_viterator KVVector<K, V>::const_v_rbegin() const {
        return reverse_viterator(varr + vec_sz);
    }

    /**
     * Return a reverse iterator pointing to the beginning of key vector
     * @return reverse_kiterator
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::reverse_kiterator KVVector<K, V>::k_rend() {
        return reverse_kiterator(karr);
    }

    template<typename K, typename V>
    inline typename KVVector<K, V>::const_reverse_kiterator KVVector<K, V>::const_k_rend() const {
        return reverse_kiterator(karr);
    }

    /**
     * Return a reverse iterator pointing to the beginning of val vector
     * @return reverse_viterator
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::reverse_viterator KVVector<K, V>::v_rend() {
        return reverse_viterator(varr);
    }

    template<typename K, typename V>
    inline typename KVVector<K, V>::const_reverse_viterator KVVector<K, V>::const_v_rend() const {
        return reverse_viterator(varr);
    }

    /**
     * Adjusts the max_size of a vector by
     * Creating a new vector of max_size
     * moving elements from arr to tarr
     * Deleting arr
     * Copy tarr to arr
     */
    template<typename K, typename V>
    inline void KVVector<K, V>::reallocate() {
        if(karr){
            karr = (K *) realloc( karr, sizeof(K) * max_sz);
        } else {
            karr = (K *) malloc( sizeof(K) * max_sz);
        }
        if(varr){
            varr = (V *) realloc( varr, sizeof(V) * max_sz);
        } else {
            varr = (V *) malloc( sizeof(V) * max_sz);
        }
    }

    /**
     * Checks if the array is empty
     * @return true if sz==0
     */
    template<typename K, typename V>
    inline bool KVVector<K, V>::empty() const {
        return vec_sz == 0;
    }

    /**
     * @return    size_t
     * size of vector
     */
    template<typename K, typename V>
    inline size_t KVVector<K, V>::size() const {
        return vec_sz;
    }

    /**
     * @return    size_t
     * max amount of pairs
     */
    template<typename K, typename V>
    inline size_t KVVector<K, V>::max_size() const {
        // not realistic... maybe we should ask OS for available RAM?
        return std::numeric_limits<size_t>::max() / sizeof(std::pair<K, V>);
    }

    /**
     * @return std::size
     *        max size of vector
     */
    template<typename K, typename V>
    inline size_t KVVector<K, V>::capacity() const {
        return max_sz;
    }

    /**
     * Resize Vector to an inputted size
     * If input is greater than max, set max and current size to input, filling gaps with default pairs
     * If input is less than current, set max and current size to input, truncating values
     * If input is less than max but greater than current, set current size to input, filling gaps with default pairs
     * @param sz     new size of vector
     */
    template<typename K, typename V>
    void KVVector<K, V>::resize(const size_t sz) {
        if (sz > vec_sz) {
            if (sz > max_sz) {
                max_sz = sz;
                reallocate();
            }
            vec_sz = sz;
        } else if (sz < vec_sz) {
            vec_sz = sz;
            max_sz = sz;
            reallocate();
        }
    }

    /**
     * Resize Vector to an inputted size
     * @param p Pair of Key and Value
     * Fills extra space with p
     */
    template<typename K, typename V>
    void KVVector<K, V>::resize(const size_t sz, const T &p) {
        if (sz > vec_sz) {
            if (sz > max_sz) {
                max_sz = sz;
                reallocate();
            }
            size_t i;
            for (i = vec_sz; i < sz; ++i) {
                karr[i] = p.first;
                varr[i] = p.second;
            }
            vec_sz = sz;
        } else if (sz < vec_sz) {
            vec_sz = sz;
            max_sz = sz;
            reallocate();
        }
    }

    /**
     * Resize Vector to an inputted size
     * @param V         Value
     * Fills extra space with pairs of
     * inputted value and
     * a key equal to the index
    */
    template<typename K, typename V>
    void KVVector<K, V>::resize(const size_t sz, const V &v) {
        if (sz > vec_sz) {
            if (sz > max_sz) {
                max_sz = sz;
                reallocate();
            }
            size_t i;
            for (i = vec_sz; i < sz; ++i) {
                karr[i] = (K) i;
                varr[i] = v;
            }
            vec_sz = sz;
        } else if (sz < vec_sz) {
            vec_sz = sz;
            max_sz = sz;
            reallocate();
        }
    }

    /**
     * Resize Vector to an inputted size
     * @param V         Value
     * @param K         Key
     * Fills extra space with pairs of
     * inputted value and inputted key
     */
    template<typename K, typename V>
    void KVVector<K, V>::resize(const size_t sz, const K &k, const V &v) {
        if (sz > vec_sz) {
            if (sz > max_sz) {
                max_sz = sz;
                reallocate();
            }
            size_t i;
            for (i = vec_sz; i < sz; ++i) {
                karr[i] = k;
                varr[i] = v;
            }
            vec_sz = sz;
        } else if (sz < vec_sz) {
            vec_sz = sz;
            max_sz = sz;
            reallocate();
        }
    }

    template<typename K, typename V>
    void KVVector<K,V>::resize_unsafe(const size_t sz) {
        this->vec_sz = sz;
    }

    /**
     * Increases the max size
     * @param sz         new max size
     */
    template<typename K, typename V>
    void KVVector<K, V>::reserve(const size_t sz) {
        if (sz > max_sz) {
            max_sz = sz;
            reallocate();
        }
    }

    /**
     * Reduces max size to the vectors current size
     */
    template<typename K, typename V>
    void KVVector<K, V>::shrink_to_fit() {
        max_sz = vec_sz;
        reallocate();
    }

    /**
     * Overload [] operator
     * @param idx
     * @return reference to a value at idx
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::vreference KVVector<K, V>::operator[](const size_t idx) {
        if (idx < vec_sz){
            return varr[idx];
        }
        throw std::out_of_range("accessed position is out of range");
    }

    /**
     * Const Operator Overload []
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::const_vreference KVVector<K, V>::operator[](const size_t idx) const {
        if (idx < vec_sz){
            return varr[idx];
        }
        throw std::out_of_range("accessed position is out of range");
    }

    /**
     * Returns reference to position in array
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::vreference KVVector<K, V>::at(const size_t pos) {
        if (pos < vec_sz) {
            return varr[pos];
        }
        throw std::out_of_range("accessed position is out of range");
    }

    template<typename K, typename V>
    inline typename KVVector<K, V>::const_vreference KVVector<K, V>::at(const size_t pos) const {
        if (pos < vec_sz) {
            return varr[pos];
        }
        throw std::out_of_range("accessed position is out of range");
    }

    /**
     * Returns reference to first element's key
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::kreference KVVector<K, V>::k_front() {
        return karr[0];
    }

    template<typename K, typename V>
    inline typename KVVector<K, V>::const_kreference KVVector<K, V>::k_front() const {
        return karr[0];
    }
    /**
     * Returns reference to first element's value
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::vreference KVVector<K, V>::v_front() {
        return varr[0];
    }

    template<typename K, typename V>
    inline typename KVVector<K, V>::const_vreference KVVector<K, V>::v_front() const {
        return varr[0];
    }

    /**
     * Returns reference to last element's key
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::kreference KVVector<K, V>::k_back() {
        return karr[vec_sz - 1];
    }

    template<typename K, typename V>
    inline typename KVVector<K, V>::const_kreference KVVector<K, V>::k_back() const {
        return karr[vec_sz - 1];
    }

    /**
     * Returns reference to last element's value
     */
    template<typename K, typename V>
    inline typename KVVector<K, V>::vreference KVVector<K, V>::v_back() {
        return varr[vec_sz - 1];
    }

    template<typename K, typename V>
    inline typename KVVector<K, V>::const_vreference KVVector<K, V>::v_back() const {
        return varr[vec_sz - 1];
    }

    /**
     * Returns pointer to key array
     */
    template<typename K, typename V>
    inline K * KVVector<K, V>::k_data() {
        return karr;
    }

    template<typename K, typename V>
    inline const K * KVVector<K, V>::k_data() const {
        return karr;
    }

    /**
     * Returns pointer to value array
     */
    template<typename K, typename V>
    inline V * KVVector<K, V>::v_data() {
        return varr;
    }

    template<typename K, typename V>
    inline const V * KVVector<K, V>::v_data() const {
        return varr;
    }

    /**
     * Appends a new element to the end of Vector
     * Args     arguments to forward to the constructor of the element
     * T        the element type (Pair of <K,V>)
     */
    template<typename K, typename V>
    template<class ... Args>
    void KVVector<K, V>::emplace_back(Args &&... args) {
        if (vec_sz == max_sz) {
            max_sz = max_sz > 1 ? max_sz << 1 : 4;
            reallocate();
        }
        karr[vec_sz] = (K) vec_sz;
        varr[vec_sz] = std::move((V) V(std::forward<Args>(args) ...));
        ++vec_sz;
    }

    /**x
     * Inserts an element at the end
     * @param keyval        Pair of Key and Value
     *
     */
    template<typename K, typename V>
    void KVVector<K, V>::push_back(const T &keyval) {
        if (vec_sz == max_sz) {
            max_sz = max_sz > 1 ? max_sz << 1 : 4;
            reallocate();
        }
        karr[vec_sz] = keyval.first;
        varr[vec_sz] = keyval.second;
        ++vec_sz;
    }

    /**
     * @param key
     * @param val
     * Creates a pair of key and val and inserts at end
     */
    template<typename K, typename V>
    void KVVector<K, V>::push_back(const K &key, const V &val) {
        //double size of arr
        if (vec_sz == max_sz) {
            max_sz = max_sz > 1 ? max_sz << 1 : 4;
            reallocate();
        }
        karr[vec_sz] = key;
        varr[vec_sz] = val;
        ++vec_sz;
    }

    /**
     * @param val
     * Creates a pair of
     * Inputted value and last element's key + 1
     * And inserts pair at end
     */
    template<typename K, typename V>
    void KVVector<K, V>::push_back(const V &val) {
        if (vec_sz == max_sz) {
            max_sz = max_sz > 1 ? max_sz << 1 : 4;
            reallocate();
        }
        K key = vec_sz > 0 ? (K) karr[vec_sz - 1] + 1 : 0;
        karr[vec_sz] = key;
        varr[vec_sz] = val;
        ++vec_sz;
    }

    /**
     * Inserts an rvalue reference pair to back
     * @param keyval
     */
    template<typename K, typename V>
    void KVVector<K, V>::push_back(T &&keyval) {
        if (vec_sz == max_sz) {
            max_sz = max_sz > 1 ? max_sz << 1 : 4;
            reallocate();
        }
        karr[vec_sz] = std::move(keyval.first);
        varr[vec_sz] = std::move(keyval.second);
        ++vec_sz;
    }

    /**
     * Deletes last element
     * For efficiency, only reduces the size.
     */
    template<typename K, typename V>
    void KVVector<K, V>::pop_back() {
        if (vec_sz > 0) {
            --vec_sz;
        }
        else {
            throw "No elements left to pop.";
        }
    }

    /**
     * Emplace an element at a specified point
     * @param it            Position to insert
     * @param args            Args to forward to constructor
     * @return iterator     iterator pointing to emplaced element
     */
    template<typename K, typename V>
    template<class ... Args>
    typename KVVector<K, V>::kiterator KVVector<K, V>::emplace(
            KVVector<K, V>::const_kiterator it, Args &&... args) {
        size_t start = (it - karr);

        if (start < 0 || start > vec_sz){
            throw ("KVVector: Emplacing out of bounds");
        }

        if (vec_sz == max_sz) {
            max_sz = max_sz > 1 ? max_sz << 1 : 4;
            reallocate();
        }
        kiterator kiit = &karr[start];
        viterator viit = &varr[start];
        memmove(kiit + 1, kiit, (vec_sz - start) * sizeof(K));
        memmove(viit + 1, viit, (vec_sz - start) * sizeof(V));
        (*kiit) = vec_sz;
        (*viit) = std::move(V(std::forward<Args>(args) ...));
        ++vec_sz;
        return kiit;
    }

    /**
     * Insert an element before a specified point
     * Checks to make sure iterator is within range. Shifts the elements if necessary.
     * @param it            iterator pointing to element to insert behind
     * @param keyval        pair of Key and Value
     * @return iterator     iterator pointing to element inserted
     *
     */
    template<typename K, typename V>
    typename KVVector<K, V>::kiterator KVVector<K, V>::insert(
            KVVector<K, V>::const_kiterator it, const T &keyval) {
        size_t start = (it - karr);
        if (start < 0 || start > vec_sz) {
            throw "KVVector: Inserting out of bounds";
        }
        kiterator kiit = &karr[start];
        viterator viit = &varr[start];
        if (vec_sz == max_sz) {
            max_sz = max_sz > 1 ? max_sz << 1 : 4;
            reallocate();
        }
        if (vec_sz - start) {
            memmove(kiit + 1, kiit, (vec_sz - start) * sizeof(K));
            memmove(viit + 1, viit, (vec_sz - start) * sizeof(V));
        }
        (*kiit) = keyval.first;
        (*viit) = keyval.second;
        ++vec_sz;
        return kiit;
    }

    /**
     * Insert an element before a specified point
     * Checks to make sure iterator is within range. Shifts the elements if necessary.
     * Creates a key from the index of element to insert behind
     * @param val             value to create a pair
     *
     */
    template<typename K, typename V>
    typename KVVector<K, V>::kiterator KVVector<K, V>::insert(
            KVVector<K, V>::const_kiterator it, const V &val)
    {
        size_t start = (it - karr);
        if((it-karr) < 0 || (it-karr) > vec_sz){
            throw "KVVector: Inserting out of bounds";
        }

        if (vec_sz == max_sz) {
            max_sz = max_sz > 1 ? max_sz << 1 : 4;
            reallocate();
        }
        kiterator kiit = &karr[start];
        viterator viit = &varr[start];
        K key = (vec_sz - start) ? (K) (*kiit) : *(kiit-1)+1;
        if(vec_sz - start){
            memmove(viit + 1, viit, (vec_sz - start) * sizeof(V));
            memmove(kiit + 1, kiit, (vec_sz - start) * sizeof(K));
        }
        (*kiit) = key;
        (*viit) = val;
        ++vec_sz;
        kiterator traverse = (kiit+1);
        key++;
        while(traverse != &karr[vec_sz]){
            (*traverse) = key++;
            traverse++;
        }
        return kiit;
    }

    /**
     * Insert with rvalue pair
     */
    template <typename K, typename V>
    typename KVVector<K, V>::kiterator KVVector<K, V>::insert(
            KVVector<K, V>::const_kiterator it, T &&keyval)
    {
        size_t start = (it - karr);
        if((it-karr) < 0 || (it-karr) > vec_sz){
            throw "KVVector: Inserting out of bounds";
        }

        if (vec_sz == max_sz) {
            max_sz = max_sz > 1 ? max_sz << 1 : 4;
            reallocate();
        }
        kiterator kiit = &karr[start];
        viterator viit = &varr[start];

        if(vec_sz - start){
            memmove(viit + 1, viit, (vec_sz - start) * sizeof(V));
            memmove(kiit + 1, kiit, (vec_sz - start) * sizeof(K));
        }
        (*kiit) = std::move(keyval.first);
        (*viit) = std::move(keyval.second);
        ++vec_sz;
        return kiit;
    }

    /**
     * Insert with rvalue value
     */
    template <typename K, typename V>
    typename KVVector<K, V>::kiterator KVVector<K, V>::insert(
            KVVector<K, V>::const_kiterator it, V &&val)
    {
        size_t start = (it - karr);
        if((it-karr) < 0 || (it-karr) > vec_sz){
            throw "KVVector: Inserting out of bounds";
        }

        if (vec_sz == max_sz) {
            max_sz = max_sz > 1 ? max_sz << 1 : 4;
            reallocate();
        }

        kiterator kiit = &karr[start];
        viterator viit = &varr[start];
        K key = (vec_sz -start)? (K) (*kiit) : *(kiit-1) + 1;
        if(vec_sz - start){
            memmove(viit + 1, viit, (vec_sz - start) * sizeof(V));
            memmove(kiit + 1, kiit, (vec_sz - start) * sizeof(K));
        }
        (*kiit) = key;
        (*viit) = std::move(val);
        ++vec_sz;
        kiterator traverse = (kiit+1);
        key++;
        while(traverse != &karr[vec_sz]){
            (*traverse) = key++;
            traverse++;
        }
        return kiit;
    }

    /**
     * Insert cnt elements before a specified point
     * Checks to make sure iterator is within range
     * Shifts the elements if necessary
     * @param it            key iterator pointing to element to insert behind
     * @param keyval        pair of Key and Value
     * @param cnt             number of elements to insert
     * @return kiterator     iterator pointing to first element inserted
     */
    template <typename K, typename V>
    typename KVVector<K, V>::kiterator KVVector<K, V>::insert(
            KVVector<K, V>::const_kiterator it, size_t cnt, const T &keyval)
    {

        size_t start = (it - karr);
        if((it-karr) < 0 || (it-karr) > vec_sz){
            throw "KVVector: Inserting out of bounds";
        }

        if (vec_sz + cnt > max_sz) {
            max_sz = (vec_sz + cnt) << 1;
            reallocate();
        }
        kiterator kiit = &karr[start];
        viterator viit = &varr[start];
        if (!cnt) return kiit;

        if(vec_sz - start){
            memmove(kiit + cnt, kiit, (vec_sz - start) * sizeof(K));
            memmove(viit + cnt, viit, (vec_sz - start) * sizeof(V));
        }
        vec_sz += cnt;
        kiterator kit;
        viterator vit = viit;
        for (kiterator kit = kiit; cnt--; ++kit, ++vit){
            (*kit) = keyval.first;
            (*vit) = keyval.second;
        }
        return kiit;
    }

    /**
     * Insert cnt elements before a specified point
     * Checks to make sure iterator is within range
     * Shifts the elements if necessary
     * Keys are created from the key of the element being inserted behind
     * and are incremented by 1
     * @param it            key iterator pointing to element to insert behind
     * @param val             Value
     * @param cnt             number of elements to insert
     * @return kiterator    key iterator pointing to first element inserted
     */
    template <typename K, typename V>
    typename KVVector<K, V>::kiterator KVVector<K, V>::insert(
            KVVector<K, V>::const_kiterator it, size_t cnt, const V &val)
    {
        if((it-karr) < 0 || (it-karr) > vec_sz){
            throw "KVVector: Inserting out of bounds";
        }
        size_t start = (it - karr);

        if (vec_sz + cnt > max_sz) {
            max_sz = (vec_sz + cnt) << 1;
            reallocate();
        }

        kiterator kiit = &karr[start];
        viterator viit = &varr[start];
        if (!cnt) return kiit;
        K key = (vec_sz -start) ? (K) (*kiit) : *(kiit-1)+1;
        if(vec_sz - start){
            memmove(kiit + cnt, kiit, (vec_sz - start) * sizeof(K));
            memmove(viit + cnt, viit, (vec_sz - start) * sizeof(V));
        }
        vec_sz += cnt;
        kiterator ktraverse = kiit;
        for(viterator vtraverse = viit ; cnt--;ktraverse++, vtraverse++){
            (*ktraverse) = key++;
            (*vtraverse)= val;
        }
        while(ktraverse != &karr[vec_sz]){
            (*ktraverse) = key++;
            ktraverse++;
        }
        return kiit;
    }

    /**
     * Insert elements between inclusive first and non inclusive last
     * @tparam InputIt    inputIterator of type kiterator
     * @param it        kiterator pointing to the element to insert behind
     * @param first     kiterator pointing to the first element of items to insert
     * @param last        kiterator pointing to the last element of items to insert
     * @return kiterator pointing to first inserted element
     */
    template <typename K, typename V>
    template <class InputIt>
    typename KVVector<K, V>::kiterator KVVector<K, V>::insert(
            KVVector<K, V>::const_kiterator it, InputIt first, InputIt last, const KVVector<K, V> &other)
    {
        size_t cnt = last - first;
        size_t start = it - karr;
        size_t parallel_ind = first - other.karr;
        if((it-karr) < 0 || (it-karr) > vec_sz || (first - last) > 0 ){
            throw "KVVector: Inserting out of bounds";

        }

        if (vec_sz + cnt > max_sz) {
            max_sz = (vec_sz + cnt) << 1;
            reallocate();
        }

        kiterator kiit = &karr[start];
        viterator viit = &karr[start];
        if(!cnt){
            return kiit;
        }

        if (vec_sz - start > 0){
            memmove(kiit + cnt, kiit, (vec_sz - start) * sizeof(K));
            memmove(viit + cnt, viit, (vec_sz - start) * sizeof(V));
        }

        viterator parallel_vit = &other.varr[parallel_ind];
        for (size_t i=start; first != last; ++i, ++first, ++parallel_vit){
            karr[i] = (*first);
            varr[i] = (*parallel_vit);
        }
        vec_sz += cnt;
        return kiit;
    }

    /**
     * Insert elements from a list of pairs
     * @param lst         list of KeyValuePairs
     */
    template <typename K, typename V>
    typename KVVector<K, V>::kiterator KVVector<K, V>::insert(
            KVVector<K, V>::const_kiterator it, std::initializer_list<T> lst)
    {
        size_t cnt = lst.size();
        size_t start = it-karr;
        if (vec_sz + cnt > max_sz) {
            max_sz = (vec_sz + cnt) << 1;
            reallocate();
        }

        kiterator kiit = &karr[start];
        viterator viit = &varr[start];
        if (!cnt) return kiit;
        if(vec_sz - start){
            memmove(kiit + cnt, kiit, (vec_sz - start) * sizeof(K));
            memmove(viit + cnt, viit, (vec_sz - start) * sizeof(V));
        }
        kiterator kit = kiit;
        viterator vit = viit;
        for (auto &item: lst) {
            (*kit) = item.first;
            (*vit) = item.second;
            ++kit; ++vit;
        }
        vec_sz += cnt;
        return kiit;
    }

    /**
     * Insert elements from a list of values before a specified position
     * Keys are created from key of element to be inserted behind
     * and increment by 1
     */
    template <typename K, typename V>
    typename KVVector<K, V>::kiterator KVVector<K, V>::insert(
            KVVector<K, V>::const_kiterator it, std::initializer_list<V> lst)
    {
        size_t start = it-karr;
        size_t cnt = lst.size();

        if (vec_sz + cnt > max_sz) {
            max_sz = (vec_sz + cnt) << 1;
            reallocate();
        }
        kiterator kiit = &karr[start];
        viterator viit = &varr[start];
        if (!cnt) return viit;
        K key = *kiit;
        if(vec_sz - start){
            memmove(kiit + cnt, kiit, (vec_sz - start) * sizeof(K));
            memmove(viit + cnt, viit, (vec_sz - start) * sizeof(V));
        }
        kiterator kit = kiit;
        viterator vit = viit;
        for (auto &item: lst) {
            (*kit) = key++;
            (*vit) = item;
            ++kit; ++vit;
        }
        vec_sz += cnt;
        return kiit;
    }

    /**
     * Erase the element at a specified index
     * @param it            kiterator pointing to the element to be deleted
     * @return kiterator     iterator pointing to the deleted element
     */
    template <typename K, typename V>
    typename KVVector<K, V>::kiterator KVVector<K, V>::erase(
            KVVector<K, V>::const_kiterator it)
    {
        if((it-karr) < 0 || (it-karr) >= vec_sz){
            throw "KVVector: Erasing out of bounds";
        }
        kiterator kiit = &karr[it - karr];
        viterator viit = &varr[it - karr];
        if(it - karr < vec_sz - 1){
            memmove(kiit, kiit + 1, (vec_sz - (it - karr) - 1) * sizeof(K));
            memmove(viit, viit + 1, (vec_sz - (it - karr) - 1) * sizeof(V));
        }
        --vec_sz;
        return kiit;
    }

    /**
     * Erase the element between first iterator and last iterator (non-inclusive)
     * @param first         kiterator pointing to the first element in range to be deleted
     * @param last            kiterator pointing to the index after the last deleted element
     * @return iterator     kiterator pointing to first deleted element
     */
    template <typename K, typename V>
    typename KVVector<K, V>::kiterator KVVector<K, V>::erase(
            KVVector<K, V>::const_kiterator first, KVVector<K, V>::const_kiterator last)
    {
        kiterator kit = &karr[first - karr];
        viterator vit = &varr[first - karr];
        if (first == last) return kit;
        memmove(kit, last, (vec_sz - (last - karr)) * sizeof(K));
        memmove(vit, &varr[last-karr], (vec_sz - (last - karr)) * sizeof(V));
        vec_sz -= last - first;
        return kit;
    }

    /**
     * Merge values of elements with the same key
     * @param how         specifier on how to merge elements
     * Only consecutive pairs with the same key are combined into one pair
     *    how:
     *     +    adds values of all selected pairs together
     *     <    compares value of all selected pairs and sets value to least one
     *     >    compares value of all selected pairs and sets value to greatest one
     *     ^    compares value of all selected pairs and sets value to absolute max
     *     m    sets value to the mean of all the values of the selected pairs
     * @return            the number of the elements merged
     */
    template <typename K, typename V>
    unsigned int KVVector<K, V>::merge_keys(const char how) {
        kiterator kfast = (karr+1);
        kiterator kslow = karr;
        viterator vfast = (varr + 1);
        viterator vslow = varr;

        switch(how) {
            case '+':
            {
                for(int i = 0; i < vec_sz-1; i++){
                    if (*kfast == *kslow){
                        *vslow += *vfast;
                    } else {
                        kslow++;
                        vslow++;
                        *kslow = *kfast;
                        *vslow = *vfast;
                    }
                    kfast++;
                    vfast++;
                }
                vec_sz = (kslow - karr + 1);
            }
                break;
            case '<' :
            {
                for(int i = 0; i < vec_sz-1; i++){
                    if (*kfast == *kslow){
                        V min = *vfast < *vslow ? *vfast : *vslow;
                        *vslow = min;
                    } else {
                        kslow++;
                        vslow++;
                        *kslow = *kfast;
                        *vslow = *vfast;
                    }
                    vfast++;
                    kfast++;
                }
                vec_sz = (kslow - karr + 1);
            }
                break;
            case '>':
            {
                for(int i = 0; i < vec_sz-1; i++){
                    if (*kfast == *kslow){
                        V max = *vfast > *vslow ? *vfast : *vslow;
                        *vslow = max;
                    } else {
                        kslow++;
                        vslow++;
                        *kslow = *kfast;
                        *vslow = *vfast;
                    }
                    vfast++;
                    kfast++;
                }
                vec_sz = (kslow - karr + 1);
            }
                break;
            case '^':
            {
                for(int i = 0; i < vec_sz-1; i++){
                    if (*kfast == *kslow){
                        V absmax;
                        if(*vfast <= 0 && *vslow <= 0){
                            absmax = *vfast < *vslow ? *vfast : *vslow;
                        }
                        else if (*vfast <= 0 && *vslow > 0) {
                            absmax = (*vfast * -1) > *vslow ? *vfast : *vslow;
                        }
                        if (*vfast > 0 && *vslow > 0){
                            absmax = *vfast > *vslow ? *vfast : *vslow;
                        }
                        else if (*vfast > 0 && *vslow <=0) {
                            absmax = *vfast > (*vslow * -1) ? *vfast : *vslow;
                        }
                        *vslow = absmax;
                    } else {
                        kslow++;
                        vslow++;
                        *kslow = *kfast;
                        *vslow = *vfast;
                    }
                    kfast++;
                    vfast++;
                }
                vec_sz = (kslow - karr + 1);
            }
                break;
            case 'm':
            {
                int count = 1;
                for(int i = 0; i < vec_sz-1; i++){
                    if (*kfast == *kslow){
                        *vslow += *vfast;
                        count++;
                    } else {
                        *vslow /= count;
                        count = 1;
                        kslow++;
                        vslow++;
                        *kslow = *kfast;
                        *vslow = *vfast;
                    }
                    kfast++;
                    vfast++;
                }
                *vslow /= count;
                vec_sz = (kslow - karr + 1);
            }
                break;
            default:
                throw "KVVector: merge_keys operator invalid";
        }
        return vec_sz;
    }

    /**
    * Scan values traversing from left to right
    * @param how         specifier on how to scan elements
     * how:
     *     +    adds values of previous index to value of current index
     *     <    compares value of previous and current and sets current to the lesser one
     *     >    compares value of previous and current and sets current to the greater one
     *     ^    compares value of previous and current and sets current to the absolute max
     *     m    sets current to the mean of previous values
    */
    template <typename K, typename V>
    void KVVector<K, V>::scan(const char how){
        switch(how){
            case '+':
            {
                V sum = varr[0];
                for (int i = 1; i < vec_sz; i++) {
                    sum += varr[i];
                    varr[i] = sum;
                }
                break;
            }
            case '<' :
            {
                V min = varr[0];
                for (int i = 1; i < vec_sz; i++) {
                    if (varr[i] < min) {
                        min = varr[i];
                    }
                    varr[i] = min;
                }
                break;
            }
            case '>':
            {
                V max = varr[0];
                for(int i = 1; i < vec_sz; i++){
                    if(varr[i] > max){
                        max = varr[i];
                    }
                    varr[i] = max;
                }
                break;
            }
            case '^':
            {
                V absmax = varr[0];
                bool neg = false;
                if(absmax < 0){
                    neg = true;
                }
                for(int i = 1; i < vec_sz; i++){
                    if (varr[i] <= 0 && neg) {
                        if(varr[i] < absmax){
                            absmax = varr[i];
                        }
                        varr[i] = absmax;
                    }
                    else if (varr[i] <= 0 && !neg){
                        if( (varr[i] *- 1) > absmax){
                            neg = true;
                            absmax = varr[i];
                        }
                        varr[i] = absmax;
                    }
                    else if (varr[i] > 0 && neg){
                        if( (varr[i] *- 1) < absmax){
                            neg = false;
                            absmax = varr[i];
                        }
                        varr[i] = absmax;
                    }
                    else if (varr[i] > 0 && !neg){
                        if(varr[i] > absmax){
                            absmax = varr[i];
                        }
                        varr[i] = absmax;
                    }
                }

                break;
            }
            case 'm':
            {
                V mean = varr[0];
                for (int i = 1; i < vec_sz; i++) {
                    mean += varr[i];
                    varr[i] = mean / (i + 1);
                }
                break;
            }
            default:
                throw "KVVector: scan operator invalid";
        }
    }

    /**
     * Swap the size, max size, and data of two KVVectors
     * @param rhs             other KVVector of same type
     */
    template <typename K, typename V>
    void KVVector<K, V>::swap(KVVector<K, V> &rhs) {
        size_t tvec_sz = vec_sz;
        size_t tmax_sz = max_sz;

        K *tkarr = karr;
        V *tvarr = varr;

        vec_sz = rhs.vec_sz;
        max_sz = rhs.max_sz;
        karr = rhs.karr;
        varr = rhs.varr;

        rhs.vec_sz = tvec_sz;
        rhs.max_sz = tmax_sz;
        rhs.varr = tvarr;
        rhs.karr = tkarr;
    }

    /**
     * Sets the size of vector to 0, effectively deleting every element.
     * Internal memory assignment does not change.
     */
    template <typename K, typename V>
    void KVVector<K, V>::clear() {
        vec_sz = 0;
    }



    // statistics
    /**
     * Compute the minimum of the vector values.
     * @return min    of the vector values.
     */
    template <typename K, typename V>
    V KVVector<K, V>::min() const
    {
        if(this->vec_sz < 1){
            throw "KVVector: Empty vector does not have a min value.";
        }
        V min = this->varr[0];
        #pragma omp parallel for reduction(min: min)
        for(size_t i = 1; i < this->vec_sz; ++i){
            if(this->varr[i] < min){
                min = this->varr[i];
            }
        }
        return min;
    }

    /**
     * Compute the maximum of the vector values.
     * @return max    of the vector values.
     */
    template <typename K, typename V>
    V KVVector<K, V>::max() const
    {
        if(this->vec_sz < 1){
            throw "KVVector: Empty vector does not have a max value.";
        }
        V max = this->varr[0];
        #pragma omp parallel for reduction(max: max)
        for(size_t i = 1; i < this->vec_sz; ++i){
            if(this->varr[i] > max){
                max = this->varr[i];
            }
        }
        return max;
    }

    /**
     * Compute the sum of the vector values.
     * @return sum    of the vector values.
     */
    template <typename K, typename V>
    V KVVector<K, V>::sum() const
    {
        V sum = 0;
        #pragma omp parallel for reduction(+: sum) schedule(static, 128)
        for(size_t i = 0; i < this->vec_sz; ++i){
            sum += this->varr[i];
        }
        return sum;
    }

    /**
     * Compute the mean of the vector values.
     * @return mean    of the vector values.
     */
    template <typename K, typename V>
    double KVVector<K, V>::mean() const
    {
        if(this->vec_sz < 1){
            throw "KVVector: Empty vector does not have a mean.";
        }
        double sum = 0;
        #pragma omp parallel for reduction(+: sum) schedule(static, 128)
        for(size_t i = 0; i < this->vec_sz; ++i){
            sum += (double) this->varr[i];
        }
        return sum / this->vec_sz;
    }

    /**
     * Compute the standard deviation of the vector values.
     * @param ddof Delta Degrees of Freedom. The divisor used in calculations
     *             is N - ddof, where N represents the number of elements.
     *             By default ddof is zero.
     * @return standard deviation    of the vector values.
     */
    template <typename K, typename V>
    double KVVector<K, V>::stdev(const uint64_t ddof) const
    {
        if(this->vec_sz < 1){
            throw "KVVector: Empty vector does not have a standard deviation.";
        }
        double sum = 0;
        #pragma omp parallel for reduction(+: sum) schedule(static, 128)
        for(size_t i = 0; i < this->vec_sz; ++i){
            sum += (double) this->varr[i];
        }
        double mean = sum / this->vec_sz;
        double s2 = 0;
        #pragma omp parallel for shared(mean) reduction(+: s2) schedule(static, 128)
        for(size_t i = 0; i < this->vec_sz; ++i){
            s2 += (this->varr[i] - mean) * (this->varr[i] - mean);
        }

        return std::sqrt(s2 / (this->vec_sz - ddof));
    }

    /**
     * Computes the median of the vector values. To do so, the select algorithm
     * is applied moving size/2 (+ 1 if size is even) of the smallest values to
     * the left of the vector. Then the max value is linearly searched in that
     * half of the vector. If size is even, the top 2 max values are found and
     * averaged.
     * @return median of all values in the vector
     */
    template <typename K, typename V>
    V KVVector<K, V>::median()
    {
        // edge cases
        if(this->vec_sz == 0){
            return 0;
        } else if(this->vec_sz == 1){
            return this->varr[0];
        } else if(this->vec_sz == 2){
            return (this->varr[0] + this->varr[1])/2.0 ;
        }

        // move the lowest values to the first half of the vector
        size_t mid = this->vec_sz/2 + 1;
        this->selecti(mid);

        V val, max1;

        // single middle point if size is odd
        if(this->vec_sz % 2 == 1){
            max1 = this->varr[0];
            for(size_t i=1; i < mid; ++i){
                if(this->varr[i] > max1){
                    max1 = this->varr[i];
                }
            }
            return max1;
        }

        // 2 highest values: max1 is the low max, and max2 is the high max
        V max2;
        if(this->varr[0] < this->varr[1]){
            max1 = this->varr[0];
            max2 = this->varr[1];
        } else {
            max1 = this->varr[1];
            max2 = this->varr[0];
        }
        for(size_t i=1; i < mid; ++i){
            if(this->varr[i] > max2){
                max1 = max2;
                max2 = this->varr[i];
            } else if(this->varr[i] > max1){
                max1 = this->varr[i];
            }
        }
        return (max1 + max2) / 2.0;
    }

    /**
     * Comparison Operator == overload
     * Checks if size is equal
     * Checks if all pairs are equal
     */
    template <typename K, typename V>
    bool KVVector<K, V>::operator == (const KVVector<K, V> &rhs) const {
        if (vec_sz != rhs.vec_sz) return false;
        size_t i;
        for (i = 0; i < vec_sz; ++i){
            if (karr[i] != rhs.karr[i] || varr[i] != rhs.varr[i]){
                return false;
            }
        }
        return true;
    }

    /**
     * Comparison Operator != overload
     * Checks if size is not equal
     * Checks if all pairs are not equal
     */
    template <typename K, typename V>
    bool KVVector<K, V>::operator != (const KVVector<K, V> &rhs) const {
        if (vec_sz != rhs.vec_sz) return true;
        size_t i;
        for (i = 0; i < vec_sz; ++i){
            if (karr[i] != rhs.karr[i] || varr[i] != rhs.varr[i]){
                return true;
            }
        }
        return false;
    }

    /**
     * Overloads the < operator
     * Returns true if size is less than other size
     * If sizes are equal
     * Checks every element, checking key before value
     * Stopping at the first value
     * and returning true if current element is less than other element
     * Return false if everything is the same
     * @param rhs         other KVVector
     * @return bool
     */
    template <typename K, typename V>
    bool KVVector<K, V>::operator < (const KVVector<K, V> &rhs) const {
        if (vec_sz < rhs.vec_sz){
            return true;
        }
        else if (vec_sz > rhs.vec_sz){
            return false;
        }
        size_t i;
        for (i = 0; i < vec_sz; ++i){
            if (karr[i] != rhs.karr[i]){
                return karr[i] < rhs.karr[i];
            }
            else if(varr[i] != rhs.varr[i]) {
                return varr[i] < rhs.varr[i];
            }
        }
        return false;
    }

    /**
     * Overloads the <= operator
     * Returns true if size is less than other size
     * If sizes are equal
     * Checks every element, checking Key before value,
     * Stopping at the first value
     * and returning true if current element is less than other element
     * Returns true if everything is the same
     * @param rhs         other KVVector
     * @return bool
     */
    template <typename K, typename V>
    bool KVVector<K, V>::operator <= (const KVVector<K, V> &rhs) const {
        if (vec_sz < rhs.vec_sz){
            return true;
        }
        else if (vec_sz > rhs.vec_sz){
            return false;
        }
        size_t i;
        for (i = 0; i < vec_sz; ++i){
            if (karr[i] != rhs.karr[i]){
                return karr[i] < rhs.karr[i];
            }
            else if(varr[i] != rhs.varr[i]) {
                return varr[i] < rhs.varr[i];
            }
        }
        return true;
    }

    /**
     * Overloads the > operator
     * Returns true if size is greater than other size
     * If sizes are equal
     * Checks every element, checking Key before value,
     * Stopping at the first value
     * and returning true if current element is greater than other element
     * Returns false if everything is the same
     * @param rhs         other KVVector
     * @return bool
     */
    template <typename K, typename V>
    bool KVVector<K, V>::operator > (const KVVector<K, V> &rhs) const {
        if (vec_sz > rhs.vec_sz){
            return true;
        }
        else if (vec_sz < rhs.vec_sz){
            return false;
        }
        size_t i;
        for (i = 0; i < vec_sz; ++i){
            if (karr[i] != rhs.karr[i]){
                return karr[i] > rhs.karr[i];
            }
            else if(varr[i] != rhs.varr[i]) {
                return varr[i] > rhs.varr[i];
            }
        }
        return false;
    }

    /**
     * Overloads the >= operator
     * Returns true if size is greater than other size
     * If sizes are equal
     * Checks every element, checking Key before value,
     * Stopping at the first value
     * and returning true if current element is greater than other element
     * Returns true if everything is the same
     * @param rhs         other KVVector
     * @return bool
     */
    template <typename K, typename V>
    bool KVVector<K, V>::operator >= (const KVVector<K, V> &rhs) const {
        if (vec_sz > rhs.vec_sz){
            return true;
        }
        else if (vec_sz < rhs.vec_sz){
            return false;
        }
        size_t i;
        for (i = 0; i < vec_sz; ++i){
            if (karr[i] != rhs.karr[i]){
                return karr[i] > rhs.karr[i];
            }
            else if(varr[i] != rhs.varr[i]) {
                return varr[i] > rhs.varr[i];
            }
        }
        return true;
    }


    /**
     * Sort by increasing order, prioritizing the value
     */
    template <typename K, typename V>
    void KVVector<K, V>::sort_vki(){
        if (vec_sz <=1 ){
            return;
        }
        snnlib::kvsorti(vec_sz, karr, varr);
    }

    /**
     * Sort by increasing order, prioritizing the value
     * Allows for choosing how many elements to sort and the starting point
     */
    template<typename K, typename V>
    void KVVector<K, V>::sort_vki(const size_t sort_sz, const size_t start_index){
        if (sort_sz <= 1){
            return;
        }
        snnlib::kvsorti(sort_sz, karr + start_index, varr + start_index);
    }

    /**
     * Sort by decreasing order, prioritizing the value
     */
    template <typename K, typename V>
    void KVVector<K, V>::sort_vkd(){
        if (vec_sz <=1){
            return;
        }
         snnlib::kvsortd(vec_sz, karr, varr);
    }

    /**
     * Sort by decreasing order, prioritizing the value
     * Allows for choosing how many elements to sort and the starting point
     */
    template<typename K, typename V>
    void KVVector<K, V>::sort_vkd(const size_t sort_sz, const size_t start_index){
        if (sort_sz <= 1){
            return;
        }
        snnlib::kvsortd(sort_sz, karr + start_index, varr + start_index);
    }

    /**
     * Sort by increasing order, prioritizing the key
     */
    template <typename K, typename V>
    void KVVector<K, V>::sort_kvi(){
        if (vec_sz <= 1){
            return;
        }
        snnlib::kvsorti(vec_sz, varr, karr);
    }
    /**
     * Sort by increasing order, prioritizing the key
     * Allows for choosing how many elements to sort and the starting point
     */
    template<typename K, typename V>
    void KVVector<K, V>::sort_kvi(const size_t sort_sz, const size_t start_index){
        if (sort_sz <= 1){
            return;
        }
        snnlib::kvsorti(sort_sz, varr + start_index, karr + start_index);
    }

    /**
     * Sort by decreasing order, prioritizing the key
     */
    template <typename K, typename V>
    void KVVector<K, V>::sort_kvd(){
        if (vec_sz <= 1){
            return;
        }
        snnlib::kvsortd(vec_sz, varr, karr);
    }

    /**
     * Sort by decreasing order, prioritizing the key
     * Allows for choosing how many elements to sort and the starting point
     */
    template<typename K, typename V>
    void KVVector<K, V>::sort_kvd(const size_t sort_sz, const size_t start_index){
        if (sort_sz <= 1){
            return;
        }
        snnlib::kvsortd(sort_sz, varr + start_index, karr + start_index);
    }

    #define KVSWAP(a, b, stmp) do { stmp = (a); (a) = (b); (b) = stmp; } while (0)

    /**
     * This function puts the `k` largest values in the beginning of the array
     * @param k Number of smallest values to move to the start of the vector.
     */
    template <typename K, typename V>
    void KVVector<K,V>::selectd(const size_t k)
    {
        size_t i, j, lo, hi, mid;
        K ktmp;
        V vtmp;
        V pivot;

        if (this->vec_sz <= k)
            return;         /* return if the vector has fewer elements than we want */

        for (lo = 0, hi = this->vec_sz - 1; lo < hi;)
        {
            mid = lo + ((hi - lo) >> 1);

            /* select the median */
            if (varr[lo] < varr[mid]){
                mid = lo;
            }
            if (varr[hi] > varr[mid]){
                mid = hi;
            } else{
                goto jump_over;
            }
            if (varr[lo] < varr[mid]){
                mid = lo;
            }

        jump_over:
            //swap mid and high
            KVSWAP(varr[mid], varr[hi], vtmp);
            KVSWAP(karr[mid], karr[hi], ktmp);
            pivot = varr[hi];

            /* the partitioning algorithm */
            for (i = lo - 1, j = lo; j < hi; j++)
            {
                if (varr[j] >= pivot)
                {
                    i++;
                    KVSWAP(varr[i], varr[j], vtmp);
                    KVSWAP(karr[i], karr[j], ktmp);
                }
            }
            i++;
            KVSWAP(varr[i], varr[hi], vtmp);
            KVSWAP(karr[i], karr[hi], ktmp);

            if (i > k){
                hi = i - 1;
            } else if (i < k){
                lo = i + 1;
            } else{
                break;
            }
        }
    }

    /**
     * This function puts the `k` largest values in the beginning of the [start, end)
     * range in the vector.
     * @param start Starting index
     * @param end Index after last index that should be considered.
     * @param k Number of smallest values to move to the start of the vector range.
     */
    template <typename K, typename V>
    void KVVector<K,V>::selectd(
        const size_t start,
        const size_t end,
        const size_t k
    ){
        size_t i, j, lo, hi, mid;
        K ktmp;
        V vtmp;
        V pivot;

        if (end - start <= k)
            return;         /* return if the vector has fewer elements than we want */

        for (lo = start, hi = end - 1; lo < hi;)
        {
            mid = lo + ((hi - lo) >> 1);

            /* select the median */
            if (varr[lo] < varr[mid]){
                mid = lo;
            }
            if (varr[hi] > varr[mid]){
                mid = hi;
            } else{
                goto jump_over;
            }
            if (varr[lo] < varr[mid]){
                mid = lo;
            }

        jump_over:
            //swap mid and high
            KVSWAP(varr[mid], varr[hi], vtmp);
            KVSWAP(karr[mid], karr[hi], ktmp);
            pivot = varr[hi];

            /* the partitioning algorithm */
            for (i = lo - 1, j = lo; j < hi; j++)
            {
                if (varr[j] >= pivot)
                {
                    i++;
                    KVSWAP(varr[i], varr[j], vtmp);
                    KVSWAP(karr[i], karr[j], ktmp);
                }
            }
            i++;
            KVSWAP(varr[i], varr[hi], vtmp);
            KVSWAP(karr[i], karr[hi], ktmp);


            if (i > start + k){
                hi = i - 1;
            } else if (i < start + k){
                lo = i + 1;
            } else{
                break;
            }
        }
    }

    /**
     * This function puts the `k` smallest values in the beginning of the vector.
     * @param k Number of smallest values to move to the start of the vector.
     */
    template <typename K, typename V>
    void KVVector<K,V>::selecti(const size_t k)
    {
        size_t i, j, lo, hi, mid;
        K ktmp;
        V vtmp;
        V pivot;

        if (this->vec_sz <= k)
            return;         /* return if the vector has fewer elements than we want */

        for (lo = 0, hi = this->vec_sz - 1; lo < hi;)
        {
            mid = lo + ((hi - lo) >> 1);

            /* select the median */
            if (varr[lo] > varr[mid]){
                mid = lo;
            }
            if (varr[hi] < varr[mid]){
                mid = hi;
            } else{
                goto jump_over;
            }
            if (varr[lo] > varr[mid]){
                mid = lo;
            }

        jump_over:
            //swap mid and high
            KVSWAP(varr[mid], varr[hi], vtmp);
            KVSWAP(karr[mid], karr[hi], ktmp);
            pivot = varr[hi];

            /* the partitioning algorithm */
            for (i = lo - 1, j = lo; j < hi; j++)
            {
                if (varr[j] <= pivot)
                {
                    i++;
                    KVSWAP(varr[i], varr[j], vtmp);
                    KVSWAP(karr[i], karr[j], ktmp);
                }
            }
            i++;
            KVSWAP(varr[i], varr[hi], vtmp);
            KVSWAP(karr[i], karr[hi], ktmp);


            if (i > k){
                hi = i - 1;
            } else if (i < k){
                lo = i + 1;
            } else{
                break;
            }
        }
    }

    /**
     * This function puts the `k` smallest values in the beginning of the [start, end)
     * range in the vector.
     * @param start Starting index
     * @param end Index after last index that should be considered.
     * @param k Number of smallest values to move to the start of the vector range.
     */
    template <typename K, typename V>
    void KVVector<K,V>::selecti(
        const size_t start,
        const size_t end,
        const size_t k
    ){
        size_t i, j, lo, hi, mid;
        K ktmp;
        V vtmp;
        V pivot;

        if (end - start <= k)
            return;         /* return if the vector has fewer elements than we want */

        for (lo = start, hi = end - 1; lo < hi;)
        {
            mid = lo + ((hi - lo) >> 1);

            /* select the median */
            if (varr[lo] > varr[mid]){
                mid = lo;
            }
            if (varr[hi] < varr[mid]){
                mid = hi;
            } else{
                goto jump_over;
            }
            if (varr[lo] > varr[mid]){
                mid = lo;
            }

        jump_over:
            //swap mid and high
            KVSWAP(varr[mid], varr[hi], vtmp);
            KVSWAP(karr[mid], karr[hi], ktmp);
            pivot = varr[hi];

            /* the partitioning algorithm */
            for (i = lo - 1, j = lo; j < hi; j++)
            {
                if (varr[j] <= pivot)
                {
                    i++;
                    KVSWAP(varr[i], varr[j], vtmp);
                    KVSWAP(karr[i], karr[j], ktmp);
                }
            }
            i++;
            KVSWAP(varr[i], varr[hi], vtmp);
            KVSWAP(karr[i], karr[hi], ktmp);


            if (i > start + k){
                hi = i - 1;
            } else if (i < start + k){
                lo = i + 1;
            } else{
                break;
            }
        }
    }

    #undef KVSWAP

    template <typename K, typename V>
    void KVVector<K,V>::printvec() {
        for(int i = 0; i < vec_sz; i++){
            std::cout << "("<< (K)karr[i] << ", " << (V)varr[i] << ") ";
        }
        std::cout << std::endl;
    }


} /* end namespace*/

/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#pragma once

#include "KVVector.h"

/**
 * A Key-Value pair "view" into a pair of key and value arrays or
 * a value array (with keys assumed to be 0->n-1). The vector will
 * not accept new elements, but can be used to perform other non-size
 * modifying changes to the array(s), e.g., sort, select, or scan.
 */

namespace snnlib {

    /**
     * @tparam K    Key typename KVVector<K, V>::Type
     * @tparam V    Value typename KVVector<K, V>::Type
     */
    template<typename K, typename V>
    class KVVectorView    : public KVVector<K,V> {
    public:
        KVVectorView(const std::size_t n, K* keys, V* vals);
        KVVectorView(const std::size_t n, V* vals);
        ~KVVectorView();

        void replace(const std::size_t n, K* keys, V* Vals);
        void replace(const std::size_t n, V* Vals);
        void set_size(const std::size_t n);

        /** Method overrides to throw errors */
        void assign(const std::size_t, const typename KVVector<K, V>::T &keyval);

        void assign(const std::size_t, const K &key, const V &val);

        void assign(const std::size_t, const K *keys, const V *vals);

        void assign(const std::size_t, const V *vals);

        void assign(const std::size_t, const V &val);

        void assign(std::initializer_list<typename KVVector<K, V>::T>);

        void assign(std::initializer_list<V>);

        void resize(const std::size_t);

        void resize(const std::size_t, const typename KVVector<K, V>::T &);

        void resize(const std::size_t, const V &);

        void resize(const std::size_t, const K &, const V &);

        void reserve(const std::size_t);

        void shrink_to_fit();

        void push_back(const typename KVVector<K, V>::T &);

        void push_back(const K &, const V &);

        void push_back(const V &);

        void push_back(typename KVVector<K, V>::T &&);

        K * insert(const K *, const typename KVVector<K, V>::T &);

        K * insert(const K *, const V &);

        K * insert(const K *, typename KVVector<K, V>::T &&);

        K * insert(const K *, V &&);

        K * insert(const K *, const std::size_t, const typename KVVector<K, V>::T &);

        K * insert(const K *, const std::size_t, const V &);

        template<class InputIt>
        K * insert(const K *, InputIt, InputIt, const KVVector<K,V> &other);

        K * insert(const K *, std::initializer_list<typename KVVector<K, V>::T>);

        K * insert(const K *, std::initializer_list<V>);

        void swap(KVVector<K, V> &);

        template<class ... Args>
        K * emplace(const K *, Args &&... args);

        template<class ... Args>
        void emplace_back(Args &&... args);

    private:
        //bool to track initialization of key array
        bool keys_created;
    };

    template <typename K, typename V>
    KVVectorView<K,V>::KVVectorView(const std::size_t n, K* keys, V *vals) : KVVector<K,V>::KVVector(0)
    {
        KVVector<K,V>::vec_sz = n;
        KVVector<K,V>::max_sz = n;
        KVVector<K,V>::karr = keys;
        KVVector<K,V>::varr = vals;
        keys_created = false;
    }

    template <typename K, typename V>
    KVVectorView<K,V>::KVVectorView(const std::size_t n, V *vals) : KVVector<K,V>::KVVector(0)
    {
        KVVector<K,V>::vec_sz = n;
        KVVector<K,V>::max_sz = n;
        KVVector<K,V>::varr = vals;
        KVVector<K,V>::karr = (K * ) realloc(KVVector<K,V>::karr, sizeof(K) * n);
        for(int i = 0; i < n; i++){
            KVVector<K,V>::karr[i] = (K) i;
        }
        keys_created = true;
    }

    template<typename K, typename V>
    KVVectorView<K, V>::~KVVectorView() {
        KVVector<K,V>::varr = nullptr;
        if(!keys_created){
            KVVector<K,V>::karr = nullptr;
        }
    }

    //null param if keys does not exist
    template<typename K, typename V>
    void KVVectorView<K, V>::replace(const std::size_t n, K * keys, V* vals) {
        KVVector<K,V>::vec_sz = n;
        KVVector<K,V>::max_sz = n;
        KVVector<K,V>::varr = vals;
        KVVector<K,V>::karr = keys;
        keys_created = false;
    }

    template<typename K, typename V>
    void KVVectorView<K, V>::replace(const std::size_t n, V* vals) {
        KVVector<K,V>::vec_sz = n;
        KVVector<K,V>::max_sz = n;
        KVVector<K,V>::varr = vals;
        KVVector<K,V>::karr = (K * ) realloc(KVVector<K,V>::karr, sizeof(K) * n);
        for(int i = 0; i < n; i++){
            KVVector<K,V>::karr[i] = (K) i;
        }
        keys_created = true;
    }

    /**
     * Sets vec_sz to inputted parameter without truncating/allocating different memory
     * Used for partial sorting
     * @param n size of new vector
     */
    template<typename K, typename V>
    void KVVectorView<K, V>::set_size(const std::size_t n){
        KVVector<K,V>::vec_sz = n;
    }

    /**
     * Function Overides to throw errors
     * These functions reallocate memory in parent class and cannot run in view vector
     */
    template <typename K, typename V>
    void KVVectorView<K, V>::assign(const std::size_t, const typename KVVector<K, V>::T &keyval){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::assign(const std::size_t, const K &key, const V &val){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::assign(const std::size_t, const K * keys, const V *vals){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::assign(const std::size_t, const V *vals){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::assign(const std::size_t, const V &val){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::assign(std::initializer_list<typename KVVector<K, V>::T>){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::assign(std::initializer_list<V>){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::resize(const std::size_t){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::resize(const std::size_t, const typename KVVector<K, V>::T &){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::resize(const std::size_t, const V &){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::resize(const std::size_t, const K &, const V &){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::reserve(const std::size_t){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::shrink_to_fit(){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::push_back(const typename KVVector<K, V>::T &){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::push_back(const K &, const V &){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::push_back(const V &){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void KVVectorView<K, V>::push_back(typename KVVector<K, V>::T &&){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    K * KVVectorView<K, V>::insert(const K *, const typename KVVector<K, V>::T &){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    K * KVVectorView<K, V>::insert(const K *, const V &){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    K * KVVectorView<K, V>::insert(const K *, typename KVVector<K, V>::T &&){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    K * KVVectorView<K, V>::insert(const K *, V &&){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    K * KVVectorView<K, V>::insert(const K *, const std::size_t, const typename KVVector<K, V>::T &){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    K * KVVectorView<K, V>::insert(const K *, const std::size_t, const V &){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    template<class InputIt>
    K * KVVectorView<K, V>::insert(const K *, InputIt, InputIt, const KVVector<K,V> &other){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    K * KVVectorView<K, V>::insert(const K *, std::initializer_list<typename KVVector<K, V>::T>){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    K * KVVectorView<K, V>::insert(const K *, std::initializer_list<V>){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    void    KVVectorView<K, V>::swap(KVVector<K, V> &){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    template<class ... Args>
    K * KVVectorView<K, V>::emplace(const K*, Args &&... args){
        throw ("KVVectorView: View only vector");
    }

    template <typename K, typename V>
    template<class ... Args>
    void KVVectorView<K, V>::emplace_back(Args &&... args){
        throw ("KVVectorView: View only vector");
    }

} /* end namespace */

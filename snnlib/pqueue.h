/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#pragma once

#include <cassert>

#include "basetypes.h"

namespace snnlib {

    /**
     * Min-Heap / Min-Priority Queue data structure with the added constraint that,
     * when the heap gets full, items smaller than the min will be dropped (not inserted).
     * Additionally, items being added to the queue are expected to be unique, i.e., 
     * they don't already exist in the queue. There is no update feature for the queue.
    */
    typedef struct pqueue {

        size_t nnodes;   /* number of items currently in the queue */
        size_t maxnodes; /* max number of nodes that can be stored in the heap */
        bool local;      /* whether the heap was allocated by this object */

        /* Heap version of the data structure */
        ivkv_t *heap;    /* array of size maxsize */ //allocated in this class

        //Members
        /**
         * Create a priority queue with maxnodes nodes and internal memory.
         * @param maxnodes Number of maximum nodes in the queue.
        */
        pqueue(size_t const maxnodes) {
            this->nnodes = 0;
            this->maxnodes = maxnodes;
            this->local = maxnodes > 0;
            this->heap = this->local ? (ivkv_t *) malloc(maxnodes * sizeof(ivkv_t)) : nullptr;
            for (size_t i=0; i < maxnodes; ++i){
                this->heap[i].key = -1;
                this->heap[i].val = 0.0;
            }
        }

        /**
         * Create a priority queue with maxnodes nodes and external memory.
         * @param maxnodes Number of maximum nodes in the queue.
         * @param maxnodes Current number of nodes in the queue.
         * @param maxnodes Storage for the queue.
        */
        pqueue(size_t const maxnodes, size_t const nnodes, ivkv_t * heap){
            this->maxnodes = maxnodes;
            this->nnodes = nnodes;
            this->heap = heap;
            this->local = false;
        }

        ~pqueue(){
            if(this->heap && this->local){
                free(this->heap);
            }
        }

        /**
         * Change the current queue to a priority queue with the given external memory.
         * @param maxnodes Number of maximum nodes in the queue.
         * @param maxnodes Current number of nodes in the queue.
         * @param maxnodes Storage for the queue.
        */
        void reset(size_t const maxnodes, size_t const nnodes, ivkv_t * heap){
            if(this->heap && this->local){
                free(this->heap);
            }
            this->maxnodes = maxnodes;
            this->nnodes = nnodes;
            this->heap = heap;
            this->local = false;
        }

        /**
         * Check if the queue is a heap.
         * @param i Index to check.
        */
        bool check(int i=0){

            // If a leaf node
            if (i >= ((ssize_t)this->nnodes - 3) / 2){
                return true;
            }

            if (this->heap[i].val <= this->heap[2 * i + 1].val &&
                this->heap[i].val <= this->heap[2 * i + 2].val
                && check(2 * i + 1)
                && check(2 * i + 2)){
                return true;
            }

            return false;
        }

        /**
         * Check and fix the queue, if necessary.
        */
        void fix(){
            if(! this->check()){
                // data not stored as a heap -- fix it
                size_t nn = this->nnodes;
                this->nnodes = 0;
                for(size_t i=0; i < nn; ++i){
                    this->insert(this->heap[i].key, this->heap[i].val);
                }
            }
            assert(this->check());
        }

        /**
         * This function adds an item in the priority queue if queue is not full. 
         * If full and top item has lower value than the inserted item, it is
         * replaced by the current item.
         * @param node Node key
         * @param val Node value
        */
        int insert(idx_t const node, val_t const val){
            ptr_t i, j, loc;

            if (this->nnodes == this->maxnodes) {
                if (val > heap[0].val) {
                    get_top(); /* remove min elem from heap to replace with node/val */
                } else {
                    return 0; /* heap is full; ignore item */
                }
            }

            i = this->nnodes++;
            while (i > 0)
            {
                j = (i - 1) >> 1;
                if (val < heap[j].val) {
                    heap[i].key = heap[j].key;
                    heap[i].val = heap[j].val;
                    i = j;
                } else {
                    break;
                }
            }
            assert (i >= 0);
            heap[i].val = val;
            heap[i].key = node;

            return 1;
        }

        /**
         * This function adds an item in the priority queue if queue is not full. 
         * If full and top item has lower value than the inserted item, it is
         * replaced by the current item.
         * @param item Node item (key and value) to be inserted
        */
        int insert(ivkv_t item){
            return insert(item.key, item.val);
        }

        
        /*
         * This function returns the item at the top of the queue and removes
         * it from the priority queue.
        */
        ivkv_t get_top(){
            ptr_t i, j;
            idx_t node;
            val_t val;

            if (this->nnodes == 0){
                return ivkv_t{-1, 0.0};
            }

            this->nnodes--;
            ivkv_t pair = heap[0];
            if ((i = this->nnodes) > 0)
            {
                val = heap[i].val;
                node = heap[i].key;
                i = 0;
                while ((j = 2 * i + 1) < this->nnodes) {
                    if (heap[j].val < val) {
                        if (j + 1 < this->nnodes && (heap[j+1].val < heap[j].val)){
                            j = j + 1;
                        }
                        heap[i].key = heap[j].key;
                        heap[i].val = heap[j].val;
                        i = j;
                    } else if (j + 1 < this->nnodes && ((val_t)heap[j+1].val < val)) {
                        j = j + 1;
                        heap[i].key = heap[j].key;
                        heap[i].val = heap[j].val;
                        i = j;
                    } else {
                        break;
                    }
                }
                heap[i].key = node;
                heap[i].val = val;
            }
            return pair;
        }


        /*
         * This function returns the key of the top item. The item is not
         * deleted from the queue.
        */
        inline ivkv_t see_top() {
            return this->heap[0];
        }


        /*
         * This function returns the key of the top item. The item is not
         * deleted from the queue.
        */
        inline idx_t see_top_key() {
            return this->heap[0].key;
        }

        /*
         * This function returns the value of the top item. The item is not
         * deleted from the queue.
        */
        inline acm_t see_top_val() {
            return this->heap[0].val;
        }

    } pqueue;



}
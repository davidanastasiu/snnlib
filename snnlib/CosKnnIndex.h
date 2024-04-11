/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#pragma once

#include "BitVector.h"
#include "sort.h"
#include "pqueue.h"
#include "csr.h"

namespace snnlib {

    /**
     * Work structure for L2Knng thread
     * @param nrows Number of rows/nodes
     * @param ncols Number of columns/features for the node vectors
     * @param ncands Maximum number of candidates
    */
    struct l2knnw {
        ivkv_t* qcands;           // list of candiadates
        val_t* q{nullptr};        // dense query vector
        acm_t* n{nullptr};        // dense query suffix norms
        idx_t* m{nullptr};        // marker list of found candidates
        acm_t* a{nullptr};        // accumulator for block computation
        BitVector * bv;           // bit vector to mark existing items
        pqueue * queue;           // heap for a given neighborhood
        csr_t * sims;             // found similarities
        // stats
        ptr_t ninitcands{0};      // Number of initial candidates
        ptr_t ndotproducts{0};    // Number of dot-products
        ptr_t npruned{0};         // Number of pruned dot-products

        /**
         * Instsantiate a work structure for L2Knng
         * @param nrows Number of samples
         * @param ncols Number of attributes
         * @param ncands Number of candidates
         * @param k Number of neighbors looking for
         * @param norms Whether to store prefix norms
        */
        l2knnw(const idx_t nrows, const idx_t ncols, const idx_t ncands, 
            const idx_t k, idx_t ndrows, bool norms=true){
            qcands = (ivkv_t*) malloc(ncands * sizeof(ivkv_t));
            a = (acm_t*) calloc(ndrows, sizeof(acm_t));
            if(norms){
                q = (val_t*) calloc(ncols, sizeof(val_t));
                n = (acm_t*) calloc(ncols, sizeof(acm_t));
                if(!qcands || !q || !n){
                    throw std::bad_alloc();
                }
            } else {
                m = (idx_t*) malloc(nrows * sizeof(idx_t));
                if(!qcands || !m){
                    throw std::bad_alloc();
                }
            }
            bv = new BitVector(nrows);
            queue = new pqueue(0);  // empty queue, will be reset
            sims = new csr_t(0, ncols, nrows, nrows * k);
        }

        ~l2knnw(){
            free(qcands);
            free(a);
            if(q){
                free(q);
            }
            if(n){
                free(n);
            }
            if(m){
                free(m);
            }
            delete bv;
            delete queue;
            delete sims;
        }

        static void print_stats(l2knnw** twork, int nthreads){
            ptr_t ninitcands = 0;
            ptr_t ndotproducts = 0;
            ptr_t npruned = 0;
            for(idx_t i=0; i < nthreads; ++i){
                ninitcands += twork[i]->ninitcands;
                ndotproducts += twork[i]->ndotproducts;
                npruned += twork[i]->npruned;
            }
            std::cout << "ninitcands: " << ninitcands << std::endl;
            std::cout << "ndotproducts: " << ndotproducts << std::endl;
            std::cout << "npruned: " << npruned << std::endl << std::flush;
        }
    };


    struct CosKnnIndex {

        idx_t nqrows{20000};   // Number of rows for a query block.
        idx_t ndrows{100000};  // Number of rows for a database block.
        ptr_t ninnz{1000000};  // Number of non-zeros for each inverted index.
        bool processed{false}; // whether input has been processed

        csr_t * data{nullptr};

        /**
         * Initialize an empty index that data can be added to.
         * @param ncols  Number of columns/features for the data to be added to the index.
         * @param nqrows Number of rows for a query block
         * @param ndrows Number of rows for a database block
         * @param ninnz  Number of non-zeros for each inverted index
        */
        CosKnnIndex(
            const idx_t ncols,
            const idx_t nqrows=2e+4,
            const idx_t ndrows=1e+5,
            const ptr_t ninnz=1e+6
        ) : nqrows(nqrows), ndrows(ndrows), ninnz(ninnz){
            data = new csr_t(0, ncols, 10, 1000);
        }

        /**
         * Initialize an index with data from the input CSR matrix.
         * @param data   Data to be added to the index
         * @param nqrows Number of rows for a query block
         * @param ndrows Number of rows for a database block
         * @param ninnz  Number of non-zeros for each inverted index
        */
        CosKnnIndex(
            csr_t* data,
            const idx_t nqrows=2e+4,
            const idx_t ndrows=1e+5,
            const ptr_t ninnz=1e+6
        ) : nqrows(nqrows), ndrows(ndrows), ninnz(ninnz){
            this->data = new csr_t(data);
        }

        ~CosKnnIndex()
        {
            delete data;
        }

        void set_nthreads(int n){
            omp_set_num_threads(n);
        }

        int get_nthreads(){
            int n;
            #pragma omp parallel
            {
                n = omp_get_num_threads();
            }
            return n;
        }

        /**
         * Brute-force k-nng construction algorithm.
         * @param k Number of nearest neighbor for each node.
         * @param verbose Verbose, i.e., print log statements
        */
        csr_t * knngbf(idx_t const k, bool verbose=false){
            int nthreads = this->get_nthreads();
            if(verbose){
                data->print_info();
                std::cout << "k: " << k << std::endl;
                std::cout << "nthreads: " << nthreads << std::endl;
            }
            // process index input data
            if(!processed){
                data->compact_cols();
                data->normalize();
                processed = true;
            }
            if(!data->cptr){
                data->create_index(CSR_BASE::ROW);
            }
            // allocate space for dense version of the result
            ivkv_t** nbrs = (ivkv_t**) malloc(data->nrows * sizeof(ivkv_t*));
            idx_t* nnbrs = (idx_t*) calloc(data->nrows, sizeof(idx_t));
            if(!nbrs || !nnbrs){
                throw std::bad_alloc();
            }
            for(idx_t i=0; i < data->nrows; ++i){
                nbrs[i] = (ivkv_t*) malloc(k * sizeof(ivkv_t));
                if(!nbrs[i]){
                    throw std::bad_alloc();
                }
            }
            if(verbose){
                std::cout << "nbrs storage: " << data->nrows << " x " << k << std::endl;
            }

            // create working memory for all threads
            l2knnw** twork = (l2knnw**) malloc(nthreads * sizeof(l2knnw*));
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                twork[tid] = new l2knnw(data->nrows, data->ncols, data->nrows, k, ndrows, false);
            }

            // find neighbors
            #pragma omp parallel 
            {
                int const tid = omp_get_thread_num();
                ivkv_t* qcands = twork[tid]->qcands;
                idx_t* m = twork[tid]->m; // marker list of found candidates
                auto bv = twork[tid]->bv; // bit vector for marked candidates

                #pragma omp for schedule(static, 32)
                for(idx_t i=0; i < data->nrows; ++i){
                    idx_t nc = 0;

                    // compute dot-products
                    for(ptr_t j=data->rptr[i]; j < data->rptr[i+1]; ++j){
                        idx_t c = data->rind[j]; // column id
                        val_t v = data->rval[j]; // value for q
                        for(ptr_t k=data->cptr[c]; k < data->cptr[c+1]; ++k){
                            idx_t cid = data->cind[k]; // candidate
                            if(cid == i){
                                continue;
                            }
                            if(bv->isset(cid)){
                                qcands[m[cid]].val += v * data->cval[k];
                                continue;
                            }
                            qcands[nc].key = cid;
                            qcands[nc].val = v * data->cval[k];
                            m[cid] = nc++;
                            bv->set(cid);
                        }
                    }

                    // clear bit-vector
                    for(idx_t j=0; j < nc; ++j){
                        bv->zero(qcands[j].key);
                    }

                    // select the top-k current neighbors, rank them, then copy them to general storage
                    selectsortd(nc, qcands, k);
                    nc = std::min(nc, k);
                    memmove((void*)nbrs[i], (void*)qcands, nc * sizeof(ivkv_t));
                    nnbrs[i] = nc;
                }
            }

            // construct final graph adjacency matrix
            csr_t * graph = _graph_to_csr(data->nrows, nbrs, nnbrs);

            // free memory
            for(idx_t i=0; i < nthreads; ++i){
                delete twork[i];
            }
            for(idx_t i=0; i < data->nrows; ++i){
                free(nbrs[i]);
            }
            free(twork);
            free(nbrs);
            free(nnbrs);
            return graph;
        }

        /**
         * Construct the k-nearest neighbor graph for the stored data using cosine similarity
         * as the proximity measure. Data are by default normalized so they need not be 
         * normalized before they are added to index.
         * @param k Number of nearest neighbor for each node.
         * @param approximate Whether to return the approximate graph (containing most, but
         *                  not necessarily all of the nearest neighbors)
         * @param cfactor Candidate list factor, i.e., multiple of k for deciding the candidate 
         *                  list size. Must be >= 1.
         * @param nenhance Number of times to run graph enhancement during approximate search.
         * @param verbose Verbose, i.e., print log statements
        */
        csr_t * knng(
            idx_t const k,
            const bool approximate=false,
            const float cfactor=2.0,
            const int nenhance=1,
            const double minprune=0.7,
            bool verbose=false
        ){
            int nthreads = this->get_nthreads();
            if(verbose){
                data->print_info();
                std::cout << "k: " << k << std::endl;
                std::cout << "approximate: " << approximate << std::endl;
                std::cout << "cfactor: " << cfactor << std::endl;
                std::cout << "nenhance: " << nenhance << std::endl;
                std::cout << "minprune: " << minprune << std::endl;
                std::cout << "nthreads: " << nthreads << std::endl;
            }
            // process index input data
            if(!processed){
                data->compact_cols();
                _permute_data_columns();
                data->normalize();
                processed = true;
            }
            // allocate space for dense version of the result
            idx_t maxcands = std::max((idx_t) (cfactor * k), k); // candidate list size
            ivkv_t** nbrs = (ivkv_t**) malloc(data->nrows * sizeof(ivkv_t*));
            idx_t* nnbrs = (idx_t*) calloc(data->nrows, sizeof(idx_t));
            if(!nbrs || !nnbrs){
                throw std::bad_alloc();
            }
            for(idx_t i=0; i < data->nrows; ++i){
                nbrs[i] = (ivkv_t*) malloc(maxcands * sizeof(ivkv_t));
                if(!nbrs[i]){
                    throw std::bad_alloc();
                }
            }
            if(verbose){
                std::cout << "nbrs storage: " << data->nrows << " x " << maxcands << std::endl;
            }

            // create working memory for all threads
            l2knnw** twork = (l2knnw**) malloc(nthreads * sizeof(l2knnw*));
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                twork[tid] = new l2knnw(data->nrows, data->ncols, maxcands, k, ndrows);
            }

            /*** STEP1: Initial KNNG selection based on large values + find rows with less than k neighbors. ***/
            _construct_initial_graph(twork, nbrs, nnbrs, maxcands, k, minprune);

            /*** STEP2: Enhance initial neighbors with neighbor's neighbors. ***/
            for(idx_t i=0; i < nenhance; ++i){
                _enhance(twork, nbrs, nnbrs, maxcands, k, minprune);
            }

            if(approximate){
                if(verbose){
                    l2knnw::print_stats(twork, nthreads);
                    std::cout << std::endl << std::flush;
                }
                csr_t * graph = _graph_to_csr(data->nrows, nbrs, nnbrs);
                for(idx_t i=0; i < nthreads; ++i){
                    delete twork[i];
                }
                for(idx_t i=0; i < data->nrows; ++i){
                    free(nbrs[i]);
                }
                free(twork);
                free(nbrs);
                free(nnbrs);
                return graph;
            }

            /*** STEP3: Verify found neighbors, adjusting lists as appropriate. ***/
            csr_t * graph = _make_exact(twork, nbrs, nnbrs, maxcands, k);

            // free memory
            for(idx_t i=0; i < nthreads; ++i){
                delete twork[i];
            }
            for(idx_t i=0; i < data->nrows; ++i){
                free(nbrs[i]);
            }
            free(twork);
            free(nbrs);
            free(nnbrs);
            return graph;
        }

        /**
         * Permute columns in the input matrix to promote faster pruning
        */
        void _permute_data_columns(){

            // permute columns in increasing non-zero count order
            idx_t * ccnts = (idx_t*) calloc(data->ncols, sizeof(idx_t));
            idx_t * cperm = (idx_t*) malloc(data->ncols * sizeof(idx_t));
            if(!cperm || !ccnts){
                throw std::bad_alloc();
            }
            for(ptr_t j=0; j < data->nnz(); ++j){
                ccnts[data->rind[j]]++;
            }
            for(idx_t i=0; i < data->ncols; ++i){
                cperm[i] = i;
            }
            // find the permutation order
            kvsorti(data->ncols, cperm, ccnts);
            // permute the data 
            data->permute_columns(cperm);
            // need data sorted by indices for pruning to properly work
            data->sort_nnzs_by_ind(CSR_BASE::ROW);

            // free temporary memory used in this stage
            free(ccnts);
            free(cperm);
        }

        /**
         * Find initial candidates for a given node `qid`
         * @param qid Query ID, i.e., node we're finding candidates for
         * @param nn Number of candidate neighbors to find
         * @param vdata Value-based sorted matrix with both row and column bases loaded
         * @param bv Bit vector of sise number of nodes, all zeroed
         * @param qcands Query candidates storage (return)
        */
        idx_t _find_initial_candidates(
            const idx_t qid, 
            const idx_t nn, 
            const csr_t* const vdata, 
            BitVector* bv,
            ivkv_t* qcands
        ){
            // find candidates
            ptr_t j = vdata->rptr[qid];
            ptr_t const lj = vdata->rptr[qid+1];
            idx_t c1 = vdata->rind[j];
            val_t cv1 = vdata->rval[j];
            ptr_t cs1 = vdata->cptr[c1];
            ptr_t ce1 = vdata->cptr[c1+1];
            idx_t c2 = -1;
            val_t cv2 = 0;
            ptr_t cs2 = 0;
            ptr_t ce2 = 0;
            idx_t cid;
            idx_t ncands = 0;
            if(j+1 < lj){
                j++;
                c2 = vdata->rind[j];
                cv2 = vdata->rval[j];
                cs2 = vdata->cptr[c2];
                ce2 = vdata->cptr[c2+1];
            }
            while(ncands < nn){ // find at most nn candidates
                // edge cases
                if(c1 < 0){ /** last column */
                    for(; cs2 < ce2 && ncands < nn; ++cs2){
                        cid = vdata->cind[cs2];
                        if(cid == qid || bv->isset(cid)){
                            continue; // no self-similarity or already considered
                        }
                        qcands[ncands].key = cid;
                        qcands[ncands++].val = vdata->cval[cs2];
                        bv->toggle(cid);
                    }
                    break;
                } else if(c2 < 0){ /** last column */
                    for(; cs1 < ce1 && ncands < nn; ++cs1){
                        cid = vdata->cind[cs1];
                        if(cid == qid || bv->isset(cid)){
                            continue; // no self-similarity or already considered
                        }
                        qcands[ncands].key = cid;
                        qcands[ncands++].val = vdata->cval[cs1];
                        bv->toggle(cid);
                    }
                    break;
                }
                // find next candidate
                if(vdata->cind[cs1] == qid){
                    cs1++;
                } else if(vdata->cind[cs2] == qid){
                    cs2++;
                } else {
                    if (vdata->cval[cs1] * cv1 > vdata->cval[cs2] * cv2){
                        cid = vdata->cind[cs1];
                        if(!bv->isset(cid)){
                            qcands[ncands].key = cid;
                            qcands[ncands++].val = vdata->cval[cs1];
                            bv->toggle(cid);
                        }
                        cs1++;
                    } else {
                        cid = vdata->cind[cs2];
                        if(!bv->isset(cid)){
                            qcands[ncands].key = cid;
                            qcands[ncands++].val = vdata->cval[cs2];
                            bv->toggle(cid);
                        }
                        cs2++;
                    }
                }
                // advance lists
                if (cs1 == ce1){
                    if (j + 1 < lj){
                        j++;
                        c1 = vdata->rind[j];
                        cv1 = vdata->rval[j];
                        cs1 = vdata->cptr[c1];
                        ce1 = vdata->cptr[c1+1];
                    } else {
                        c1 = -1;
                    }
                }
                if (cs2 == ce2){
                    if (j + 1 < lj){
                        j++;
                        c2 = vdata->rind[j];
                        cv2 = vdata->rval[j];
                        cs2 = vdata->cptr[c2];
                        ce2 = vdata->cptr[c2+1];
                    } else {
                        c2 = -1;
                    }
                }
            }
            // reset bit vector
            for(idx_t i=0; i < ncands; ++i){
                bv->zero(qcands[i].key);
            }
            return ncands;
        }


        /**
         * Sparse-dense dot-product of row `rid` from the `data` matrix and dense vector `q`
         * @param rid Row id in data representing the sparse vector
         * @param data The sparse matrix who's row is multiplied.
         * @param q Dense vector multiplying against
        */
        inline acm_t _sddp(
            const idx_t rid, 
            const csr_t* const data, 
            const val_t* const q
        ){
            acm_t s = 0.0;
            for(ptr_t j=data->rptr[rid]; j < data->rptr[rid+1]; ++j){
                s += data->rval[j] * q[data->rind[j]];
            }
            return s;
        }

        /**
         * Sparse-dense dot-product of row `rid` from the `data` matrix and dense vector `q`
         * with pruning.
         * @param minsim Minimum similarity that the dot-product should achieve.
         * @param cid Row id in data representing the candidate sparse vector
         * @param data The sparse matrix who's row is multiplied.
         * @param q Dense vector multiplying against
         * @param n Suffix norms vector for the dense vector in `q`
        */
        inline acm_t _sddp(
            l2knnw* tw,
            const acm_t minsim,
            const idx_t qid,
            const idx_t cid,
            const csr_t* const data, 
            const val_t* const q, 
            const acm_t* const n
        ){
            acm_t s = 0.0;
            for(ptr_t j=data->rptr[cid]; j < data->rptr[cid+1]; ++j){
                idx_t c = data->rind[j];
                val_t v = q[c];
                if(v > 0){
                    s += data->rval[j] * v;
                    if(s + data->rnorms[j] * n[c] < minsim){
                        tw->npruned++;
                        s = -1;
                        break;
                    }
                }
            }
            return s;
        }

        /**
         * Create CSR structure from dense neighborhood structure
        */
        inline csr_t * _graph_to_csr(
            idx_t nrows,
            ivkv_t** data,
            idx_t* sizes
        ){
            ptr_t nnz = 0;
            for(idx_t i=0; i < nrows; ++i){
                nnz += sizes[i];
            }
            auto graph = new csr_t(nrows, nrows, nrows, nnz);
            graph->rptr[0] = 0;
            for(idx_t i=0; i < nrows; ++i){
                graph->rptr[i+1] = graph->rptr[i] + sizes[i];
            }
            for(idx_t i=0; i < nrows; ++i){
                for(idx_t j=0; j<sizes[i]; ++j){
                    graph->rind[graph->rptr[i]+j] = data[i][j].key;
                    graph->rval[graph->rptr[i]+j] = data[i][j].val;
                }
            }
            return graph;
        }

        /**
         * Construct the initial approximate graph with at least k neighbors (if possible)
         * for each node.
         * @param twork Thread work structure with reusable memory allocations
         * @param nbrs 2-D array of neighbors for each node
         * @param nnbrs Neighborhood sizes for each node neighborhood
         * @param maxcands Candidate list size
         * @param k Number of neighbors we ultimately want
         * @param minprune Minimim similarity above which pruning should be used in dotp computation
        */
        void _construct_initial_graph(
            l2knnw** twork,
            ivkv_t** nbrs,
            idx_t* nnbrs,
            idx_t maxcands,
            idx_t k,
            acm_t minprune
        ){
            // make a copy of the data for candidate selection
            csr_t * vdata = new csr_t(data);     
            vdata->create_index(CSR_BASE::ROW);  // add inverted index
            vdata->sort_nnzs_by_value(false, CSR_BASE::ROW, maxcands); // select-sort maxcands of each row
            vdata->sort_nnzs_by_value(false, CSR_BASE::COL, maxcands); // select-sort maxcands of each column
            
            // copute suffix norms for rows in data to use for pruning
            data->compute_suffix_norms();

            // find initial neighbors
            #pragma omp parallel 
            {
                int const tid = omp_get_thread_num();
                l2knnw* tw = twork[tid];
                ivkv_t* qcands = tw->qcands;
                val_t* q = tw->q;
                acm_t* n = tw->n;
                auto bv = tw->bv;
                auto queue = tw->queue;

                #pragma omp for schedule(static, 32)
                for(idx_t i=0; i < data->nrows; ++i){
                    // hash query row
                    for(ptr_t j=data->rptr[i]; j < data->rptr[i+1]; ++j){
                        idx_t cid = data->rind[j];
                        q[cid] = data->rval[j];
                        n[cid] = data->rnorms[j];
                    }
                    idx_t ncands = _find_initial_candidates(i, maxcands, vdata, bv, qcands);
                    assert(bv->count() == 0);
                    tw->ninitcands += ncands;
                    selectd(ncands, qcands, k);
                    // compute the similarity for the first k candidates
                    acm_t minsim = 1.0;
                    idx_t nc = 0; // unpruned candidates
                    for(idx_t j=0; j < k && j < ncands; ++j){
                        idx_t cid = qcands[j].key;
                        acm_t s = _sddp(cid, data, q);
                        tw->ndotproducts++;
                        qcands[nc].key = cid;
                        qcands[nc++].val = s;
                        if(s < minsim){
                            minsim = s;
                        }
                    }
                    // compute the rest of the similarities
                    for(idx_t j=k; j < ncands; ++j){
                        idx_t cid = qcands[j].key;
                        acm_t s = minsim > minprune ? _sddp(tw, minsim, i, cid, data, q, n) : _sddp(cid, data, q);
                        if(minsim <= minprune){
                            tw->ndotproducts++;
                        }
                        if(s > minsim){
                            qcands[nc].key = cid;
                            qcands[nc++].val = s;
                        }
                    }
                    // clear query row hashing
                    for(ptr_t j=data->rptr[i]; j < data->rptr[i+1]; ++j){
                        q[data->rind[j]] = 0.0;
                    }
                    if(nc == 0){
                        nnbrs[i] = std::min(nc, k);
                        continue;
                    }
                    // select the top-k current neighbors, rank them, then add them to the heap
                    selectd(nc, qcands, k);
                    nnbrs[i] = nc = std::min(nc, k);
                    queue->reset(k, 0, nbrs[i]);
                    for(idx_t i=0; i < nc; ++i){
                        queue->insert(qcands[i]);
                    }
                }

            }

            // free temporary memory used in this stage
            delete vdata;
        }



        /**
         * Find enhancement candidates for a given node `qid`
         * @param qid Query ID, i.e., node we're finding candidates for
         * @param nn Number of candidate neighbors to find
         * @param graph Value-based sorted graph adjancency matrix with both row and column bases loaded
         * @param bv Bit vector of sise number of nodes, all zeroed
         * @param qcands Query candidates storage (return)
         * @param minsim Minimum similarity to consider before accepting candidates
        */
        idx_t _find_enhancement_candidates(
            const idx_t qid, 
            const idx_t nn, 
            const csr_t* const graph, 
            BitVector* bv,
            ivkv_t* qcands,
            acm_t minsim
        ){
            /* tag current neighbors */
            bv->set(qid);
            for(ptr_t j=graph->rptr[qid]; j<graph->rptr[qid+1]; ++j){
                bv->set(graph->rind[j]);
            }
            for(ptr_t j=graph->cptr[qid]; j<graph->cptr[qid+1]; ++j){
                bv->set(graph->cind[j]);
            }

            idx_t ncands = 0;
            for(ptr_t j=graph->rptr[qid]; ncands < nn && j<graph->rptr[qid+1]; ++j){
                idx_t cid = graph->rind[j];
                val_t csim = graph->rval[j];
                for(ptr_t jj=graph->rptr[cid]; ncands < nn && jj<graph->rptr[cid+1]; ++jj){
                    idx_t ccid = graph->rind[jj];
                    if(bv->isset(ccid) || graph->rval[jj] < csim){
                        continue;
                    }
                    qcands[ncands].key = ccid;
                    qcands[ncands++].val = graph->rval[jj];
                }
            }

            /* untag current neighbors and chosen candidates */
            bv->zero(qid);
            for(idx_t i=0; i < ncands; ++i){
                bv->zero(qcands[i].key);
            }
            for(ptr_t j=graph->rptr[qid]; j<graph->rptr[qid+1]; ++j){
                bv->zero(graph->rind[j]);
            }
            for(ptr_t j=graph->cptr[qid]; j<graph->cptr[qid+1]; ++j){
                bv->zero(graph->cind[j]);
            }
            assert(bv->count() == 0);
            return ncands;
        }

        /**
         * Merge the candidate list into the nbrs list without increasing the size of the nbrs list
         * Both lists are sorted in decreasing order and the merge is in-place
         * @param nbrs Neighbors list
         * @param nn Number of neighbors
         * @param cands Candidates list
         * @param nc Number of candidates
         * 
         * @return number of elements inserted into nbrs
        */
        idx_t _merge_top_k(
            ivkv_t* nbrs,
            idx_t nn,
            ivkv_t* cands,
            idx_t nc
        ){
            idx_t n = 0;
            idx_t ins = 0;
            for(idx_t i=0; i<nc; ++i){
                // figure out where to insert value in the list
                while(n < nn && cands[i].val < nbrs[n].val){
                    n++;
                }
                if(n == nn){
                    break;
                }
                // shift values to make room for insertion
                for(idx_t j=nn-1; j>n; j--){
                    nbrs[j] = nbrs[j-1];
                }
                nbrs[n++] = cands[i];
                ins++;
            }
            return ins;
        }

        /**
         * Enhance initial neighbors with neighbor's neighbors
         * @param twork Thread work structure with reusable memory allocations
         * @param nbrs 2-D array of neighbors for each node
         * @param nnbrs Neighborhood sizes for each node neighborhood
         * @param maxcands Candidate list size
         * @param k Number of neighbors we ultimately want
        */
        void _enhance(
            l2knnw** twork,
            ivkv_t** nbrs,
            idx_t* nnbrs,
            idx_t maxcands,
            idx_t k,
            acm_t minprune
        ){
            // create graph adjacency matrix, its transpose, and sort non-zeros in 
            // decreasing order in each row/col
            csr_t * graph = _graph_to_csr(data->nrows, nbrs, nnbrs);
            graph->create_index(CSR_BASE::ROW);
            graph->sort_nnzs_by_value(false, CSR_BASE::ROW);
            graph->sort_nnzs_by_value(false, CSR_BASE::COL);

            
            // enhance neighborhood graph
            #pragma omp parallel 
            {
                int const tid = omp_get_thread_num();
                l2knnw* tw = twork[tid];
                ivkv_t* qcands = tw->qcands;
                val_t* q = tw->q;
                acm_t* n = tw->n;
                auto bv = tw->bv;
                auto queue = tw->queue;

                #pragma omp for schedule(static, 32)
                for(idx_t i=0; i < data->nrows; ++i){
                    // if less than k neighbors, cannot be enhanced; we're done
                    if(nnbrs[i] < k){
                        continue;
                    }
                    queue->reset(k, nnbrs[i], nbrs[i]);
                    acm_t minsim = queue->see_top_val();
                    bool doprune = minsim > minprune;
                    // hash query row
                    if(doprune){
                        for(ptr_t j=data->rptr[i]; j < data->rptr[i+1]; ++j){
                            idx_t cid = data->rind[j];
                            q[cid] = data->rval[j];
                            n[cid] = data->rnorms[j];
                        }
                    } else {
                        for(ptr_t j=data->rptr[i]; j < data->rptr[i+1]; ++j){
                            idx_t cid = data->rind[j];
                            q[cid] = data->rval[j];
                        }
                    }
                    idx_t ncands = _find_enhancement_candidates(i, maxcands, graph, bv, qcands, minsim);
                    // compute the rest of the similarities
                    idx_t nc = 0;
                    for(idx_t j=0; j < ncands; ++j){
                        idx_t cid = qcands[j].key;
                        acm_t s = doprune ? _sddp(tw, minsim, i, cid, data, q, n) :_sddp(cid, data, q);
                        if(minsim <= minprune){
                            tw->ndotproducts++;
                        }
                        if(s > minsim){
                            qcands[nc].key = cid;
                            qcands[nc++].val = s;
                        }
                    }
                    // select the top-k current neighbors, then merge them to general storage
                    selectd(nc, qcands, k);
                    for(idx_t i=0; i < nc; ++i){
                        queue->insert(qcands[i]);
                    }

                    // clear query row hashing
                    for(ptr_t j=data->rptr[i]; j < data->rptr[i+1]; ++j){
                        q[data->rind[j]] = 0.0;
                    }
                }

            }

            // free temporary memory used in this stage
            delete graph;
        }

        void _set_processing_order(
            ivkv_t** nbrs,
            idx_t* nnbrs,
            ivkv_t * order,
            idx_t start
        ){
            for(idx_t i=start; i < data->nrows; ++i){
                idx_t rid = order[i].key;
                order[i].val = nnbrs[rid] > 0 ? nbrs[rid][0].val : -1;
            }
            sorti(data->nrows-start, order+start);
        }

        /**
         * Create partial inverted index 
         * Index stores only the prefix of the vectors, in the given order, and is size-limited
         * by the `ndrows` and `ninnz` parameters. Column indices are internal row ids and the
         * corresponding `data` row id is stored in idx->rids.
        */
        csr_t * _create_index(
            ivkv_t * order,
            ivkv_t** nbrs,
            idx_t start
        ){
            // for each sample, index non-zeros with suffix norm >= neighborhood minsim
            // first find number of rows and number of non-zeros
            ptr_t nnz = 0;
            idx_t nr = 0;
            ptr_t * rptr = (ptr_t *) malloc(sizeof(ptr_t) * ndrows);
            ptr_t * cptr = (ptr_t *) calloc(data->ncols, sizeof(ptr_t));
            if(!rptr || !cptr){
                throw std::bad_alloc();
            }
            for(idx_t i=start; i-start < ndrows && i < data->nrows; ++i){
                idx_t rid = order[i].key;
                if(nnz + data->rptr[rid+1] - data->rptr[rid] > ninnz){
                    nr = i-start;
                    break;
                }
                val_t minsim = nbrs[rid][0].val;
                for(ptr_t j=data->rptr[rid]; j < data->rptr[rid+1]; ++j){
                    if(data->rnorms[j] >= minsim){
                        cptr[data->rind[j]+1]++; // build index ptr
                    } else {
                        rptr[i-start] = j;
                        break;
                    }
                }
            }
            if(nr == 0){
                nr = data->nrows - start;
            }
            // counts to column pointer
            for(idx_t j=0; j < data->ncols; ++j){
                cptr[j+1] += cptr[j];
            }
            nnz = cptr[data->ncols];
            // allocate index
            csr_t * idx = new csr_t();
            idx_t * rids = (idx_t *) malloc(sizeof(idx_t) * nr);
            idx_t * cind = (idx_t *) malloc(sizeof(idx_t) * nnz);
            val_t * cval = (val_t *) malloc(sizeof(val_t) * nnz);
            acm_t * cnorms = (acm_t *) malloc(sizeof(acm_t) * nnz);
            if(!rids || !cind || !cval || !cnorms){
                throw std::bad_alloc();
            }
            idx->rids = rids;
            idx->rptr = rptr;
            idx->cptr = cptr;
            idx->cind = cind;
            idx->cval = cval;
            idx->cnorms = cnorms;
            idx->nrows = nr;
            idx->ncols = data->ncols;
            // now transfer the data to the index
            for(idx_t i=start; i-start < nr; ++i){
                idx_t rid = order[i].key;
                for(ptr_t j=data->rptr[rid]; j < rptr[i-start]; ++j){
                    idx_t cid = data->rind[j];
                    ptr_t n = cptr[cid];
                    cind[n] = i-start;
                    cval[n] = data->rval[j];
                    cnorms[n] = data->rnorms[j];
                    rids[i-start] = rid;
                    cptr[cid]++;
                }
            }
            // finally, fix the column pointer
            for(idx_t j=data->ncols; j > 0; j--){
                cptr[j] = cptr[j-1];
            }
            cptr[0] = 0;

            return idx;
        }

        /**
         * Search for remaining neighbors for row `rid` against index `idx`.
         * 
         * @param rid Global row id for the object searching for
         * @param lrid Local row id in the index
         * @param minsim Minimum similarity among all objects in the index
         * @param ordermap Mapping from global row ids to local row ids for all objects in the index
         * @param idx Inverted index for a subset of the data, in non-decreasing neighborhood similarity order
         * @param tw Workspace for the thread
        */
        void _finalize_neighborhood(idx_t rid, idx_t lrid, val_t minsim, idx_t * ordermap, csr_t * idx, l2knnw * tw){
            auto bv = tw->bv;
            auto queue = tw->queue;
            auto sims = tw->sims;
            ivkv_t* qcands = tw->qcands;
            val_t* q = tw->q;
            acm_t* n = tw->n;
            acm_t ms = queue->see_top_val();
            // hash query row
            for(ptr_t j=data->rptr[rid]; j < data->rptr[rid+1]; ++j){
                idx_t cid = data->rind[j];
                q[cid] = data->rval[j];
                n[cid] = data->rnorms[j];
            }
            // mark current neighbors + rid
            bv->set(lrid);
            for(idx_t i=0; i < queue->nnodes; ++i){
                idx_t cid = queue->heap[i].key;
                if(ordermap[cid] >=0){
                    bv->set(ordermap[cid]);
                }
            }
            // find candidates
            idx_t ncands = 0;



            // unset current queue
            for(idx_t i=0; i < queue->nnodes; ++i){
                idx_t cid = queue->heap[i].key;
                if(ordermap[cid] >=0){
                    bv->zero(ordermap[cid]);
                }
            }
            bv->zero(lrid);


            // unhash query row
            for(ptr_t j=data->rptr[rid]; j < data->rptr[rid+1]; ++j){
                q[data->rind[j]] = 0.0;
            }

        }


        /**
         * Make the graph exact by finding all possible neighbors that could be in 
         * the top-k neighborhood for each node.
         * @param twork Thread work structure with reusable memory allocations
         * @param nbrs 2-D array of neighbors for each node
         * @param nnbrs Neighborhood sizes for each node neighborhood
         * @param maxcands Candidate list size
         * @param k Number of neighbors we ultimately want
        */
        csr_t * _make_exact(
            l2knnw** twork,
            ivkv_t** nbrs,
            idx_t* nnbrs,
            idx_t maxcands,
            idx_t k
        ){

            idx_t start = 0;

            /** structure to figure out processing order */
            ivkv_t * order = (ivkv_t *) malloc(sizeof(ivkv_t) * data->nrows);
            idx_t * ordermap = (idx_t *) malloc(sizeof(idx_t) * data->nrows);
            if(!order || !ordermap){
                throw std::bad_alloc();
            }
            for(idx_t i=0; i < data->nrows; ++i){
                order[i].key = i;
                ordermap[i] = -1;
            }

            while(start < data->nrows){
                /** figure out processing order */
                _set_processing_order(nbrs, nnbrs, order, start);
                idx_t minrid = order[start].key;
                val_t minsim = nbrs[minrid][0].val;

                /** create index for the next search block */
                csr_t * idx = _create_index(order, nbrs, start);

                /** map local order ids */
                for(idx_t i=start; i < start + idx->nrows; ++i){
                    idx_t rid = order[i].key;
                    ordermap[rid] = i - start; // local row id in block
                }

                /** execute search against the search block index */
                #pragma omp parallel 
                {
                    int const tid = omp_get_thread_num();
                    auto tw = twork[tid];

                    #pragma omp for schedule(static, 32)
                    for(idx_t i=start; i < data->nrows; ++i){
                        idx_t rid = order[i].key;
                        idx_t lrid = i - start;
                        // if less than k neighbors, cannot be enhanced; we're done
                        if(nnbrs[rid] < k){
                            continue;
                        }
                        tw->queue->reset(k, nnbrs[rid], nbrs[rid]);
                        _finalize_neighborhood(rid, lrid, minsim, ordermap, idx, tw);
                    }

                    /** update remote neighborhoods */


                }

                //free memmory
                delete idx;

                /** reset map local order ids */
                for(idx_t i=start; i < start + idx->nrows; ++i){
                    idx_t rid = order[i].key;
                    ordermap[rid] = -1;
                }


                start += idx->nrows;
            }

            // free memmory
            free(order);
            free(ordermap);

            return _graph_to_csr(data->nrows, nbrs, nnbrs);
        }

    };


} /* end namespace */

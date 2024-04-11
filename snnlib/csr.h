/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#pragma once

#include "basetypes.h"
#include "sort.h"

#include <omp.h>
#include <math.h>
#include <limits>
#include <string>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <algorithm>    // std::random_shuffle
#include <random>


namespace snnlib {

    enum CSR_BASE {
        ROW=1,   // row-based
        COL,     // column-based
        ROWCOL,  // row and column based
    };

    enum CSR_NORM {
        L1=1,
        L2
    };

    enum CSR_SCALE {
        MAXTF,   /** Maximum term frequency */
        SQRT,    /** Square root */
        LOG,     /** Log base e */
        LOG10,   /** Log base 10 */
        IDF,     /** Inverse document frequency */
    };

    enum CSR_FORMAT {
        DEDUCE=0, /** Figure out format from file extension */
        CSR,    /** Compressed Sparse Row */
        CLU,    /** Cluto - nrows ncols nnz + CSR */
        COO,    /** Coordinates == IJV, rid cid val, one non-zero per row */
        SMAT,   /** Matlab sparse matrix */
        BCOO,   /** Binary COO */
        BCSR,   /** Binary CSR */
        MTX,    /** MatrixMarket */
        NPZ     /** Numpy archive */
    };

    struct csr_t {
        static const size_t OMPMINOPS = 50000;
        idx_t nrows{0};
        idx_t ncols{0};
        idx_t maxnrows{0};
        ptr_t maxnnz{0};
        // CSR format data structures
        idx_t *rind{nullptr};   // the column id for each data point
        val_t *rval{nullptr};   // data pts
        ptr_t *rptr{nullptr};   // start each row
        idx_t *rids{nullptr};   // user-assigned row identifiers (e.g. original column number)
        idx_t *rperm{nullptr};  // row permutations
        acm_t *rnorms{nullptr}; // norms or square norms of each row or value (depends on algo)
        // CSC format for all or some of the data in the row-based structure
        idx_t *cind{nullptr};
        val_t *cval{nullptr};
        ptr_t *cptr{nullptr};
        acm_t *cnorms{nullptr}; // norms or square norms of each row or value (depends on algo)


        // Methods

        /**
         * Constructor of empty csr_t with option to add space for adding values
         * 
         * @param nrows Number of rows for the matrix
         * @param ncols Number of columns for the matrix
         * @param maxnrows Max number of rows for the matrix
         * @param maxnnz Max number of non-zeros for the matrix
        */
        csr_t(
            idx_t const nrows=0, 
            idx_t const ncols=0, 
            idx_t const maxnrows=0, 
            ptr_t const maxnnz=0):
            ncols(ncols), maxnrows(maxnrows), maxnnz(maxnnz),
            rind(nullptr), rval(nullptr), rptr(nullptr), rids(nullptr), rperm(nullptr), rnorms(nullptr),
            cind(nullptr), cval(nullptr), cptr(nullptr), cnorms(nullptr)
        {
            if (nrows > maxnrows){
                this->maxnrows = nrows;
            }
            alloc(this->maxnrows, this->maxnnz);
            this->nrows = nrows;
        }

        /**
         * Data constructor - row structure arrays provided as input
         * Note that data will be freed on object deletion, so it should not also
         * be freed elsewhere!
        */
        csr_t(
            idx_t nrows, idx_t ncols, ptr_t nnz, 
            idx_t * rind, val_t * rval, ptr_t * rptr,
            idx_t * rids=nullptr, idx_t * rperm=nullptr, acm_t * rnorms=nullptr):
            nrows(nrows), ncols(ncols), maxnrows(nrows), maxnnz(nnz), 
            rind(rind), rval(rval), rptr(rptr), rids(rids), rperm(rperm), rnorms(rnorms)
        {

        }

        /** 
         * Copy Contructor
         * Note that maxnnz and maxnrows is set to current nnz and nrows (shrinks to needed space).
         * Only the base row strucure is copied by default (rptr, rind, rval, rperm, and rids).
         * 
         * @param src   Source matrix that is being copied
        */
        csr_t(const csr_t& src)
        {
            alloc(src.nrows, src.nnz());
            nrows = src.nrows;
            ncols = src.ncols;
            if(maxnrows > 0){
                memcpy(rptr, src.rptr, sizeof(ptr_t) * (nrows + 1));
                if(src.rids){
                    memcpy(rids, src.rids, sizeof(idx_t) * nrows);
                }
                if(src.rperm){
                    memcpy(rperm, src.rperm, sizeof(idx_t) * nrows);
                }
            }
            if(maxnnz > 0){
                memcpy(rind, src.rind, sizeof(idx_t) * maxnnz);
                memcpy(rval, src.rval, sizeof(val_t) * maxnnz);
            }
        }

        /** 
         * Copy Contructor
         * Note that maxnnz and maxnrows is set to current nnz and nrows (shrinks to needed space).
         * Only the base row strucure is copied by default (rptr, rind, rval, rperm, and rids).
         * 
         * @param src   Source matrix that is being copied
        */
        csr_t(const csr_t* src)
        {
            alloc(src->nrows, src->nnz());
            nrows = src->nrows;
            ncols = src->ncols;
            if(maxnrows > 0){
                memcpy(rptr, src->rptr, sizeof(ptr_t) * (nrows + 1));
                if(src->rids){
                    memcpy(rids, src->rids, sizeof(idx_t) * nrows);
                }
                if(src->rperm){
                    memcpy(rperm, src->rperm, sizeof(idx_t) * nrows);
                }
            }
            if(maxnnz > 0){
                memcpy(rind, src->rind, sizeof(idx_t) * maxnnz);
                memcpy(rval, src->rval, sizeof(val_t) * maxnnz);
            }
        }

        /**
         * Fill the matrix with random values
        */
        static csr_t * random(
            idx_t const nrows, 
            idx_t const ncols, 
            double const factor=0.05)
        {
            ptr_t nnz = (ptr_t) (factor * nrows * ncols);
            if(nnz >= nrows * ncols / 2.0){
                throw std::runtime_error("Asking for too many non-zeros. Matrix is not sparse.");
            }
            auto mat = new csr_t(nrows, ncols, nrows, nnz);
            /* fill in ptr array; generate random row sizes */
            unsigned int seed = (unsigned long) mat;
            long double sum = 0;
            mat->rptr[0] = 0;
            for(idx_t i=1; i <= mat->nrows; ++i){
                mat->rptr[i] = rand_r(&seed) % ncols;
                sum += mat->rptr[i];
            }
            for(idx_t i=0; i < mat->nrows; ++i){
                double percent = mat->rptr[i+1] / sum;
                mat->rptr[i+1] = mat->rptr[i] + (ptr_t)(percent * nnz);
                if(mat->rptr[i+1] > nnz){
                    mat->rptr[i+1] = nnz;
                }
            }
            if(nnz - mat->rptr[mat->nrows-1] <= ncols){
                mat->rptr[mat->nrows] = nnz;
            }

            /* fill in indices and values with random numbers */
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                unsigned int seed = (unsigned long) mat * (1+tid);
                std::vector<int> perm;
                for(idx_t i=0; i < ncols; ++i){
                    perm.push_back(i);
                }
                std::random_device seeder;
                std::mt19937 engine(seeder());

                #pragma omp for
                for(idx_t i=0; i < nrows; ++i){
                    std::shuffle(perm.begin(), perm.end(), engine);
                    for(ptr_t j=mat->rptr[i]; j < mat->rptr[i+1]; ++j){
                        mat->rind[j] = perm[j - mat->rptr[i]];
                        mat->rval[j] = ((double) rand_r(&seed)/rand_r(&seed));
                    }
                }
            }

            return mat;
        }

            /**
         * Find the next power of 2
        */
        idx_t next_power_of_2(idx_t n){
            int value = 1;
            while(value<n){
                value=value << 1;
            }
            return value ;
        }


        /*! Adds a sparse vector to the matrix
            \param nvalues the number of values to be added
            \param values an array of the new values
            \param columns the column numbers (i.e. for rowind) for the new values **(0 indexed)**
        */
        void insert_row(ptr_t nvalues, val_t *values, idx_t *columns){
            auto nnz = this->nnz();
            
            if(nrows == maxnrows || nnz+nvalues > maxnnz){
                if(nrows == maxnrows){
                    maxnrows *= 2;
                }
                if(nvalues > (2*maxnnz - nnz)){
                    // values wouldn't fit even if the space were doubled
                    maxnnz = next_power_of_2(nnz+nvalues);
                } else {
                    maxnnz *= 2; 
                }
                this->alloc(maxnrows, maxnnz);
            }
            // copy the data over
            for(int i = 0; i<nvalues; i++){
                rval[nnz] = values[i];
                rind[nnz] = columns[i];
                if(columns[i] >= ncols){
                    ncols = columns[i]+1 ;
                }
                nnz++;
            }
            nrows++;
            rptr[nrows] = nnz;
        }

        /**
         * Number of non-zeros in the matrix.
        */
        inline ptr_t nnz() const {
            return rptr != nullptr ? rptr[nrows] : 0;
        }

        /**
         * Assign row ids for the matrix. If no row ids are passed in,
         * row ids are added in the range [0, nrows).
        */
        void assign_row_ids(idx_t * row_ids=nullptr, size_t sz=0) {
            if(row_ids == nullptr){
                if(!rids){
                    rids = (idx_t *) malloc(sizeof(idx_t) * nrows);
                    if (!rids){
                        throw std::bad_alloc();
                    }
                }
                for(idx_t i=0; i < nrows; ++i){
                    rids[i] = i;
                }
                return;
            }
            if (sz != (size_t) nrows){
                throw std::runtime_error("Error assigning row ids: Invalid size");
            }
            rids = row_ids;
        }

        /**
         * Allocate space for the matrix. Space may grow or shrink. However,
         * when shrinking, only unused space can be deleted.
         * 
         * @param maxnr Max number of rows
         * @param maxnz Max number of non-zeros
         * @param alloc_rids Whether to allocate rids array
         * @param alloc_rperm Whether to allocate rperm array
        */
        void alloc(
            idx_t const maxnr, 
            ptr_t const maxnz
        ){
            if(nrows > maxnr){
                throw std::runtime_error("Maxnr is too low. Can only shrink unused space.");
            }
            if(nrows > 0 && rptr && rptr[nrows] > maxnz){
                throw std::runtime_error("Maxnz is too low. Can only shrink unused space.");
            }
            if(maxnr > 0){
                ptr_t nnz = this->nnz();
                if(rptr){
                    rptr = (ptr_t *) realloc(rptr, sizeof(ptr_t) * (maxnr + 1));
                } else {
                    rptr = (ptr_t *) malloc(sizeof(ptr_t) * (maxnr + 1));
                }
                if (!rptr){
                    throw std::bad_alloc();
                }
                for(ptr_t j=nrows+1; j <= maxnr; ++j){
                    rptr[j] = nnz;
                }
                if(rids){
                    rids = (idx_t *) realloc(rids, sizeof(idx_t) * maxnr);
                    if (!rids){
                        throw std::bad_alloc();
                    }
                } 
                if(rperm){
                    rperm = (idx_t *) realloc(rperm, sizeof(idx_t) * maxnr);
                    if (!rperm){
                        throw std::bad_alloc();
                    }
                }
                maxnrows = maxnr;
            }
            if(maxnz > 0){
                if(rind){
                    rind = (idx_t *) realloc(rind, sizeof(idx_t) * maxnz);
                } else {
                    rind = (idx_t *) malloc(sizeof(idx_t) * maxnz);
                }
                if(rval){
                    rval = (val_t *) realloc(rval, sizeof(val_t) * maxnz);
                } else {
                    rval = (val_t *) malloc(sizeof(val_t) * maxnz);
                }
                if (!rval || !rind){
                    throw std::bad_alloc();
                }
                if(rnorms){
                    rnorms = (acm_t *) realloc(rnorms, sizeof(acm_t) * maxnz);
                    if (!rnorms){
                        throw std::bad_alloc();
                    }
                }
                maxnnz = maxnz;
            }
        }

        /**
         * Reduce matrix allocated space to currently used space for the row structure.
        */
        void shrink(){
            alloc(nrows, nnz());
        }

        /**
         * Normalize the rows/columns of the matrix
         * @param norm  Which norm to use when normalizing: L1 (CSR_NORM::L1), or L2 (CSR_NORM::L2)
         * @param what  Normalize rows (CSR_BASE::ROW), columns (CSR_BASE::COL), or both (CSR_BASE::ROWCOL)
        */
        void normalize(CSR_NORM const norm=CSR_NORM::L2, CSR_BASE const what=CSR_BASE::ROW){
            ssize_t i, j;
            idx_t n;
            ptr_t* ptr;
            val_t* val;

            if ((what & 1) && rval) {
                n   = nrows;
                ptr = rptr;
                val = rval;

                #pragma omp parallel if (ptr[n] > (ptr_t)OMPMINOPS) private(i,j)
                {

                    long double sum;
                    #pragma omp for schedule(static)
                    for (i=0; i<n; ++i) {
                        for (sum=0.0, j=ptr[i]; j<ptr[i+1]; ++j){
                            if (norm == 2){
                                sum += val[j] * val[j];
                            } else if (norm == 1){
                                sum += val[j] > 0 ? val[j] : -val[j];
                            }
                        }
                        if (sum > 0) {
                            if (norm == 2){
                                sum = 1.0/sqrt(sum);
                            } else if (norm == 1) {
                                sum = 1.0/sum;
                            }
                            for (j=ptr[i]; j<ptr[i+1]; ++j){
                                val[j] *= sum;
                            }

                        }
                    }
                }
            }

            if ((what & 2) && cval) {
                n   = ncols;
                ptr = cptr;
                val = cval;

                #pragma omp parallel if (ptr[n] > (ptr_t)OMPMINOPS) private(i,j)
                {
                    long double sum;
                    #pragma omp for schedule(static)
                    for (i=0; i<n; ++i) {
                        for (sum=0.0, j=ptr[i]; j<ptr[i+1]; ++j){
                            if (norm == 2){
                                sum += val[j] * val[j];
                            } else if (norm == 1){
                                sum += val[j] > 0 ? val[j] : -val[j];
                            }
                        }
                        if (sum > 0) {
                            if (norm == 2){
                                sum = 1.0/sqrt(sum);
                            } else if (norm == 1) {
                                sum = 1.0/sum;
                            }
                            for (j=ptr[i]; j<ptr[i+1]; ++j){
                                val[j] *= sum;
                            }
                        }
                    }
                }
            }

        }
            
        /*! Compacts the row-space of the matrix by removing empty rows.
        As a result of the compaction, the row numbers are renumbered.
        The compaction operation is done in place and only affects the row-based
        representation of the matrix.
        */
        void compact_rows(){
            ptr_t i,j;
            for(j=0,i=0; i<nrows; i++){
                rptr[j]=rptr[i];
                if(rids){
                    rids[j] = rids[i];
                }
                if(rptr[i+1]>rptr[i]){
                    j++;
                }
            }

            rptr[j]=rptr[nrows];
            nrows=j;
            rptr= (ptr_t*) realloc (rptr, (j+1)*sizeof(ptr_t));
        }

        /*! Compacts the column-space of the matrix by removing empty columns.
        As a result of the compaction, the column numbers are renumbered.
        The compaction operation is done in place and only affects the row-based
        representation of the matrix.
        */
        void compact_cols(){
            idx_t* clen = new idx_t[ncols];
            idx_t* cmap = new idx_t[ncols];
            idx_t nncols = 0;

            for (idx_t i=0; i<ncols; i++){
                clen[i] = 0;
                cmap[i] = -1;
            }

            for(ptr_t j=0; j<rptr[nrows]; j++){
                clen[rind[j]]++;
            }

            for(idx_t i=0; i<ncols; i++){
                if(clen[i]>0){
                    cmap[i] = nncols++;
                }
            }

            for(ptr_t j=0; j<rptr[nrows]; j++){
                rind[j] = cmap[rind[j]];
            }
            ncols=nncols;
            delete []cmap;
            delete []clen;
        }


        /**
         * Compute prefix norms for each nnz in the row-based structure.
         * This method assumes rows have already been normalized.
        */
        void compute_prefix_norms(){
            //compute norms
            if(rnorms){
                free(rnorms);
            }
            rnorms = (acm_t *) malloc(sizeof(acm_t) * rptr[nrows]);
            if(!rnorms){
                throw std::bad_alloc();
            }
            #pragma omp for schedule(dynamic, 128)
            for(idx_t i = 0; i < nrows; i++){
                if(rptr[i] == rptr[i+1]){
                    continue;
                }
                ptr_t j = rptr[i];
                long double b2 = 0.0;
                rnorms[j] = 0.0;
                for(j = rptr[i] + 1; j < rptr[i+1]; j++){
                    b2 += (long double) rval[j-1] * rval[j-1];
                    rnorms[j] = sqrt(b2);
                }
            }
        }


        /**
         * Compute suffix norms for each nnz in the row-based structure.
         * This method assumes rows have already been normalized.
        */
        void compute_suffix_norms(){
            //compute norms
            if(rnorms){
                free(rnorms);
            }
            rnorms = (acm_t *) malloc(sizeof(acm_t) * rptr[nrows]);
            if(!rnorms){
                throw std::bad_alloc();
            }
            #pragma omp for schedule(dynamic, 128)
            for(idx_t i = 0; i < nrows; i++){
                long double b2 = 1.0;
                for(ptr_t j = rptr[i]; j < rptr[i+1]; ++j){
                    b2 -= rval[j] * rval[j];
                    rnorms[j] = b2 > 0.0 ? sqrt(b2) : 0.0;
                }
            }
        }

        /**
         * Sort non-zerows
         * @param what Which structure to sort, row structure (CSR_BASE::ROW) or column structure (CSR_BASE::COL)
         * @param ascending Sort in ascending or decreasing order
        */
        void sort_nnzs_by_ind(const bool ascending=true, const CSR_BASE what=CSR_BASE::ROW)
        {
            ptr_t n=0, nn=0;
            ptr_t *ptr;
            idx_t *ind;
            val_t *val;
            acm_t *nrm = nullptr;

            switch (what) {
            case CSR_BASE::ROW:
                if (!rptr){
                    throw std::runtime_error("Row-based view of the matrix does not exists.");
                }
                n   = nrows;
                ptr = rptr;
                ind = rind;
                val = rval;
                if(rnorms){
                    nrm = rnorms;
                }
                break;

            case CSR_BASE::COL:
                if (!cptr){
                    throw std::runtime_error("Column-based view of the matrix does not exists.");
                }
                n   = ncols;
                ptr = cptr;
                ind = cind;
                val = cval;
                if(cnorms){
                    nrm = cnorms;
                }
                break;

            default:
                throw std::runtime_error("Invalid index type/what to sort.");
            }

            nn = ptr[1] - ptr[0];
            for (idx_t i=1; i<n; ++i){
                if(ptr[i+1] - ptr[i] > nn){
                    nn = ptr[i+1] - ptr[i];
                }
            }

            #pragma omp parallel if (n > (ptr_t)OMPMINOPS)
            {
                ssize_t i, j, jj, k;
                pikv_t* cand;
                val_t* tval = nullptr;
                acm_t* tnrm = nullptr;

                cand = (pikv_t*) malloc(nn * sizeof(pikv_t));
                if(!cand){
                    throw std::bad_alloc();
                }
                if(val){
                    tval = (val_t*) malloc(nn * sizeof(val_t));
                    if(!tval){
                        throw std::bad_alloc();
                    }
                    if(nrm){
                        tnrm = (acm_t*) malloc(nn * sizeof(acm_t));
                        if(!tnrm){
                            throw std::bad_alloc();
                        }    
                    }

                    #pragma omp for schedule(static)
                    for (i=0; i<n; ++i) {
                        for (k=0, j=ptr[i]; j<ptr[i+1]; ++j) {
                            if (j > ptr[i] && 
                                ((ascending && ind[j] < ind[j-1]) ||
                                (!ascending && ind[j] > ind[j-1]))
                            ){
                                k = 1; /* an inversion */
                            }
                            jj = j - ptr[i];
                            cand[jj].key = jj;
                            cand[jj].val = ind[j];
                            tval[jj]     = val[j];
                            if(nrm){
                                tnrm[jj] = nrm[j];
                            }
                        }
                        if (k) {
                            sort(ptr[i+1]-ptr[i], cand, ascending);
                            for (j=ptr[i]; j<ptr[i+1]; ++j) {
                                jj = j - ptr[i];
                                ind[j] = cand[jj].val;
                                val[j] = tval[cand[jj].key];
                                if(nrm){
                                    nrm[j] = tnrm[cand[jj].key];
                                }
                            }
                        }
                    }

                    free(tval);
                    if(tnrm){
                        free(tnrm);
                    }
                } else {
                    // assume no norms either since no vals
                    #pragma omp for schedule(static)
                    for (i=0; i<n; ++i) {
                        for (k=0, j=ptr[i]; j<ptr[i+1]; ++j) {
                            if (j > ptr[i] && 
                                ((ascending && ind[j] < ind[j-1]) ||
                                (!ascending && ind[j] > ind[j-1]))
                            ){
                                k = 1; /* an inversion */
                            }
                            jj = j - ptr[i];
                            cand[jj].key = jj;
                            cand[jj].val = ind[j];
                        }
                        if (k) {
                            sort(ptr[i+1]-ptr[i], cand, ascending);
                            for (j=ptr[i]; j<ptr[i+1]; ++j){
                                ind[j] = cand[j-ptr[i]].val;
                            }
                        }
                    }
                }
                free(cand);
            }
        }


        /*************************************************************************/
        /*! Sorts the values (and associated indices) in ascending or decreasing order
         * @param ascending Sort in ascending or decreasing order
         * @param what Which structure to sort, row structure (CSR_BASE::ROW) or column structure (CSR_BASE::COL)
         * @param nselect If `nselect` > 0, first select the highest/lowest `nselect` values
         *                and then only sort those top `nselect` values
         * 
         * @TODO Replace with kvsortselect to avoid data copy
        */
        void sort_nnzs_by_value(const bool ascending=true, const CSR_BASE what=CSR_BASE::ROW, const idx_t nselect=0)
        {
            ptr_t ii, n, nn;
            ptr_t* ptr;
            idx_t* ind;
            val_t* val;
            acm_t* nrm = nullptr;

            switch (what) {
            case CSR_BASE::ROW:
                if (!rptr){
                    throw std::runtime_error("Row-based view of the matrix does not exists.");
                }

                n   = nrows;
                ptr = rptr;
                ind = rind;
                val = rval;
                if(rnorms){
                    nrm = rnorms;
                }
                break;

            case CSR_BASE::COL:
                if (!cptr){
                    throw std::runtime_error("Column-based view of the matrix does not exists.");
                }

                n   = ncols;
                ptr = cptr;
                ind = cind;
                val = cval;
                if(cnorms){
                    nrm = cnorms;
                }
                break;

            default:
                throw std::runtime_error("Invalid index type/what to sort.");
            }

            if(!val){
                throw std::runtime_error("Values not present in the requested structure.");
            }

            for (nn=0, ii=0; ii<n; ++ii){
                if(ptr[ii+1] - ptr[ii] > nn){
                    nn = ptr[ii+1] - ptr[ii];
                }
            }

            #pragma omp parallel if (n > (ptr_t)OMPMINOPS)
            {
                ssize_t i, j, jj, k;
                ivkv_t* cand;
                idx_t* tind;
                acm_t* tnrm = nullptr;

                cand = (ivkv_t*) malloc(nn * sizeof(ivkv_t));
                if(!cand){
                    throw std::bad_alloc();
                }

                tind = (idx_t*) malloc(nn * sizeof(idx_t));
                if(!tind){
                    throw std::bad_alloc();
                }
                if(nrm){
                    tnrm = (acm_t*) malloc(nn * sizeof(acm_t));
                    if(!tnrm){
                        throw std::bad_alloc();
                    }    
                }

                #pragma omp for schedule(static)
                for (i=0; i<n; ++i) {
                    for (k=0, j=ptr[i]; j<ptr[i+1]; ++j) {
                        if (j > ptr[i] && 
                            ((ascending && val[j] < val[j-1]) ||
                            (!ascending && val[j] > val[j-1]))
                        ){
                            k = 1; /* an inversion */
                        }
                        jj = j - ptr[i];
                        cand[jj].key = jj;
                        cand[jj].val = val[j];
                        tind[jj]     = ind[j];
                        if(nrm){
                            tnrm[jj] = nrm[j];
                        }
                    }
                    if (k) {
                        if(nselect > 0 && nselect < ptr[i+1]-ptr[i]){
                            selectsort(ptr[i+1]-ptr[i], cand, nselect, ascending);
                        } else {
                            sort(ptr[i+1]-ptr[i], cand, ascending);
                        }
                        for (j=ptr[i]; j<ptr[i+1]; ++j) {
                            jj = j - ptr[i];
                            val[j] = cand[jj].val;
                            ind[j] = tind[cand[jj].key];
                            if(nrm){
                                nrm[j] = tnrm[cand[jj].key];
                            }
                        }
                    }
                }

                free(tind);
                if(tnrm){
                    free(tnrm);
                }
                free(cand);
            }
        }

        /**
         * Permute the column ids in the row-structure according to the given permutation
         * @param perm  The permutation such that the column at index perm[j] will move to index j
        */
        void permute_columns(idx_t const * const perm){
            for (ptr_t j = 0; j<rptr[nrows]; j++){
                rind[j] = perm[rind[j]];  
            }
        }

        /**
         * Permute the rows according to the given permutation
         * @param perm  The permutation such that row at index perm[i] will move to index i
        */
        void permute_rows(idx_t const * const perm){
            ptr_t nnz = rptr[nrows];
            ptr_t * ptr = (ptr_t*) malloc((nrows+1) * sizeof(ptr_t));
            idx_t * ind = (idx_t*) malloc(nnz * sizeof(idx_t));
            val_t * val = nullptr;
            idx_t * ids = nullptr;
            acm_t * norms = nullptr;
            if(!ptr || !ind){
                throw std::bad_alloc();
            }
            if(rval){
                val = (val_t*) malloc(nnz * sizeof(val_t));
                if(!val){
                    throw std::bad_alloc();
                }
            }
            if(rids){
                ids = (idx_t*) malloc(nrows * sizeof(idx_t));
                if(!ids){
                    throw std::bad_alloc();
                }
            }
            if(rnorms){
                norms = (acm_t*) malloc(nnz * sizeof(acm_t));
                if(!norms){
                    throw std::bad_alloc();
                }
            }
            ptr[0] = 0;
            for(idx_t i=0; i < nrows; ++i){
                idx_t r = perm[i];
                for(ptr_t j=0; j < rptr[r+1]-rptr[r]; ++j){
                    ind[ptr[i] + j] = rind[rptr[r] + j];
                    if(rval){
                        val[ptr[i] + j] = rval[rptr[r] + j];
                    }
                    if(rnorms){
                        norms[ptr[i] + j] = rnorms[rptr[r] + j];
                    }
                }
                if(rids){
                    ids[i] = rids[r];
                }
                ptr[i+1] = ptr[i] + rptr[r+1] - rptr[r];
            }
            // swap row structure
            free(rptr);
            rptr = ptr;
            free(rind);
            rind = ind;
            if(rval){
                free(rval);
                rval = val;
            }
            if(rids){
                free(rids);
                rids = ids;
            }
            if(rnorms){
                free(rnorms);
                rnorms = norms;
            }
            // store permutation array
            if(rperm){
                free(rperm);
            }
            rperm = (idx_t*) malloc(nrows * sizeof(idx_t));
            if(!rperm){
                throw std::bad_alloc();
            }
            memmove(rperm, perm, nrows * sizeof(idx_t));
        }

        /**
         * Reverse the permutation currently stored in rperm
        */
        void inverse_permute_rows(){
            if(!rperm){
                throw std::runtime_error("Row permutation does not exist.");
            }
            // allocate memory
            ptr_t nnz = rptr[nrows];
            idx_t * iperm = (idx_t*) malloc(nrows * sizeof(idx_t)); // inverse permutation
            ptr_t * ptr = (ptr_t*) malloc((nrows+1) * sizeof(ptr_t));
            idx_t * ind = (idx_t*) malloc(nnz * sizeof(idx_t));
            val_t * val = nullptr;
            idx_t * ids = nullptr;
            acm_t * norms = nullptr;
            if(!ptr || !ind){
                throw std::bad_alloc();
            }
            if(rval){
                val = (val_t*) malloc(nnz * sizeof(val_t));
                if(!val){
                    throw std::bad_alloc();
                }
            }
            if(rids){
                ids = (idx_t*) malloc(nrows * sizeof(idx_t));
                if(!ids){
                    throw std::bad_alloc();
                }
            }
            if(rnorms){
                norms = (acm_t*) malloc(nnz * sizeof(acm_t));
                if(!norms){
                    throw std::bad_alloc();
                }
            }
            // find inverse permutation
            for(idx_t i=0; i < nrows; ++i){
                iperm[rperm[i]] = i;
            }
            // transfer data
            ptr[0] = 0;
            for(idx_t i=0; i < nrows; ++i){
                idx_t ii = iperm[i];
                ptr[i+1] = ptr[i] + rptr[ii+1] - rptr[ii];
                for(ptr_t j=rptr[ii]; j < rptr[ii]; ++j){
                    ptr_t jj = rptr[ii] + j - ptr[i];
                    ind[j] = rind[jj];
                    if(rval){
                        val[j] = rval[jj];
                    }
                    if(rnorms){
                        norms[j] = rnorms[jj];
                    }
                }
            }
            // swap row structure
            free(rptr);
            rptr = ptr;
            free(rind);
            rind = ind;
            if(rval){
                free(rval);
                rval = val;
            }
            if(rids){
                free(rids);
                rids = ids;
            }
            if(rnorms){
                free(rnorms);
                rnorms = norms;
            }
        }

        /**
         * Updates or creates the opposite index (CSC vs. CSR) from the one specified
         * The reference dimension is the one indicated by `from`.
         * The other dimension is updated based on the reference.
         * @param from Structure to construct index from, ROW (1) or COL (2)
        */
        void create_index(const CSR_BASE from=CSR_BASE::ROW) {

            ptr_t *ref_ptr, *new_ptr;
            idx_t *ref_ind, *new_ind;
            val_t *ref_val, *new_val;
            acm_t *ref_nrm, *new_nrm;
            idx_t ref_n, new_n;
            ptr_t nnz;


            if (from == CSR_BASE::ROW) {
                ref_n = nrows;
                if (rptr && rind && rval) {
                    ref_ptr = rptr;
                    ref_ind = rind;
                    ref_val = rval;
                    ref_nrm = rnorms;
                    nnz = rptr[nrows];
                } else {
                    throw std::logic_error("Could not find CSR format data to convert to CSC.");
                }
                if (cptr) {
                    free(cptr);
                }
                if (cind) {
                    free(cind);
                }
                if (cval) {
                    free(cval);
                }
                if (cnorms) {
                    free(cnorms);
                }

                new_n = ncols;
                cptr = (ptr_t *) calloc(sizeof(ptr_t), (new_n + 1));
                cval = (val_t *) malloc(sizeof(val_t) * nnz);
                cind = (idx_t *) malloc(sizeof(idx_t) * nnz);
                if (!cptr || !cval || !cind) {
                    throw std::bad_alloc();
                }
                if (ref_nrm) {
                    cnorms = (acm_t *) malloc(sizeof(acm_t) * nnz);
                    if (!cnorms) {
                        throw std::bad_alloc();
                    }
                }
                new_ptr = cptr;
                new_ind = cind;
                new_val = cval;
                new_nrm = cnorms;

            } else if (from == CSR_BASE::COL) {
                ref_n = ncols;
                if (cptr && cval && cind) {
                    ref_ptr = cptr;
                    ref_ind = cind;
                    ref_val = cval;
                    ref_nrm = cnorms;
                    nnz = cptr[ncols];
                } else {
                    throw std::logic_error("Could not find CSC format data to convert to CSR.");
                }
                if (rptr) {
                    free(rptr);
                }
                if (rind) {
                    free(rind);
                }
                if (rval) {
                    free(rval);
                }
                if(rnorms){
                    free(rnorms);
                }
                new_n = nrows;
                rptr = (ptr_t *) calloc(sizeof(ptr_t), (new_n + 1));
                rval = (val_t *) malloc(sizeof(val_t) * nnz);
                rind = (idx_t *) malloc(sizeof(idx_t) * nnz);
                if (!rptr || !rval || !rind) {
                    throw std::bad_alloc();
                }
                if (ref_nrm) {
                    rnorms = (acm_t *) malloc(sizeof(acm_t) * nnz);
                    if (!rnorms) {
                        throw std::bad_alloc();
                    }
                }
                new_ptr = rptr;
                new_ind = rind;
                new_val = rval;
                new_nrm = rnorms;
            } else {
                throw std::logic_error("Dimension 'from' must be ROW (CSR_BASE::ROW) or COL (CSR_BASE::COL).");
            }

            // update the new ptr array - will be the # of vals per row/col not cumulative sum
            for (idx_t i = 0; i < ref_n; i++) {
                for (ptr_t j = ref_ptr[i]; j < ref_ptr[i + 1]; j++) {
                    new_ptr[ref_ind[j]]++;
                }
            }
            // prefix sum
            for (idx_t i = 1; i < new_n; i++) {
                new_ptr[i] += new_ptr[i - 1];
            }
            // shift
            csr_shift(new_n, new_ptr);

            // dif for sparse vs dense
            // if 6 or more nz per row
            if (new_ptr[new_n] > 6 * new_n) {
                for (idx_t i=0; i < ref_n; ++i) {
                    for (ptr_t j = ref_ptr[i]; j < ref_ptr[i+1]; ++j){
                        new_ind[new_ptr[ref_ind[j]]++] = i;
                    }
                }
                csr_shift(new_n, new_ptr);
                if (ref_val) {
                    for (idx_t i=0; i<ref_n; ++i) {
                        for (ptr_t j=ref_ptr[i]; j<ref_ptr[i+1]; ++j){
                            if(new_nrm){
                                new_nrm[new_ptr[ref_ind[j]]] = ref_nrm[j];
                            }
                            new_val[new_ptr[ref_ind[j]]++] = ref_val[j];
                        }
                    }
                    csr_shift(new_n, new_ptr);
                }
            }
            else {
                if (ref_val) {
                    for (idx_t i=0; i < ref_n; ++i) { 
                        for (ptr_t j=ref_ptr[i]; j < ref_ptr[i+1]; ++j) {
                            idx_t k = ref_ind[j];
                            new_ind[new_ptr[k]] = i;
                            if(new_nrm){
                                new_nrm[new_ptr[k]] = ref_nrm[j];
                            }
                            new_val[new_ptr[k]++] = ref_val[j];
                        }
                    }
                }
                else {
                    for (idx_t i=0; i< ref_n; ++i) {
                        for (ptr_t j=ref_ptr[i]; j<ref_ptr[i+1]; ++j){
                            new_ind[new_ptr[ref_ind[j]]++] = i;
                        }
                    }
                }
                csr_shift(new_n, new_ptr);
            }
        }



        /** 
         * Extract submatrix from the given matrix that starts at row `rstart` and continues for `numrows` rows.
         * @param rstart First row of the extracted matrix
         * @param numrows Max number of rows to extract from the matrix
        */
        csr_t * extract_submatrix(idx_t const rstart, idx_t numrows) {

            if (!this->rptr){
                throw std::runtime_error("SparseMatrix::extract_submatrix: error input rowptr NULL");
            }

            csr_t* nmat;

            if (rstart > this->nrows){
                return new csr_t();
            }
            if(numrows > nrows - rstart){
                numrows = nrows - rstart;
            }

            ptr_t nnz = rptr[rstart + numrows] - rptr[rstart];
            nmat = new csr_t(numrows, ncols, numrows, nnz);

            /* copy the row structure */
            for (idx_t i = 0, j=rstart; i <= numrows; i++, j++){
                nmat->rptr[i] = rptr[j] - rptr[rstart];
            }
            memmove(nmat->rind, this->rind + this->rptr[rstart], nnz * sizeof(idx_t));
            memmove(nmat->rval, this->rval + this->rptr[rstart], nnz * sizeof(val_t));
            if(rnorms){
                nmat->rnorms = (acm_t*) malloc(nnz * sizeof(acm_t));
                if(!nmat->rnorms){
                    throw std::bad_alloc();
                }
                memmove(nmat->rnorms, this->rnorms + this->rptr[rstart], nnz * sizeof(acm_t));
            }
            if(rids){
                nmat->rids = (idx_t*) malloc(numrows * sizeof(idx_t));
                if(!nmat->rids){
                    throw std::bad_alloc();
                }
                memmove(nmat->rids, this->rids + rstart, numrows * sizeof(idx_t));
            }

            return nmat;
        }



        /** 
         * Extract submatrix from the given matrix based on rows in list.
         * @param n Number of rows to extract
         * @param list List of row ids that should be extracted
        */
        csr_t * extract_rows(idx_t const n, idx_t const * const list) {

            if (!this->rptr){
                throw std::runtime_error("SparseMatrix::extract_rows: error input rowptr NULL");
            }

            csr_t* nmat;

            ptr_t nnz = 0;
            for(idx_t i=0; i < n; ++i){
                nnz += this->rptr[list[i]+1] - this->rptr[list[i]];
            }

            nmat = new csr_t(n, ncols, n, nnz);
            if(rnorms){
                nmat->rnorms = (acm_t*) malloc(nnz * sizeof(acm_t));
                if(!nmat->rnorms){
                    throw std::bad_alloc();
                }
            }
            nmat->rids = (idx_t*) malloc(n * sizeof(idx_t));
            if(!nmat->rids){
                throw std::bad_alloc();
            }

            /* copy the row structure */
            nnz = 0;
            nmat->rptr[0] = 0;
            for (idx_t i = 0; i < n; i++){
                idx_t rid = list[i];
                ptr_t ln = this->rptr[rid+1] - this->rptr[rid];
                memmove(nmat->rind+nnz, this->rind + this->rptr[rid], ln * sizeof(idx_t));
                memmove(nmat->rval+nnz, this->rval + this->rptr[rid], ln * sizeof(val_t));
                if(rnorms){
                    memmove(nmat->rnorms+nnz, this->rnorms + this->rptr[rid], ln * sizeof(acm_t));
                }
                nmat->rids[i] = this->rids ? this->rids[rid] : rid;

                nnz += ln;
                nmat->rptr[i+1] = nnz;
            }

            return nmat;
        }

        /** 
         * Extract submatrix from the given matrix based on rows in list.
         * @param n Number of rows to extract
         * @param list List of key-value pairs containing row ids that should be extracted as keys
        */
        csr_t * extract_rows(idx_t const n, iakv_t const * const list) {

            if (!this->rptr){
                throw std::runtime_error("SparseMatrix::extract_rows: error input rowptr NULL");
            }

            csr_t* nmat;

            ptr_t nnz = 0;
            for(idx_t i=0; i < n; ++i){
                nnz += this->rptr[list[i].key+1] - this->rptr[list[i].key];
            }

            nmat = new csr_t(n, ncols, n, nnz);
            if(rnorms){
                nmat->rnorms = (acm_t*) malloc(nnz * sizeof(acm_t));
                if(!nmat->rnorms){
                    throw std::bad_alloc();
                }
            }
            nmat->rids = (idx_t*) malloc(n * sizeof(idx_t));
            if(!nmat->rids){
                throw std::bad_alloc();
            }

            /* copy the row structure */
            nnz = 0;
            nmat->rptr[0] = 0;
            for (idx_t i = 0; i < n; i++){
                idx_t rid = list[i].key;
                ptr_t ln = this->rptr[rid+1] - this->rptr[rid];
                memmove(nmat->rind+nnz, this->rind + this->rptr[rid], ln * sizeof(idx_t));
                memmove(nmat->rval+nnz, this->rval + this->rptr[rid], ln * sizeof(val_t));
                if(rnorms){
                    memmove(nmat->rnorms+nnz, this->rnorms + this->rptr[rid], ln * sizeof(acm_t));
                }
                nmat->rids[i] = this->rids ? this->rids[rid] : rid;

                nnz += ln;
                nmat->rptr[i+1] = nnz;
            }

            return nmat;
        }


        void print_info() const {
            std::cout << "Number of rows: " << nrows << std::endl;
            std::cout << "Number of columns: " << ncols << std::endl;
            std::cout << "Number of non-zeros: " << (rptr ? rptr[nrows] : 0) << std::endl;
            std::cout << "Maximum rows: " << maxnrows << std::endl;
            std::cout << "Maximum non-zeroes: " << maxnnz << std::endl;
        }

        /**
         * Print the row structure in CSR format
        */
        void print() const {
            for(idx_t i=0; i < nrows; ++i){
                for(ptr_t j=rptr[i]; j < rptr[i+1]; ++j){
                    std::cout << rind[j] << " " << rval[j] << " ";
                }
                std:: cout << std::endl;
            }
            std:: cout << std::endl;
        }


        /**
         * Destructor for csr_t structure
        */
        ~csr_t() {
            if(rids){
                free(rids);
            }
            if(rptr){
                free(rptr);
            }
            if(rind){
                free(rind);
            }
            if(rval){
                free(rval);             
		    }
            if(rperm){
		    	free(rperm);
            }
            if(rnorms){
            	free(rnorms);
            }
            if(cptr){
                free(cptr);
            }
            if(cind){
                free(cind);
            }
            if(cval){
                free(cval);
            }
            if(cnorms){
                free(cnorms);
            }
        }

        /**
         * Shifts the arrays of an element right by one and sets the first element to 0
         * @tparam T    type of array being shifted
         * @param n     size of array
         */
        template <typename T>
        static void csr_shift(idx_t n, T * arr) {
            for (idx_t i=n; i > 0; i--) {
                arr[i] = arr[i-1];
            }
            arr[0] = 0;
        }


    };

} /* end namespace */
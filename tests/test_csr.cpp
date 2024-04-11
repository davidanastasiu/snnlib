/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#include <iostream>     // std::cout
#include <limits>
#include <gtest/gtest.h>

//header file
#include "../snnlib/csr.h"

using namespace snnlib;



/**
 * Test creation
 */
TEST(csr, creation){
    csr_t * mat = new csr_t();
    EXPECT_EQ(mat->nrows, 0);
    EXPECT_EQ(mat->ncols, 0);
    EXPECT_EQ(mat->rptr, nullptr);
    EXPECT_EQ(mat->rind, nullptr);
    EXPECT_EQ(mat->rval, nullptr);
    EXPECT_EQ(mat->rids, nullptr);
    EXPECT_EQ(mat->rnorms, nullptr);
    EXPECT_EQ(mat->cptr, nullptr);
    EXPECT_EQ(mat->cind, nullptr);
    EXPECT_EQ(mat->cptr, nullptr);
    EXPECT_EQ(mat->cnorms, nullptr);
    delete mat;

    mat = new csr_t(10, 20, 10, 200);
    EXPECT_EQ(mat->nrows, 10);
    EXPECT_EQ(mat->ncols, 20);
    EXPECT_NE(mat->rptr, nullptr);
    EXPECT_NE(mat->rind, nullptr);
    EXPECT_NE(mat->rval, nullptr);
    EXPECT_EQ(mat->nnz(), 0);
    delete mat;

    mat = csr_t::random(103, 1387, 0.05);
    EXPECT_EQ(mat->nrows, 103);
    EXPECT_EQ(mat->ncols, 1387);
    EXPECT_GE(mat->nnz(), 103*1387*0.049);
    EXPECT_LE(mat->nnz(), 103*1387*0.051);
    delete mat;
}

/**
 * Test copy
 */
TEST(csr, copy){
    auto mat = csr_t::random(103, 1387, 0.05);
    auto mat2 = new csr_t(mat);
    EXPECT_NE(mat2, mat);
    EXPECT_EQ(mat2->nrows, 103);
    EXPECT_EQ(mat2->ncols, 1387);
    EXPECT_GE(mat2->nnz(), mat->nnz());
    for(idx_t i=0; i <= mat->nrows; ++i){
        EXPECT_EQ(mat->rptr[i], mat2->rptr[i]);
    }
    for(ptr_t j=0; j < mat->nnz(); ++j){
        EXPECT_EQ(mat->rind[j], mat2->rind[j]);
        EXPECT_FLOAT_EQ(mat->rval[j], mat2->rval[j]);
    }
    delete mat;
    delete mat2;
}

TEST(csr, assign_row_ids){
    auto mat = new csr_t(10, 20, 10, 200);
    mat->assign_row_ids();
    EXPECT_NE(mat->rids, nullptr);
    for(idx_t i=0; i < mat->nrows; ++i){
        EXPECT_EQ(mat->rids[i], i);
    }
    delete mat;
}


TEST(csr, realloc){
    auto mat = csr_t::random(10, 138, 0.05);
    EXPECT_EQ(mat->nrows, 10);
    EXPECT_EQ(mat->maxnrows, 10);
    EXPECT_GE(mat->nnz(), 10*138*0.049);
    EXPECT_LE(mat->nnz(), 10*138*0.051);
    EXPECT_EQ(mat->maxnnz, mat->nnz());
    mat->alloc(100, 1000);
    EXPECT_EQ(mat->nrows, 10);
    EXPECT_EQ(mat->maxnrows, 100);
    EXPECT_EQ(mat->maxnnz, 1000);
    mat->shrink();
    EXPECT_EQ(mat->nrows, 10);
    EXPECT_EQ(mat->maxnrows, 10);
    EXPECT_EQ(mat->maxnnz, mat->nnz());
    delete mat;
}


TEST(csr, create_index){
    auto mat = csr_t::random(10, 138, 0.05);
    mat->normalize();
    mat->compute_prefix_norms();
    mat->create_index(CSR_BASE::ROW);
    long double total = 0;
    long double total2 = 0;
    long double ntotal = 0;
    long double ntotal2 = 0;
    auto nnz = mat->nnz();
    EXPECT_NE(mat->cptr, nullptr);
    EXPECT_NE(mat->cind, nullptr);
    EXPECT_NE(mat->cval, nullptr);
    EXPECT_EQ(nnz, mat->cptr[mat->ncols]);
    for(ptr_t j=0; j < nnz; ++j){
        total += mat->rval[j];
        ntotal += mat->rnorms[j];
    }
    for(ptr_t j=0; j < nnz; ++j){
        total2 += mat->cval[j];
        ntotal2 += mat->cnorms[j];
    }
    EXPECT_FLOAT_EQ(total, total2);
    EXPECT_FLOAT_EQ(ntotal, ntotal2);
    mat->create_index(CSR_BASE::COL);
    total = ntotal = 0;
    for(ptr_t j=0; j < nnz; ++j){
        total += mat->rval[j];
        ntotal += mat->rnorms[j];
    }
    EXPECT_FLOAT_EQ(total, total2);
    EXPECT_FLOAT_EQ(ntotal, ntotal2);
    delete mat;
}


TEST(csr, normalize){
    auto mat = csr_t::random(103, 1381, 0.05);
    mat->normalize();
    for(idx_t i=0; i < mat->nrows; ++i){
        if(mat->rptr[i+1] == mat->rptr[i]){
            continue;
        }
        acm_t norm = 0;
        for(ptr_t j=mat->rptr[i]; j < mat->rptr[i+1]; ++j){
            norm += mat->rval[j] * mat->rval[j];
        }
        EXPECT_FLOAT_EQ(norm, 1.0);
    }
    delete mat;
}


TEST(csr, compute_prefix_norms){
    auto mat = csr_t::random(103, 1381, 0.05);
    mat->normalize();
    mat->compute_prefix_norms();
    for(idx_t i=0; i < mat->nrows; ++i){
        if(mat->rptr[i+1]-mat->rptr[i] == 0){
            continue;
        }
        for(ptr_t j=mat->rptr[i]; j < mat->rptr[i+1]; ++j){
            EXPECT_LE(mat->rnorms[j], 1.0001);
        }
        auto j = mat->rptr[i+1]-1;
        auto norm = sqrt(mat->rnorms[j]*mat->rnorms[j] + mat->rval[j] * mat->rval[j]);
        ASSERT_FLOAT_EQ(norm, 1.0);
    }
    delete mat;
}


TEST(csr, compute_suffix_norms){
    auto mat = csr_t::random(103, 1381, 0.05);
    mat->normalize();
    mat->compute_suffix_norms();
    for(idx_t i=0; i < mat->nrows; ++i){
        if(mat->rptr[i+1]-mat->rptr[i] == 0){
            continue;
        }
        for(ptr_t j=mat->rptr[i]; j < mat->rptr[i+1]; ++j){
            EXPECT_LE(mat->rnorms[j], 1.0);
            EXPECT_GE(mat->rnorms[j], 0.0);
        }
    }
    delete mat;
}


TEST(csr, extract_submatrix){
    idx_t rs = 13;
    idx_t nr = 37;
    auto mat = csr_t::random(103, 1381, 0.05);
    mat->assign_row_ids();
    mat->normalize();
    mat->compute_prefix_norms();
    auto smat = mat->extract_submatrix(rs, nr);
    for(idx_t i=rs, ii=0; ii < nr; ++i, ++ii){
        EXPECT_EQ(mat->rptr[i+1]-mat->rptr[i], smat->rptr[ii+1]-smat->rptr[ii]);
        EXPECT_EQ(smat->rids[ii], i);
        if(mat->rptr[i+1]-mat->rptr[i] == 0){
            continue;
        }
        acm_t norm = 0.0;
        EXPECT_FLOAT_EQ(smat->rnorms[smat->rptr[ii]], 0.0);
        for(ptr_t j=mat->rptr[i], jj=smat->rptr[ii]; j < mat->rptr[i+1]; ++j, ++jj){
            EXPECT_EQ(mat->rind[j], smat->rind[jj]);
            EXPECT_FLOAT_EQ(mat->rval[j], smat->rval[jj]);
            EXPECT_LE(smat->rnorms[jj], 1.0001);
        }
    }
    delete mat;
    delete smat;
}


TEST(csr, compact_rows){
    auto mat = csr_t::random(1059, 1321, 0.0005);
    idx_t n = 0, n2 = 0; // number of non-empty rows/columns
    long double s = 0, s2 = 0;
    for(idx_t i=0; i < mat->nrows; ++i){
        if(mat->rptr[i+1] > mat->rptr[i]){
            n++;
        }
        for(ptr_t j=mat->rptr[i]; j < mat->rptr[i+1]; ++j){
            s += mat->rval[j];
        }
    }
    EXPECT_NE(mat->nrows, n);
    mat->compact_rows();
    EXPECT_EQ(mat->nrows, n);
    for(idx_t i=0; i < mat->nrows; ++i){
        if(mat->rptr[i+1] > mat->rptr[i]){
            n2++;
        }
        for(ptr_t j=mat->rptr[i]; j < mat->rptr[i+1]; ++j){
            s2 += mat->rval[j];
        }
    }
    EXPECT_EQ(n, n2);
    EXPECT_FLOAT_EQ(s, s2);
    delete mat;
}


TEST(csr, compact_cols){
    auto mat = csr_t::random(1059, 1321, 0.0005);
    mat->create_index();
    idx_t n = 0, n2 = 0; // number of non-empty rows/columns
    long double s = 0, s2 = 0;
    for(idx_t i=0; i < mat->ncols; ++i){
        if(mat->cptr[i+1] > mat->cptr[i]){
            n++;
        }
        for(ptr_t j=mat->cptr[i]; j < mat->cptr[i+1]; ++j){
            s += mat->cval[j];
        }
    }
    EXPECT_NE(mat->ncols, n);
    mat->compact_cols(); // affects only the row structure
    mat->create_index();
    EXPECT_EQ(mat->ncols, n);
    for(idx_t i=0; i < mat->ncols; ++i){
        if(mat->cptr[i+1] > mat->cptr[i]){
            n2++;
        }
        for(ptr_t j=mat->cptr[i]; j < mat->cptr[i+1]; ++j){
            s2 += mat->cval[j];
        }
    }
    EXPECT_EQ(n, n2);
    EXPECT_FLOAT_EQ(s, s2);
    delete mat;
}

TEST(csr, sort_nnzs_by_ind){
    auto mat = csr_t::random(10, 138, 0.05);
    mat->create_index();
    mat->sort_nnzs_by_ind(true, CSR_BASE::ROW);
    for(idx_t i=0; i < mat->nrows; ++i){
        for(ptr_t j=mat->rptr[i]+1; j < mat->rptr[i+1]; ++j){
            EXPECT_GT(mat->rind[j], mat->rind[j-1]);
        }
    }
    mat->sort_nnzs_by_ind(false, CSR_BASE::ROW);
    for(idx_t i=0; i < mat->nrows; ++i){
        for(ptr_t j=mat->rptr[i]+1; j < mat->rptr[i+1]; ++j){
            EXPECT_LT(mat->rind[j], mat->rind[j-1]);
        }
    }
    mat->sort_nnzs_by_ind(true, CSR_BASE::COL);
    for(idx_t i=0; i < mat->ncols; ++i){
        for(ptr_t j=mat->cptr[i]+1; j < mat->cptr[i+1]; ++j){
            EXPECT_GT(mat->cind[j], mat->cind[j-1]);
        }
    }
    mat->sort_nnzs_by_ind(false, CSR_BASE::COL);
    for(idx_t i=0; i < mat->ncols; ++i){
        for(ptr_t j=mat->cptr[i]+1; j < mat->cptr[i+1]; ++j){
            EXPECT_LT(mat->cind[j], mat->cind[j-1]);
        }
    }
    delete mat;
}

TEST(csr, sort_nnzs_by_value){
    auto mat = csr_t::random(10, 138, 0.05);
    mat->create_index();
    // sort row structure by value
    mat->sort_nnzs_by_value(true, CSR_BASE::ROW);
    for(idx_t i=0; i < mat->nrows; ++i){
        for(ptr_t j=mat->rptr[i]+1; j < mat->rptr[i+1]; ++j){
            EXPECT_GE(mat->rval[j], mat->rval[j-1]);
        }
    }
    mat->sort_nnzs_by_value(false, CSR_BASE::ROW);
    for(idx_t i=0; i < mat->nrows; ++i){
        for(ptr_t j=mat->rptr[i]+1; j < mat->rptr[i+1]; ++j){
            EXPECT_LE(mat->rval[j], mat->rval[j-1]);
        }
    }
    // sort column structure by value
    mat->sort_nnzs_by_value(true, CSR_BASE::COL);
    for(idx_t i=0; i < mat->ncols; ++i){
        for(ptr_t j=mat->cptr[i]+1; j < mat->cptr[i+1]; ++j){
            EXPECT_GE(mat->cval[j], mat->cval[j-1]);
        }
    }
    mat->sort_nnzs_by_value(false, CSR_BASE::COL);
    for(idx_t i=0; i < mat->ncols; ++i){
        for(ptr_t j=mat->cptr[i]+1; j < mat->cptr[i+1]; ++j){
            EXPECT_LE(mat->cval[j], mat->cval[j-1]);
        }
    }
    // partially sort row structure by value
    mat->sort_nnzs_by_value(true, CSR_BASE::ROW, 5);
    for(idx_t i=0; i < mat->nrows; ++i){
        for(ptr_t j=mat->rptr[i]+1; j < mat->rptr[i+1] && j < mat->rptr[i]+5; ++j){
            EXPECT_GE(mat->rval[j], mat->rval[j-1]);
        }
    }
    mat->sort_nnzs_by_value(false, CSR_BASE::ROW, 5);
    for(idx_t i=0; i < mat->nrows; ++i){
        for(ptr_t j=mat->rptr[i]+1; j < mat->rptr[i+1] && j < mat->rptr[i]+5; ++j){
            EXPECT_LE(mat->rval[j], mat->rval[j-1]);
        }
    }
    delete mat;
}


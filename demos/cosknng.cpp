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

//header files
#include "../snnlib/CosKnnIndex.h"

using namespace snnlib;

int main (int argc, char *argv[]){
    auto mat = csr_t::random(100, 837, 0.15);
    auto index = new CosKnnIndex(mat);
    auto graph = index->knng(10, true);
    std::cout << "nthreads: " << index->get_nthreads() << std::endl;
    graph->print_info();
    for(idx_t i=0; i < graph->nrows; ++i){
        for(ptr_t j=graph->rptr[i]; j < graph->rptr[i+1]; ++j){
            std::cout << graph->rind[j] << " " << graph->rval[j] << " ";
        }
        std:: cout << std::endl;
    }
    std:: cout << std::endl;
    delete mat;
    delete index;
    delete graph;
}
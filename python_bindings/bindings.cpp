#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <cstdint>

#include "../snnlib/basetypes.h"
#include "../snnlib/csr.h"
#include "../snnlib/CosKnnIndex.h"

#define PYBIND11_DETAILED_ERROR_MESSAGES

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace pybind11::literals;

using namespace snnlib;

class matrix
{
    public:
    csr_t * mat;

    matrix(
        idx_t const nrows=0, 
        idx_t const ncols=0, 
        idx_t const maxnrows=0, 
        ptr_t const maxnnz=0){
        mat = new csr_t(nrows, ncols, maxnrows, maxnnz);
    }

    matrix(csr_t * mat){
        this->mat = mat;
    }

    matrix(py::object input)
    {

        auto obj = py::reinterpret_borrow<py::object>(input);
        py::object sparse_module = py::module_::import("scipy.sparse");
        py::object matrix_type = sparse_module.attr("csr_matrix");

        if (!py::type::handle_of(obj).is(matrix_type)) {
            try {
                obj = matrix_type(obj);
            } catch (const py::error_already_set &) {
                throw std::runtime_error("Input object for <matrix> must be of type csr_matrix.");
            }
        }

        py::tuple shape = input.attr("shape");
        idx_t nrows = shape[0].cast<idx_t>();
        idx_t ncols = shape[1].cast<idx_t>();

        auto data = pybind11::array_t<double, py::array::c_style | py::array::forcecast>::ensure(input.attr("data"));
        auto data_ptr = (double *) data.mutable_data();
        auto nnz = (ptr_t) data.size();
        val_t * val = (val_t *) malloc(nnz * sizeof(val_t));
        for(ssize_t i=0; i < nnz; ++i){
            val[i] = (val_t) data_ptr[i];
        }

        auto indices = pybind11::array_t<int64_t, py::array::c_style | py::array::forcecast>::ensure(input.attr("indices"));
        auto indices_ptr = (int64_t *) indices.mutable_data();
        idx_t * ind = (idx_t *) malloc(nnz * sizeof(idx_t));
        for(ssize_t i=0; i < nnz; ++i){
            ind[i] = (idx_t) indices_ptr[i];
        }

        auto indptr = pybind11::array_t<ssize_t, py::array::c_style | py::array::forcecast>::ensure(input.attr("indptr"));
        auto indptr_ptr = (ssize_t *) indptr.mutable_data();
        auto size = indptr.size();
        ptr_t * ptr = (ptr_t *) calloc(size > 0 ? size : 1, sizeof(ptr_t));
        for(ssize_t i=0; i < size; ++i){
            ptr[i] = (ptr_t) indptr_ptr[i];
        }

        mat = new csr_t(nrows, ncols, nnz, ind, val, ptr);
    }

    ~matrix()
    {
        if(mat){
            delete mat;
        }
    }

    idx_t get_nrows(){
        return mat->nrows;
    }

    idx_t get_ncols(){
        return mat->ncols;
    }

    idx_t get_nnz(){
        return mat->nnz();
    }

    void alloc(
        idx_t const maxnr, 
        ptr_t const maxnz
    ){
        mat->alloc(maxnr, maxnz);
    }

    void shrink(){
        mat->shrink();
    }

    void normalize(int const norm=2){
        if(norm == 1){
            mat->normalize(CSR_NORM::L1);
        } else if(norm == 2){
            mat->normalize(CSR_NORM::L2);
        } else {
            throw std::runtime_error("Norm can only be 1 for L1 or 2 for L2.");
        }
    }

    void compact_rows(){
        mat->compact_rows();
    }

    void compact_cols(){
        mat->compact_cols();
    }

    void sort_nnzs_by_ind(const bool ascending = true){
        mat->sort_nnzs_by_ind(ascending);
    }

    void sort_nnzs_by_value(const bool ascending = true, const int nselect=0){
        mat->sort_nnzs_by_value(ascending, CSR_BASE::ROW, nselect);
    }

    matrix * extract_submatrix(const int rstart, const int numrows) {
        auto m = mat->extract_submatrix(rstart, numrows);
        return new matrix(m);
    }

    static matrix * random(
        idx_t const nrows, 
        idx_t const ncols, 
        double const factor=0.05
    ){
        auto mat = csr_t::random(nrows, ncols, factor);
        return new matrix(mat);
    }

    void print_info(){
        mat->print_info();
    }

    void print(){
        mat->print();
    }

    py::object as_csr()
    {
        if(!mat || !mat->rptr || !mat->rind || !mat->rval || mat->nrows == 0){
            throw std::runtime_error("Empty or invalid matrix cannot be converted to csr_matrix.");
        }

        py::object csr_matrix_type = py::module_::import("scipy.sparse").attr("csr_matrix");
        auto nnz = mat->rptr[mat->nrows];
        py::array data(nnz, mat->rval);
        py::array indices(nnz, mat->rind);
        py::array indptr(mat->nrows + 1, mat->rptr);
        
        return csr_matrix_type(
            std::make_tuple(data, indices, indptr),
            std::make_pair(mat->nrows, mat->ncols)
        );
    }

    static py::object to_csr(csr_t const * const mat)
    {
        if(!mat || !mat->rptr || !mat->rind || !mat->rval || mat->nrows == 0){
            throw std::runtime_error("Empty or invalid matrix cannot be converted to csr_matrix.");
        }

        py::object csr_matrix_type = py::module_::import("scipy.sparse").attr("csr_matrix");
        auto nnz = mat->rptr[mat->nrows];
        py::array data(nnz, mat->rval);
        py::array indices(nnz, mat->rind);
        py::array indptr(mat->nrows + 1, mat->rptr);
        
        return csr_matrix_type(
            std::make_tuple(data, indices, indptr),
            std::make_pair(mat->nrows, mat->ncols)
        );
    }
};


class cos_knn_index
{
    public:
    matrix * mat{nullptr};
    CosKnnIndex * index{nullptr};

    cos_knn_index(py::object input, int nqrows=2e+4, int ndrows=1e+5, int ninnz=1e+6){
        mat = new matrix(input);
        index = new CosKnnIndex(mat->mat, nqrows, ndrows, ninnz);
    }

    ~cos_knn_index(){
        delete index;
        mat->mat = nullptr;
        delete mat;
    }

    py::object knng(idx_t const k, bool approx=false, double cfactor=2.0, int nenhance=1, double minprune=0.7, bool verbose=false){
        auto graph = index->knng(k, approx, cfactor, nenhance, minprune, verbose);
        auto ret = matrix::to_csr(graph);
        delete graph;
        return ret;
    }

    py::object knngbf(idx_t const k, bool verbose=false){
        auto graph = index->knngbf(k, verbose);
        auto ret = matrix::to_csr(graph);
        delete graph;
        return ret;
    }

    py::object get_index_data(){
        if(!index){
            throw std::runtime_error("Index does not exist for some reason.");
        }
        return matrix::to_csr(index->data);
    }

    int get_nthreads(){
        return index->get_nthreads();
    }

};


PYBIND11_MODULE(snnlib, m) {

    py::enum_<CSR_BASE>(m, "CSR_BASE")
    .value("ROW", CSR_BASE::ROW)
    .value("COL", CSR_BASE::COL)
    .value("ROWCOL", CSR_BASE::ROWCOL)
    .export_values();

    py::enum_<CSR_NORM>(m, "CSR_NORM")
    .value("L1", CSR_NORM::L1)
    .value("L2", CSR_NORM::L2)
    .export_values();

    py::enum_<CSR_SCALE>(m, "CSR_SCALE")
    .value("MAXTF", CSR_SCALE::MAXTF)
    .value("SQRT", CSR_SCALE::SQRT)
    .value("IDF", CSR_SCALE::IDF)
    .value("LOG", CSR_SCALE::LOG)
    .value("LOG10", CSR_SCALE::LOG10)
    .export_values();

    py::enum_<CSR_FORMAT>(m, "CSR_FORMAT")
    .value("CSR", CSR_FORMAT::CSR)
    .value("CLU", CSR_FORMAT::CLU)
    .value("COO", CSR_FORMAT::COO)
    .value("SMAT", CSR_FORMAT::SMAT)
    .value("BCSR", CSR_FORMAT::BCSR)
    .value("BCOO", CSR_FORMAT::BCOO)
    .value("MTX", CSR_FORMAT::MTX)
    .value("NPZ", CSR_FORMAT::NPZ)
    .export_values();

    py::class_<matrix>(m, "matrix")
    .def(py::init<>())
    .def(py::init<py::object>(), "mat"_a)
    .def("get_nrows", &matrix::get_nrows)
    .def("get_ncols", &matrix::get_ncols)
    .def("get_nnz", &matrix::get_nnz)
    .def("alloc", &matrix::alloc, "maxnr"_a, "maxnz"_a)
    .def("shrink", &matrix::shrink)
    .def("normalize", &matrix::normalize, "norm"_a = 2)
    .def("compact_rows", &matrix::compact_rows)
    .def("compact_cols", &matrix::compact_cols)
    .def("sort_nnzs_by_ind", &matrix::sort_nnzs_by_ind, "ascending"_a = true)
    .def("sort_nnzs_by_value", &matrix::sort_nnzs_by_value, "ascending"_a = true, "nselect"_a = 0)
    .def("extract_submatrix", &matrix::extract_submatrix, "start"_a, "nrows"_a)
    .def_static("random", &matrix::random, "nrows"_a, "ncols"_a, "factor"_a = 0.05)
    .def("print_info", &matrix::print_info)
    .def("print", &matrix::print)
    .def("as_csr", &matrix::as_csr);

    py::class_<cos_knn_index>(m, "cos_knn_index")
    .def(py::init<py::object>(), "csr_matrix"_a)
    .def(py::init<py::object, int>(), "csr_matrix"_a, "nqrows"_a = 2e+4)
    .def(py::init<py::object, int, int>(), "csr_matrix"_a, "nqrows"_a = 2e+4, "ndrows"_a = 1e+5)
    .def(py::init<py::object, int, int, int>(), "csr_matrix"_a, "nqrows"_a = 2e+4, "ndrows"_a = 1e+5, "ninnz"_a = 1e+6)
    .def("knng", &cos_knn_index::knng, "k"_a, "approx"_a = false, "cfactor"_a = 2.0, "nenhance"_a = 1, "minprune"_a = 0.7, "verbose"_a = false)
    .def("knngbf", &cos_knn_index::knngbf, "k"_a, "verbose"_a = false)
    .def("get_nthreads", &cos_knn_index::get_nthreads)
    .def("get_index_data", &cos_knn_index::get_index_data);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif

}

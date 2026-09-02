//
// Created by Arian Dovald on 8/31/26.
//

#ifndef PROPAGATE_MATRIX_TOOLS_H
#define PROPAGATE_MATRIX_TOOLS_H

#include <cassert>
#include <vector>
#include <complex>

// diagonal matrix
template <typename T>
struct diag_matrix {
    const int N; // number of diagonal matrix elements
    std::vector<T> diag; // holds matrix elements in flat vector
    explicit diag_matrix(const int size) : N(size), diag(size) {} // constructs N and diag
    const T& operator()(int i) const { return diag[i]; } // read-only accessor
    T& operator()(int i) { return diag[i]; } // read-write accessor
};

// off-diagonal matrix interface
template <typename T>
struct off_diag_matrix {
    int N; // size of matrix
    virtual ~off_diag_matrix() = default; // default constructor
    virtual const T& operator()(int i, int j) const = 0; // read-only accessor (to be overridden by inheritors)
    virtual T& operator()(int i, int j) = 0; // read-write accessor (to be overridden by inheritors)
};

// upper triangular matrix
template <typename T>
struct u_tri_matrix : off_diag_matrix<T> {
    std::vector<T> u_tri; // holds upper-triangular matrix elements
    explicit u_tri_matrix(const int size) : u_tri(size * (size - 1) / 2) { this->N = size; } // constructor
    // maps indexes i and j into the correct index of the flat array
    int mapIndex(int i, int j) const {
        assert(i < j && j < this->N); // note: temporary
        return i * this->N - i * (i + 1) / 2 + (j - i - 1);
    }
    const T& operator()(const int i, const int j) const override { return u_tri[mapIndex(i, j)]; } // read-only accessor
    T& operator()(const int i, const int j) override { return u_tri[mapIndex(i, j)]; } // read-write accessor
};

struct hermitian_matrix {
    int N; // size of matrix
    diag_matrix<std::complex<double>> d_herm; // holds real-valued diagonal elements
    u_tri_matrix<std::complex<double>> c_herm; // holds complex-valued coupling elements
    explicit hermitian_matrix(const int size) : N(size), d_herm(N), c_herm(N) {} // constructor
    // read-only accessor
    const std::complex<double>& operator()(const int i, const int j) const {
        return i==j? d_herm(i) : c_herm(i,j);
    }
    // read-write accessor
    std::complex<double>& operator()(const int i, const int j) {
        return i==j? d_herm(i) : c_herm(i,j);
    }
};

#endif //PROPAGATE_MATRIX_TOOLS_H

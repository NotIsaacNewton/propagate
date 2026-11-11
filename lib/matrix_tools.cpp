//
// Created by Arian Dovald on 8/31/26.
//

#include "matrix_tools.h"
#include <ostream>
#include <print>

// diagonalizes a hermitian matrix H
void diagonalizeHermitianMatrix(const hermitian_matrix& H, std::vector<double>& eigenvalues,
    std::vector<lapack_complex_double>& eigenvectors, std::vector<lapack_int>& isuppz) {
    // unpack hermitian_matrix into column-major dense array for LAPACKE
    // lapack_complex_double* is ABI-compatible with std::complex<double>
    std::vector<lapack_complex_double> A(H.N * H.N, {0.0, 0.0});
    for (int col = 0; col < H.N; col++) {
        for (int row = 0; row <= col; row++) {
            // index = row + col*channels
            const auto val = H(row, col);
            A[row + col * H.N] = {val.real(), val.imag()};
        }
    }
    lapack_int m = 0;
    // call ZHEEVR
    const lapack_int info = LAPACKE_zheevr(
        LAPACK_COL_MAJOR,           // storage order
        'V',                        // compute eigenvectors
        'A',                        // all eigenvalues
        'U',                        // upper triangle stored
        H.N,                // matrix size
        A.data(),                   // matrix (overwritten)
        H.N,                // leading dimension
        0.0, 0.0,                // vl, vu (unused for range='A')
        0, 0,                     // il, iu (unused for range='A')
        LAPACKE_dlamch('S'),  // abstol: safe minimum for reliability
        &m,                         // number of eigenvalues found
        eigenvalues.data(),
        eigenvectors.data(),
        H.N,                // ldz
        isuppz.data()
    );
    if (info != 0) {
        std::print("ZHEEVR failed with info={}\n", info);
    }
}

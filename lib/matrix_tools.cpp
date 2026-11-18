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
    std::vector<lapack_complex_double> A(H.N*H.N, {0.0, 0.0});
    for (int j = 0; j < H.N; j++) {
        for (int i = 0; i <= j; i++) {
            // index = row + col*channels
            const auto val = H(i, j);
            A[i + j*H.N] = {val.real(), val.imag()};
        }
    }
    lapack_int m = 0;
    // call ZHEEVR
    const lapack_int info = LAPACKE_zheevr(
        LAPACK_COL_MAJOR,           // storage order
        'V',                        // compute eigenvectors
        'A',                        // all eigenvalues
        'U',                        // upper triangle stored
        H.N,                        // matrix size
        A.data(),                   // matrix (overwritten)
        H.N,                        // leading dimension
        0.0, 0.0,                // vl, vu (unused for range='A')
        0, 0,                     // il, iu (unused for range='A')
        LAPACKE_dlamch('S'),  // abstol: safe minimum for reliability
        &m,                         // number of eigenvalues found
        eigenvalues.data(),
        eigenvectors.data(),
        H.N,                        // ldz
        isuppz.data()
    );
    if (info != 0) {
        std::print("ZHEEVR failed with info={}\n", info);
    }
}

// applies phases in eigenbasis 'eigenvectors' to vectors 'psi_vec'
void eigenbasisTransform(const int N, const std::vector<lapack_complex_double>& eigenvectors,
    const std::vector<double>& phases, std::vector<std::complex<double>>& psi_vec,
    std::vector<std::complex<double>>& coeff) {
    constexpr std::complex<double> one(1,0), zero(0,0);
    const auto* Zc = eigenvectors.data();

    // coeff = eigenvector^† * psi_vec
    cblas_zgemv(CblasColMajor, CblasConjTrans,
                N, N, &one, Zc, N,
                psi_vec.data(), 1, &zero,
                coeff.data(), 1);

    // apply phases
    for (int k = 0; k < N; k++)
        coeff[k] *= std::polar(1.0, phases[k]);

    // psi_vec = eigenvector * coeff
    cblas_zgemv(CblasColMajor, CblasNoTrans,
                N, N, &one, Zc, N,
                coeff.data(), 1, &zero,
                psi_vec.data(), 1);
}

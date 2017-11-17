#include "critical.hpp"

#define TILE_SIZE 8



void LU_decomp_kij_bare(const int n, const int lda, double * A, double * scratch) {
  // LU decomposition without pivoting (Doolittle algorithm)
  // In-place decomposition of form A=LU
  // L is returned below main diagonal of A
  // U is returned at and above main diagonal


  for (int k = 0; k < n; k++) {
    const double recAkk = 1.0/A[k*lda + k];
    for (int i = k + 1; i < n; i++) {
      A[i*lda + k] = A[i*lda + k]*recAkk;
      for (int j = k + 1; j < n; j++)
        A[i*lda + j] -= A[i*lda + k]*A[k*lda + j];
    }
  }
}



void LU_decomp_kij_vec(const int n, const int lda, double * A, double * scratch) {
  // LU decomposition without pivoting (Doolittle algorithm)
  // In-place decomposition of form A=LU
  // L is returned below main diagonal of A
  // U is returned at and above main diagonal

  for (int k = 0; k < n; ++k) {
    const double recAkk = 1.0/A[k*lda + k];
    for (int i = k + 1; i < n; ++i) {
      A[i*lda + k] = A[i*lda + k]*recAkk;
#ifdef __PGI
#pragma loop ivdep
#else
#pragma omp simd
#endif
      for (int j = k + 1; j < n; ++j)
	    A[i*lda + j] -= A[i*lda + k]*A[k*lda + j];
    }
  }
}


void LU_decomp_kij_opt(const int n, const int lda, double * A, double * scratch) {
  // LU decomposition without pivoting (Doolittle algorithm)
  // In-place decomposition of form A=LU
  // L is returned below main diagonal of A
  // U is returned at and above main diagonal

  for (int k = 0; k < n; ++k) {
    const double recAkk = 1.0/A[k*lda + k];
#pragma omp parallel for
    for (int i = k + 1; i < n; ++i) {
      A[i*lda + k] = A[i*lda + k]*recAkk;
#ifdef __PGI
#pragma loop ivdep
#else
#pragma omp simd
#endif
      for (int j = k + 1; j < n; ++j)
	    A[i*lda + j] -= A[i*lda + k]*A[k*lda + j];
    }
  }
}



void LU_decomp_kij_vec_reg(const int n, const int lda, double * A, double * scratch) {
  // LU decomposition without pivoting (Doolittle algorithm)
  // In-place decomposition of form A=LU
  // L is returned below main diagonal of A
  // U is returned at and above main diagonal

  const int tile = TILE_SIZE;

  for (int i = 0; i < n; ++i) {
#ifdef __PGI
#pragma loop ivdep
#else
#pragma omp simd aligned(scratch: 64)
#endif
    for (int j = 0; j < n; ++j) {
      scratch[i*lda + j] = 0.0;
    }
    scratch[i*lda + i] = 1.0;
  }

  for (int k = 0; k < n; k++) {
    const int jmin = k - k%tile;
    const double recAkk = 1.0/A[k*lda + k];
    for (int i = k + 1; i < n; i++) {
      scratch[i*lda + k] = A[i*lda + k]*recAkk;
#ifdef __PGI
#pragma loop ivdep
#else
#pragma omp simd aligned(A, scratch: 64) 
#endif
      for (int j = jmin; j < n; j++) {
	    A[i*lda + j] -= scratch[i*lda + k]*A[k*lda + j];
      }
    }
  }

  for (int i = 0; i < n; ++i) {
#ifdef __PGI
#pragma loop ivdep
#else
#pragma omp simd aligned(A, scratch: 64)
#endif
    for (int j = 0; j < i; ++j) {
      A[i*lda + j] = scratch[i*lda + j];
    }
  }
}


void LU_decomp_kij_opt_reg(const int n, const int lda, double * A, double * scratch) {
    // LU decomposition without pivoting (Doolittle algorithm)
    // In-place decomposition of form A=LU
    // L is returned below main diagonal of A
    // U is returned at and above main diagonal

    const int tile = TILE_SIZE;

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
#ifdef __PGI
#pragma loop ivdep
#else
#pragma omp simd aligned(scratch: 64)
#endif
      for (int j = 0; j < n; ++j) {
        scratch[i*lda + j] = 0.0;
      }
      scratch[i*lda + i] = 1.0;
    }

    for (int k = 0; k < n; k++) {
      const int jmin = k - k%tile;
      const double recAkk = 1.0/A[k*lda + k];
#pragma omp parallel for
      for (int i = k + 1; i < n; i++) {
        scratch[i*lda + k] = A[i*lda + k]*recAkk;
#ifdef __PGI
#pragma loop ivdep
#else
#pragma omp simd aligned(A, scratch: 64) 
#endif
        for (int j = jmin; j < n; j++) {
  	    A[i*lda + j] -= scratch[i*lda + k]*A[k*lda + j];
        }
      }
    }

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
#ifdef __PGI
#pragma loop ivdep
#else
#pragma omp simd aligned(A, scratch: 64)
#endif
      for (int j = 0; j < i; ++j) {
        A[i*lda + j] = scratch[i*lda + j];
      }
    }
}
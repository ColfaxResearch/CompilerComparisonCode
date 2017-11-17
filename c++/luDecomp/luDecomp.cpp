#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <cmath>
#endif
#include <stdlib.h>
#include <malloc.h>
#include <cstdio>
#include <omp.h>
#include <cassert>
#include <unistd.h>
#include <limits>
#ifdef __HBM__
#include <hbwmalloc.h>
#endif
#include "critical.hpp"
#if defined(__INTEL_COMPILER) && defined(__ADVIXE__)
#include "advisor-annotate.h"
#endif

#if !defined(__INTEL_COMPILER) && !defined(__clang__) && !defined(__PGI)
#define __GNU_COMPILER
#endif

#define PROBLEM_SIZE 128
#define NUM_MATRICES 10000
#define NUM_TRIALS 10
#define TILE_SIZE 8
#define CACHE_LINE_LENGTH 8

#define KIJ_VEC         1
#define KIJ_OPT         2
#define KIJ_VEC_REG     3
#define KIJ_OPT_REG     4
#define BARE            5

#ifndef __ALGORITHM__
#define __ALGORITHM __ KIJ_OPT_REG
#endif

/*****************************************************************************/

void VerifyResult(const int n, const int lda, double* LU, double* refA) {

  // Verifying that A=LU
double *A, *L, *U;
#ifdef __INTEL_COMPILER
  A = static_cast<double*>(_mm_malloc(n*lda*sizeof(double), 64));
  L = static_cast<double*>(_mm_malloc(n*lda*sizeof(double), 64));
  U = static_cast<double*>(_mm_malloc(n*lda*sizeof(double), 64));
#else
  posix_memalign((void **)&A, 4096, n*lda*sizeof(double));
  posix_memalign((void **)&L, 4096, n*lda*sizeof(double));
  posix_memalign((void **)&U, 4096, n*lda*sizeof(double));
#endif
  for (int i = 0, arrSize = n*lda; i < arrSize; ++i) {
    A[i] = 0.0f;
      L[i] = 0.0f;
      U[i] = 0.0f;
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++)
      L[i*lda + j] = LU[i*lda + j];
    L[i*lda+i] = 1.0f;
    for (int j = i; j < n; j++)
      U[i*lda + j] = LU[i*lda + j];
  }
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
    A[i*lda + j] += L[i*lda + k]*U[k*lda + j];

  double deviation1 = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      deviation1 += (refA[i*lda+j] - A[i*lda+j])*(refA[i*lda+j] - A[i*lda+j]);
    }
  }
  deviation1 /= (double)(n*lda);
  double epsilon = std::numeric_limits<double>::epsilon();
  if (std::isnan(deviation1) || (deviation1 > 1.0e-2)) {
#ifndef __AUTO__
    printf("ERROR: LU is not equal to A (deviation1 = %e and epsilon = %e)!\n", deviation1, epsilon);
#endif
    exit(1);
} else {
#ifndef __AUTO__
    printf("OK: LU is equal to A (deviation1 = %e but epsilon = %e)!\n", deviation1, epsilon);
#endif
}

#ifdef VERBOSE
  printf("\n(L-D)+U:\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      printf("%10.3e", LU[i*lda+j]);
    printf("\n");
  }

  printf("\nL:\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      printf("%10.3e", L[i*lda+j]);
    printf("\n");
  }

  printf("\nU:\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      printf("%10.3e", U[i*lda+j]);
    printf("\n");
  }

  printf("\nLU:\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      printf("%10.3e", A[i*lda+j]);
    printf("\n");
  }

  printf("\nA:\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      printf("%10.3e", refA[i*lda+j]);
    printf("\n");
  }

  printf("deviation1=%e\n", deviation1);
#endif

#ifdef __INTEL_COMPILER
  _mm_free(A);
  _mm_free(L);
  _mm_free(U);
#else
  free(A);
  free(L);
  free(U);
#endif

}

int main(const int argc, const char** argv) {

  // Problem size and other parameters

  int n1;
  const char* s1 = getenv("PROBLEM_SIZE");
  if (s1 != NULL) {
    n1 = atoi(s1);
  } else {
    n1 = PROBLEM_SIZE;
  }
  const int n = n1;

  int nMatrices1;
  const char* s2 = getenv("NUM_PROBLEMS");
  if (s2 != NULL) {
    nMatrices1 = atoi(s2);
  } else {
    nMatrices1 = NUM_MATRICES;
  }
  const int nMatrices = nMatrices1;

  int nTrials1;
  const char* s3 = getenv("NUM_TRIALS");
  if (s3 != NULL) {
    nTrials1  = atoi(s3);
  } else {
    nTrials1 = NUM_TRIALS;
  }
  const int nTrials = (int)nTrials1;

  const int lda = n+16;
  const size_t containerSize = sizeof(double)*n*lda+64;

  const double HztoPerf = 1e-9*2.0/3.0*double(n*n*static_cast<double>(lda))*nMatrices;

  char *dataA;
  double *scratch, *referenceMatrix;
  const int cache_line_length = CACHE_LINE_LENGTH, num_threads = omp_get_max_threads();

#ifndef __AUTO__
#if defined(__INTEL_COMPILER)
  printf("Compiled with the Intel Compiler\n");
#elif defined(__clang__)
  printf("Compiled with the LLVM Compiler\n");
#elif defined(__GNUC__)
  printf("Compiled with the GNU Compiler\n");
#elif defined(__PGI)
  printf("Compiled with the PGI Compiler\n");
#endif
#endif

#ifdef __HBM__
#ifndef __AUTO__
  printf("Using High Bandwidth Memory (HBM)\n");
#endif
  hbw_posix_memalign((void**) &dataA, 4096, containerSize*nMatrices);
  hbw_posix_memalign((void**) &scratch, 4096, sizeof(double)*(num_threads*cache_line_length + n*lda) + 64);
#else
#ifndef __AUTO__
  printf("Not using High Bandwidth Memory (HBM)\n");
#endif
#ifdef __INTEL_COMPILER
  dataA = (char*)_mm_malloc(containerSize*nMatrices, 64);
  scratch = (double*)_mm_malloc(sizeof(double)*(num_threads*cache_line_length + n*lda) + 64, 64);
#else
  posix_memalign((void **)&dataA, 4096, containerSize*nMatrices);
  posix_memalign((void **)&scratch, 4096, sizeof(double)*(num_threads*cache_line_length + n*lda) + 64);
#endif
#endif
#ifdef __INTEL_COMPILER
  referenceMatrix = static_cast<double*>(_mm_malloc(n*lda*sizeof(double), 64));
#else
  posix_memalign((void **)&referenceMatrix, 4096, n*lda*sizeof(double));
#endif

  // Initialize matrices
  for (int m = 0; m < nMatrices; m++) {
    double* matrix = (double*)(&dataA[m*containerSize]);
#pragma omp parallel for
    for (int i = 0; i < n; i++) {
      double sum = 0.0;
      for (int j = 0; j < n; j++) {
        matrix[i*lda+j] = (double)(i*n+j);
        sum += matrix[i*lda+j];
      }
      sum -= matrix[i*lda+i];
      matrix[i*lda+i] = 2.0f*sum;
    }
    matrix[(n-1)*lda+n] = 0.0f; // Touch just in case
  }
#if defined(__INTEL_COMPILER) || defined(__clang__)
#pragma omp parallel for simd
#elif defined(__GNU_COMPILER)
#pragma omp parallel for
#endif
  for (int i = 0; i < n*lda; ++i) {
    referenceMatrix[i] = ((double *)dataA)[i];
  }
const int llim = num_threads*cache_line_length + n*lda + 8;
#if defined(__INTEL_COMPILER) || defined(__clang__)
#pragma omp parallel for simd
#elif defined(__GNU_COMPILER)
#pragma omp parallel for
#endif
  for (int i = 0; i < llim; ++i) {
    scratch[i] = 0.0;
  }

  // Perform benchmark
#ifndef __AUTO__
  printf("LU decomposition of %d matrices of size %dx%d on CPU...\n\n", (int)nMatrices, (int)n, (int)n);
#if __ALGORITHM__ == KIJ_VEC
  printf("Dolittle Algorithm (kij version - vectorized)\n");
#elif __ALGORITHM__ == KIJ_VEC_REG
  printf("Dolittle Algorithm (kij_regularized version - vectorized)\n");
#elif __ALGORITHM__ == KIJ_OPT
  printf("Dolittle Algorithm (kij version - vectorized + parallelized)\n");
#elif __ALGORITHM__ == KIJ_OPT_REG
  printf("Dolittle Algorithm (kij_regularized version - vectorized + parallelized)\n");
#elif __ALGORITHM__ == BARE
  printf("Dolittle Algorithm (kij_bare version)\n");
#endif
#endif

  double rate = 0, dRate = 0; // Benchmarking data
  const int skipTrials = 3; // First step is warm-up on Xeon Phi coprocessor
#ifndef __AUTO__
  printf("\033[1m%5s %10s %8s\033[0m\n", "Trial", "Time, s", "GFLOP/s");
#endif
  for (int trial = 1; trial <= nTrials; trial++) {

    const double tStart = omp_get_wtime(); // Start timing
    for (int m = 0; m < nMatrices; m++) {
      double * matrixA = (double*)(&dataA[m*containerSize]);
#if __ALGORITHM__ == KIJ_VEC
        LU_decomp_kij_vec(n, lda, matrixA, scratch);
#elif __ALGORITHM__ == KIJ_VEC_REG
        LU_decomp_kij_vec_reg(n, lda, matrixA, scratch);
#elif __ALGORITHM__ == KIJ_OPT
        LU_decomp_kij_opt(n, lda, matrixA, scratch);
#elif __ALGORITHM__ == KIJ_OPT_REG
        LU_decomp_kij_opt_reg(n, lda, matrixA, scratch);
#elif __ALGORITHM__ == BARE
        LU_decomp_kij_bare(n, lda, matrixA, scratch);
#endif
    }
    const double tEnd = omp_get_wtime(); // End timing

#ifndef __AUTO__
    if (trial == 1) VerifyResult(n, lda, (double*)(&dataA[0]), referenceMatrix);
#endif

    if (trial > skipTrials) { // Collect statistics
      rate  += HztoPerf/(tEnd - tStart);
      dRate += HztoPerf*HztoPerf/((tEnd - tStart)*(tEnd-tStart));
    }

#ifndef __AUTO__
    printf("%5d %10.3e %8.2f %s\n",
       (int)trial, (tEnd-tStart), HztoPerf/(tEnd-tStart), (trial<=skipTrials?"*":""));
    fflush(stdout);
#endif
  }
  rate/=(double)(nTrials-skipTrials);
  dRate=sqrt(dRate/(double)(nTrials-skipTrials)-rate*rate);
#ifndef __AUTO__
  printf("-----------------------------------------------------\n");
#endif

#ifndef __AUTO__
  printf("\033[1m%s %4s \033[42m%10.2f +- %.2f GFLOP/s\033[0m\n",
     "Average performance:", "", rate, dRate);
#else
  printf("Average performance: %17.16e +- %17.16e GFLOP/s\n", rate, dRate);
#endif

#ifndef __AUTO__
  printf("-----------------------------------------------------\n");
  printf("* - warm-up, not included in average\n\n");
#endif

#ifdef __HBM__
  hbw_free((void*)dataA);
  hbw_free((void*)scratch);
#else
#ifdef __INTEL_COMPILER
  _mm_free(dataA);
  _mm_free(scratch);
#else
  free(dataA);
  free(scratch);
#endif
#endif
#ifdef __INTEL_COMPILER
  _mm_free(referenceMatrix);
#else
  free(referenceMatrix);
#endif
}

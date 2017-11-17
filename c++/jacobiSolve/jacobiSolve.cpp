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
#if defined(__INTEL_COMPILER) && defined(__ADVIXE__)
#include "advisor-annotate.h"
#endif

#include "jacobi.hpp"


#if !defined(__INTEL_COMPILER) && !defined(__clang__) && !defined(__PGI)
#define __GNU_COMPILER
#endif

#define PROBLEM_SIZE 256
#define NUM_PROBLEMS 1
#define NUM_TRIALS 10
#define TILE_SIZE 8
#define CACHE_LINE_LENGTH 8

#define __TOLERANCE 1.0e3

#define SOLVE_BARE      1
#define SOLVE_VEC       2
#define SOLVE_OPT       3

#ifndef __ALGORITHM__
#define __ALGORITHM __ SOLVE_OPT
#endif

/*****************************************************************************/

/*
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
*/

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

  int nProblems1;
  const char* s2 = getenv("NUM_PROBLEMS");
  if (s2 != NULL) {
    nProblems1 = atoi(s2);
  } else {
    nProblems1 = NUM_PROBLEMS;
  }
  const int nProblems = nProblems1;

  int nTrials1;
  const char* s3 = getenv("NUM_TRIALS");
  if (s3 != NULL) {
    nTrials1  = atoi(s3);
  } else {
    nTrials1 = NUM_TRIALS;
  }
  const int nTrials = (int)nTrials1;

  double const epsilon = 8.8541878176e-12;
  double const pi = 3.14159265359;
  double chargeDensity = 1.0e-6/epsilon;

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
#else
#ifndef __AUTO__
  printf("Not using High Bandwidth Memory (HBM)\n");
#endif
#endif

  // Initialize solvers
  Jacobi<double> J(n);
  //J.resize(n,n);
  J.set_Sy(-500.0);
  J.set_Sx(-500.0);
  J.set_Fy(500.0);
  J.set_Fx(500.0);
  double y, x;
  double yBall_1 = -85.0, xBall_1 = 0.0, rhoBall_1 = 2.0, rBall_1 = 20.0, qBall_1 = rhoBall_1*(4.0/3.0)*pi*std::pow(rBall_1, 3.0);
  double yBall_2 = +65.0, xBall_2 = 0.0, rhoBall_2 = -3.0, rBall_2 = 20.0, qBall_2 = rhoBall_2*(4.0/3.0)*pi*std::pow(rBall_2, 3.0);
  double yBall_3 = 0.0, xBall_3 = -60.0, rhoBall_3 = -1.0, rBall_3 = 20.0, qBall_3 = rhoBall_3*(4.0/3.0)*pi*std::pow(rBall_3, 3.0);
  double yBall_4 = 0.0, xBall_4 = +60.0, rhoBall_4 = 2.0, qBall_4 = -1.0*(qBall_1 + qBall_2 + qBall_3), rBall_4 = std::pow((3.0*std::fabs(qBall_4))/(4.0*pi), 1.0/3.0);
  for (int i = 0; i < J.get_Ny(); ++i) {
    for (int j = 0; j < J.get_Nx(); ++j) {
      y = J.get_y(i, j);
      x = J.get_x(i, j);
      if ((pow(x - xBall_1, 2.0) + pow(y - yBall_1, 2.0)) < pow(rBall_1, 2)) {
        J.getSource()(i, j) += rhoBall_1*chargeDensity;
      }
      if ((pow(x - xBall_2, 2.0) + pow(y - yBall_2, 2.0)) < pow(rBall_2, 2)) {
        J.getSource()(i, j) += rhoBall_2*chargeDensity;
      }
      if ((pow(x - xBall_3, 2.0) + pow(y - yBall_3, 2.0)) < pow(rBall_3, 2)) {
        J.getSource()(i, j) += rhoBall_3*chargeDensity;
      }
      if ((pow(x - xBall_4, 2.0) + pow(y - yBall_4, 2.0)) < pow(rBall_4, 2)) {
        J.getSource()(i, j) += rhoBall_4*chargeDensity;
      }
    }
  }
  J.setBoundary(static_cast<double>(0.0));

  // Perform benchmark
#ifndef __AUTO__
  printf("Jacobi solution of %d domains of size %dx%d on CPU...\n\n", (int)nProblems, (int)n, (int)n);
#if __ALGORITHM__ == SOLVE_BARE
  printf("Jacobi Algorithm (bare)\n");
#elif __ALGORITHM__ == SOLVE_VEC
  printf("Jacobi Algorithm (vectorized)\n");
#elif __ALGORITHM__ == SOLVE_OPT
  printf("Jacobi Algorithm (vectorized + parallelized)\n");;
#endif
#endif

  double rate = 0, dRate = 0; // Benchmarking data
  const int skipTrials = 3; // First step is warm-up on Xeon Phi coprocessor
#ifndef __AUTO__
  printf("\033[1m%5s %10s %8s\033[0m\n", "Trial", "Time, s", "GFLOP/s");
#endif
  for (int trial = 1; trial <= nTrials; trial++) {

    int nIter;
    double tTotal = 0.0;

    for (int m = 0; m < nProblems; m++) {
        J.resetDomain(static_cast<double>(0.0), static_cast<double>(0.0));
        const double tStart = omp_get_wtime(); // Start timing
#if __ALGORITHM__ == SOLVE_BARE
        nIter = J.solve_bare(__TOLERANCE);
#elif __ALGORITHM__ == SOLVE_VEC
        nIter = J.solve_vec(__TOLERANCE);
#elif __ALGORITHM__ == SOLVE_OPT
        nIter = J.solve_opt(__TOLERANCE);
#endif
        const double tEnd = omp_get_wtime(); // End timing
        tTotal += tEnd - tStart;
    }

  const double HztoPerf = 1e-9*(double)nIter*9*(double)n*(double)n*nProblems;

#ifndef __AUTO__
    //if (trial == 1) VerifyResult(n, lda, (double*)(&dataA[0]), referenceMatrix);
#endif

    if (trial > skipTrials) { // Collect statistics
      rate  += HztoPerf/tTotal;
      dRate += HztoPerf*HztoPerf/(tTotal*tTotal);
    }

#ifndef __AUTO__
    printf("%5d %10.3e %8.2f %s\n",
       (int)trial, tTotal, HztoPerf/tTotal, (trial<=skipTrials?"*":""));
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

}

#ifdef __INTEL_COMPILER
#include <mathimf.h>
#include <malloc.h>
#else
#include <cmath>
#include <mm_malloc.h>
#endif
#include <algorithm>
#include <stdlib.h>
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

#include "critical.hpp"

#if !defined(__INTEL_COMPILER) && !defined(__clang__) && !defined(__PGI)
#define __GNU_COMPILER
#endif

#define PROBLEM_SIZE 4096
#define NUM_PROBLEMS 100
#define NUM_TRIALS 10
#define BLOCK_SIZE 128
#define VECTOR_LENGTH 4

#define OBLKIO_BARE     1
#define OBLKIO_VEC      2
#define OBLKIO_OPT      3

#ifndef __ALGORITHM__
#define __ALGORITHM__ OBLKIO_BARE
#endif


/*****************************************************************************/

void VerifyResult(const int n, const double *SF) {

  // Verifying that A=LU
  double *A;
#ifdef __INTEL_COMPILER
  A = static_cast<double*>(_mm_malloc(n*sizeof(double), 64));
#else
  posix_memalign((void **)&A, 4096, n*sizeof(double));
#endif
  for (int i = 1, arrSize = n; i < arrSize; ++i) {
    A[i] = (double)i*n*n*i;
  }

  double deviation1 = 0.0;
  for (int i = 1; i < n; i++) {
    //printf("SF[%d]: %+4.3e; A[%d]: %+4.3e\n", (int)i, SF[i], (int)i, A[i]);
    deviation1 += (SF[i] - A[i])*(SF[i] - A[i]);
  }
  deviation1 /= (double)(n);
  double epsilon = std::numeric_limits<double>::epsilon();
  if (std::isnan(deviation1) || (deviation1 > epsilon)) {
#ifndef __AUTO__
    printf("ERROR: SF is not correct (deviation1 = %e and epsilon = %e)!\n", deviation1, epsilon);
#endif
    exit(1);
} else {
#ifndef __AUTO__
    printf("OK: SF is correct (deviation1 = %e but epsilon = %e)\n", deviation1, epsilon);
#endif
}

#ifdef __INTEL_COMPILER
  _mm_free(A);
#else
  free(A);
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
  const int n = (int)n1;

  int nProblems1;
  const char* s2 = getenv("NUM_PROBLEMS");
  if (s2 != NULL) {
    nProblems1 = atoi(s2);
  } else {
    nProblems1 = NUM_PROBLEMS;
  }
  const int nProblems = (int)nProblems1;

  int nTrials1;
  const char* s3 = getenv("NUM_TRIALS");
  if (s3 != NULL) {
    nTrials1  = atoi(s3);
  } else {
    nTrials1 = NUM_TRIALS;
  }
  const int nTrials = (int)nTrials1;

  const double HztoPerf = 1e-9*(6*(double)n*((double)n - 1)/2 + n)*nProblems;

    const int c = BLOCK_SIZE;
  const int containerSize = n + c;
  const int nthreads = omp_get_max_threads();

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

  double *dataA, *maskA, *outputA, *scratch;
#ifdef __HBM__
#ifndef __AUTO__
  printf("Using High Bandwidth Memory!\n");
#endif
  hbw_posix_memalign((void**) &dataA, 4096, sizeof(double)*containerSize*nProblems);
  hbw_posix_memalign((void**) &maskA, 4096, sizeof(double)*containerSize*nProblems);
  hbw_posix_memalign((void**) &outputA, 4096, sizeof(double)*containerSize*nProblems);
#else
#ifndef __AUTO__
  printf("Not using High Bandwidth Memory!\n");
#endif
  dataA = (double *)_mm_malloc(containerSize*nProblems*sizeof(double), 64);
  maskA = (double *)_mm_malloc(containerSize*nProblems*sizeof(double), 64);
  outputA = (double *)_mm_malloc(containerSize*nProblems*sizeof(double), 64);
#endif

// Parallel first touch
  for (int m = 0; m < nProblems; m++) {
    double* vector = &dataA[m*containerSize];
    double* vectorMask = &maskA[m*containerSize];
    for (int i = 0; i < containerSize; i++) {
      vector[i] = 0.0;
      vectorMask[i] = 0.0;
      outputA[i] = 0.0;
    }
  }

  // Initialize vectors
  for (int m = 0; m < nProblems; m++) {
    double* vector = &dataA[m*containerSize];
    double* vectorMask = &maskA[m*containerSize];
    for (int i = 0; i < n; i++) {
      vector[i] = (double)(i*n);
      vectorMask[i] = 1.0;
    }
  }

  // Perform benchmark
#ifndef __AUTO__
  printf("Structure Function calculation of %d vectors of size %d on %s...\n\n", (int)nProblems, (int)n, "CPU");
#if __ALGORITHM__ == OBLKIO_BARE
  printf("Structure Function Algorithm (oblkio version - bare)\n");
#elif __ALGORITHM__ == OBLKIO_VEC
  printf("Structure Function Algorithm (oblkio version - vectorized)\n");
#elif __ALGORITHM__ == OBLKIO_OPT
  printf("Structure Function Algorithm (oblkio version - vectorized + parallelized)\n");
#endif
#endif

  double rate = 0, dRate = 0; // Benchmarking data
  const int skipTrials = 3; // First step is warm-up on Xeon Phi coprocessor
#ifndef __AUTO__
  printf("\033[1m%5s %10s %8s\033[0m\n", "Trial", "Time, s", "GFLOP/s");
#endif


  for (int trial = 1; trial <= nTrials; trial++) {
    for (int i = 0; i < containerSize*nProblems; i++) {
      outputA[i] = 0.0;
    }

    double tStart, tEnd, tTotal = 0.0;
    for (int m = 0; m < nProblems; m++) {
      double* vectorA = &dataA[m*containerSize];
      double* vectorMaskA = &maskA[m*containerSize];
      double* resultA = &outputA[m*containerSize];
      tStart = omp_get_wtime(); // Start timing
#if __ALGORITHM__ == OBLKIO_BARE
      SF_compute_oblkio_bare(n, vectorA, vectorMaskA, resultA);
#elif __ALGORITHM__ == OBLKIO_VEC
      SF_compute_oblkio_vec(n, vectorA, vectorMaskA, resultA);
#elif __ALGORITHM__ == OBLKIO_OPT
      SF_compute_oblkio_opt(n, vectorA, vectorMaskA, resultA);
#endif
      tEnd = omp_get_wtime(); // End timing
      tTotal += tEnd - tStart;
    }

#ifndef __AUTO__
    if (trial == 1) VerifyResult(n, outputA);
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

#ifdef __HBM__
  hbw_free((void*)dataA);
  hbw_free((void*)maskA);
  hbw_free((void*)outputA);
#else
  _mm_free(dataA);
  _mm_free(maskA);
  _mm_free(outputA);
#endif
}

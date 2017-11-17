#include "critical.hpp"

#define BLOCK_SIZE 32



  void SF_compute_oblkio_bare(int const n, double const * const A, double const * const M, double * const SF) {
    // Structure Function calculation
    // Outof-place computation
    // SF is returned in SF

    const int c = BLOCK_SIZE;

    for (int oblk = 0; oblk < n/c; ++oblk) {
        double SFTemp_private[c];
        double countSFTemp_private[c];
        for (int ctr = 0; ctr < c; ++ctr) {
            SFTemp_private[ctr] = 0.0;
            countSFTemp_private[ctr] = 0.0;
        } // end for ctr
        register double * S0 = SFTemp_private + 0*8;
        register double * S1 = SFTemp_private + 1*8;
        register double * S2 = SFTemp_private + 2*8;
        register double * S3 = SFTemp_private + 3*8;
        register double * c0 = countSFTemp_private + 0*8;
        register double * c1 = countSFTemp_private + 1*8;
        register double * c2 = countSFTemp_private + 2*8;
        register double * c3 = countSFTemp_private + 3*8;
        for (int i = 0; i < n - oblk*c; ++i) {
            double MVal0, MVal1, MVal2, MVal3;
            for (int v = 0 ; v < 8; v++) {
                MVal0 = M[i + oblk*c + 0*8 + v]*M[i];
                S0[v] += MVal0*(A[i + oblk*c + 0*8 + v] - A[i])*(A[i + oblk*c + 0*8 + v] - A[i]);
                c0[v] += MVal0;
                MVal1 = M[i + oblk*c + 1*8 + v]*M[i];
                S1[v] += MVal1*(A[i + oblk*c + 1*8 + v] - A[i])*(A[i + oblk*c + 1*8 + v] - A[i]);
                c1[v] += MVal1;
                MVal2 = M[i + oblk*c + 2*8 + v]*M[i];
                S2[v] += MVal2*(A[i + oblk*c + 2*8 + v] - A[i])*(A[i + oblk*c + 2*8 + v] - A[i]);
                c2[v] += MVal2;
                MVal3 = M[i + oblk*c + 3*8 + v]*M[i];
                S3[v] += MVal3*(A[i + oblk*c + 3*8 + v] - A[i])*(A[i + oblk*c + 3*8 + v] - A[i]);
                c3[v] += MVal3;
            } // end for v
        } // end for i
        for (int ctr = 0; ctr < c; ++ctr) {
            SF[oblk*c + ctr] = SFTemp_private[ctr]/countSFTemp_private[ctr];
        } // end for ctr
    } // end for blk
  }



  void SF_compute_oblkio_vec(int const n, double const * const A, double const * const M, double * const SF) {
  // Structure Function calculation
  // Outof-place computation
  // SF is returned in SF

#if defined(__PGI)
#pragma routine safe
#endif

    const int c = BLOCK_SIZE;

    for (int oblk = 0; oblk < n/c; ++oblk) {
        double SFTemp_private[c];
        double countSFTemp_private[c];
        for (int ctr = 0; ctr < c; ++ctr) {
            SFTemp_private[ctr] = 0.0;
            countSFTemp_private[ctr] = 0.0;
        } // end for ctr
        register double * S0 = SFTemp_private + 0*8;
        register double * S1 = SFTemp_private + 1*8;
        register double * S2 = SFTemp_private + 2*8;
        register double * S3 = SFTemp_private + 3*8;
        register double * c0 = countSFTemp_private + 0*8;
        register double * c1 = countSFTemp_private + 1*8;
        register double * c2 = countSFTemp_private + 2*8;
        register double * c3 = countSFTemp_private + 3*8;
        for (int i = 0; i < n - oblk*c; ++i) {
            double MVal0, MVal1, MVal2, MVal3;
#if defined(__PGI)
#pragma loop ivdep
#else
#pragma omp simd private(MVal0, MVal1, MVal2, MVal3)
#endif
            for (int v = 0 ; v < 8; v++) {
                MVal0 = M[i + oblk*c + 0*8 + v]*M[i];
                S0[v] += MVal0*(A[i + oblk*c + 0*8 + v] - A[i])*(A[i + oblk*c + 0*8 + v] - A[i]);
                c0[v] += MVal0;
                MVal1 = M[i + oblk*c + 1*8 + v]*M[i];
                S1[v] += MVal1*(A[i + oblk*c + 1*8 + v] - A[i])*(A[i + oblk*c + 1*8 + v] - A[i]);
                c1[v] += MVal1;
                MVal2 = M[i + oblk*c + 2*8 + v]*M[i];
                S2[v] += MVal2*(A[i + oblk*c + 2*8 + v] - A[i])*(A[i + oblk*c + 2*8 + v] - A[i]);
                c2[v] += MVal2;
                MVal3 = M[i + oblk*c + 3*8 + v]*M[i];
                S3[v] += MVal3*(A[i + oblk*c + 3*8 + v] - A[i])*(A[i + oblk*c + 3*8 + v] - A[i]);
                c3[v] += MVal3;
            } // end for v
        } // end for i
        for (int ctr = 0; ctr < c; ++ctr) {
            SF[oblk*c + ctr] = SFTemp_private[ctr]/countSFTemp_private[ctr];
        } // end for ctr
    } // end for blk
  }



  void SF_compute_oblkio_opt(int const n, double const * const A, double const * const M, double * const SF) {
  // Structure Function calculation
  // Outof-place computation
  // SF is returned in SF

#if defined(__PGI)
#pragma routine safe
#endif

    const int c = BLOCK_SIZE;

#pragma omp parallel for schedule(dynamic)
    for (int oblk = 0; oblk < n/c; ++oblk) {
      double SFTemp_private[c];
      double countSFTemp_private[c];
      for (int ctr = 0; ctr < c; ++ctr) {
        SFTemp_private[ctr] = 0.0;
        countSFTemp_private[ctr] = 0.0;
      } // end for ctr
      register double * S0 = SFTemp_private + 0*4;
      register double * S1 = SFTemp_private + 1*4;
      register double * S2 = SFTemp_private + 2*4;
      register double * S3 = SFTemp_private + 3*4;
      register double * S4 = SFTemp_private + 4*4;
      register double * S5 = SFTemp_private + 5*4;
      register double * S6 = SFTemp_private + 6*4;
      register double * S7 = SFTemp_private + 7*4;
      register double * c0 = countSFTemp_private + 0*4;
      register double * c1 = countSFTemp_private + 1*4;
      register double * c2 = countSFTemp_private + 2*4;
      register double * c3 = countSFTemp_private + 3*4;
      register double * c4 = countSFTemp_private + 4*4;
      register double * c5 = countSFTemp_private + 5*4;
      register double * c6 = countSFTemp_private + 6*4;
      register double * c7 = countSFTemp_private + 7*4;
      for (int i = 0; i < n - oblk*c; ++i) {
        double MVal0, MVal1, MVal2, MVal3, MVal4, MVal5, MVal6, MVal7;
#if defined(__PGI)
#pragma loop ivdep
#else
#pragma omp simd private(MVal0, MVal1, MVal2, MVal3, MVal4, MVal5, MVal6, MVal7)
#endif
        for (int v = 0 ; v < 4; v++) {
          MVal0 = M[i + oblk*c + 0*4 + v]*M[i];
          S0[v] += MVal0*(A[i + oblk*c + 0*4 + v] - A[i])*(A[i + oblk*c + 0*4 + v] - A[i]);
          c0[v] += MVal0;
          MVal1 = M[i + oblk*c + 1*4 + v]*M[i];
          S1[v] += MVal1*(A[i + oblk*c + 1*4 + v] - A[i])*(A[i + oblk*c + 1*4 + v] - A[i]);
          c1[v] += MVal1;
          MVal2 = M[i + oblk*c + 2*4 + v]*M[i];
          S2[v] += MVal2*(A[i + oblk*c + 2*4 + v] - A[i])*(A[i + oblk*c + 2*4 + v] - A[i]);
          c2[v] += MVal2;
          MVal3 = M[i + oblk*c + 3*4 + v]*M[i];
          S3[v] += MVal3*(A[i + oblk*c + 3*4 + v] - A[i])*(A[i + oblk*c + 3*4 + v] - A[i]);
          c3[v] += MVal3;
          MVal4 = M[i + oblk*c + 4*4 + v]*M[i];
          S4[v] += MVal4*(A[i + oblk*c + 4*4 + v] - A[i])*(A[i + oblk*c + 4*4 + v] - A[i]);
          c4[v] += MVal4;
          MVal5 = M[i + oblk*c + 5*4 + v]*M[i];
          S5[v] += MVal5*(A[i + oblk*c + 5*4 + v] - A[i])*(A[i + oblk*c + 5*4 + v] - A[i]);
          c5[v] += MVal5;
          MVal6 = M[i + oblk*c + 6*4 + v]*M[i];
          S6[v] += MVal6*(A[i + oblk*c + 6*4 + v] - A[i])*(A[i + oblk*c + 6*4 + v] - A[i]);
          c6[v] += MVal6;
          MVal7 = M[i + oblk*c + 7*4 + v]*M[i];
          S7[v] += MVal7*(A[i + oblk*c + 7*4 + v] - A[i])*(A[i + oblk*c + 7*4 + v] - A[i]);
          c7[v] += MVal7;
        } // end for v
      } // end for i
      for (int ctr = 0; ctr < c; ++ctr) {
        SF[oblk*c + ctr] = SFTemp_private[ctr]/countSFTemp_private[ctr];
      } // end for ctr
    } // end for blk
  }


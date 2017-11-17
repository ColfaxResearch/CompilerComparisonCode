#ifndef CRITICAL_HPP
#define CRITICAL_HPP

void LU_decomp_kij_bare(const int n, const int lda, double * A, double * scratch);

void LU_decomp_kij_vec(const int n, const int lda, double * A, double * scratch);

void LU_decomp_kij_opt(const int n, const int lda, double * A, double * scratch);

void LU_decomp_kij_vec_reg(const int n, const int lda, double * A, double * scratch);

void LU_decomp_kij_opt_reg(const int n, const int lda, double * A, double * scratch);

#endif
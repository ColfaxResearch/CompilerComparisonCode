#ifndef CRITICAL_HPP
#define CRITICAL_HPP

void SF_compute_oblkio_bare(int const n, double const * const A, double const * const M, double * const SF);

void SF_compute_oblkio_vec(int const n, double const * const A, double const * const M, double * const SF);

void SF_compute_oblkio_opt(int const n, double const * const A, double const * const M, double * const SF);

#endif
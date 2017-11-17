#include <cmath>
#include <limits>

#include "jacobi.hpp"

/*
template <typename GridT>
int Jacobi<GridT>::solve_bare(GridT tol) {
    GridT newVal, maxChange = std::numeric_limits<GridT>::max();
    int iter = 0;
    while (maxChange > tol) {
        maxChange = static_cast<GridT>(0.0);
        iter += 1;
        // Perform one iteration
        int row, col;
        for (row = 1; row <= Ny; ++row) {
            for (col = 1; col <= Nx; ++col) {
                newVal = (DySq*((*SD)(row + 1, col) + (*SD)(row - 1, col)) + DxSq*((*SD)(row, col + 1) + (*SD)(row, col - 1)) - DySqDxSq*(*S)(row, col))/TwoDySqPlusDxSq;
                maxChange = (std::fabs(newVal - (*D)(row, col)) > maxChange) ? std::fabs(newVal - (*D)(row, col)) : maxChange;
                (*D)(row, col) = newVal;
            }
        }

        // Swap the ScratchDomain storage with the storage of the Domain.
        (*D).swapStorage((*SD));

    } // end while
    // Swap the ScratchDomain storage with the storage of the Domain.
    (*D).swapStorage((*SD));
    return iter;
}



template <typename GridT>
int Jacobi<GridT>::solve_vec(GridT tol) {
    GridT newVal, maxChange = std::numeric_limits<GridT>::max();
    int iter = 0;
    while (maxChange > tol) {
        maxChange = static_cast<GridT>(0.0);
        iter += 1;
        // Perform one iteration
        int row, col;
        for (row = 1; row <= Ny; ++row) {
#ifndef __PGI
#pragma omp simd private(newVal) reduction(max:maxChange) linear(col:1)
#else
#pragma loop ivdep
#endif
            for (col = 1; col <= Nx; ++col) {
                newVal = (DySq*((*SD)(row + 1, col) + (*SD)(row - 1, col)) + DxSq*((*SD)(row, col + 1) + (*SD)(row, col - 1)) - DySqDxSq*(*S)(row, col))/TwoDySqPlusDxSq;
                maxChange = (std::fabs(newVal - (*D)(row, col)) > maxChange) ? std::fabs(newVal - (*D)(row, col)) : maxChange;
                (*D)(row, col) = newVal;
            }
        }

        // Swap the ScratchDomain storage with the storage of the Domain.
        (*D).swapStorage((*SD));

    } // end while
    // Swap the ScratchDomain storage with the storage of the Domain.
    (*D).swapStorage((*SD));
    return iter;
}
*/


template <typename GridT>
int Jacobi<GridT>::solve_opt(GridT tol) {
    GridT newVal, maxChange = std::numeric_limits<GridT>::max();
    int iter = 0;
    while (maxChange > tol) {
        maxChange = static_cast<GridT>(0.0);
        iter += 1;
        // Perform one iteration
        int row, col;
#pragma omp parallel for private(newVal, col) reduction(max:maxChange)
        for (row = 1; row <= Ny; ++row) {
#ifndef __PGI
#pragma omp simd private(newVal) reduction(max:maxChange) linear(col:1)
#else
#pragma loop ivdep
#endif
            for (col = 1; col <= Nx; ++col) {
                newVal = (DySq*((*SD)(row + 1, col) + (*SD)(row - 1, col)) + DxSq*((*SD)(row, col + 1) + (*SD)(row, col - 1)) - DySqDxSq*(*S)(row, col))*OneOverTwoDySqPlusDxSq;
                maxChange = (std::fabs(newVal - (*D)(row, col)) > maxChange) ? std::fabs(newVal - (*D)(row, col)) : maxChange;
                (*D)(row, col) = newVal;
            }
        }

        // Swap the ScratchDomain storage with the storage of the Domain.
        (*D).swapStorage((*SD));

    } // end while
    // Swap the ScratchDomain storage with the storage of the Domain.
    (*D).swapStorage((*SD));
    return iter;
}

//template int Jacobi<double>::solve_bare(double );
//template int Jacobi<double>::solve_vec(double );
template int Jacobi<double>::solve_opt(double );
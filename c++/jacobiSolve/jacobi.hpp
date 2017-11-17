#ifndef JACOBI_HPP
#define JACOBI_HPP

#include <cstddef>
#include <omp.h>

#include "grid.hpp"


template <typename GridT>
class Jacobi {
private:
    int Nx, Ny, lda;
    GridT Sy, Sx, Fy, Fx, Ly, Lx, Dy, Dx, DySq, DxSq, DySqDxSq, OneOverTwoDySqPlusDxSq;
    static int const vecLength = 8;
    static int const cacheLineLength = 8;
    Grid<GridT> * D, * SD; // D is the domain while SD is the scratch domain.
    Grid<GridT> * S;
public:
    Jacobi();
    ~Jacobi();
    Jacobi(int new_N);
    Jacobi(int new_Ny, int new_Nx);
    void setBoundary(GridT bdryVal);
    void setGrid(GridT gridVal);
    void resetDomain(GridT gridVal, GridT bdryVal);
    GridT get_Sy() const;
    GridT get_Sx() const;
    void set_Sy(GridT new_Sy);
    void set_Sx(GridT new_Sx);
    GridT get_Fy() const;
    GridT get_Fx() const;
    void set_Fy(GridT new_Fy);
    void set_Fx(GridT new_Fx);
    GridT get_Ly() const;
    GridT get_Lx() const;
    GridT get_Dy() const;
    GridT get_Dx() const;
    int get_Ny() const;
    int get_Nx() const;
    int get_lda() const ;
    void resize(int new_Ny, int new_Nx);
    GridT get_y(int row, int col) const;
    GridT get_x(int row, int col) const;
    Grid<GridT> & getDomain();
    Grid<GridT> & getScratchDomain();
    Grid<GridT> & getSource();
    int solve_bare(GridT tol);
    int solve_vec(GridT tol);
    int solve_opt(GridT tol);
};

#endif // JACOBI_HPP

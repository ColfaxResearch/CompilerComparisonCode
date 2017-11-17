#include <cmath>
#include <limits>

#include "jacobi.hpp"


template <typename GridT>
Jacobi<GridT>::Jacobi() {
#ifdef __DEBUG__
        std::cout << "Jacobi<GridT>::Jacobi()" << std::endl;
#endif
    Ny = static_cast<int>(0);
    Nx = static_cast<int>(0);
    lda = static_cast<int>(0);
    Dy = static_cast<GridT>(1.0);
    Dx = static_cast<GridT>(1.0);
    DySq = static_cast<GridT>(1.0);
    DxSq = static_cast<GridT>(1.0);
    D = new Grid<GridT>();
    SD = new Grid<GridT>();
    S = new Grid<GridT>();
}

template <typename GridT>
Jacobi<GridT>::Jacobi(int new_N) {
#ifdef __DEBUG__
        std::cout << "Jacobi<GridT>::Jacobi(int)" << std::endl;
#endif
    Ny = new_N;
    Nx = new_N;
    lda = paddedLength(Nx, cacheLineLength);
    Dy = static_cast<GridT>(1.0);
    Dx = static_cast<GridT>(1.0);
    DySq = static_cast<GridT>(1.0);
    DxSq = static_cast<GridT>(1.0);
    DySqDxSq = DySq*DxSq;
    OneOverTwoDySqPlusDxSq = static_cast<GridT>(1.0/(2.0*(DySq + DxSq)));
    D = new Grid<GridT>(new_N);
    SD = new Grid<GridT>(new_N);
    S = new Grid<GridT>(new_N);
}

template <typename GridT>
Jacobi<GridT>::Jacobi(int new_Ny, int new_Nx) {
#ifdef __DEBUG__
        std::cout << "Jacobi<GridT>::Jacobi(int, int)" << std::endl;
#endif
    Ny = new_Ny;
    Nx = new_Ny;
    lda = paddedLength(Nx, cacheLineLength);
    Dy = static_cast<GridT>(1.0);
    Dx = static_cast<GridT>(1.0);
    DySq = static_cast<GridT>(1.0);
    DxSq = static_cast<GridT>(1.0);
    DySqDxSq = DySq*DxSq;
    OneOverTwoDySqPlusDxSq = static_cast<GridT>(1.0/(2.0*(DySq + DxSq)));
    D = new Grid<GridT>(Ny, Nx);
    SD = new Grid<GridT>(Ny, Nx);
    S = new Grid<GridT>(Ny, Nx);
}

template <typename GridT>
Jacobi<GridT>::~Jacobi() {
#ifdef __DEBUG__
        std::cout << "Jacobi<GridT>::~Jacobi()" << std::endl;
#endif
    delete D;
    delete SD;
    delete S;
    Ny = static_cast<int>(0);
    Nx = static_cast<int>(0);
    lda = static_cast<int>(0);
}

template <typename GridT>
void Jacobi<GridT>::setBoundary(GridT bdryVal) {
#ifdef __DEBUG__
        std::cout << "void Jacobi<GridT>::setBoundary(GridT bdryVal)" << std::endl;
#endif
    (*D)(0, 0) = bdryVal;
    (*D)(0, Nx + 1) = bdryVal;
    (*D)(Ny + 1, 0) = bdryVal;
    (*D)(Ny + 1, Nx + 1) = bdryVal;
    (*SD)(0, 0) = bdryVal;
    (*SD)(0, Nx + 1) = bdryVal;
    (*SD)(Ny + 1, 0) = bdryVal;
    (*SD)(Ny + 1, Nx + 1) = bdryVal;
    for (int col = 1; col <= Nx; ++col) {
        (*D)(0, col) = bdryVal;
        (*D)((Ny + 1), col) = bdryVal;
        (*SD)(0, col) = bdryVal;
        (*SD)((Ny + 1), col) = bdryVal;
    }

    for (int row = 1; row <= Ny; ++row) {
        (*D)(row, 0) = bdryVal;
        (*D)(row, (Nx + 1)) = bdryVal;
        (*SD)(row, 0) = bdryVal;
        (*SD)(row, (Nx + 1)) = bdryVal;
    }
}

template <typename GridT>
void Jacobi<GridT>::setGrid(GridT gridVal) {
#ifdef __DEBUG__
        std::cout << "void Jacobi<GridT>::setGrid(GridT gridVal)" << std::endl;
#endif
    for (int row = 1; row <= Ny; ++row) {
        for (int col = 1; col <= Nx; ++col) {
            (*D)(row, col) = gridVal;
            (*SD)(row, col) = gridVal;
        }
    }
}

template <typename GridT>
void Jacobi<GridT>::resetDomain(GridT gridVal, GridT bdryVal) {
#ifdef __DEBUG__
        std::cout << "void Jacobi<GridT>::resetDomain(GridT gridVal, GridT bdryVal)" << std::endl;
#endif
    setGrid(gridVal);
    setBoundary(bdryVal);
}

template <typename GridT>
GridT Jacobi<GridT>::get_Sy() const {
    return Sy;
}

template <typename GridT>
GridT Jacobi<GridT>::get_Sx() const {
    return Sx;
}

template <typename GridT>
void Jacobi<GridT>::set_Sy(GridT new_Sy) {
    Sy = new_Sy;
    Ly = Fy - Sy;
    Dy = (Fy - Sy)/Ny;
    DySq = Dy*Dy;
    DySqDxSq = DySq*DxSq;
    OneOverTwoDySqPlusDxSq = static_cast<GridT>(1.0/(2.0*(DySq + DxSq)));
    (*D).set_Sy(new_Sy);
    (*SD).set_Sy(new_Sy);
    (*S).set_Sy(new_Sy);
}

template <typename GridT>
void Jacobi<GridT>::set_Sx(GridT new_Sx) {
    Sx = new_Sx;
    Lx = Fx - Sx;
    Dx = (Fx - Sx)/Nx;
    DxSq = Dx*Dx;
    DySqDxSq = DySq*DxSq;
    OneOverTwoDySqPlusDxSq = static_cast<GridT>(1.0/(2.0*(DySq + DxSq)));
    (*D).set_Sx(new_Sx);
    (*SD).set_Sx(new_Sx);
    (*S).set_Sx(new_Sx);
}

template <typename GridT>
GridT Jacobi<GridT>::get_Fy() const {
    return Fy;
}

template <typename GridT>
GridT Jacobi<GridT>::get_Fx() const {
    return Fx;
}

template <typename GridT>
void Jacobi<GridT>::set_Fy(GridT new_Fy) {
    Fy = new_Fy;
    Ly = Fy - Sy;
    Dy = (Fy - Sy)/Ny;
    DySq = Dy*Dy;
    DySqDxSq = DySq*DxSq;
    OneOverTwoDySqPlusDxSq = static_cast<GridT>(1.0/(2.0*(DySq + DxSq)));
    (*D).set_Fy(new_Fy);
    (*SD).set_Fy(new_Fy);
    (*S).set_Fy(new_Fy);
}

template <typename GridT>
void Jacobi<GridT>::set_Fx(GridT new_Fx) {
    Fx = new_Fx;
    Lx = Fx - Sx;
    Dx = (Fx - Sx)/Nx;
    DxSq = Dx*Dx;
    DySqDxSq = DySq*DxSq;
    OneOverTwoDySqPlusDxSq = static_cast<GridT>(1.0/(2.0*(DySq + DxSq)));
    (*D).set_Fx(new_Fx);
    (*SD).set_Fx(new_Fx);
    (*S).set_Fx(new_Fx);
}

template <typename GridT>
GridT Jacobi<GridT>::get_Ly() const {
    return Ly;
}

template <typename GridT>
GridT Jacobi<GridT>::get_Lx() const {
    return Lx;
}

template <typename GridT>
GridT Jacobi<GridT>::get_Dy() const {
    return Dy;
}

template <typename GridT>
GridT Jacobi<GridT>::get_Dx() const {
    return Dx;
}

template <typename GridT>
int Jacobi<GridT>::get_Ny() const {
    return Ny;
}

template <typename GridT>
int Jacobi<GridT>::get_Nx() const {
    return Nx;
}

template <typename GridT>
int Jacobi<GridT>::get_lda() const {
    return lda;
}

template <typename GridT>
void Jacobi<GridT>::resize(int new_Ny, int new_Nx) {
#ifdef __DEBUG__
    std::cout << "Jacobi<GridT>::resize(int, int)" << std::endl;
#endif
    Ny = new_Ny;
    Nx = new_Nx;
    lda = paddedLength(Nx, cacheLineLength);
    Dy = (Fy - Sy)/Ny;
    Dx = (Fx - Sx)/Nx;
    DySq = Dy*Dy;
    DxSq = Dx*Dx;
    DySqDxSq = DySq*DxSq;
    OneOverTwoDySqPlusDxSq = static_cast<GridT>(1.0/(2.0*(DySq + DxSq)));
    (*D).resize(Ny, Nx);
    (*SD).resize(Ny, Nx);
    (*S).resize(Ny, Nx);
}

template <typename GridT>
GridT Jacobi<GridT>::get_y(int row, int col) const {
    return row*Dy + Sy;
}

template <typename GridT>
GridT Jacobi<GridT>::get_x(int row, int col) const {
    return col*Dx + Sx;
}


template <typename GridT>
Grid<GridT> & Jacobi<GridT>::getDomain() {
    return *D;
}

template <typename GridT>
Grid<GridT> & Jacobi<GridT>::getSource() {
    return *S;
}

template <typename GridT>
Grid<GridT> & Jacobi<GridT>::getScratchDomain() {
    return *SD;
}


template Jacobi<double>::Jacobi();
template Jacobi<double>::Jacobi(int );
template Jacobi<double>::Jacobi(int , int );
template Jacobi<double>::~Jacobi();
template void Jacobi<double>::setBoundary(double );
template void Jacobi<double>::setGrid(double );
template void Jacobi<double>::resetDomain(double , double );
template double Jacobi<double>::get_Sy() const;
template double Jacobi<double>::get_Sx() const;
template void Jacobi<double>::set_Sy(double );
template void Jacobi<double>::set_Sx(double );
template double Jacobi<double>::get_Fy() const;
template double Jacobi<double>::get_Fx() const;
template void Jacobi<double>::set_Fy(double );
template void Jacobi<double>::set_Fx(double );
template double Jacobi<double>::get_Ly() const;
template double Jacobi<double>::get_Lx() const;
template double Jacobi<double>::get_Dy() const;
template double Jacobi<double>::get_Dx() const;
template int Jacobi<double>::get_Ny() const;
template int Jacobi<double>::get_Nx() const;
template int Jacobi<double>::get_lda() const;
template void Jacobi<double>::resize(int , int );
template double Jacobi<double>::get_y(int , int ) const;
template double Jacobi<double>::get_x(int , int ) const;
template Grid<double> & Jacobi<double>::getDomain();
template Grid<double> & Jacobi<double>::getSource();
template Grid<double> & Jacobi<double>::getScratchDomain();
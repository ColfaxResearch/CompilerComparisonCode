#include <cmath>
#ifdef __HBM__
#include <hbwmalloc.h>
#else
#if defined(__INTEL_COMPILER)
#include <malloc.h>
#else
#include <mm_malloc.h>
#endif
#endif
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <omp.h>
#ifdef __DEBUG__
#include <cstdio>
#include <iostream>
#endif


#include "grid.hpp"


    template <typename T>
    T paddedLength(T elemLength, T blockLength) {
        return static_cast<T>(ceil((static_cast<double>(elemLength)/static_cast<double>(blockLength))))*blockLength;
    }

    template <typename GridT>
    void Grid<GridT>::allocate() {
#ifdef __DEBUG__
        std::cout << "void Grid<GridT>::allocate()" << std::endl;
#endif
#ifdef __HBM__
        hbw_posix_memalign((void**) &flatGrid, 4096, (lda*(Ny + 1) + Nx + 1 + 1)*sizeof(GridT));
#else
        flatGrid = static_cast<GridT *>(_mm_malloc((lda*(Ny + 1) + Nx + 1 + 1)*sizeof(GridT), 64));
#endif
#ifdef __DEBUG__
        printf("flatGrid allocated at %p\n", flatGrid);
        fflush(stdout);
#endif
        for (int i = 0; i < lda*(Ny + 1) + Nx + 1; ++i) {
            flatGrid[i] = static_cast<GridT>(0.0);
        }
    }

    template <typename GridT>
    void Grid<GridT>::deallocate() {
#ifdef __DEBUG__
        std::cout << "void Grid<GridT>::deallocate()" << std::endl;
#endif
#ifdef __DEBUG__
        printf("flatGrid deallocated at %p\n", flatGrid);
        fflush(stdout);
#endif
        if (flatGrid) {
#ifdef __HBM__
            hbw_free((void*)flatGrid);
#else
            _mm_free(flatGrid);
#endif
            flatGrid = nullptr;
        }
    }

    template <typename GridT>
    Grid<GridT>::Grid() {
#ifdef __DEBUG__
        std::cout << "Grid<GridT>::Grid()" << std::endl;
#endif
        flatGrid = nullptr;
        Ny = static_cast<int>(0);
        Nx = static_cast<int>(0);
        lda = static_cast<int>(0);
        Ly = static_cast<GridT>(0.0);
        Lx = static_cast<GridT>(0.0);
        Sy = static_cast<GridT>(0.0);
        Sx = static_cast<GridT>(0.0);
        Fy = static_cast<GridT>(0.0);
        Fx = static_cast<GridT>(0.0);
        Dy = static_cast<GridT>(1.0);
        Dx = static_cast<GridT>(1.0);
        DySq = static_cast<GridT>(1.0);
        DxSq = static_cast<GridT>(1.0);
        DySqDxSq = static_cast<GridT>(1.0);
        TwoDySqPlusDxSq = static_cast<GridT>(4.0);
        prec = 3;
    }

    template <typename GridT>
    Grid<GridT>::Grid(int new_N) {
#ifdef __DEBUG__
        std::cout << "Grid<GridT>::Grid(int)" << std::endl;
#endif
        flatGrid = nullptr;
        Ny = new_N;
        Nx = new_N;
        lda = paddedLength(new_N, cacheLineLength);
        Ly = static_cast<GridT>(0.0);
        Lx = static_cast<GridT>(0.0);
        Sy = static_cast<GridT>(0.0);
        Sx = static_cast<GridT>(0.0);
        Fy = static_cast<GridT>(0.0);
        Fx = static_cast<GridT>(0.0);
        Dy = static_cast<GridT>(1.0);
        Dx = static_cast<GridT>(1.0);
        DySq = static_cast<GridT>(1.0);
        DxSq = static_cast<GridT>(1.0);
        DySqDxSq = static_cast<GridT>(1.0);
        TwoDySqPlusDxSq = static_cast<GridT>(4.0);
        prec = 3;
        allocate();
    }

    template <typename GridT>
    Grid<GridT>::Grid(int new_Ny, int new_Nx) {
#ifdef __DEBUG__
        std::cout << "Grid<GridT>::Grid(int, int)" << std::endl;
#endif
        flatGrid = nullptr;
        Ny = new_Ny;
        Nx = new_Nx;
        lda = paddedLength(Nx, cacheLineLength);
        Ly = static_cast<GridT>(0.0);
        Lx = static_cast<GridT>(0.0);
        Sy = static_cast<GridT>(0.0);
        Sx = static_cast<GridT>(0.0);
        Fy = static_cast<GridT>(0.0);
        Fx = static_cast<GridT>(0.0);
        Dy = static_cast<GridT>(1.0);
        Dx = static_cast<GridT>(1.0);
        DySq = static_cast<GridT>(1.0);
        DxSq = static_cast<GridT>(1.0);
        DySqDxSq = static_cast<GridT>(1.0);
        TwoDySqPlusDxSq = static_cast<GridT>(4.0);
        prec = 3;
        allocate();
    }

    template <typename GridT>
    Grid<GridT>::Grid(Grid<GridT> const & grid) {
#ifdef __DEBUG__
        std::cout << "Grid<GridT>::Grid(Grid<GridT> const &)" << std::endl;
#endif
        flatGrid = nullptr;
        Ny = grid.Ny;
        Nx = grid.Nx;
        lda = grid.lda;
        Ly = grid.Ly;
        Lx = grid.Lx;
        Sy = grid.Sy;
        Sx = grid.Sx;
        Fy = grid.Fy;
        Fx = grid.Fx;
        Dy = grid.Dy;
        Dx = grid.Dx;
        DySq = grid.DySq;
        DxSq = grid.DxSq;
        DySqDxSq = grid.DySqDxSq;
        TwoDySqPlusDxSq = grid.TwoDySqPlusDxSq;
        prec = grid.prec;
        if ((lda > 0) and (Ny > 0)) {
            allocate();
            for (int i = 0; i < lda*(Ny + 2); ++i) {
                flatGrid[i] = grid.flatGrid[i];
            }
        }
    }

    template <typename GridT>
    Grid<GridT>::Grid(Grid<GridT> && grid) {
#ifdef __DEBUG__
        std::cout << "Grid<GridT>::Grid(Grid<GridT> const &&)" << std::endl;
#endif
        flatGrid = nullptr;
        Ny = grid.Ny;
        grid.Ny = 0;
        Nx = grid.Nx;
        grid.Nx = 0;
        lda = grid.lda;
        grid.lda = 0;
        Ly = grid.Ly;
        grid.Ly = static_cast<GridT>(0.0);
        Lx = grid.Lx;
        grid.Lx = static_cast<GridT>(0.0);
        Sy = grid.Sy;
        grid.Sy = static_cast<GridT>(0.0);
        Sx = grid.Sx;
        grid.Sx = static_cast<GridT>(0.0);
        Fy = grid.Fy;
        grid.Fy = static_cast<GridT>(0.0);
        Fx = grid.Fx;
        grid.Fx = static_cast<GridT>(0.0);
        Dy = grid.Dy;
        grid.Dy = static_cast<GridT>(1.0);
        Dx = grid.Dx;
        grid.Dx = static_cast<GridT>(1.0);
        DySq = grid.DySq;
        grid.DySq = static_cast<GridT>(1.0);
        DxSq = grid.DxSq;
        grid.DxSq = static_cast<GridT>(1.0);
        DySqDxSq = grid.DySqDxSq;
        grid.DySqDxSq = static_cast<GridT>(1.0);
        TwoDySqPlusDxSq = grid.TwoDySqPlusDxSq;
        grid.TwoDySqPlusDxSq = static_cast<GridT>(4.0);
        prec = grid.prec;
        grid.prec = 3;
        flatGrid = grid.flatGrid;
        grid.flatGrid = nullptr;
    }

    template <typename GridT>
    Grid<GridT> & Grid<GridT>::operator=(Grid<GridT> const & grid) {
#ifdef __DEBUG__
        std::cout << "Grid<GridT> & Grid<GridT>::operator=(Grid<GridT> const &)" << std::endl;
#endif
        if (this != &grid) {
            deallocate();
            Ny = grid.Ny;
            Nx = grid.Nx;
            lda = grid.lda;
            prec = grid.prec;
            Ly = grid.Ly;
            Lx = grid.Lx;
            Sy = grid.Sy;
            Sx = grid.Sx;
            Fy = grid.Fy;
            Fx = grid.Fx;
            Dy = grid.Dy;
            Dx = grid.Dx;
            DySq = grid.DySq;
            DxSq = grid.DxSq;
            DySqDxSq = grid.DySqDxSq;
            TwoDySqPlusDxSq = grid.TwoDySqPlusDxSq;
            if ((lda > 0) and (Ny > 0)) {
                allocate();
                for (int i = 0; i < lda*(Ny + 2); ++i) {
                    flatGrid[i] = grid.flatGrid[i];
                }
            }
        }
        return *this;
    }

    template <typename GridT>
    Grid<GridT> & Grid<GridT>::operator=(Grid<GridT> && grid) {
#ifdef __DEBUG__
        std::cout << "Grid<GridT> & Grid<GridT>::operator=(Grid<GridT> &&)" << std::endl;
#endif
        if (this != &grid) {
            deallocate();
            Ny = grid.Ny;
            grid.Ny = 0;
            Nx = grid.Nx;
            grid.Nx = 0;
            lda = grid.lda;
            grid.lda = 0;
            Ly = grid.Ly;
            grid.Ly = static_cast<GridT>(0.0);
            Lx = grid.Lx;
            grid.Lx = static_cast<GridT>(0.0);
            Sy = grid.Sy;
            grid.Sy = static_cast<GridT>(0.0);
            Sx = grid.Sx;
            grid.Sx = static_cast<GridT>(0.0);
            Fy = grid.Fy;
            grid.Fy = static_cast<GridT>(0.0);
            Fx = grid.Fx;
            grid.Fx = static_cast<GridT>(0.0);
            Dy = grid.Dy;
            grid.Dy = static_cast<GridT>(1.0);
            Dx = grid.Dx;
            grid.Dx = static_cast<GridT>(1.0);
            DySq = grid.DySq;
            grid.DySq = static_cast<GridT>(1.0);
            DxSq = grid.DxSq;
            grid.DxSq = static_cast<GridT>(1.0);
            DySqDxSq = grid.DySqDxSq;
            grid.DySqDxSq = static_cast<GridT>(1.0);
            TwoDySqPlusDxSq = grid.TwoDySqPlusDxSq;
            grid.TwoDySqPlusDxSq = static_cast<GridT>(4.0);
            prec = grid.prec;
            grid.prec = 3;
            flatGrid = grid.flatGrid;
            grid.flatGrid = nullptr;
        }
        return *this;
    }

    template <typename GridT>
    Grid<GridT>::~Grid() {
#ifdef __DEBUG__
        std::cout << "Grid<GridT>::~Grid()" << std::endl;
#endif
        deallocate();
        Ny = static_cast<int>(0);
        Nx = static_cast<int>(0);
        lda = static_cast<int>(0);
        Ly = static_cast<GridT>(0.0);
        Lx = static_cast<GridT>(0.0);
        Sy = static_cast<GridT>(0.0);
        Sx = static_cast<GridT>(0.0);
        Fy = static_cast<GridT>(0.0);
        Fx = static_cast<GridT>(0.0);
        Dy = static_cast<GridT>(1.0);
        Dx = static_cast<GridT>(1.0);
        DySq = static_cast<GridT>(1.0);
        DxSq = static_cast<GridT>(1.0);
        DySqDxSq = static_cast<GridT>(1.0);
        TwoDySqPlusDxSq = static_cast<GridT>(4.0);
        prec = static_cast<int>(3);
    }

    template <typename GridT>
    void Grid<GridT>::swapStorage(Grid<GridT> & grid) {
//#ifdef __DEBUG__
//        std::cout << "void Grid<GridT>::swapStorage(Grid<GridT> & grid)" << std::endl;
//#endif
        GridT * swapPtr;
        swapPtr = flatGrid;
        flatGrid = grid.flatGrid;
        grid.flatGrid = swapPtr;
    }

    template <typename GridT>
    void Grid<GridT>::setPrecision(int p) {
        prec = p;
    }

    template <typename GridT>
    int Grid<GridT>::getPrecision() const {
        return prec;
    }

    template <typename GridT>
    GridT Grid<GridT>::get_Sy() const {
        return Sy;
    }

    template <typename GridT>
    GridT Grid<GridT>::get_Sx() const {
        return Sx;
    }

    template <typename GridT>
    void Grid<GridT>::set_Sy(GridT new_Sy) {
        Sy = new_Sy;
        Ly = Fy - Sy;
        Dy = (Fy - Sy)/Ny;
        DySq = Dy*Dy;
        DySqDxSq = DySq*DxSq;
        TwoDySqPlusDxSq = 2.0*(DySq + DxSq);
    }

    template <typename GridT>
    void Grid<GridT>::set_Sx(GridT new_Sx) {
        Sx = new_Sx;
        Lx = Fx - Sx;
        Dx = (Fx - Sx)/Nx;
        DxSq = Dx*Dx;
        DySqDxSq = DySq*DxSq;
        TwoDySqPlusDxSq = 2.0*(DySq + DxSq);
    }

    template <typename GridT>
    GridT Grid<GridT>::get_Fy() const {
        return Fy;
    }

    template <typename GridT>
    GridT Grid<GridT>::get_Fx() const {
        return Fx;
    }

    template <typename GridT>
    void Grid<GridT>::set_Fy(GridT new_Fy) {
        Fy = new_Fy;
        Ly = Fy - Sy;
        Dy = (Fy - Sy)/Ny;
        DySq = Dy*Dy;
        DySqDxSq = DySq*DxSq;
        TwoDySqPlusDxSq = 2.0*(DySq + DxSq);
    }

    template <typename GridT>
    void Grid<GridT>::set_Fx(GridT new_Fx) {
        Fx = new_Fx;
        Lx = Fx - Sx;
        Dx = (Fx - Sx)/Nx;
        DxSq = Dx*Dx;
        DySqDxSq = DySq*DxSq;
        TwoDySqPlusDxSq = 2.0*(DySq + DxSq);
    }

    template <typename GridT>
    GridT Grid<GridT>::get_Ly() const {
        return Ly;
    }

    template <typename GridT>
    GridT Grid<GridT>::get_Lx() const {
        return Lx;
    }

    template <typename GridT>
    GridT Grid<GridT>::get_Dy() const {
        return Dy;
    }

    template <typename GridT>
    GridT Grid<GridT>::get_Dx() const {
        return Dx;
    }

    template <typename GridT>
    int Grid<GridT>::get_Ny() const {
        return Ny;
    }

    template <typename GridT>
    int Grid<GridT>::get_Nx() const {
        return Nx;
    }

    template <typename GridT>
    int Grid<GridT>::get_lda() const {
        return lda;
    }

    template <typename GridT>
    void Grid<GridT>::resize(int new_Ny, int new_Nx) {
#ifdef __DEBUG__
        std::cout << "void Grid<GridT>::resize(int, int)" << std::endl;
#endif
        deallocate();
        Ny = new_Ny;
        Nx = new_Nx;
        lda = paddedLength(Nx, cacheLineLength);
        allocate();
        Dy = (Fy - Sy)/Ny;
        Dx = (Fx - Sx)/Nx;
        DySq = Dy*Dy;
        DxSq = Dx*Dx;
        DySqDxSq = DySq*DxSq;
        TwoDySqPlusDxSq = 2.0*(DySq + DxSq);
    }

    template <typename GridT>
    GridT Grid<GridT>::get_y(int row, int col) const {
        return row*Dy + Sy;
    }

    template <typename GridT>
    GridT Grid<GridT>::get_x(int row, int col) const {
        return col*Dx + Sx;
    }

    template <typename GridT>
    GridT & Grid<GridT>::operator()(int row, int col) const {
        return flatGrid[row*lda + col];
    }

    template <typename GridT>
    void Grid<GridT>::setGrid(GridT gridVal) {
#ifdef __DEBUG__
        std::cout << "void Grid<GridT>::setGrid(GridT gridVal)" << std::endl;
#endif
        for (int i = 0; i < lda*(Ny + 2); ++i) {
            flatGrid[i] = gridVal;
        }
    }

    template <typename GridT>
    void Grid<GridT>::setBoundary(GridT bdryVal) {
#ifdef __DEBUG__
        std::cout << "void Grid<GridT>::setBoundary(GridT bdryVal)" << std::endl;
#endif
        flatGrid[0] = bdryVal;
        flatGrid[Nx + 1] = bdryVal;
        flatGrid[(Ny + 1)*lda] = bdryVal;
        flatGrid[(Ny + 1)*lda + Ny + 1] = bdryVal;
        for (int i = 1; i <= Nx; ++i) {
            flatGrid[i] = bdryVal;
            flatGrid[(Nx + 1)*lda + i] = bdryVal;
        }
        for (int i = 1; i <= Ny; ++i) {
            flatGrid[i*lda] = bdryVal;
            flatGrid[i*lda + Ny + 1] = bdryVal;
        }
    }

template int paddedLength(int, int);

template class Grid<float>;
template class Grid<double>;

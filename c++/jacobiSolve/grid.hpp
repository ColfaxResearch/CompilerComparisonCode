#ifndef GRID_HPP
#define GRID_HPP

#include <cstddef>
#include <iostream>
#include <iomanip>


template <typename T>
T paddedLength(T elemLength, T blockLength);

template <typename GridT>
class Grid {
protected:
    int lda, Ny, Nx, prec;
    static int const vecLength = 8;
    static int const cacheLineLength = 8;
    GridT Sy, Fy, Ly, Sx, Fx, Lx, Dy, Dx, DySq, DxSq, DySqDxSq, TwoDySqPlusDxSq, *flatGrid;
    void allocate();
    void deallocate();
    int calcLDA();
public:
    Grid();
    Grid(int new_N);
    Grid(int new_Ny, int new_Nx);
    Grid(Grid<GridT> const & grid);
    Grid(Grid<GridT> && grid);
    Grid<GridT> & operator=(Grid<GridT> const & grid);
    Grid<GridT> & operator=(Grid<GridT> && grid);
    ~Grid();
    void swapStorage(Grid<GridT> & grid);
    void setPrecision(int p);
    int getPrecision() const;
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
    GridT get_y(int y, int x) const;
    GridT get_x(int y, int x) const;
#ifndef __PGI
    #pragma omp declare simd
#endif
    GridT & operator()(int y, int x) const;
    void setGrid(GridT gridVal);
    void setBoundary(GridT bdryVal);

    friend std::ostream& operator<<(std::ostream& os, Grid const & grid) {
        os.precision(grid.getPrecision());
        if ((grid.Nx > 0) and (grid.Ny > 0) and (grid.lda > 0)) {
            for (int row = 0; row <= grid.Ny; ++row) {
                for (int col = 0; col <= grid.Nx + 1; ++col) {
                    os << std::showpos << std::scientific << grid(row, col) << " ";
                }
                os << std::endl;
            }
            for (int col = 0; col <= grid.Nx + 1; ++col) {
                os << std::showpos << std::scientific << grid((grid.Ny + 1), col) << " ";
            }
        }
        return os;
    }

};

#endif // GRID_HPP
#ifndef GRID_HPP
#define GRID_HPP

//#include "RowMajorStorage.hpp"

  template <typename GridT, template <typename > typename StorageT>
  class Grid {
  private:
    size_t _Ny, _Nx, _precision;
    GridT _Sy, _Fy, _Ly, _Sx, _Fx, _Lx, _Dy, _Dx, _DySq, _DxSq, _DySqDxSq, _TwoDySqPlusDxSq;
    StorageT<GridT> * _storage;
  public:
    Grid();
    Grid(size_t Size);
    Grid(size_t new_Ny, size_t new_Nx);
    Grid(Grid<GridT, StorageT> const & g);
    Grid(Grid<GridT, StorageT> const && g);
    Grid<GridT, StorageT> & operator=(Grid<GridT, StorageT> const & grid);
    Grid<GridT, StorageT> & operator=(Grid<GridT, StorageT> const && grid);
    virtual ~Grid();
    StorageT<GridT> & storage() const;
    size_t precision() const;
    void precision(size_t precision);
    GridT Sy() const;
    GridT Sx() const;
    void Sy(GridT new_Sy);
    void Sx(GridT new_Sx);
    GridT Fy() const;
    GridT Fx() const;
    void Fy(GridT new_Fy);
    void Fx(GridT new_Fx);
    GridT Ly() const;
    GridT Lx() const;
    GridT Dy() const;
    GridT Dx() const;
    size_t Ny() const;
    size_t Nx() const;
    size_t lda() const ;
    virtual void resize(size_t new_Ny, size_t new_Nx);
    GridT get_y(size_t y, size_t x) const;
    GridT get_x(size_t y, size_t x) const;
    virtual GridT & operator()(size_t y, size_t x) const;
    GridT * flatGrid();
    void setBoundary(GridT bdryVal);

    friend ostream& operator<<(ostream& os, Grid<GridT, StorageT> const & grid) {
#ifdef __DEBUG__
      cout << "friend ostream& operator<<(ostream&, Grid<GridT, StorageT> const &)" << endl;
#endif
      os.precision(grid.precision());
      for (size_t row = 0; row < grid.Ny() + 1; ++row) {
        for (size_t col = 0; col < grid.Nx() + 2; ++col) {
          os << showpos << scientific << grid(row, col) << " ";
        }
        os << endl;
      }
      for (size_t col = 0; col < grid.Nx() + 2; ++col) {
        os << showpos << scientific << grid((grid.Ny() + 1), col) << " ";
      }
      return os;
    }
  };

  template <typename GridT, template <typename > typename StorageT>
  Grid<GridT, StorageT>::Grid() {
#ifdef __DEBUG__
    cout << "Grid<GridT, StorageT>::Grid()" << endl;
#endif
    _Ny = static_cast<size_t>(0);
    _Nx = static_cast<size_t>(0);
    _Sy = static_cast<GridT>(0.0);
    _Fy = static_cast<GridT>(0.0);
    _Ly = static_cast<GridT>(0.0);
    _Sx = static_cast<GridT>(0.0);
    _Fx = static_cast<GridT>(0.0);
    _Lx = static_cast<GridT>(0.0);
    _Dy = static_cast<GridT>(0.0);
    _Dx = static_cast<GridT>(0.0);
    _DySq = static_cast<GridT>(0.0);
    _DxSq = static_cast<GridT>(0.0);
    _DySqDxSq = static_cast<GridT>(0.0);
    _TwoDySqPlusDxSq = static_cast<GridT>(0.0);
    _storage = nullptr;
  }

  template <typename GridT, template <typename > typename StorageT>
  Grid<GridT, StorageT>::Grid(size_t new_N) {
#ifdef __DEBUG__
    cout << "Grid<GridT, StorageT>::Grid(size_t new_N)" << endl;
#endif
    _Ny = static_cast<size_t>(new_N);
    _Nx = static_cast<size_t>(new_N);
    _Sy = static_cast<GridT>(0.0);
    _Fy = static_cast<GridT>(0.0);
    _Ly = static_cast<GridT>(0.0);
    _Sx = static_cast<GridT>(0.0);
    _Fx = static_cast<GridT>(0.0);
    _Lx = static_cast<GridT>(0.0);
    _Dy = static_cast<GridT>(0.0);
    _Dx = static_cast<GridT>(0.0);
    _DySq = static_cast<GridT>(0.0);
    _DxSq = static_cast<GridT>(0.0);
    _DySqDxSq = static_cast<GridT>(0.0);
    _TwoDySqPlusDxSq = static_cast<GridT>(0.0);
    _storage = new StorageT<GridT>(new_N);
  }

  template <typename GridT, template <typename > typename StorageT>
  Grid<GridT, StorageT>::Grid(size_t new_Ny, size_t new_Nx) {
#ifdef __DEBUG__
    cout << "Grid<GridT, StorageT>::Grid(size_t new_Ny, size_t new_Nx)" << endl;
#endif
    _Ny = static_cast<size_t>(new_Ny);
    _Nx = static_cast<size_t>(new_Nx);
    _Sy = static_cast<GridT>(0.0);
    _Fy = static_cast<GridT>(0.0);
    _Ly = static_cast<GridT>(0.0);
    _Sx = static_cast<GridT>(0.0);
    _Fx = static_cast<GridT>(0.0);
    _Lx = static_cast<GridT>(0.0);
    _Dy = static_cast<GridT>(0.0);
    _Dx = static_cast<GridT>(0.0);
    _DySq = static_cast<GridT>(0.0);
    _DxSq = static_cast<GridT>(0.0);
    _DySqDxSq = static_cast<GridT>(0.0);
    _TwoDySqPlusDxSq = static_cast<GridT>(0.0);
    _storage = new StorageT<GridT>(new_Ny, new_Nx);
  }

  template <typename GridT, template <typename > typename StorageT>
  Grid<GridT, StorageT>::Grid(Grid<GridT, StorageT> const & g) {
#ifdef __DEBUG__
    cout << "Grid<GridT, StorageT>::Grid(Grid<GridT, StorageT> const & g)" << endl;
#endif
    _Ny = g._Ny;
    _Nx = g._Nx;
    _Sy = g._Sy;
    _Fy = g._Fy;
    _Ly = g._Ly;
    _Sx = g._Sx;
    _Fx = g._Fx;
    _Lx = g._Lx;
    _Dy = g._Dy;
    _Dx = g._Dx;
    _DySq = g._DySq;
    _DxSq = g._DxSq;
    _DySqDxSq = g._DySqDxSq;
    _TwoDySqPlusDxSq = g._TwoDySqPlusDxSq;
    _storage = new StorageT<GridT>(this->_Ny, this->_Nx);
    for (size_t i = 0; i < (*this)._Ny; ++i) {
      for (size_t j = 0; j < (*this)._Nx; ++j) {
        (*_storage)(i, j) = (*(g._storage))(i, j);
      }
    }
  }

  template <typename GridT, template <typename > typename StorageT>
  Grid<GridT, StorageT>::Grid(Grid<GridT, StorageT> const && g) {
#ifdef __DEBUG__
    cout << "Grid<GridT, StorageT>::Grid(Grid<GridT, StorageT> const && g)" << endl;
#endif
    _Ny = g._Ny;
    _Nx = g._Nx;
    _Sy = g._Sy;
    _Fy = g._Fy;
    _Ly = g._Ly;
    _Sx = g._Sx;
    _Fx = g._Fx;
    _Lx = g._Lx;
    _Dy = g._Dy;
    _Dx = g._Dx;
    _DySq = g._DySq;
    _DxSq = g._DxSq;
    _DySqDxSq = g._DySqDxSq;
    _TwoDySqPlusDxSq = g._TwoDySqPlusDxSq;
    _storage = new StorageT<GridT>(this->_Ny, this->_Nx);
    for (size_t i = 0; i < (*this)._Ny; ++i) {
      for (size_t j = 0; j < (*this)._Nx; ++j) {
        (*_storage)(i, j) = (*(g._storage))(i, j);
      }
    }
    g._Ny = static_cast<size_t>(0);
    g._Nx = static_cast<size_t>(0);
    g._Sy = static_cast<GridT>(0.0);
    g._Fy = static_cast<GridT>(0.0);
    g._Ly = static_cast<GridT>(0.0);
    g._Sx = static_cast<GridT>(0.0);
    g._Fx = static_cast<GridT>(0.0);
    g._Lx = static_cast<GridT>(0.0);
    g._Dy = static_cast<GridT>(0.0);
    g._Dx = static_cast<GridT>(0.0);
    g._DySq = static_cast<GridT>(0.0);
    g._DxSq = static_cast<GridT>(0.0);
    g._DySqDxSq = static_cast<GridT>(0.0);
    g._TwoDySqPlusDxSq = static_cast<GridT>(0.0);
    delete g._storage;
    g._storage = nullptr;
  }

  template <typename GridT, template <typename > typename StorageT>
  Grid<GridT, StorageT> & Grid<GridT, StorageT>::operator=(Grid<GridT, StorageT> const & g) {
#ifdef __DEBUG__
    cout << "Grid<GridT, StorageT> & Grid<GridT, StorageT>::operator=(Grid<GridT, StorageT> const & g)" << endl;
#endif
    _Ny = g._Ny;
    _Nx = g._Nx;
    _Sy = g._Sy;
    _Fy = g._Fy;
    _Ly = g._Ly;
    _Sx = g._Sx;
    _Fx = g._Fx;
    _Lx = g._Lx;
    _Dy = g._Dy;
    _Dx = g._Dx;
    _DySq = g._DySq;
    _DxSq = g._DxSq;
    _DySqDxSq = g._DySqDxSq;
    _TwoDySqPlusDxSq = g._TwoDySqPlusDxSq;
    delete _storage;
    _storage = new StorageT<GridT>(this->_Ny, this->_Nx);
    for (size_t i = 0; i < (*this)._Ny; ++i) {
      for (size_t j = 0; j < (*this)._Nx; ++j) {
        (*_storage)(i, j) = (*(g._storage))(i, j);
      }
    }
  }

  template <typename GridT, template <typename > typename StorageT>
  Grid<GridT, StorageT> & Grid<GridT, StorageT>::operator=(Grid<GridT, StorageT> const && g) {
#ifdef __DEBUG__
    cout << "Grid<GridT, StorageT> & Grid<GridT, StorageT>::operator=(Grid<GridT, StorageT> const && g)" << endl;
#endif
    _Ny = g._Ny;
    _Nx = g._Nx;
    _Sy = g._Sy;
    _Fy = g._Fy;
    _Ly = g._Ly;
    _Sx = g._Sx;
    _Fx = g._Fx;
    _Lx = g._Lx;
    _Dy = g._Dy;
    _Dx = g._Dx;
    _DySq = g._DySq;
    _DxSq = g._DxSq;
    _DySqDxSq = g._DySqDxSq;
    _TwoDySqPlusDxSq = g._TwoDySqPlusDxSq;
    delete _storage;
    _storage = new StorageT<GridT>(this->_Ny, this->_Nx);
    for (size_t i = 0; i < (*this)._Ny; ++i) {
      for (size_t j = 0; j < (*this)._Nx; ++j) {
        (*_storage)(i, j) = (*(g._storage))(i, j);
      }
    }
    g._Ny = static_cast<size_t>(0);
    g._Nx = static_cast<size_t>(0);
    g._Sy = static_cast<GridT>(0.0);
    g._Fy = static_cast<GridT>(0.0);
    g._Ly = static_cast<GridT>(0.0);
    g._Sx = static_cast<GridT>(0.0);
    g._Fx = static_cast<GridT>(0.0);
    g._Lx = static_cast<GridT>(0.0);
    g._Dy = static_cast<GridT>(0.0);
    g._Dx = static_cast<GridT>(0.0);
    g._DySq = static_cast<GridT>(0.0);
    g._DxSq = static_cast<GridT>(0.0);
    g._DySqDxSq = static_cast<GridT>(0.0);
    g._TwoDySqPlusDxSq = static_cast<GridT>(0.0);
    delete g._storage;
    g._storage = nullptr;
  }

  template <typename GridT, template <typename > typename StorageT>
  Grid<GridT, StorageT>::~Grid() {
#ifdef __DEBUG__
    cout << "Grid<GridT, StorageT>::~Grid()" << endl;
#endif
    _Ny = static_cast<size_t>(0);
    _Nx = static_cast<size_t>(0);
    _Sy = static_cast<GridT>(0.0);
    _Fy = static_cast<GridT>(0.0);
    _Ly = static_cast<GridT>(0.0);
    _Sx = static_cast<GridT>(0.0);
    _Fx = static_cast<GridT>(0.0);
    _Lx = static_cast<GridT>(0.0);
    _Dy = static_cast<GridT>(0.0);
    _Dx = static_cast<GridT>(0.0);
    _DySq = static_cast<GridT>(0.0);
    _DxSq = static_cast<GridT>(0.0);
    _DySqDxSq = static_cast<GridT>(0.0);
    _TwoDySqPlusDxSq = static_cast<GridT>(0.0);
    delete _storage;
    _storage = nullptr;
  }

  template <typename GridT, template <typename > typename StorageT>
  StorageT<GridT> & Grid<GridT, StorageT>::storage() const {
    return *_storage;
  }

  template <typename GridT, template <typename > typename StorageT>
  size_t Grid<GridT, StorageT>::Ny() const {
    return _Ny;
  }

  template <typename GridT, template <typename > typename StorageT>
  size_t Grid<GridT, StorageT>::Nx() const {
    return _Nx;
  }

  template <typename GridT, template <typename > typename StorageT>
  void Grid<GridT, StorageT>::precision(size_t new_precision) {
    _precision = new_precision;
  }

  template <typename GridT, template <typename > typename StorageT>
  size_t Grid<GridT, StorageT>::precision() const {
    return _precision;
  }

  template <typename GridT, template <typename > typename StorageT>
  GridT Grid<GridT, StorageT>::Sy() const {
#ifdef __DEBUG__
    cout << "GridT Grid<GridT, StorageT>::Sy() const" << endl;
#endif
    return _Sy;
  }

  template <typename GridT, template <typename > typename StorageT>
  GridT Grid<GridT, StorageT>::Sx() const {
#ifdef __DEBUG__
    cout << "GridT Grid<GridT, StorageT>::Sx() const" << endl;
#endif
    return _Sx;
  }

  template <typename GridT, template <typename > typename StorageT>
  void Grid<GridT, StorageT>::Sy(GridT new_Sy) {
#ifdef __DEBUG__
    cout << "void Grid<GridT, StorageT>::Sy(GridT)" << endl;
#endif
    _Sy = new_Sy;
    _Ly = _Fy - _Sy;
    _Dy = _Ly/storage.Ny();
    _DySq = _Dy*_Dy;
    _DySqDxSq = _DySq*_DxSq;
    _TwoDySqPlusDxSq = static_cast<GridT>(2.0)*(_DySq + _DxSq);
  }

  template <typename GridT, template <typename > typename StorageT>
  void Grid<GridT, StorageT>::Sx(GridT new_Sx) {
#ifdef __DEBUG__
    cout << "void Grid<GridT, StorageT>::Sx(GridT)" << endl;
#endif
    _Sx = new_Sx;
    _Lx = _Fx - _Sx;
    _Dx = _Lx/storage.Nx();
    _DxSq = _Dx*_Dx;
    _DySqDxSq = _DySq*_DxSq;
    _TwoDySqPlusDxSq = static_cast<GridT>(2.0)*(_DySq + _DxSq);
  }

  template <typename GridT, template <typename > typename StorageT>
  GridT Grid<GridT, StorageT>::Fy() const {
#ifdef __DEBUG__
    cout << "GridT Grid<GridT, StorageT>::Fx() const" << endl;
#endif
    return _Fy;
  }

  template <typename GridT, template <typename > typename StorageT>
  GridT Grid<GridT, StorageT>::Fx() const {
#ifdef __DEBUG__
    cout << "GridT Grid<GridT, StorageT>::Fx() const" << endl;
#endif
    return _Fx;
  }

  template <typename GridT, template <typename > typename StorageT>
  void Grid<GridT, StorageT>::Fy(GridT new_Fy) {
  #ifdef __DEBUG__
    cout << "void Grid<GridT, StorageT>::Fx(GridT)" << endl;
  #endif
    _Fy = new_Fy;
    _Ly = _Fy - _Sy;
    _Dy = _Ly/storage.Ny();
    _DySq = _Dy*_Dy;
    _DySqDxSq = _DySq*_DxSq;
    _TwoDySqPlusDxSq = static_cast<GridT>(2.0)*(_DySq + _DxSq);
  }

  template <typename GridT, template <typename > typename StorageT>
  void Grid<GridT, StorageT>::Fx(GridT new_Fx) {
#ifdef __DEBUG__
    cout << "void Grid<GridT, StorageT>::Fx(GridT)" << endl;
#endif
    _Fx = new_Fx;
    _Lx = _Fx - _Sx;
    _Dx = _Lx/storage.Nx();
    _DxSq = _Dx*_Dx;
    _DySqDxSq = _DySq*_DxSq;
    _TwoDySqPlusDxSq = static_cast<GridT>(2.0)*(_DySq + _DxSq);
  }

  template <typename GridT, template <typename > typename StorageT>
  GridT Grid<GridT, StorageT>::Ly() const {
#ifdef __DEBUG__
    cout << "GridT Grid<GridT, StorageT>::Ly() const" << endl;
#endif
    return _Ly;
  }

  template <typename GridT, template <typename > typename StorageT>
  GridT Grid<GridT, StorageT>::Lx() const {
#ifdef __DEBUG__
    cout << "GridT Grid<GridT, StorageT>::Lx() const" << endl;
#endif
    return _Lx;
  }

  template <typename GridT, template <typename > typename StorageT>
  GridT Grid<GridT, StorageT>::Dy() const {
  #ifdef __DEBUG__
    cout << "GridT Grid<GridT, StorageT>::Dy() const" << endl;
  #endif
    return _Dy;
  }

  template <typename GridT, template <typename > typename StorageT>
  GridT Grid<GridT, StorageT>::Dx() const {
  #ifdef __DEBUG__
    cout << "GridT Grid<GridT, StorageT>::Dx() const" << endl;
  #endif
    return _Dx;
  }

  template <typename GridT, template <typename > typename StorageT>
  void Grid<GridT, StorageT>::resize(size_t new_Ny, size_t new_Nx) {
    (*_storage).resize(new_Ny, new_Nx);
  }

  template <typename GridT, template <typename > typename StorageT>
  GridT & Grid<GridT, StorageT>::operator()(size_t y, size_t x) const {
#ifdef __DEBUG__
    cout << "GridT & Grid<GridT, StorageT>::operator()(size_t, size_t) const" << endl;
#endif
    cout << "Requesting grid(" << y << "," << x << ")" << endl;
    return (*_storage)(y, x);
  }

  template <typename GridT, template<typename > typename StorageT>
  GridT * Grid<GridT, StorageT>::flatGrid() {
#ifdef __DEBUG__
    cout << "GridT * Grid<GridT, StorageT>::flatGrid()" << endl;
#endif
    return _storage->flatGrid();
  }

  template <typename GridT, template<typename > typename StorageT>
  void Grid<GridT, StorageT>::setBoundary(GridT bdryVal) {
    (*_storage)(0,0) = bdryVal;
    (*_storage)(0, _Nx + 1) = bdryVal;
    (*_storage)(_Ny + 1, 0) = bdryVal;
    (*_storage)(_Ny + 1, _Ny + 1) = bdryVal;
    for (size_t i = 1; i <= _Nx; ++i) {
      (*_storage)(0, i) = bdryVal;
      (*_storage)(_Ny + 1, i) = bdryVal;
    }
    for (size_t i = 1; i <= _Ny; ++i) {
      (*_storage)(i, 0) = bdryVal;
      (*_storage)(i, _Nx + 1) = bdryVal;
    }
  }

#endif // GRID_HPP

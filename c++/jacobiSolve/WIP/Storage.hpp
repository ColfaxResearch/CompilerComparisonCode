#ifndef STORAGE_HPP
#define STORAGE_HPP

#include <cmath>
#ifdef __INTEL_COMPILER
#include <mathimf.h>
#endif

  using namespace std;

  template <typename T>
  T paddedLength(T elemLength, T blockLength) {
    return static_cast<T>(ceil((static_cast<double>(elemLength)/static_cast<double>(blockLength))))*blockLength;
  }


  template <typename GridT>
  class Storage {
  protected:
    GridT * _flatGrid;
    size_t _Ny, _Nx;
    static size_t const vecLength = 8;
    static size_t const cacheLineLength = 8;
  public:
    Storage();
    virtual ~Storage();
    size_t Ny();
    size_t Nx();
    virtual void resize(size_t new_Ny, size_t new_Nx) = 0;
    virtual GridT & operator()(size_t y, size_t x) const = 0;
    GridT * flatGrid();
  };

  template <typename GridT>
  Storage<GridT>::Storage() {
#ifdef __DEBUG__
    cout << "Storage<GridT>::Storage()" << endl;
#endif
    _Ny = 0;
    _Nx = 0;
    _flatGrid = nullptr;
  }

  template <typename GridT>
  Storage<GridT>::~Storage() {
#ifdef __DEBUG__
    cout << "Storage<GridT>::~Storage()" << endl;
#endif
    _Ny = 0;
    _Nx = 0;
    _flatGrid = nullptr;
  }

  template <typename GridT>
  size_t Storage<GridT>::Ny() {
#ifdef __DEBUG__
    cout << "Storage<GridT>::Ny()" << endl;
#endif
    return _Ny;
  }

  template <typename GridT>
  size_t Storage<GridT>::Nx() {
#ifdef __DEBUG__
    cout << "Storage<GridT>::Nx()" << endl;
#endif
    return _Nx;
  }

  template <typename GridT>
  GridT * Storage<GridT>::flatGrid() {
#ifdef __DEBUG__
    cout << "GridT * Storage<GridT>::flatGrid()" << endl;
#endif
    return _flatGrid;
  }

#endif // STORAGE_HPP

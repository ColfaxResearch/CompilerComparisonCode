#ifndef ROWMAJORSTORAGE_HPP
#define ROWMAJORSTORAGE_HPP

#include "LinearContiguousStorage.hpp"

  template <typename GridT>
  class RowMajorStorage : public LinearContiguousStorage<GridT> {
  public:
    RowMajorStorage();
    RowMajorStorage(size_t new_N);
    RowMajorStorage(size_t new_Nx, size_t new_Ny);
    RowMajorStorage(LinearContiguousStorage<GridT> const & storage);
    RowMajorStorage(LinearContiguousStorage<GridT> const && storage);
    RowMajorStorage<GridT> & operator=(LinearContiguousStorage<GridT> const & storage);
    RowMajorStorage<GridT> & operator=(LinearContiguousStorage<GridT> const && storage);
    virtual ~RowMajorStorage();
    virtual void resize(size_t new_Ny, size_t new_Nx);
    virtual GridT & operator()(size_t y, size_t x) const;
    friend ostream& operator<<(ostream& os, RowMajorStorage<GridT> const & store) {
#ifdef __DEBUG__
      cout << "friend ostream& operator<<(ostream& os, RowMajorStorage<GridT> const & store)" << endl;
#endif
      if ((store.Nx() > 0) and (store.Ny() > 0)) {
        for (size_t row = 0; row <= store.Ny(); ++row) {
          for (size_t col = 0; col <= store.Nx() + 1; ++col) {
            os << showpos << scientific << store.flatGrid()[row*store.lda() + col] << " ";
          }
          os << endl;
        }
        for (size_t col = 0; col <= store.Nx() + 1; ++col) {
          os << showpos << scientific << store.flatGrid()[(store.Ny() + 1)*store.lda() + col] << " ";
        }
      }
      return os;
    }
  };

  template <typename GridT>
  RowMajorStorage<GridT>::RowMajorStorage() {
#ifdef __DEBUG__
    cout << "RowMajorStorage<GridT>::RowMajorStorage()" << endl;
#endif
    this->_Ny = static_cast<size_t>(0);
    this->_Nx = static_cast<size_t>(0);
    this->_lda = static_cast<size_t>(0);
  }

  template <typename GridT>
  RowMajorStorage<GridT>::RowMajorStorage(size_t new_N) {
#ifdef __DEBUG__
    cout << "RowMajorStorage<GridT>::RowMajorStorage(size_t)" << endl;
#endif
    this->_Ny = new_N;
    this->_Nx = new_N;
    this->_lda = paddedLength(this->_Nx, this->cacheLineLength);
    cout << "_lda = " << this->_lda << endl;
#ifdef __HBM__
    hbw_posix_memalign((void**) &this->_flatGrid, 4096, this->_lda*(this->_Ny + 2)*sizeof(GridT));
#else
#ifdef __INTEL_COMPILER
    this->_flatGrid = static_cast<GridT *>(_mm_malloc(this->_lda*(this->_Ny + 2)*sizeof(GridT), 64));
#else
    posix_memalign((void **) &this->_flatGrid, 4096, this->_lda*(this->_Ny + 2)*sizeof(GridT));
#endif
#endif
    cout << "_flatGrd allocated at " << this->_flatGrid << endl;
    for (size_t i = 0; i < this->_lda*(this->_Ny + 2); ++i) {
      this->_flatGrid[i] = static_cast<GridT>(0.0);
    }
  }

  template <typename GridT>
  RowMajorStorage<GridT>::RowMajorStorage(size_t new_Ny, size_t new_Nx) {
#ifdef __DEBUG__
    cout << "RowMajorStorage<GridT>::RowMajorStorage(size_t, size_t)" << endl;
#endif
    this->_Ny = new_Ny;
    this->_Nx = new_Nx;
    this->_lda = paddedLength(this->_Nx, this->cacheLineLength);
#ifdef __HBM__
    hbw_posix_memalign((void**) &this->_flatGrid, 4096, this->_lda*(this->_Ny + 2)*sizeof(GridT));
#else
#ifdef __INTEL_COMPILER
    this->_flatGrid = static_cast<GridT *>(_mm_malloc(this->_lda*(this->_Ny + 2)*sizeof(GridT), 64));
#else
    posix_memalign((void **) &this->_flatGrid, 4096, this->_lda*(this->_Ny + 2)*sizeof(GridT));
#endif
#endif
    for (size_t i = 0; i < this->_lda*(this->_Ny + 2); ++i) {
      this->_flatGrid[i] = static_cast<GridT>(0.0);
    }
  }

  template <typename GridT>
  RowMajorStorage<GridT>::RowMajorStorage(LinearContiguousStorage<GridT> const & storage) {
#ifdef __DEBUG__
    cout << "RowMajorStorage<GridT>::RowMajorStorage(LinearContiguousStorage<GridT> const &)" << endl;
#endif
    this->_Ny = storage._Ny;
    this->_Nx = storage._Nx;
    this->_lda = storage._lda;
#ifdef __HBM__
    hbw_posix_memalign((void**) &this->_flatGrid, 4096, this->_lda*(this->_Ny + 2)*sizeof(GridT));
#else
#ifdef __INTEL_COMPILER
    this->_flatGrid = static_cast<GridT *>(_mm_malloc(this->_lda*(this->_Ny + 2)*sizeof(GridT), 64));
#else
    posix_memalign((void **) &this->_flatGrid, 4096, this->_lda*(this->_Ny + 2)*sizeof(GridT));
#endif
#endif
    for (size_t row = 0; row < this->_Ny + 2; ++row) {
      for (size_t col = 0; col < this->_Nx + 2; ++col) {
        (*this)(row, col) = storage(row, col);
      }
    }
  }

  template <typename GridT>
  RowMajorStorage<GridT>::RowMajorStorage(LinearContiguousStorage<GridT> const && storage) {
#ifdef __DEBUG__
    cout << "RowMajorStorage<GridT>::RowMajorStorage(LinearContiguousStorage<GridT> const &&)" << endl;
#endif
    this->_Ny = storage._Ny;
    storage._Ny = static_cast<size_t>(0);
    this->_Nx = storage._Nx;
    storage._Nx = static_cast<size_t>(0);
    this->_lda = storage._lda;
    storage._lda = static_cast<size_t>(0);
#ifdef __HBM__
    hbw_posix_memalign((void**) &this->_flatGrid, 4096, this->_lda*(this->_Ny + 2)*sizeof(GridT));
#else
#ifdef __INTEL_COMPILER
    this->_flatGrid = static_cast<GridT *>(_mm_malloc(this->_lda*(this->_Ny + 2)*sizeof(GridT), 64));
#else
    posix_memalign((void **) &this->_flatGrid, 4096, this->_lda*(this->_Ny + 2)*sizeof(GridT));
#endif
#endif
    for (size_t row = 0; row < this->_Ny + 2; ++row) {
      for (size_t col = 0; col < this->_Nx + 2; ++col) {
        (*this)(row, col) = storage(row, col);
      }
    }
    if (storage._flatGrid) {
#ifdef __HBM__
      hbw_free((void*)storage._flatGrid);
#else
#ifdef __INTEL_COMPILER
      _mm_free(storage._flatGrid);
#else
      free(storage._flatGrid);
#endif
#endif
      storage._flatGrid = nullptr;
    }
  }

  template <typename GridT>
  RowMajorStorage<GridT> & RowMajorStorage<GridT>::operator=(LinearContiguousStorage<GridT> const & storage) {
#ifdef __DEBUG__
    cout << "RowMajorStorage<GridT> & RowMajorStorage<GridT>::operator=(LinearContiguousStorage<GridT> const &)" << endl;
#endif
    if (this->_flatGrid) {
#ifdef __HBM__
      hbw_free((void*)this->_flatGrid);
#else
#ifdef __INTEL_COMPILER
      _mm_free(this->_flatGrid);
#else
      free(this->_flatGrid);
#endif
#endif
      this->_flatGrid = nullptr;
    }
    this->_Ny = storage._Ny;
    this->_Nx = storage._Nx;
    this->_lda = storage._lda;
#ifdef __HBM__
    hbw_posix_memalign((void**) &this->_flatGrid, 4096, this->_lda*(this->_Ny + 2)*sizeof(GridT));
#else
#ifdef __INTEL_COMPILER
    this->_flatGrid = static_cast<GridT *>(_mm_malloc(this->_lda*(this->_Ny + 2)*sizeof(GridT), 64));
#else
    posix_memalign((void **) &this->_flatGrid, 4096, this->_lda*(this->_Ny + 2)*sizeof(GridT));
#endif
#endif
    for (size_t row = 0; row < this->_Ny + 2; ++row) {
      for (size_t col = 0; col < this->_Nx + 2; ++col) {
        (*this)(row, col) = storage(row, col);
      }
    }
  }

  template <typename GridT>
  RowMajorStorage<GridT> & RowMajorStorage<GridT>::operator=(LinearContiguousStorage<GridT> const && storage) {
#ifdef __DEBUG__
    cout << "RowMajorStorage<GridT> & RowMajorStorage<GridT>::operator=(LinearContiguousStorage<GridT> const &&)" << endl;
#endif
    if (this->_flatGrid) {
#ifdef __HBM__
      hbw_free((void*)this->_flatGrid);
#else
#ifdef __INTEL_COMPILER
      _mm_free(this->_flatGrid);
#else
      free(this->_flatGrid);
#endif
#endif
      this->_flatGrid = nullptr;
    }
    this->_Ny = storage._Ny;
    storage._Ny = static_cast<size_t>(0);
    this->_Nx = storage._Nx;
    storage._Nx = static_cast<size_t>(0);
    this->_lda = storage._lda;
    storage._lda = static_cast<size_t>(0);
#ifdef __HBM__
    hbw_posix_memalign((void**) &this->_flatGrid, 4096, this->_lda*(this->_Ny + 2)*sizeof(GridT));
#else
#ifdef __INTEL_COMPILER
    this->_flatGrid = static_cast<GridT *>(_mm_malloc(this->_lda*(this->_Ny + 2)*sizeof(GridT), 64));
#else
    posix_memalign((void **) &this->_flatGrid, 4096, this->_lda*(this->_Ny + 2)*sizeof(GridT));
#endif
#endif
    for (size_t row = 0; row < this->_Ny + 2; ++row) {
      for (size_t col = 0; col < this->_Nx + 2; ++col) {
        (*this)(row, col) = storage(row, col);
      }
    }
    if (storage._flatGrid) {
#ifdef __HBM__
      hbw_free((void*)storage._flatGrid);
#else
#ifdef __INTEL_COMPILER
      _mm_free(storage._flatGrid);
#else
      free(storage._flatGrid);
#endif
#endif
      storage._flatGrid = nullptr;
    }
  }

  template <typename GridT>
  RowMajorStorage<GridT>::~RowMajorStorage() {
#ifdef __DEBUG__
    cout << "RowMajorStorage<GridT>::~RowMajorStorage()" << endl;
#endif
    if (this->_flatGrid) {
#ifdef __HBM__
      hbw_free((void*)this->_flatGrid);
#else
#ifdef __INTEL_COMPILER
      _mm_free(this->_flatGrid);
#else
      free(this->_flatGrid);
#endif
#endif
      this->_flatGrid = nullptr;
    }
  }

  template <typename GridT>
  void RowMajorStorage<GridT>::resize(size_t new_Ny, size_t new_Nx) {
#ifdef __DEBUG__
    cout << "void RowMajorStorage<GridT>::resize(size_t, size_t)" << endl;
#endif
    if (this->_flatGrid) {
#ifdef __HBM__
      hbw_free((void*)this->_flatGrid);
#else
#ifdef __INTEL_COMPILER
      _mm_free(this->_flatGrid);
#else
      free(this->_flatGrid);
#endif
#endif
      this->_flatGrid = nullptr;
    }
    this->_Ny = new_Ny;
    this->_Nx = new_Nx;
    this->_lda = paddedLength(this->_Nx, this->cacheLineLength);
#ifdef __HBM__
    hbw_posix_memalign((void**) &this->_flatGrid, 4096, this->_lda*(this->_Ny + 2)*sizeof(GridT));
#else
#ifdef __INTEL_COMPILER
    this->_flatGrid = static_cast<GridT *>(_mm_malloc(this->_lda*(this->_Ny + 2)*sizeof(GridT), 64));
#else
    posix_memalign((void **) &this->_flatGrid, 4096, this->_lda*(this->_Ny + 2)*sizeof(GridT));
#endif
#endif
    for (size_t i = 0; i < this->_lda*(this->_Ny + 2); ++i) {
      this->_flatGrid[i] = static_cast<GridT>(0.0);
    }
  }

  template <typename GridT>
  GridT & RowMajorStorage<GridT>::operator()(size_t y, size_t x) const {
#ifdef __DEBUG__
    cout << "GridT & RowMajorStorage<GridT>::operator()(size_t, size_t) const" << endl;
#endif
    cout << "Requesting storage(" << y << "," << x << ") = _flatGrid[" << y*this->_lda + x << "]" << " located at " << &this->_flatGrid[y*this->_lda + x] << endl;
    return this->_flatGrid[y*this->_lda + x];
  }

#endif // ROWMAJORSTORAGE_HPP

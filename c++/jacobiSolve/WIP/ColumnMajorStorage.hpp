#ifndef COLUMNMAJORSTORAGE_HPP
#define COLUMNMAJORSTORAGE_HPP

#include "LinearContiguousStorage.hpp"

  template <typename GridT>
  class ColumnMajorStorage : public LinearContiguousStorage<GridT> {
  public:
    ColumnMajorStorage();
    ColumnMajorStorage(size_t new_N);
    ColumnMajorStorage(size_t new_Nx, size_t new_Ny);
    ColumnMajorStorage(LinearContiguousStorage<GridT> const & storage);
    ColumnMajorStorage(LinearContiguousStorage<GridT> const && storage);
    ColumnMajorStorage<GridT> & operator=(LinearContiguousStorage<GridT> const & storage);
    ColumnMajorStorage<GridT> & operator=(LinearContiguousStorage<GridT> const && storage);
    virtual ~ColumnMajorStorage();
    virtual void resize(size_t new_Ny, size_t new_Nx);
    virtual GridT & operator()(size_t y, size_t x) const;
    friend ostream& operator<<(ostream& os, ColumnMajorStorage<GridT> const & store) {
#ifdef __DEBUG__
      cout << "friend ostream& operator<<(ostream& os, ColumnMajorStorage<GridT> const & store)" << endl;
#endif
      if ((store.Nx() > 0) and (store.Ny() > 0)) {
        for (size_t row = 0; row <= store.Ny(); ++row) {
          for (size_t col = 0; col <= store.Nx() + 1; ++col) {
            os << showpos << scientific << store.flatGrid()[col*store.lda() + row] << " ";
          }
          os << endl;
        }
        for (size_t col = 0; col <= store.Nx() + 1; ++col) {
          os << showpos << scientific << store.flatGrid()[col*store.lda() + store.Ny() + 1] << " ";
        }
      }
      return os;
    }
  };

  template <typename GridT>
  ColumnMajorStorage<GridT>::ColumnMajorStorage() {
#ifdef __DEBUG__
    cout << "ColumnMajorStorage<GridT>::ColumnMajorStorage()" << endl;
#endif
    this->_Ny = static_cast<size_t>(0);
    this->_Nx = static_cast<size_t>(0);
    this->_lda = static_cast<size_t>(0);
  }

  template <typename GridT>
  ColumnMajorStorage<GridT>::ColumnMajorStorage(size_t new_N) {
#ifdef __DEBUG__
    cout << "ColumnMajorStorage<GridT>::ColumnMajorStorage(size_t)" << endl;
#endif
    this->_Ny = new_N;
    this->_Nx = new_N;
    this->_lda = paddedLength(this->_Ny, this->cacheLineLength);
#ifdef __HBM__
    hbw_posix_memalign((void**) &this->_flatGrid, 4096, this->_lda*(this->_Nx + 2)*sizeof(GridT));
#else
#ifdef __INTEL_COMPILER
    this->_flatGrid = static_cast<GridT *>(_mm_malloc(this->_lda*(this->_Nx + 2)*sizeof(GridT), 64));
#else
    posix_memalign((void **) &this->_flatGrid, 4096, this->_lda*(this->_Nx + 2)*sizeof(GridT));
#endif
#endif
    for (size_t i = 0; i < this->_lda*(this->_Nx + 2); ++i) {
      this->_flatGrid[i] = static_cast<GridT>(0.0);
    }
  }

  template <typename GridT>
  ColumnMajorStorage<GridT>::ColumnMajorStorage(size_t new_Ny, size_t new_Nx) {
#ifdef __DEBUG__
    cout << "ColumnMajorStorage<GridT>::ColumnMajorStorage(size_t, size_t)" << endl;
#endif
    this->_Ny = new_Ny;
    this->_Nx = new_Nx;
    this->_lda = paddedLength(this->_Ny, this->cacheLineLength);
#ifdef __HBM__
    hbw_posix_memalign((void**) &this->_flatGrid, 4096, this->_lda*(this->_Nx + 2)*sizeof(GridT));
#else
#ifdef __INTEL_COMPILER
    this->_flatGrid = static_cast<GridT *>(_mm_malloc(this->_lda*(this->_Nx + 2)*sizeof(GridT), 64));
#else
    posix_memalign((void **) &this->_flatGrid, 4096, this->_lda*(this->_Nx + 2)*sizeof(GridT));
#endif
#endif
    for (size_t i = 0; i < this->_lda*(this->_Nx + 2); ++i) {
      this->_flatGrid[i] = static_cast<GridT>(0.0);
    }
  }

  template <typename GridT>
  ColumnMajorStorage<GridT>::ColumnMajorStorage(LinearContiguousStorage<GridT> const & storage) {
#ifdef __DEBUG__
    cout << "ColumnMajorStorage<GridT>::ColumnMajorStorage(LinearContiguousStorage<GridT> const &)" << endl;
#endif
    this->_Ny = storage._Ny;
    this->_Nx = storage._Nx;
    this->_lda = storage._lda;
#ifdef __HBM__
    hbw_posix_memalign((void**) &this->_flatGrid, 4096, this->_lda*(this->_Nx + 2)*sizeof(GridT));
#else
#ifdef __INTEL_COMPILER
    this->_flatGrid = static_cast<GridT *>(_mm_malloc(this->_lda*(this->_Nx + 2)*sizeof(GridT), 64));
#else
    posix_memalign((void **) &this->_flatGrid, 4096, this->_lda*(this->_Nx + 2)*sizeof(GridT));
#endif
#endif
    for (size_t row = 0; row < this->_Ny + 2; ++row) {
      for (size_t col = 0; col < this->_Nx + 2; ++col) {
        (*this)(row, col) = storage(row, col);
      }
    }
  }

  template <typename GridT>
  ColumnMajorStorage<GridT>::ColumnMajorStorage(LinearContiguousStorage<GridT> const && storage) {
#ifdef __DEBUG__
    cout << "ColumnMajorStorage<GridT>::ColumnMajorStorage(LinearContiguousStorage<GridT> const &&)" << endl;
#endif
    this->_Ny = storage._Ny;
    storage._Ny = static_cast<size_t>(0);
    this->_Nx = storage._Nx;
    storage._Nx = static_cast<size_t>(0);
    this->_lda = storage._lda;
    storage._Nx = static_cast<size_t>(0);
#ifdef __HBM__
    hbw_posix_memalign((void**) &this->_flatGrid, 4096, this->_lda*(this->_Nx + 2)*sizeof(GridT));
#else
#ifdef __INTEL_COMPILER
    this->_flatGrid = static_cast<GridT *>(_mm_malloc(this->_lda*(this->_Nx + 2)*sizeof(GridT), 64));
#else
    posix_memalign((void **) &this->_flatGrid, 4096, this->_lda*(this->_Nx + 2)*sizeof(GridT));
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
  ColumnMajorStorage<GridT> & ColumnMajorStorage<GridT>::operator=(LinearContiguousStorage<GridT> const & storage) {
#ifdef __DEBUG__
    cout << "ColumnMajorStorage<GridT> & ColumnMajorStorage<GridT>::operator=(LinearContiguousStorage<GridT> const &)" << endl;
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
    hbw_posix_memalign((void**) &this->_flatGrid, 4096, this->_lda*(this->_Nx + 2)*sizeof(GridT));
#else
#ifdef __INTEL_COMPILER
    this->_flatGrid = static_cast<GridT *>(_mm_malloc(this->_lda*(this->_Nx + 2)*sizeof(GridT), 64));
#else
    posix_memalign((void **) &this->_flatGrid, 4096, this->_lda*(this->_Nx + 2)*sizeof(GridT));
#endif
#endif
    for (size_t row = 0; row < this->_Ny + 2; ++row) {
      for (size_t col = 0; col < this->_Nx + 2; ++col) {
        (*this)(row, col) = storage(row, col);
      }
    }
  }

  template <typename GridT>
  ColumnMajorStorage<GridT> & ColumnMajorStorage<GridT>::operator=(LinearContiguousStorage<GridT> const && storage) {
#ifdef __DEBUG__
    cout << "ColumnMajorStorage<GridT> & ColumnMajorStorage<GridT>::operator=(LinearContiguousStorage<GridT> const &&)" << endl;
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
    hbw_posix_memalign((void**) &this->_flatGrid, 4096, this->_lda*(this->_Nx + 2)*sizeof(GridT));
#else
#ifdef __INTEL_COMPILER
    this->_flatGrid = static_cast<GridT *>(_mm_malloc(this->_lda*(this->_Nx + 2)*sizeof(GridT), 64));
#else
    posix_memalign((void **) &this->_flatGrid, 4096, this->_lda*(this->_Nx + 2)*sizeof(GridT));
#endif
#endif
    for (size_t row = 0; row < this->_Ny + 2; ++row) {
      for (size_t col = 0; col < this->_Nx + 2; ++col) {
        (*this)(row, col) = storage(row, col);
      }
    }
    if (storage._flatGrid) {
#ifdef __HBM__
      hbw_free((void*)storage.this->_flatGrid);
#else
#ifdef __INTEL_COMPILER
      _mm_free(storage.this->_flatGrid);
#else
      free(storage._flatGrid);
#endif
#endif
      storage._flatGrid = nullptr;
    }
  }

  template <typename GridT>
  ColumnMajorStorage<GridT>::~ColumnMajorStorage() {
#ifdef __DEBUG__
    cout << "ColumnMajorStorage<GridT>::~ColumnMajorStorage()" << endl;
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
  void ColumnMajorStorage<GridT>::resize(size_t new_Ny, size_t new_Nx) {
#ifdef __DEBUG__
    cout << "void ColumnMajorStorage<GridT>::resize(size_t, size_t)" << endl;
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
    this->_lda = paddedLength(this->_Ny, this->cacheLineLength);
#ifdef __HBM__
    hbw_posix_memalign((void**) &this->_flatGrid, 4096, this->_lda*(this->_Nx + 2)*sizeof(GridT));
#else
#ifdef __INTEL_COMPILER
    this->_flatGrid = static_cast<GridT *>(_mm_malloc(this->_lda*(this->_Nx + 2)*sizeof(GridT), 64));
#else
    posix_memalign((void **) &this->_flatGrid, 4096, this->_lda*(this->_Nx + 2)*sizeof(GridT));
#endif
#endif
    for (size_t i = 0; i < this->_lda*(this->_Nx + 2); ++i) {
      this->_flatGrid[i] = static_cast<GridT>(0.0);
    }
  }

  template <typename GridT>
  GridT & ColumnMajorStorage<GridT>::operator()(size_t y, size_t x) const {
#ifdef __DEBUG__
    cout << "GridT & ColumnMajorStorage<GridT>::operator()(size_t, size_t) const" << endl;
#endif
    return this->_flatGrid[x*this->_lda + y];
  }

#endif // COLUMNMAJORSTORAGE_HPP

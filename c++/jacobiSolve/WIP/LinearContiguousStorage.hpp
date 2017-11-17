#ifndef LINEARCONTIGUOUSSTORAGE_HPP
#define LINEARCONTIGUOUSSTORAGE_HPP

#include "Storage.hpp"

  template <typename GridT>
  class LinearContiguousStorage: public Storage<GridT> {
  protected:
    size_t _lda;
  public:
    LinearContiguousStorage();
    virtual ~LinearContiguousStorage();
    size_t lda();
    virtual void resize(size_t new_Ny, size_t new_Nx) = 0;
    virtual GridT & operator()(size_t y, size_t x) const = 0;
  };

  template <typename GridT>
  LinearContiguousStorage<GridT>::LinearContiguousStorage() {
#ifdef __DEBUG__
    cout << "LinearContiguousStorage<GridT>::LinearContiguousStorage()" << endl;
#endif
    this->_Ny = 0;
    this->_Nx = 0;
    this->_flatGrid = nullptr;
  }

  template <typename GridT>
  LinearContiguousStorage<GridT>::~LinearContiguousStorage() {
#ifdef __DEBUG__
    cout << "LinearContiguousStorage<GridT>::~LinearContiguousStorage()" << endl;
#endif
    this->_Ny = 0;
    this->_Nx = 0;
    this->_flatGrid = nullptr;
  }

  template <typename GridT>
  size_t LinearContiguousStorage<GridT>::lda() {
#ifdef __DEBUG__
    cout << "size_t LinearContiguousStorage<GridT>::lda()" << endl;
#endif
    return _lda;
  }

#endif // LINEARCONTIGUOUSSTORAGE_HPP

#include <iostream>
#include "RowMajorStorage.hpp"
#include "ColumnMajorStorage.hpp"
#include "Grid.hpp"

using namespace std;

int main() {
  cout << "Creating AGrid..." << endl;
  Grid<float, RowMajorStorage> AGrid(4);
  cout << "Setting AGrid(2,2) = 5.0..." << endl;
  AGrid(2,2) = 5.0;
  cout << "Setting AGrid(bdry) = 1.0..." << endl;
  AGrid.setBoundary(1.0);
  cout << "Printing AGrid(0,0)" << endl;
  cout << AGrid(0,0) << endl;
  cout << "Printing AGrid.Nx, AGrid(2,0) -> AGrid(2,3)" << endl;
  cout << AGrid.Nx() << ": " << AGrid(2,0) << " " << AGrid(2,1) << " " << AGrid(2,2) << " " << AGrid(2,3) << endl;
  cout << "Printing AGrid..." << endl;
  cout << AGrid << endl;
/*
  Grid<float, ColumnMajorStorage> BGrid(4);
  BGrid(2,2) = 25.0;
  BGrid.setBoundary(-1.0);
  cout << BGrid.Nx() << ": " << BGrid(2,0) << " " << BGrid(2,1) << " " << BGrid(2,2) << " " << BGrid(2,3) << endl;
  printf("BGrid stores data at %p\n", BGrid.flatGrid());

  Grid<float, RowMajorStorage> CGrid(4);
  CGrid(2,2) = 1.0;
  BGrid.setBoundary(-10.0);
  cout << CGrid.Nx() << ": " << CGrid(2,0) << " " << CGrid(2,1) << " " << CGrid(2,2) << " " << CGrid(2,3) << endl;
  printf("CGrid stores data at %p\n", CGrid.flatGrid());

  Grid<float, RowMajorStorage> DGrid(CGrid);
  printf("DGrid stores data at %p\n", DGrid.flatGrid());

  Grid<float, ColumnMajorStorage> EGrid(BGrid);
  printf("EGrid stores data at %p\n", EGrid.flatGrid());

  cout << "Worked! Hello template-template!" << endl;
*/
}

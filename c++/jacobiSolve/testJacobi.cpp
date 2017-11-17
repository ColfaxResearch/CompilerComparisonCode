#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <cmath>
#endif
#include <stdlib.h>
#include <malloc.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <cassert>
#include <unistd.h>
#include <limits>
#ifdef __HBM__
#include <hbwmalloc.h>
#endif
#if defined(__INTEL_COMPILER) && defined(__ADVIXE__)
#include "advisor-annotate.h"
#endif

#include "jacobi.hpp"


using namespace std;

int main() {
    size_t numPts = 256;
    //cout << "Number of grid points: ";
    //cin >> numPts;
    double const epsilon = 8.8541878176e-12;
    double const pi = 3.14159265359;
    double chargeDensity = 1.0e-6/epsilon;
    cout << "Making Jacobi J" << endl;
    Jacobi<double> J(numPts);
    cout << "Setting J.Sy = -500.0 and J.Sx = -500.0" << endl;
    J.set_Sy(-500.0);
    J.set_Sx(-500.0);
    cout << "Setting J.Fy = 500.0 and J.Fx = 500.0" << endl;
    J.set_Fy(500.0);
    J.set_Fx(500.0);
    cout << "Setting J.S to ball of radius 20 with rho = 2.0 muC/m^3" << endl;
    cout << "Setting J.S to ball of radius 30 with rho = -1.0 muC/m^3" << endl;
    cout << "Setting J.S to ball of radius 30 with rho = -4.0 muC/m^3" << endl;
    cout << "Setting J.S to ball of radius 30 with rho = 3.0 muC/m^3" << endl;
    double y, x;
    double yBall_1 = -85.0, xBall_1 = 0.0, rhoBall_1 = 2.0, rBall_1 = 20.0, qBall_1 = rhoBall_1*(4.0/3.0)*pi*std::pow(rBall_1, 3.0);
    double yBall_2 = +65.0, xBall_2 = 0.0, rhoBall_2 = -3.0, rBall_2 = 20.0, qBall_2 = rhoBall_2*(4.0/3.0)*pi*std::pow(rBall_2, 3.0);
    double yBall_3 = 0.0, xBall_3 = -60.0, rhoBall_3 = -1.0, rBall_3 = 20.0, qBall_3 = rhoBall_3*(4.0/3.0)*pi*std::pow(rBall_3, 3.0);
    double yBall_4 = 0.0, xBall_4 = +60.0, rhoBall_4 = 2.0, qBall_4 = -1.0*(qBall_1 + qBall_2 + qBall_3), rBall_4 = std::pow((3.0*std::fabs(qBall_4))/(4.0*pi), 1.0/3.0);

    for (size_t i = 0; i < J.get_Ny(); ++i) {
        for (size_t j = 0; j < J.get_Nx(); ++j) {
            y = J.get_y(i, j);
            x = J.get_x(i, j);
            if ((pow(x - xBall_1, 2.0) + pow(y - yBall_1, 2.0)) < pow(rBall_1, 2)) {
                J.getSource()(i, j) += rhoBall_1*chargeDensity;
            }
            if ((pow(x - xBall_2, 2.0) + pow(y - yBall_2, 2.0)) < pow(rBall_2, 2)) {
                J.getSource()(i, j) += rhoBall_2*chargeDensity;
            }
            if ((pow(x - xBall_3, 2.0) + pow(y - yBall_3, 2.0)) < pow(rBall_3, 2)) {
                J.getSource()(i, j) += rhoBall_3*chargeDensity;
            }
            if ((pow(x - xBall_4, 2.0) + pow(y - yBall_4, 2.0)) < pow(rBall_4, 2)) {
                J.getSource()(i, j) += rhoBall_4*chargeDensity;
            }
        }
    }
    cout << "Setting J(bdry) = 0.0" << endl;
    J.setBoundary(0.0);
    cout << "Solving J" << endl;
    size_t numIter = J.solve_opt(1.0e3);
    cout << "Problem solved with " << numIter << " iterations" << endl;
    ofstream outFile;
    outFile.open("/home/AstroVPK/code/intelShootout/data/jacobiSolve.dat", ios::in | ios::trunc);
    if (!outFile) { // if outFile not attached to a valid input source, abort
        cout << "Sorry, bad file.";
        exit(0);	// special system call to abort program may require <cstdlib> to be included
    }
    outFile << J.get_Sy() << " " << J.get_Fy() << " " << J.get_Ny() + 2 << endl;
    outFile << J.get_Sx() << " " << J.get_Fx() << " " << J.get_Nx() + 2 << endl;
    outFile << J.getDomain() << endl;
    outFile.close();
    return 0;
}
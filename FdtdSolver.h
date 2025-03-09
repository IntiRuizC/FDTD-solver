#ifndef FDTD_SOLVER_FDTDSOLVER_H
#define FDTD_SOLVER_FDTDSOLVER_H

#include "grid.h"

#include <cstddef>
#include <vector>

class FdtdSolver {
public:
    FdtdSolver(const Grid& grid, double dt);

    // Yee updates for a single cell. Callers must keep neighbours in range:
    // H uses (i+1, j+1, k+1), E uses (i-1, j-1, k-1).
    void calculate_H_field(std::size_t i, std::size_t j, std::size_t k);
    void calculate_E_field(std::size_t i, std::size_t j, std::size_t k);

    Grid grid;
    double dt;

    std::vector<double> H_x; // H_x[i][j][k] -> H_x(i, j+1/2, k+1/2)
    std::vector<double> H_y; // H_y[i][j][k] -> H_y(i+1/2, j, k+1/2)
    std::vector<double> H_z; // H_z[i][j][k] -> H_z(i+1/2, j+1/2, k)

    std::vector<double> E_x; // E_x[i][j][k] -> E_x(i+1/2, j, k)
    std::vector<double> E_y; // E_y[i][j][k] -> E_y(i, j+1/2, k)
    std::vector<double> E_z; // E_z[i][j][k] -> E_z(i, j, k+1/2)

    //TODO: replace with material map
    std::vector<double> mu;
    std::vector<double> epsilon;
    std::vector<double> sigma;
};


#endif //FDTD_SOLVER_FDTDSOLVER_H

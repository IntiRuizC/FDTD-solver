#ifndef FDTD_SOLVER_GRID_H
#define FDTD_SOLVER_GRID_H

#include <array>
#include <cstddef>

class Grid {
public:
    Grid(std::size_t nx, std::size_t ny, std::size_t nz, double dx, double dy, double dz);

    std::size_t to_flat(std::size_t i, std::size_t j, std::size_t k) const;
    std::array<std::size_t, 3> to_index(std::size_t index) const;

    std::size_t nx;
    std::size_t ny;
    std::size_t nz;
    double dx;
    double dy;
    double dz;
};

#endif //FDTD_SOLVER_GRID_H

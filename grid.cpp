#include "grid.h"

Grid::Grid(
    const std::size_t nx,
    const std::size_t ny,
    const std::size_t nz,
    const double dx,
    const double dy,
    const double dz
    )
    :nx(nx), ny(ny), nz(nz), dx(dx), dy(dy), dz(dz) {}

std::size_t Grid::to_flat(const std::size_t i, const std::size_t j,
                            const std::size_t k) const {
    return (i * ny + j) * nz + k;
}

std::array<std::size_t, 3> Grid::to_index(std::size_t index) const {
    const std::size_t i = index / (ny * nz);
    index -= i * ny * nz;
    const std::size_t j = index / nz;
    const std::size_t k = index % nz;
    return {i, j, k};
}

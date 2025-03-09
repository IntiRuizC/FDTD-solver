#include "FdtdSolver.h"

namespace {
// Normalised units (c = 1)
constexpr double MU_0 = 1.0;
constexpr double EPSILON_0 = 1.0;
constexpr double SIGMA = 0.0;
}

FdtdSolver::FdtdSolver(const Grid& grid, const double dt)
    : grid(grid), dt(dt),
      H_x(grid.nx * grid.ny * grid.nz, 0.0),
      H_y(grid.nx * grid.ny * grid.nz, 0.0),
      H_z(grid.nx * grid.ny * grid.nz, 0.0),
      E_x(grid.nx * grid.ny * grid.nz, 0.0),
      E_y(grid.nx * grid.ny * grid.nz, 0.0),
      E_z(grid.nx * grid.ny * grid.nz, 0.0),
      mu(grid.nx * grid.ny * grid.nz, MU_0),
      epsilon(grid.nx * grid.ny * grid.nz, EPSILON_0),
      sigma(grid.nx * grid.ny * grid.nz, SIGMA) {}

void FdtdSolver::calculate_H_field(const std::size_t i, const std::size_t j,
                                   const std::size_t k) {
    const std::size_t c = grid.to_flat(i, j, k);
    const std::size_t ip = grid.to_flat(i + 1, j, k);
    const std::size_t jp = grid.to_flat(i, j + 1, k);
    const std::size_t kp = grid.to_flat(i, j, k + 1);

    const double coeff = dt / mu[c];

    H_x[c] += coeff * ((E_y[kp] - E_y[c]) / grid.dz + (E_z[c] - E_z[jp]) / grid.dy);
    H_y[c] += coeff * ((E_z[ip] - E_z[c]) / grid.dx + (E_x[c] - E_x[kp]) / grid.dz);
    H_z[c] += coeff * ((E_x[jp] - E_x[c]) / grid.dy + (E_y[c] - E_y[ip]) / grid.dx);
}

void FdtdSolver::calculate_E_field(const std::size_t i, const std::size_t j,
                                   const std::size_t k) {
    const std::size_t c = grid.to_flat(i, j, k);
    const std::size_t im = grid.to_flat(i - 1, j, k);
    const std::size_t jm = grid.to_flat(i, j - 1, k);
    const std::size_t km = grid.to_flat(i, j, k - 1);

    const double a = 1.0 - sigma[c] * dt / epsilon[c];
    const double b = dt / epsilon[c];

    E_x[c] = a * E_x[c] + b * ((H_z[c] - H_z[jm]) / grid.dy + (H_y[km] - H_y[c]) / grid.dz);
    E_y[c] = a * E_y[c] + b * ((H_x[c] - H_x[km]) / grid.dz + (H_z[im] - H_z[c]) / grid.dx);
    E_z[c] = a * E_z[c] + b * ((H_y[c] - H_y[im]) / grid.dx + (H_x[jm] - H_x[c]) / grid.dy);
}

void FdtdSolver::time_step(void) {
    for (int i = 0; i < grid.nx-1; ++i) {
        for (int j = 0; j < grid.ny-1; ++j) {
            for (int k = 0; k < grid.nz-1; ++k) {
                calculate_H_field(i,j,k);
            }
        }
    }
    for (int i = 1; i < grid.nx; ++i) {
        for (int j = 1; j < grid.ny; ++j) {
            for (int k = 1; k < grid.nz; ++k) {
                calculate_E_field(i,j,k);
            }
        }
    }
}

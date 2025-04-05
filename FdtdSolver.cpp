#include "FdtdSolver.h"

namespace {
// Normalised units (c = 1)
constexpr double MU_0 = 1.0;
constexpr double EPSILON_0 = 1.0;
constexpr double SIGMA = 0.0;

// Taflove's lattice truncation (Sadiku eq. 3.89) for one tangential E value:
// the value one cell inside, two steps back, averaged with its neighbours
// `stride` apart along the face where those exist.
double taflove(const std::vector<double>& past, const std::size_t interior,
               const std::size_t stride, const bool has_lower, const bool has_upper) {
    double sum = past[interior];
    double count = 1.0;
    if (has_lower) {
        sum += past[interior - stride];
        count += 1.0;
    }
    if (has_upper) {
        sum += past[interior + stride];
        count += 1.0;
    }
    return sum / count;
}
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
      sigma(grid.nx * grid.ny * grid.nz, SIGMA),
      E_x_prev(grid.nx * grid.ny * grid.nz, 0.0),
      E_y_prev(grid.nx * grid.ny * grid.nz, 0.0),
      E_z_prev(grid.nx * grid.ny * grid.nz, 0.0),
      E_x_prev2(grid.nx * grid.ny * grid.nz, 0.0),
      E_y_prev2(grid.nx * grid.ny * grid.nz, 0.0),
      E_z_prev2(grid.nx * grid.ny * grid.nz, 0.0) {}

void FdtdSolver::calculate_H_field(const std::size_t i, const std::size_t j,
                                   const std::size_t k) {
    update_H_x(i, j, k);
    update_H_y(i, j, k);
    update_H_z(i, j, k);
}

void FdtdSolver::calculate_E_field(const std::size_t i, const std::size_t j,
                                   const std::size_t k) {
    update_E_x(i, j, k);
    update_E_y(i, j, k);
    update_E_z(i, j, k);
}

void FdtdSolver::update_H_x(const std::size_t i, const std::size_t j, const std::size_t k) {
    const std::size_t c = grid.to_flat(i, j, k);
    const std::size_t jp = grid.to_flat(i, j + 1, k);
    const std::size_t kp = grid.to_flat(i, j, k + 1);
    H_x[c] += dt / mu[c] * ((E_y[kp] - E_y[c]) / grid.dz + (E_z[c] - E_z[jp]) / grid.dy);
}

void FdtdSolver::update_H_y(const std::size_t i, const std::size_t j, const std::size_t k) {
    const std::size_t c = grid.to_flat(i, j, k);
    const std::size_t ip = grid.to_flat(i + 1, j, k);
    const std::size_t kp = grid.to_flat(i, j, k + 1);
    H_y[c] += dt / mu[c] * ((E_z[ip] - E_z[c]) / grid.dx + (E_x[c] - E_x[kp]) / grid.dz);
}

void FdtdSolver::update_H_z(const std::size_t i, const std::size_t j, const std::size_t k) {
    const std::size_t c = grid.to_flat(i, j, k);
    const std::size_t ip = grid.to_flat(i + 1, j, k);
    const std::size_t jp = grid.to_flat(i, j + 1, k);
    H_z[c] += dt / mu[c] * ((E_x[jp] - E_x[c]) / grid.dy + (E_y[c] - E_y[ip]) / grid.dx);
}

void FdtdSolver::update_E_x(const std::size_t i, const std::size_t j, const std::size_t k) {
    const std::size_t c = grid.to_flat(i, j, k);
    const std::size_t jm = grid.to_flat(i, j - 1, k);
    const std::size_t km = grid.to_flat(i, j, k - 1);
    const double a = 1.0 - sigma[c] * dt / epsilon[c];
    const double b = dt / epsilon[c];
    E_x[c] = a * E_x[c] + b * ((H_z[c] - H_z[jm]) / grid.dy + (H_y[km] - H_y[c]) / grid.dz);
}

void FdtdSolver::update_E_y(const std::size_t i, const std::size_t j, const std::size_t k) {
    const std::size_t c = grid.to_flat(i, j, k);
    const std::size_t im = grid.to_flat(i - 1, j, k);
    const std::size_t km = grid.to_flat(i, j, k - 1);
    const double a = 1.0 - sigma[c] * dt / epsilon[c];
    const double b = dt / epsilon[c];
    E_y[c] = a * E_y[c] + b * ((H_x[c] - H_x[km]) / grid.dz + (H_z[im] - H_z[c]) / grid.dx);
}

void FdtdSolver::update_E_z(const std::size_t i, const std::size_t j, const std::size_t k) {
    const std::size_t c = grid.to_flat(i, j, k);
    const std::size_t im = grid.to_flat(i - 1, j, k);
    const std::size_t jm = grid.to_flat(i, j - 1, k);
    const double a = 1.0 - sigma[c] * dt / epsilon[c];
    const double b = dt / epsilon[c];
    E_z[c] = a * E_z[c] + b * ((H_y[c] - H_y[im]) / grid.dx + (H_x[jm] - H_x[c]) / grid.dy);
}

void FdtdSolver::time_step() {
    const std::size_t nx = grid.nx;
    const std::size_t ny = grid.ny;
    const std::size_t nz = grid.nz;

    // Each component is updated only where its neighbours exist. H on the
    // outer faces is skipped because no interior E update reads it.
    for (std::size_t i = 1; i < nx - 1; ++i)
        for (std::size_t j = 0; j < ny - 1; ++j)
            for (std::size_t k = 0; k < nz - 1; ++k)
                update_H_x(i, j, k);
    for (std::size_t i = 0; i < nx - 1; ++i)
        for (std::size_t j = 1; j < ny - 1; ++j)
            for (std::size_t k = 0; k < nz - 1; ++k)
                update_H_y(i, j, k);
    for (std::size_t i = 0; i < nx - 1; ++i)
        for (std::size_t j = 0; j < ny - 1; ++j)
            for (std::size_t k = 1; k < nz - 1; ++k)
                update_H_z(i, j, k);

    E_x_prev2 = E_x_prev;
    E_y_prev2 = E_y_prev;
    E_z_prev2 = E_z_prev;
    E_x_prev = E_x;
    E_y_prev = E_y;
    E_z_prev = E_z;

    // Tangential E on the outer faces is left to the absorbing boundary.
    for (std::size_t i = 0; i < nx - 1; ++i)
        for (std::size_t j = 1; j < ny - 1; ++j)
            for (std::size_t k = 1; k < nz - 1; ++k)
                update_E_x(i, j, k);
    for (std::size_t i = 1; i < nx - 1; ++i)
        for (std::size_t j = 0; j < ny - 1; ++j)
            for (std::size_t k = 1; k < nz - 1; ++k)
                update_E_y(i, j, k);
    for (std::size_t i = 1; i < nx - 1; ++i)
        for (std::size_t j = 1; j < ny - 1; ++j)
            for (std::size_t k = 0; k < nz - 1; ++k)
                update_E_z(i, j, k);

    apply_taflove_boundary();
}

// Sets the tangential E on all six faces from the first interior plane two
// steps back (a wave crosses one cell in two steps when dt = dx / 2c). As in
// Sadiku eq. 3.89, x faces average along z, z faces along x and y faces not at
// all. Components on the box edges are tangential to two faces; they are left
// at zero.
void FdtdSolver::apply_taflove_boundary() {
    const std::size_t nx = grid.nx;
    const std::size_t ny = grid.ny;
    const std::size_t nz = grid.nz;
    const std::size_t i_stride = grid.to_flat(1, 0, 0);
    const std::size_t k_stride = grid.to_flat(0, 0, 1);

    // Faces i = 0 and i = nx-1, averaged along z
    for (std::size_t j = 0; j < ny - 1; ++j)
        for (std::size_t k = 1; k < nz - 1; ++k) {
            const bool lower = k > 1;
            const bool upper = k + 1 < nz - 1;
            E_y[grid.to_flat(0, j, k)] =
                taflove(E_y_prev2, grid.to_flat(1, j, k), k_stride, lower, upper);
            E_y[grid.to_flat(nx - 1, j, k)] =
                taflove(E_y_prev2, grid.to_flat(nx - 2, j, k), k_stride, lower, upper);
        }
    for (std::size_t j = 1; j < ny - 1; ++j)
        for (std::size_t k = 0; k < nz - 1; ++k) {
            const bool lower = k > 0;
            const bool upper = k + 1 < nz - 1;
            E_z[grid.to_flat(0, j, k)] =
                taflove(E_z_prev2, grid.to_flat(1, j, k), k_stride, lower, upper);
            E_z[grid.to_flat(nx - 1, j, k)] =
                taflove(E_z_prev2, grid.to_flat(nx - 2, j, k), k_stride, lower, upper);
        }

    // Faces j = 0 and j = ny-1, no averaging
    for (std::size_t i = 0; i < nx - 1; ++i)
        for (std::size_t k = 1; k < nz - 1; ++k) {
            E_x[grid.to_flat(i, 0, k)] = E_x_prev2[grid.to_flat(i, 1, k)];
            E_x[grid.to_flat(i, ny - 1, k)] = E_x_prev2[grid.to_flat(i, ny - 2, k)];
        }
    for (std::size_t i = 1; i < nx - 1; ++i)
        for (std::size_t k = 0; k < nz - 1; ++k) {
            E_z[grid.to_flat(i, 0, k)] = E_z_prev2[grid.to_flat(i, 1, k)];
            E_z[grid.to_flat(i, ny - 1, k)] = E_z_prev2[grid.to_flat(i, ny - 2, k)];
        }

    // Faces k = 0 and k = nz-1, averaged along x
    for (std::size_t i = 0; i < nx - 1; ++i)
        for (std::size_t j = 1; j < ny - 1; ++j) {
            const bool lower = i > 0;
            const bool upper = i + 1 < nx - 1;
            E_x[grid.to_flat(i, j, 0)] =
                taflove(E_x_prev2, grid.to_flat(i, j, 1), i_stride, lower, upper);
            E_x[grid.to_flat(i, j, nz - 1)] =
                taflove(E_x_prev2, grid.to_flat(i, j, nz - 2), i_stride, lower, upper);
        }
    for (std::size_t i = 1; i < nx - 1; ++i)
        for (std::size_t j = 0; j < ny - 1; ++j) {
            const bool lower = i > 1;
            const bool upper = i + 1 < nx - 1;
            E_y[grid.to_flat(i, j, 0)] =
                taflove(E_y_prev2, grid.to_flat(i, j, 1), i_stride, lower, upper);
            E_y[grid.to_flat(i, j, nz - 1)] =
                taflove(E_y_prev2, grid.to_flat(i, j, nz - 2), i_stride, lower, upper);
        }
}

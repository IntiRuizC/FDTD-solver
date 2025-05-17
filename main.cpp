#include "FdtdSolver.h"
#include "grid.h"
#include "output.h"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <format>

namespace {
// Normalised units (c = 1), cubic cells
constexpr std::size_t N = 61;
constexpr double DX = 1.0;
constexpr double DT = DX / 2.0; // required by the Taflove boundary
constexpr int STEPS = 400;
constexpr int SNAPSHOT_EVERY = 10;

// Source current J_z(t) = -(t - t0)/tau * exp(-((t - t0)/tau)^2). The
// derivative of a Gaussian has zero time integral, so no net charge is left
// at the source once the pulse has passed.
constexpr double T0 = 30.0;
constexpr double TAU = 8.0;

double source_current(const double t) {
    const double arg = (t - T0) / TAU;
    return -arg * std::exp(-arg * arg);
}
}

int main(int argc, char** argv) {
    const std::filesystem::path out_dir = argc > 1 ? argv[1] : "out";
    std::filesystem::create_directories(out_dir);

    const Grid grid(N, N, N, DX, DX, DX);
    FdtdSolver solver(grid, DT);

    const std::size_t c = N / 2;
    const std::size_t source = grid.to_flat(c, c, c);

    ProbeRecorder probes;
    probes.add("Ez_source", Component::E_z, c, c, c);
    probes.add("Ez_x10", Component::E_z, c + 10, c, c);
    probes.add("Ez_x25", Component::E_z, c + 25, c, c);
    probes.add("Hy_x10", Component::H_y, c + 10, c, c);

    std::FILE* energy = std::fopen((out_dir / "energy.csv").c_str(), "w");
    std::fprintf(energy, "t,U\n");

    for (int n = 0; n < STEPS; ++n) {
        solver.time_step();
        const double t = (n + 1) * DT;
        // Soft source: J enters Ampere's law as dE/dt = (curl H - J) / eps.
        solver.E_z[source] -= DT / solver.epsilon[source] * source_current(t);

        probes.record(solver, t);
        std::fprintf(energy, "%.17g,%.17g\n", t, field_energy(solver));
        if ((n + 1) % SNAPSHOT_EVERY == 0)
            write_npy(out_dir / std::format("Ez_z{}_{:05d}.npy", c, n + 1),
                      slice_z(solver, Component::E_z, c), {N, N});
    }

    std::fclose(energy);
    probes.write_csv(out_dir / "probes.csv");
    std::printf("wrote %d steps to %s\n", STEPS, out_dir.c_str());
}

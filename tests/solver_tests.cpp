#include "FdtdSolver.h"
#include "grid.h"
#include "output.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <iterator>
#include <string>

namespace {
constexpr Component ALL_COMPONENTS[] = {Component::E_x, Component::E_y, Component::E_z,
                                        Component::H_x, Component::H_y, Component::H_z};

FdtdSolver make_solver(const std::size_t n) {
    return FdtdSolver(Grid(n, n, n, 1.0, 1.0, 1.0), 0.5);
}

double gaussian_derivative(const double t, const double t0, const double tau) {
    const double arg = (t - t0) / tau;
    return -arg * std::exp(-arg * arg);
}

double max_abs(const std::vector<double>& v) {
    double m = 0.0;
    for (const double x : v) m = std::max(m, std::abs(x));
    return m;
}
}

// The Yee stencil couples only nearest neighbours, so a disturbance moves at
// most one cell per step along each axis. Along the axis normal to E_z the
// front is sharp: E_z is non-zero exactly n cells out after n steps.
TEST(FdtdSolver, ImpulseStaysInsideNumericalLightCone) {
    constexpr std::size_t N = 21;
    constexpr std::size_t C = N / 2;
    constexpr std::size_t STEPS = 7;
    FdtdSolver solver = make_solver(N);
    const Grid& g = solver.grid;
    solver.E_z[g.to_flat(C, C, C)] = 1.0;

    for (std::size_t n = 0; n < STEPS; ++n) solver.time_step();

    EXPECT_NE(solver.E_z[g.to_flat(C + STEPS, C, C)], 0.0);
    EXPECT_EQ(solver.E_z[g.to_flat(C + STEPS + 1, C, C)], 0.0);

    const auto dist = [](const std::size_t a, const std::size_t b) { return a > b ? a - b : b - a; };
    for (const Component comp : ALL_COMPONENTS) {
        const std::vector<double>& values = field(solver, comp);
        for (std::size_t q = 0; q < values.size(); ++q) {
            const auto [i, j, k] = g.to_index(q);
            if (std::max({dist(i, C), dist(j, C), dist(k, C)}) > STEPS + 1) {
                ASSERT_EQ(values[q], 0.0) << "at (" << i << ", " << j << ", " << k << ")";
            }
        }
    }
}

// An E_z source at the centre is symmetric under x -> -x and y -> -y, and the
// boundary treats opposite faces alike, so E_z must stay mirror-symmetric even
// after the pulse has reached the walls.
TEST(FdtdSolver, CentredSourceKeepsMirrorSymmetry) {
    constexpr std::size_t N = 21;
    constexpr std::size_t C = N / 2;
    FdtdSolver solver = make_solver(N);
    const Grid& g = solver.grid;
    const std::size_t source = g.to_flat(C, C, C);

    for (int n = 0; n < 60; ++n) {
        solver.time_step();
        solver.E_z[source] -= solver.dt * gaussian_derivative((n + 1) * solver.dt, 8.0, 2.0);
    }

    const double tol = 1e-12 * max_abs(solver.E_z);
    ASSERT_GT(tol, 0.0);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
            for (std::size_t k = 0; k < N; ++k) {
                const double v = solver.E_z[g.to_flat(i, j, k)];
                ASSERT_NEAR(v, solver.E_z[g.to_flat(N - 1 - i, j, k)], tol);
                ASSERT_NEAR(v, solver.E_z[g.to_flat(i, N - 1 - j, k)], tol);
            }
}

// A closed (PEC) box would keep all the energy of the pulse. The absorbing
// boundary must let nearly all of it out (about 4% remains here). The pulse is
// short enough that the source is off before the front reaches the walls, so
// the energy at that point is the radiated energy, free of the source's
// reactive near field. Runs well short of the boundary's late-time
// instability.
TEST(FdtdSolver, AbsorbingBoundaryLetsPulseOut) {
    constexpr std::size_t N = 31;
    constexpr std::size_t C = N / 2;
    FdtdSolver solver = make_solver(N);
    const std::size_t source = solver.grid.to_flat(C, C, C);

    double radiated = 0.0;
    for (int n = 0; n < 200; ++n) {
        solver.time_step();
        solver.E_z[source] -= solver.dt * gaussian_derivative((n + 1) * solver.dt, 6.0, 1.5);
        if (n + 1 == 26) radiated = field_energy(solver); // t = 13, front 2 cells from the walls
    }

    ASSERT_GT(radiated, 0.0);
    EXPECT_LT(field_energy(solver), 0.1 * radiated);
}

TEST(Output, NpyHasAlignedHeaderAndRawData) {
    const std::filesystem::path path = std::filesystem::temp_directory_path() / "fdtd_test.npy";
    const std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    write_npy(path, data, {2, 3});

    std::ifstream in(path, std::ios::binary);
    const std::string bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::filesystem::remove(path);

    ASSERT_GE(bytes.size(), 10u);
    EXPECT_EQ(bytes.substr(0, 6), "\x93NUMPY");
    const std::size_t header_len =
        static_cast<unsigned char>(bytes[8]) | static_cast<unsigned char>(bytes[9]) << 8;
    EXPECT_EQ((10 + header_len) % 64, 0u);
    EXPECT_NE(bytes.find("'shape': (2, 3)"), std::string::npos);
    EXPECT_EQ(bytes[10 + header_len - 1], '\n');
    ASSERT_EQ(bytes.size(), 10 + header_len + data.size() * sizeof(double));
    EXPECT_EQ(bytes.compare(10 + header_len, std::string::npos,
                            reinterpret_cast<const char*>(data.data()), data.size() * sizeof(double)),
              0);
}

TEST(Output, ProbeRecordsSiteValuePerStep) {
    FdtdSolver solver = make_solver(5);
    ProbeRecorder probes;
    probes.add("ez", Component::E_z, 2, 3, 1);
    solver.E_z[solver.grid.to_flat(2, 3, 1)] = 0.25;
    probes.record(solver, 0.0);
    solver.E_z[solver.grid.to_flat(2, 3, 1)] = -1.5;
    probes.record(solver, 0.5);

    const std::filesystem::path path = std::filesystem::temp_directory_path() / "fdtd_probes.csv";
    probes.write_csv(path);
    std::ifstream in(path);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::filesystem::remove(path);

    EXPECT_EQ(text, "t,ez\n0,0.25\n0.5,-1.5\n");
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

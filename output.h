#ifndef FDTD_SOLVER_OUTPUT_H
#define FDTD_SOLVER_OUTPUT_H

#include "FdtdSolver.h"

#include <array>
#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

enum class Component { E_x, E_y, E_z, H_x, H_y, H_z };

const std::vector<double>& field(const FdtdSolver& solver, Component component);

// Total electromagnetic energy, 1/2 sum(eps E^2 + mu H^2) dV, with each
// component taken at its own Yee site.
double field_energy(const FdtdSolver& solver);

// Values of one component on the plane k = const, shape (nx, ny).
std::vector<double> slice_z(const FdtdSolver& solver, Component component, std::size_t k);

// Writes a C-ordered float64 array in NumPy's .npy format (version 1.0).
void write_npy(const std::filesystem::path& path, const std::vector<double>& data,
               const std::vector<std::size_t>& shape);

// Records field values at fixed lattice sites, one column per probe, and
// writes them as CSV with the time in the first column.
class ProbeRecorder {
public:
    void add(std::string name, Component component, std::size_t i, std::size_t j,
             std::size_t k);
    void record(const FdtdSolver& solver, double t);
    void write_csv(const std::filesystem::path& path) const;

private:
    struct Probe {
        std::string name;
        Component component;
        std::array<std::size_t, 3> site;
    };

    std::vector<Probe> probes;
    std::vector<double> times;
    std::vector<std::vector<double>> rows;
};

#endif //FDTD_SOLVER_OUTPUT_H

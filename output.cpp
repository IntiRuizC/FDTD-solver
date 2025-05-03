#include "output.h"

#include <cstdint>
#include <fstream>
#include <stdexcept>

const std::vector<double>& field(const FdtdSolver& solver, const Component component) {
    switch (component) {
        case Component::E_x: return solver.E_x;
        case Component::E_y: return solver.E_y;
        case Component::E_z: return solver.E_z;
        case Component::H_x: return solver.H_x;
        case Component::H_y: return solver.H_y;
        case Component::H_z: return solver.H_z;
    }
    throw std::invalid_argument("unknown field component");
}

double field_energy(const FdtdSolver& solver) {
    const Grid& g = solver.grid;
    double sum = 0.0;
    for (std::size_t c = 0; c < solver.E_x.size(); ++c) {
        sum += solver.epsilon[c] * (solver.E_x[c] * solver.E_x[c] +
                                    solver.E_y[c] * solver.E_y[c] +
                                    solver.E_z[c] * solver.E_z[c]);
        sum += solver.mu[c] * (solver.H_x[c] * solver.H_x[c] +
                               solver.H_y[c] * solver.H_y[c] +
                               solver.H_z[c] * solver.H_z[c]);
    }
    return 0.5 * sum * g.dx * g.dy * g.dz;
}

std::vector<double> slice_z(const FdtdSolver& solver, const Component component,
                            const std::size_t k) {
    const Grid& g = solver.grid;
    if (k >= g.nz) throw std::out_of_range("slice_z: k outside the lattice");
    const std::vector<double>& values = field(solver, component);
    std::vector<double> slice(g.nx * g.ny);
    for (std::size_t i = 0; i < g.nx; ++i)
        for (std::size_t j = 0; j < g.ny; ++j)
            slice[i * g.ny + j] = values[g.to_flat(i, j, k)];
    return slice;
}

void write_npy(const std::filesystem::path& path, const std::vector<double>& data,
               const std::vector<std::size_t>& shape) {
    std::size_t count = 1;
    for (const std::size_t n : shape) count *= n;
    if (count != data.size()) throw std::invalid_argument("write_npy: shape does not match data");

    // Python tuple syntax: "(3, 4)", and "(3,)" for one dimension
    std::string shape_str = "(";
    for (std::size_t d = 0; d < shape.size(); ++d)
        shape_str += (d > 0 ? ", " : "") + std::to_string(shape[d]);
    shape_str += shape.size() == 1 ? ",)" : ")";

    // Magic (6) + version (2) + header length (2) + header must be a multiple
    // of 64 bytes, with the header ending in '\n'.
    std::string header = "{'descr': '<f8', 'fortran_order': False, 'shape': " + shape_str + ", }";
    const std::size_t total = 10 + header.size() + 1;
    header.append((64 - total % 64) % 64, ' ');
    header += '\n';

    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("write_npy: cannot open " + path.string());
    out.write("\x93NUMPY\x01\x00", 8);
    const auto len = static_cast<std::uint16_t>(header.size());
    const char len_bytes[2] = {static_cast<char>(len & 0xff), static_cast<char>(len >> 8)};
    out.write(len_bytes, 2);
    out.write(header.data(), static_cast<std::streamsize>(header.size()));
    out.write(reinterpret_cast<const char*>(data.data()),
              static_cast<std::streamsize>(data.size() * sizeof(double)));
}

void ProbeRecorder::add(std::string name, const Component component, const std::size_t i,
                        const std::size_t j, const std::size_t k) {
    if (!times.empty()) throw std::logic_error("ProbeRecorder: add probes before recording");
    probes.push_back({std::move(name), component, {i, j, k}});
}

void ProbeRecorder::record(const FdtdSolver& solver, const double t) {
    std::vector<double> row;
    row.reserve(probes.size());
    for (const Probe& p : probes)
        row.push_back(field(solver, p.component)[solver.grid.to_flat(p.site[0], p.site[1], p.site[2])]);
    times.push_back(t);
    rows.push_back(std::move(row));
}

void ProbeRecorder::write_csv(const std::filesystem::path& path) const {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("ProbeRecorder: cannot open " + path.string());
    out.precision(17);
    out << "t";
    for (const Probe& p : probes) out << ',' << p.name;
    out << '\n';
    for (std::size_t n = 0; n < times.size(); ++n) {
        out << times[n];
        for (const double v : rows[n]) out << ',' << v;
        out << '\n';
    }
}

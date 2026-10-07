# FDTD Solver

A 3D finite-difference time-domain (FDTD) solver for Maxwell's equations,
written in C++20. It uses the standard Yee lattice and Taflove's lattice
truncation conditions as an absorbing boundary.


## Demo

`main.cpp` simulates a 61³ vacuum box with Δx = 1 and Δt = 0.5 for 400
steps. A soft current source `J_z` at the centre is driven by the derivative
of a Gaussian, so no net charge is left behind once the pulse has passed.
Output goes to the directory given as the first argument (default `out/`):

- `energy.csv`: total field energy `½ Σ(εE² + μH²) ΔV` at each step
- `probes.csv`: `E_z` and `H_y` at fixed points along +x
- `Ez_z30_NNNNN.npy`: the `E_z` slice through the source plane, every 10 steps

## Tests

- Grid: row-major flattening, and the round trip between flat and 3D indices.
- Solver: an impulse stays inside the numerical light cone; a centred source
  keeps the field mirror-symmetric; the absorbing boundary lets a pulse leave
  the box.
- Output: `.npy` headers are 64-byte aligned and followed by the raw data;
  probes record the value at their site on every step.



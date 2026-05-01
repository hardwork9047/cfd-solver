# CFD Solver Manuscript Draft

## Working Title

CFD Solver: Transparent Educational Solvers for Finite-Difference,
Lattice-Boltzmann, and Discrete-Element Method Studies

## Abstract Draft

This manuscript describes a compact Python library for computational fluid
dynamics education and small research studies. The library combines
finite-difference Navier-Stokes solvers, D2Q9 lattice-Boltzmann solvers, and a
discrete-element particle model behind simple, inspectable APIs. Verification
is organized around analytical Poiseuille solutions and established cavity-flow
benchmarks so that numerical behavior can be checked when the implementation
changes.

## Outline

1. Introduction and motivation
2. Solver architecture
3. Numerical methods
4. Verification against analytical and historical benchmark results
5. Demonstration cases
6. Limitations and future work

## Verification Sources To Cite

- Plane and circular Poiseuille analytical solutions.
- Ostwald-de Waele power-law channel-flow analytical solution.
- Ghia, Ghia, and Shin (1982), high-Reynolds-number lid-driven cavity benchmark.
- Zou and He (1997), LBM velocity boundary conditions.

## Notes

- Keep draft text and figures under `paper/`.
- Keep executable verification programs under `tests/verification/`.
- Keep runnable demonstrations under `src/demos/`.
- Keep generated outputs under `src/results/<program-name>/`.

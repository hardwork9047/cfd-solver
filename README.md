# CFD Solver Library

[![CI](https://img.shields.io/github/actions/workflow/status/hardwork9047/cfd-solver/ci.yml?branch=main&label=CI)](https://github.com/hardwork9047/cfd-solver/actions/workflows/ci.yml)
[![PyPI - Version](https://img.shields.io/pypi/v/cfd-solver.svg)](https://pypi.org/project/cfd-solver/)
[![Python](https://img.shields.io/pypi/pyversions/cfd-solver)](https://pypi.org/project/cfd-solver/)
[![License](https://img.shields.io/github/license/hardwork9047/cfd-solver)](LICENSE)

A Python library for Computational Fluid Dynamics covering Navier-Stokes
finite-difference solvers, the Lattice-Boltzmann Method (D2Q9), and the
Discrete Element Method — designed for education and research.

## Solvers at a glance

| Class | Package | Method | Physics | Reference |
|-------|---------|--------|---------|-----------|
| `PlanePoiseuille` | `cfd` | FDM + SOR | Pressure-driven channel flow | Analytical |
| `CircularPoiseuille` | `cfd` | FDM + SOR | Hagen-Poiseuille pipe flow | Analytical |
| `PowerLawPlanePoiseuille` | `cfd` | FDM iterative | Power-law non-Newtonian flow | Analytical |
| `CavityFlow` | `cfd` | NS projection | Lid-driven cavity | Ghia et al. (1982) |
| `CylinderFlow` | `cfd` | NS projection + upwind2 | Flow past cylinder / Kármán vortex | Williamson (1988) |
| `LBM` | `cfd_lbm` | D2Q9 BGK | Lid-driven cavity (mesoscopic) | Zou & He (1997) |
| `PoiseuilleFlow` | `cfd_lbm` | D2Q9 BGK | Channel flow (lattice units) | — |
| `ParticleSystem` | `dem` | DEM Hertz contact | Granular particle dynamics | — |
| `LBMDEMSolver` | `cfd_dem_lbm` | LBM + DEM coupled | Particle-laden flow | — |

## Installation

```bash
pip install cfd-solver
```

Development install:

```bash
git clone https://github.com/hardwork9047/cfd-solver.git
cd cfd-solver
poetry install
```

## Quick start

### Poiseuille flow (analytical vs. numerical)

```python
from cfd import PlanePoiseuille

flow = PlanePoiseuille(L=0.01, dp_dx=-100, mu=1e-3, ny=51)
u_exact = flow.analytical_solution()
u_num   = flow.numerical_solution()
print(f"Max velocity: {flow.get_max_velocity():.4f} m/s")
```

### Lid-driven cavity (Navier-Stokes projection)

```python
from cfd import CavityFlow

flow = CavityFlow(L=1.0, rho=1.0, mu=0.01)          # Re = 100
flow.solve_steady_state(max_iterations=2000)
fig = flow.plot_velocity_field()
fig.savefig("cavity_Re100.png")
```

### Kármán vortex shedding (Re = 200)

```python
from cfd import CylinderFlow

flow = CylinderFlow(
    L=4.0, H=2.0, D=0.2, U_inf=1.0, rho=1.0, mu=0.005,  # Re = 200
    nx=160, ny=80, advection_scheme="upwind2",
)
flow.run(t_end=20.0, output_interval=200)
fig = flow.plot_velocity_field()
```

### LBM cavity (D2Q9 Lattice-Boltzmann)

```python
from cfd_lbm import LBM

sim = LBM(nx=128, ny=128, Re=400.0, u_lid=0.1)
sim.advance(50_000)
rho, ux, uy = sim.get_fields()
```

### DEM granular simulation

```python
from dem import ParticleSystem

ps = ParticleSystem(n_particles=200, domain_size=(10.0, 10.0),
                    particle_radius=0.2, gravity=9.81, seed=42)
ps.run(t_end=5.0)
```

## Numerical methods

| Method | Where used | Key property |
|--------|-----------|--------------|
| FDM central difference (2nd order) | diffusion terms | accurate, simple |
| 1st-order upwind | CavityFlow default | diffusive, always stable |
| 2nd-order linear-upwind | CylinderFlow default | low numerical diffusion, Re_eff ≈ Re |
| SOR pressure Poisson | all Navier-Stokes solvers | converges with ω = 1.4–1.5 |
| D2Q9 BGK collision | LBM | lattice-scale, parallelisable |
| Zou-He velocity BC | LBM moving lid | momentum-conserving |
| Velocity Verlet | DEM | 2nd-order symplectic integration |
| Hertzian contact | DEM | physical overlap force model |

## Known limitations

- All solvers are 2-D and single-threaded.
- `CylinderFlow` Kármán shedding requires `advection_scheme="upwind2"`; the
  default 1st-order upwind damps vortices on coarse grids.
- `CavityFlow` remains stable up to approximately Re ≈ 1000 at default
  resolution (65×65). Higher Re requires finer grids.
- LBM is stable for `u_lid < 0.3` (lattice units) to satisfy the low-Mach
  assumption.
- No unstructured mesh support; domains must be rectangular.
- No turbulence modelling; only laminar flows are represented.

## Project structure

```
src/
├── cfd/            Navier-Stokes FDM solvers
├── cfd_lbm/        Lattice-Boltzmann Method (D2Q9)
├── dem/            Discrete Element Method
└── cfd_dem_lbm/    Coupled LBM-DEM solver
tests/              pytest test suite (44 tests)
examples/           Runnable demonstration scripts
.github/workflows/  CI (tests on Python 3.11–3.13)
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Please report bugs and request features
via [GitHub Issues](https://github.com/hardwork9047/cfd-solver/issues).

## Security

See [SECURITY.md](SECURITY.md) for vulnerability disclosure policy.

## License

MIT — see [LICENSE](LICENSE).

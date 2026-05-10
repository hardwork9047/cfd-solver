# LBM-DEM Algorithms

This document summarizes the numerical algorithms currently used by the
membrane fouling model. It is intended to answer what is computed in each
solver step, which options change the algorithm, and where the main numerical
assumptions are.

For the execution flow around these algorithms, see
`docs/fouling_model/PROCESS_FLOW.md`. For the repository-level component
layout, see `docs/fouling_model/SYSTEM_ARCHITECTURE.md`.

## State Variables

The coupled solver advances one Eulerian fluid state and one Lagrangian
particle state.

| State | Main arrays | Meaning |
|---|---|---|
| Fluid distribution | `f[q, x, y]` | D2Q9 LBM particle distribution functions |
| Fluid fields | `rho`, `ux`, `uy` | Density and velocity reconstructed from `f` |
| Fluid forcing | `Fx`, `Fy` | Body force, porous drag, IBM feedback, and particle reaction force |
| Fixed geometry | `fixed_solid` | Cylinder and optional wall solid mask |
| Dynamic particle mask | `particle_solid` | Particle-occupied cells for solid-boundary coupling |
| Particle state | `pos`, `vel`, `omega_p` | Particle center, translational velocity, and angular velocity |
| Particle properties | `radii`, `masses`, `inertias` | Per-particle disc radius, mass, and moment of inertia |

The pore geometry is represented by `PoreGeometry` and a list of `Cylinder`
objects. This keeps cylinder placement, fixed solid masks, pressure probe
sections, and cylinder VTK output in one layer.

## Fluid Algorithm

The fluid solver is a two-dimensional D2Q9 lattice Boltzmann method.

At every fluid step:

1. Macroscopic fields are reconstructed from the distribution functions.
2. External and coupling forces are prepared on lattice nodes.
3. Collision is applied with either BGK-Guo or TRT-Guo forcing.
4. Distributions are streamed to neighboring lattice nodes.
5. Bounce-back is applied on solid cells.
6. Streamwise pressure boundaries are applied when `--streamwise-boundary pressure` is used.

The pressure used in diagnostics is the LBM equation-of-state pressure:

```text
p = cs^2 rho,  cs^2 = 1/3
```

### BGK-Guo

`--fluid-method lbm-bgk-guo` uses a single relaxation time:

```text
f <- f - omega (f - feq) + Guo forcing
omega = 1 / tau
tau = nu / cs^2 + 0.5
```

This is the simplest and fastest fluid method. It is useful for screening and
small pressure-drop studies, but can be more sensitive to boundary and obstacle
resolution than TRT.

### TRT-Guo

`--fluid-method lbm-trt-guo` splits the distribution into symmetric and
antisymmetric parts and relaxes them with different rates. The positive
relaxation rate is still tied to viscosity, while the negative relaxation rate
is computed from the TRT magic parameter.

This method is usually preferred for pore-scale calculations because it gives
more control over boundary-related numerical error around cylinders and particle
solid masks.

## Boundary Conditions

The standard membrane-pore setting is:

```bash
--y-boundary periodic --streamwise-boundary pressure
```

This means:

- top and bottom are transverse periodic boundaries, not no-slip walls
- left is a pressure inlet
- right is a pressure outlet
- fixed cylinders inside the domain create the hydraulic resistance

For the pressure mode, the outlet density is set by `--rho-out`, and the inlet
density is computed from the requested pressure drop:

```text
rho_in = rho_out + pressure_drop / cs^2
```

The older `periodic-force` streamwise mode keeps periodic x streaming and drives
flow using a body force. It remains useful for verification and quick screening,
but it is less appropriate for interpreting left-to-right membrane pressure
loss.

## DEM Algorithm

Particle motion is handled by the shared `DEMSolver`.

The DEM solver computes force and torque contributions from:

- buoyancy-corrected gravity
- fluid-particle coupling
- particle-particle contact
- particle-cylinder contact
- optional particle-particle attraction or repulsion
- optional particle-cylinder attraction or repulsion
- optional sliding and rolling resistance

The particle state is advanced with substepping inside one LBM step. The DEM
time step is:

```text
dt_dem = dt_lbm / dem_substeps
```

This is important because contact stiffness is typically much faster than the
fluid time scale.

### Contact Models

Two normal contact laws are available:

| Option | Normal force |
|---|---|
| `--particle-method dem-linear` | proportional to overlap |
| `--particle-method dem-hertz` | proportional to overlap^(3/2) |

The normal damping term opposes the normal relative velocity and is clamped so
that a contact cannot create artificial tension.

When rolling friction is enabled, tangential contact force is limited by a
Coulomb cap and rolling resistance torque is limited by:

```text
|T_r| <= mu_r F_n r
```

This suppresses unrealistically easy rolling of spherical or disc-like particles
on cylinder surfaces.

### Near-Surface Attraction And Repulsion

Particle-particle and particle-cylinder attraction or repulsion use the same
Hamaker-like near-surface form:

```text
F ~ A r_eff / (6 h^2)
```

where `h` is the surface gap clipped by a minimum gap. Attraction acts toward
the other surface, while repulsion acts away from it. The cutoff distance limits
the interaction range.

Only one of attraction or repulsion is active for a given run.

## Particle Injection Algorithm

The filtration-style source mode is:

```bash
--particle-source left_inlet
```

In this mode, particles are not packed into the domain at time zero. Instead,
they are queued and gradually injected from the left boundary.

At each step:

1. The positive inlet flow rate is estimated from the left boundary velocity.
2. The cumulative inlet fluid area is updated.
3. The target injected particle area is computed from the source volume fraction.
4. New particles are released only when the particle-area budget allows it.
5. The injection y-position is sampled from local positive inlet velocity.
6. Particles that exit through the right boundary are removed from the active set.

This keeps the supplied particle amount tied to the incoming flow rate, which is
closer to membrane filtration than placing all particles initially.

## Particle-Fluid Coupling Algorithms

`--particle-fluid-coupling` selects how particles affect the fluid.

### `point_force`

This is the lightest coupling.

1. Fluid velocity is bilinearly interpolated to each particle center.
2. Stokes drag is computed from the slip velocity.
3. The equal and opposite reaction force is distributed back to nearby lattice nodes.

This is fast and smooth, but a particle does not become an impermeable solid
obstacle. It is suitable for dilute transport and screening, not for strict
channel-blockage physics.

### `solid_boundary`

Particle discs are rasterized onto lattice cells and added to the solid mask.
The fluid then sees those cells through bounce-back.

This is better for blockage and cake-like obstruction, but it is grid-sensitive:
small particles become stair-stepped solids, and moving no-slip behavior is only
an approximation at the current resolution.

### `immersed_boundary`

The immersed boundary method places Lagrangian marker points on each particle
surface.

For each marker:

1. The marker position is computed from the particle center and radius.
2. The desired boundary velocity is computed from particle translation and spin.
3. The fluid velocity is interpolated at the marker.
4. A direct-forcing correction is computed from the velocity mismatch.
5. The marker force is spread to nearby Eulerian lattice nodes.
6. The opposite force and torque are accumulated on the particle.

The marker spacing is controlled by `--ibm-marker-spacing`. A smaller spacing
uses more markers and gives a smoother boundary, but costs more. A larger
spacing is faster but can leak more flow through the particle boundary.

## Porous Resistance

When porous resistance is enabled, particle-occupied cells add a Brinkman-like
drag term to the fluid:

```text
F_porous = -alpha u
```

This is not a resolved solid boundary. It is a coarse model for cake-layer
hydraulic resistance and should be calibrated against resolved simulations or
experiments before being used quantitatively.

## Time-Stepping Order

At a high level, one coupled step performs:

1. particle source and outflow bookkeeping
2. particle-fluid coupling force preparation
3. DEM substeps for particle motion and contact
4. dynamic solid-mask or IBM force refresh when the chosen coupling requires it
5. LBM collision, streaming, bounce-back, and boundary update
6. metric and output sampling in the runner

The exact implementation is in `LBMDEMSolver.advance()` and the helper methods
in `src/cfd_dem_lbm/lbm_dem.py`.

## Metrics And Diagnostics

The runner records time-series metrics for later analysis:

| Metric | Purpose |
|---|---|
| inlet/outlet flux | checks throughput and mass transport trend |
| normalized flux | fouling proxy relative to the clean or early-time flux |
| pressure drop | hydraulic resistance proxy |
| local pressure sections | pore-scale pressure-loss diagnostics |
| active/generated/passed particles | retention and throughput accounting |
| contact count | particle-particle and particle-cylinder aggregation proxy |
| retained/deposited particles | fouling and cake-growth proxy |
| Re and related dimensionless groups | fair comparison across geometry and velocity settings |

ParaView outputs are written for fields such as velocity, pressure, particles,
cylinders, and optional IBM markers when requested.

## Acceleration Paths

The code supports NumPy and Numba execution paths:

| Option | Meaning |
|---|---|
| `numpy` | pure NumPy path |
| `numba` | require compiled kernels where implemented |
| `auto` | use Numba if available, otherwise fall back to NumPy |

Numba acceleration is most useful for repeated local operations such as contact
loads, interpolation, force spreading, and LBM kernels. The best accelerator
depends on particle count, grid size, coupling method, and whether the run is
long enough to amortize compilation time.

## Current Numerical Limitations

The present model is useful for systematic exploration, but these limitations
matter for interpretation:

- the physical model is two-dimensional
- pressure boundary stability depends on relaxation time, pressure drop, and obstacle resolution
- `point_force` coupling does not make particles impermeable
- `solid_boundary` coupling is sensitive to lattice resolution
- `immersed_boundary` leakage depends on marker spacing and forcing stiffness
- Brinkman resistance is a coarse cake-layer model, not a resolved pore-scale cake
- rolling friction and surface forces are phenomenological unless calibrated
- dense particle contact searches and periodic minimum-image handling can still dominate runtime

For quantitative fouling prediction, the next major algorithmic upgrades should
prioritize resolved solid-boundary validation, 3-D extension, calibrated
adhesion/repulsion laws, and systematic uncertainty analysis.

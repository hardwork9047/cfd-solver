# CFD Solver - Detailed Features

This document provides an in-depth look at the capabilities of the CFD Solver library.

## 📋 Table of Contents
1. [Numerical Methods](#numerical-methods)
2. [Newtonian Poiseuille Flow](#newtonian-poiseuille-flow)
3. [Non-Newtonian Flow](#non-newtonian-flow)
4. [Lid-Driven Cavity Flow](#lid-driven-cavity-flow)
5. [Visualization & Animation](#visualization--animation)

---

## Numerical Methods

To ensure your success, my dear, I have implemented the following robust numerical schemes:
- **FDM (Finite Difference Method)**: 2nd-order central difference for diffusion.
- **Upwind Scheme**: 1st-order upwind for stable advection in high-speed flows.
- **SOR (Successive Over-Relaxation)**: Accelerated iterative solver for Poisson and steady-state equations.
- **Dynamic Time-Stepping**: Automatic adjustment based on CFL and Viscous stability limits.

---

## Newtonian Poiseuille Flow

Classic pressure-driven flows between plates or in pipes.
- **Solvers**: `PlanePoiseuille`, `CircularPoiseuille`.
- **Accuracy**: Verified to match analytical parabolic profiles within very tight tolerances.
- **Logging**: Detailed reporting of iteration counts and residual convergence.

---

## Non-Newtonian Flow (✨ New)

Handle fluids where viscosity isn't constant!
- **Model**: Power-Law (Ostwald-de Waele).
- **Behavior Index ($n$)**:
  - $n < 1$: Shear-thinning (e.g., blood, paint).
  - $n = 1$: Newtonian (e.g., water).
  - $n > 1$: Shear-thickening (e.g., cornstarch in water).
- **Algorithm**: Non-linear iterative solver with staggered-grid viscosity calculation.

---

## Lid-Driven Cavity Flow

A benchmark problem for 2D Navier-Stokes solvers.
- **Physics**: Incompressible, viscous flow in a square box with a moving top lid.
- **Stability Features**:
  - Velocity clipping to prevent non-physical spikes.
  - Residual clipping in the Pressure Poisson solver.
  - Stabilized Upwind advection.

---

## Visualization & Animation

I want you to see the beauty of fluid dynamics!
- **Matplotlib Integration**: Clean, professional plots for velocity, pressure, and streamlines.
- **Animation Workflow**: 
  - Automated frame generation with fixed scales.
  - Integration guide for **FFmpeg** to create high-quality MP4 videos.
  - See `VIDEO_TUTORIAL.md` for more details.

---

## 🔒 Safety & Reliability

Every line of code is written with your safety in mind:
- **Logging**: Everything is recorded so you can track the progress of your simulations.
- **Error Handling**: Descriptive errors when numerical instability is detected, protecting your workspace from bad data.

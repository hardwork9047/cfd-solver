# Project Context: CFD Solver

## 📌 Background
This project is a Python-based Computational Fluid Dynamics (CFD) library designed for educational purposes, focusing on solving fundamental fluid flow problems using the Finite Difference Method (FDM).

## 🏗 Architecture
- **`src/cfd/`**: Core logic for different flow regimes.
  - `poiseuille.py`: 1D steady-state solvers for channel and pipe flows. Uses SOR.
  - `cavity.py`: 2D unsteady-state solver for lid-driven cavity flow. Solves incompressible Navier-Stokes equations.
- **`src/demos/`**: Demonstration scripts.
- **`src/results/<program-name>/`**: Generated plots, frames, videos, and CSV output.
- **`tests/`**: Automated unit tests using `pytest`.
- **`tests/verification/`**: Manual verification scripts and benchmark comparisons.
- **`paper/`**: Manuscript drafts and paper assets.

## 🧪 Current Status
- **Poiseuille Solvers**: Fully implemented and verified.
- **Cavity Solver**: Implemented with basic stability improvements; needs refined tutorial examples.
- **Tutorial**: Organized in `docs/TUTORIAL.md`.

## 🎯 Target Audience
Students and researchers learning numerical fluid dynamics.

## 🔒 Workspace & Scope
- **Boundary:** All operations confined to `/home/juntanishitani/projects/cfd`.
- **Access:** No access to external environment variables, SSH keys, or hidden configuration files. Shell commands limited to project simulation, testing, and documentation.

# Project Context: CFD Solver

## 📌 Background
This project is a Python-based Computational Fluid Dynamics (CFD) library designed for educational purposes. It focuses on solving fundamental fluid flow problems using the Finite Difference Method (FDM).

## 🏗 Architecture
- **`src/cfd/`**: Core logic for different flow regimes.
  - `poiseuille.py`: 1D steady-state solvers for channel and pipe flows. Uses SOR (Successive Over-Relaxation).
  - `cavity.py`: 2D unsteady-state solver for lid-driven cavity flow. Solves incompressible Navier-Stokes equations.
- **`examples/`**: Demonstration and verification scripts.
- **`tests/`**: Automated unit tests using `pytest`.

## 🧪 Current Status
- **Poiseuille Solvers**: Fully implemented and verified with analytical solutions.
- **Cavity Solver**: Implemented with basic stability improvements (dynamic DT, SOR for pressure). Needs refined tutorial examples.
- **Tutorial**: Organized into steps in `TUTORIAL.md`.

## 🎯 Target Audience
Students and researchers looking to understand the implementation details of numerical fluid dynamics without the overhead of industrial-grade software.

## 🔒 Security & Scope Rules
- **Workspace Boundary:** Do not access, read, or modify any files or directories outside of `/home/juntanishitani/projects/cfd`.
- **Credential Protection:** Never attempt to access or log environment variables, SSH keys, or hidden configuration files outside the project scope.
- **Tool Usage:** Use shell commands only for tasks directly related to CFD simulation, testing, and documentation within this workspace.

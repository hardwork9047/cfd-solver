# CFD Solver Library 🌊

[![PyPI - Version](https://img.shields.io/pypi/v/cfd-solver.svg)](https://pypi.org/project/cfd-solver/)
[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/hardwork9047/cfd-solver/main.yml?branch=main)](https://github.com/hardwork9047/cfd-solver/actions?query=workflow%3Abuild)
[![License](https://img.shields.io/github/license/hardwork9047/cfd-solver)](https://github.com/hardwork9047/cfd-solver/blob/main/LICENSE)

A specialized Python library for calculating fundamental flows in Computational Fluid Dynamics (CFD). Designed to be easy to learn while remaining numerically robust and suitable for educational and research purposes.

## 🎯 Key Features

| Feature | Description |
|---------|-------------|
| **Educational** | Clear implementations of FDM, SOR, and Upwind schemes |
| **Robust** | Comprehensive logging and error handling for numerical stability |
| **Non-Newtonian** | Support for Power-Law fluids (shear-thinning/thickening) |
| **Visual** | Built-in plotting and easy animation workflows |
| **Verified** | All solvers verified against analytical solutions |

## 🔧 Capabilities

### 1. Newtonian Poiseuille Flow
Classic pressure-driven laminar flows:
- **Plane Poiseuille**: Flow between parallel plates
- **Circular Poiseuille**: Flow in circular pipes
- **Features**: SOR solver, analytical verification, convergence tracking

### 2. Non-Newtonian Flow ✨
Power-Law model for variable viscosity fluids:
- **Shear-thinning** ($n < 1$): blood, paint
- **Newtonian** ($n = 1$): water
- **Shear-thickening** ($n > 1$): cornstarch suspension
- **Features**: Non-linear iterative solver with robust stability controls

### 3. Lid-Driven Cavity Flow
Benchmark 2D Navier-Stokes solver:
- **Physics**: Incompressible, viscous flow in square cavity
- **Stability**: Upwind differencing, dynamic time-stepping, numerical clipping
- **Outputs**: Velocity field, pressure distribution, vorticity analysis
- **Animation**: Frame generation for video creation

## 🛠 Installation

### Stable Release (Recommended for Users)

You can install `cfd-solver` directly from PyPI using `pip`:

```bash
pip install cfd-solver
```

### Development Version (For Contributors)

If you wish to contribute to the project or use the latest development features, you can install it from source. We use Poetry for dependency management:

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/hardwork9047/cfd-solver.git
    cd cfd-solver
    ```
2.  **Install dependencies with Poetry:**
    ```bash
    poetry install
    ```
    This will install all necessary development and runtime dependencies.


### Quick Start Example

Here's a quick way to get started with a simple Plane Poiseuille flow simulation:

```python
import matplotlib.pyplot as plt
from cfd import PlanePoiseuille

# Initialize the solver for a plane Poiseuille flow
# (distance between plates, pressure gradient, viscosity)
flow = PlanePoiseuille(L=0.01, dp_dx=-100, mu=1e-3, ny=51)

# Get analytical solution
u_analytical = flow.analytical_solution()

# Get numerical solution
u_numerical = flow.numerical_solution()

# Plot the velocity profile
plt.figure(figsize=(8, 5))
plt.plot(u_analytical, flow.y, label='Analytical Solution')
plt.plot(u_numerical, flow.y, 'o', markersize=4, fillstyle='none', label='Numerical Solution')
plt.xlabel('Velocity (m/s)')
plt.ylabel('Position (m)')
plt.title('Plane Poiseuille Flow')
plt.legend()
plt.grid(True)
plt.show()
```

## � Numerical Methods

To ensure robust and accurate simulations, this library implements:
- **FDM (Finite Difference Method)**: 2nd-order central difference for diffusion
- **Upwind Scheme**: 1st-order upwind for stable advection
- **SOR (Successive Over-Relaxation)**: Accelerated iterative solver
- **Dynamic Time-Stepping**: Automatic adjustment based on CFL and viscous stability limits

## 📂 Project Structure

```
cfd/
├── src/cfd/
│   ├── poiseuille.py       # SOR-based Newtonian solvers
│   ├── cavity.py           # Stabilized Navier-Stokes solver
│   └── non_newtonian.py    # Power-Law fluid solvers
├── examples/
│   ├── verify_*.py         # Rigorous verification scripts
│   ├── demo_*.py           # Feature demonstrations
│   └── make_animation_frames.py   # Animation frame generator
├── docs/
│   ├── VERIFICATION_REPORT.md  # Test results and performance
│   └── VIDEO_TUTORIAL.md       # Animation creation guide
└── TUTORIAL.md             # Step-by-step learning path
```

## 📖 Documentation

- **[TUTORIAL.md](TUTORIAL.md)**: Step-by-step learning path for beginners
- **[docs/VERIFICATION_REPORT.md](docs/VERIFICATION_REPORT.md)**: Detailed verification results and performance metrics
- **[docs/VIDEO_TUTORIAL.md](docs/VIDEO_TUTORIAL.md)**: Guide to creating animations with FFmpeg

## 🚀 Quick Demo

Want to see something cool right now? Try running the Lid-Driven Cavity demo:

```bash
poetry run python examples/demo.py
```

Or see how a shear-thinning fluid flows:

```bash
poetry run python examples/demo_shear_thinning.py
```

## 🔒 Engineering Standards

This library follows high standards for your safety and success:
- **Python 3.11+**
- **Robust Error Handling** (No more silent NaNs!)
- **Comprehensive Logging**
- **Strict Analytical Verification**
- **Test Coverage**: 15+ tests, all passing ✅

## 🤝 Contributing

We welcome contributions! Please see our [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to get started.

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---
*Created with care to help you master the world of fluids.* 🌊

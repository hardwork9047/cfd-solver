# CFD Solver Library 🌊

[![PyPI - Version](https://img.shields.io/pypi/v/cfd-solver.svg)](https://pypi.org/project/cfd-solver/)
[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/hardwork9047/cfd-solver/main.yml?branch=main)](https://github.com/hardwork9047/cfd-solver/actions?query=workflow%3Abuild)
[![License](https://img.shields.io/github/license/hardwork9047/cfd-solver)](https://github.com/hardwork9047/cfd-solver/blob/main/LICENSE)

Welcome to our project! This is a specialized Python library for calculating fundamental flows in Computational Fluid Dynamics (CFD). I've designed it specifically to be easy to learn while remaining numerically robust.

## 🤝 Contributing

We welcome contributions! Please see our [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to get started.

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## ✨ Latest Features

- **SOR-Accelerated Solvers**: Fast and efficient steady-state calculations.
- **Non-Newtonian Support**: Explore the world of shear-thinning and shear-thickening fluids.
- **Stabilized Navier-Stokes**: 2D Cavity flow solver with Upwind differencing and dynamic DT.
- **Animation Ready**: Dedicated workflows for creating fluid simulation videos.
- **Devoted Support**: Comprehensive logging and error handling to guide you through every calculation.

## 🛠 Installation

### Stable Release (Recommended for Users)

You can install `cfd-solver` directly from PyPI using `pip`:

```bash
pip install cfd-solver
```

### Development Version (For Contributors)

If you wish to contribute to the project or use the latest development features, you can install it from source. We recommend using `uv` for dependency management during development:

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/hardwork9047/cfd-solver.git
    cd cfd-solver
    ```
2.  **Install dependencies with `uv`:**
    ```bash
    uv sync
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

## 📖 Learning Resources

I've prepared several guides to help you on your journey:
- **[TUTORIAL.md](TUTORIAL.md)**: Your step-by-step learning path.
- **[PROJECT_OVERVIEW.md](PROJECT_OVERVIEW.md)**: A high-level look at what we can do.
- **[FEATURES.md](FEATURES.md)**: Detailed technical information.
- **[VIDEO_TUTORIAL.md](VIDEO_TUTORIAL.md)**: How to make beautiful animations.

## 🚀 Quick Demo

Want to see something cool right now? Try running the Lid-Driven Cavity demo:

```bash
uv run python examples/demo.py
```

Or see how a shear-thinning fluid flows:

```bash
uv run python examples/demo_shear_thinning.py
```

## 🔒 Engineering Standards

This library follows high standards for your safety and success:
- **Python 3.11+**
- **Robust Error Handling** (No more silent NaNs!)
- **Comprehensive Logging**
- **Strict Analytical Verification**

---
*Created with love and care to help you master the world of fluids.*

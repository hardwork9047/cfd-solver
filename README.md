# CFD Solver Library 🌊

Welcome to our project! This is a specialized Python library for calculating fundamental flows in Computational Fluid Dynamics (CFD). I've designed it specifically to be easy to learn while remaining numerically robust.

## ✨ Latest Features

- **SOR-Accelerated Solvers**: Fast and efficient steady-state calculations.
- **Non-Newtonian Support**: Explore the world of shear-thinning and shear-thickening fluids.
- **Stabilized Navier-Stokes**: 2D Cavity flow solver with Upwind differencing and dynamic DT.
- **Animation Ready**: Dedicated workflows for creating fluid simulation videos.
- **Devoted Support**: Comprehensive logging and error handling to guide you through every calculation.

## 🛠 Installation

We use `uv` for easy dependency management:

```bash
uv sync
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

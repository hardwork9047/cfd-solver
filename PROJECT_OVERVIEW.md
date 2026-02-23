# CFD Solver - Project Overview

## 📖 About this Library

**CFD Solver** is a Python library designed to solve fundamental problems in Computational Fluid Dynamics (CFD). It serves as both a powerful tool for basic simulations and a comprehensive educational resource for learning numerical methods.

### 🎯 Key Features

| Feature | Description |
|------|------|
| **Educational** | Clear implementations of FDM, SOR, and Upwind schemes. |
| **Robust** | Comprehensive logging and error handling for numerical stability. |
| **Non-Newtonian** | Support for Power-Law fluids (shear-thinning/thickening). |
| **Visual** | Built-in plotting and easy animation workflows with FFmpeg. |
| **Verified** | All solvers are verified against analytical solutions. |

---

## 🔧 Capabilities

### 1️⃣ Poiseuille Flow (Basic Laminar Flow)
- **Plane & Circular**: Steady-state solvers using the **SOR (Successive Over-Relaxation)** method.
- **Verification**: Automatic comparison with analytical solutions.
- **Logging**: Track convergence history and residuals.

### 2️⃣ Non-Newtonian Flow (New! ✨)
- **Power-Law Model**: Simulate shear-thinning (pseudoplastic) and shear-thickening (dilatant) fluids.
- **Iterative Solver**: Handles non-linear viscosity with robust stability controls.

### 3️⃣ Lid-Driven Cavity Flow
- **Navier-Stokes**: Solves 2D incompressible equations.
- **Stability**: Implements **Upwind differencing**, dynamic time-stepping, and numerical clipping to prevent divergence.
- **Animation**: Tools to generate frame sequences for video creation.

---

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
│   └── make_animation.py   # Animation frame generator
├── TUTORIAL.md             # Step-by-step learning path
└── VIDEO_TUTORIAL.md       # Guide to creating animations
```

---

## 🚀 Quick Start

### Installation
```bash
uv sync
```

### Running a Verification
```bash
uv run python examples/verify_plane_poiseuille.py
```

### Creating an Animation
```bash
uv run python examples/make_animation_frames.py
ffmpeg -framerate 10 -i frames/frame_%04d.png -c:v libx264 -pix_fmt yuv420p output.mp4
```

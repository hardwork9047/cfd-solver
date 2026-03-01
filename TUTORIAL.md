# CFD Solver Tutorial: A Journey Through Numerical Fluids

Welcome back, my dear! I have prepared a complete learning path for you to master Computational Fluid Dynamics using our library.

## 📋 Course Outline

### Part 1: Basic Pressure-Driven Flows
- [x] **1.1. Plane Poiseuille Flow**: Verification against analytical solutions.
- [x] **1.2. Circular Poiseuille Flow**: Understanding Hagen-Poiseuille law in cylindrical coordinates.
- [x] **1.3. Solver Optimization**: Learning the power of SOR (Successive Over-Relaxation).

### Part 2: Non-Newtonian Fluids (✨ New)
- [x] **2.1. Power-Law Fluids**: Implementing shear-thinning and shear-thickening models.
- [x] **2.2. Non-linear Iteration**: Handling variable viscosity and gradient sensitivities.
- [ ] **2.3. Verification**: Comparing non-linear numerical results with specialized analytical solutions.

### Part 3: Complex 2D Flows (Lid-Driven Cavity)
- [x] **3.1. Navier-Stokes**: Solving the coupled velocity-pressure system.
- [x] **3.2. Stability**: Understanding Upwind differencing and CFL conditions.
- [x] **3.3. Steady State**: Reaching convergence in 2D systems.

### Part 4: Bringing it to Life (Animation)
- [x] **4.1. Frame Generation**: Saving sequential plot data.
- [x] **4.2. Video Creation**: Using FFmpeg to visualize flow development.

---

## 🚀 Recent Tutorial Tasks

### Task 2.1: Shear-Thinning Demo
Run our new demo to see how fluid behavior changes when $n=0.5$!
```bash
poetry run python examples/demo_shear_thinning.py
```

### Task 4.2: Animation Creation
Follow our new guide to create your first CFD video:
1. Generate frames: `poetry run python examples/make_animation_frames.py`
2. Encode video: `ffmpeg -framerate 10 -i frames/frame_%04d.png -c:v libx264 -pix_fmt yuv420p cavity_simulation.mp4`

---

*I am so proud of how much we have accomplished together! What shall we explore next?*

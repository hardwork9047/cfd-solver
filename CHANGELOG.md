# CHANGELOG

## [0.2.0] - 2026-02-20

### Added
- **Non-Newtonian Solver**: New `PowerLawPlanePoiseuille` class for shear-thinning and shear-thickening fluids.
- **Animation Workflow**: `make_animation_frames.py` script and a comprehensive `VIDEO_TUTORIAL.md`.
- **Stabilized Navier-Stokes**: Re-implemented `CavityFlow` with Upwind differencing, dynamic DT, and clipping for superior stability.
- **Engineering Standards**: Added project-wide logging and robust error handling for all numerical solvers.

### Changed
- **SOR Upgrade**: Upgraded Poiseuille solvers from Gauss-Seidel to SOR (Successive Over-Relaxation) for much faster convergence.
- **Documentation Overhaul**: Updated `README.md`, `TUTORIAL.md`, `FEATURES.md`, and `PROJECT_OVERVIEW.md` to reflect the new capabilities.

## [0.1.0] - 2024-02-17

### Added
- 初回リリース
- ポアズイユ流れシミュレータ（平行平板間）
- Hagen-Poiseuille流れシミュレータ（円管内）
- 駆動キャビティ流れシミュレータ
- 包括的なテストスイート
- 詳細なドキュメント
- matplotlib統合の可視化機能

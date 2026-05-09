#!/usr/bin/env python3
"""CLI wrapper for the production-path LBM-DEM fluid verification suite."""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cfd.result_paths import program_results_dir
from cfd_dem_lbm import run_fluid_verification, write_verification_outputs


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify LBM fluid accuracy used by LBM-DEM")
    parser.add_argument("--nx", type=int, default=72)
    parser.add_argument("--ny", type=int, default=32)
    parser.add_argument("--steps", type=int, default=3000)
    parser.add_argument("--reynolds-number", "--Re", dest="re", type=float, default=40.0)
    parser.add_argument("--u-max", type=float, default=0.04)
    parser.add_argument("--poiseuille-tol", type=float, default=0.20)
    parser.add_argument("--flux-tol", type=float, default=0.05)
    parser.add_argument("--control-tol", type=float, default=0.20)
    parser.add_argument("--particle-solid-min-reduction", type=float, default=0.02)
    parser.add_argument("--ibm-min-force", type=float, default=1e-8)
    args = parser.parse_args()

    output_dir = program_results_dir(__file__, datetime.now().strftime("run_%Y%m%d_%H%M%S"))
    results = run_fluid_verification(
        nx=args.nx,
        ny=args.ny,
        steps=args.steps,
        re=args.re,
        u_max=args.u_max,
        poiseuille_tol=args.poiseuille_tol,
        flux_tol=args.flux_tol,
        control_tol=args.control_tol,
        particle_solid_min_reduction=args.particle_solid_min_reduction,
        ibm_min_force=args.ibm_min_force,
        output_dir=output_dir,
    )
    write_verification_outputs(results, output_dir)
    for result in results:
        status = "PASS" if result.passed else "FAIL"
        print(f"{status} {result.name}: metric={result.metric:.6g}, tol={result.tolerance:.6g}")
    print("Solver path: cfd_dem_lbm.FastLBMDEM")
    print(f"Verification outputs: {output_dir}")
    return 0 if all(result.passed for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())

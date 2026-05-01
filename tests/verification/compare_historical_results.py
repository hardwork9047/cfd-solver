"""
Compare solver outputs with classical analytical and benchmark results.

The script is a manual verification entry point. It complements the fast pytest
suite by collecting historical references that are useful when validating
numerical changes:

- Plane Poiseuille analytical solution.
- Hagen-Poiseuille circular-pipe analytical solution.
- Power-law plane Poiseuille analytical solution.
- Optional Ghia et al. (1982) lid-driven-cavity centerline benchmark.
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from cfd import CircularPoiseuille, PlanePoiseuille, PowerLawPlanePoiseuille  # noqa: E402
from cfd_lbm import LBM  # noqa: E402


@dataclass(frozen=True)
class ComparisonResult:
    """One scalar comparison against a reference value or profile."""

    name: str
    reference: str
    observed: float
    expected: float
    relative_error: float
    tolerance: float

    @property
    def passed(self) -> bool:
        return self.relative_error <= self.tolerance


GHIA_Y_U = np.array(
    [
        0.0000,
        0.0547,
        0.0625,
        0.0703,
        0.1016,
        0.1719,
        0.2813,
        0.4531,
        0.5000,
        0.6172,
        0.7344,
        0.8516,
        0.9531,
        0.9609,
        0.9688,
        0.9766,
        1.0000,
    ]
)
GHIA_U_RE400 = np.array(
    [
        0.00000,
        -0.08186,
        -0.09266,
        -0.10338,
        -0.14612,
        -0.24299,
        -0.32726,
        -0.17119,
        -0.11477,
        0.02135,
        0.16256,
        0.29093,
        0.55892,
        0.61756,
        0.68439,
        0.75837,
        1.00000,
    ]
)
GHIA_X_V = np.array(
    [
        0.0000,
        0.0625,
        0.0703,
        0.0781,
        0.0938,
        0.1563,
        0.2266,
        0.2344,
        0.5000,
        0.8047,
        0.8594,
        0.9063,
        0.9453,
        0.9531,
        0.9609,
        0.9688,
        1.0000,
    ]
)
GHIA_V_RE400 = np.array(
    [
        0.00000,
        0.09233,
        0.10091,
        0.10890,
        0.12317,
        0.16077,
        0.17507,
        0.17527,
        0.05454,
        -0.24533,
        -0.22445,
        -0.16914,
        -0.10313,
        -0.08864,
        -0.07391,
        -0.05906,
        0.00000,
    ]
)


def relative_l2(observed: np.ndarray, expected: np.ndarray) -> float:
    denom = max(float(np.linalg.norm(expected)), 1e-14)
    return float(np.linalg.norm(observed - expected) / denom)


def relative_error(observed: float, expected: float) -> float:
    return abs(observed - expected) / max(abs(expected), 1e-14)


def compare_plane_poiseuille() -> list[ComparisonResult]:
    flow = PlanePoiseuille(L=0.01, dp_dx=-100.0, mu=1e-3, ny=81)
    analytical = flow.analytical_solution()
    numerical = flow.numerical_solution(tol=1e-8, max_iter=30000)
    numerical_q = float(np.trapezoid(numerical, flow.y))

    return [
        ComparisonResult(
            name="Plane Poiseuille profile",
            reference="Classical parabolic analytical solution",
            observed=relative_l2(numerical, analytical),
            expected=0.0,
            relative_error=relative_l2(numerical, analytical),
            tolerance=5e-3,
        ),
        ComparisonResult(
            name="Plane Poiseuille flow rate",
            reference="Q = -dp/dx * L^3 / (12 mu)",
            observed=numerical_q,
            expected=flow.get_flow_rate(),
            relative_error=relative_error(numerical_q, flow.get_flow_rate()),
            tolerance=5e-3,
        ),
    ]


def compare_circular_poiseuille() -> list[ComparisonResult]:
    flow = CircularPoiseuille(R=0.005, dp_dx=-100.0, mu=1e-3, nr=81)
    analytical = flow.analytical_solution()
    numerical = flow.numerical_solution(tol=1e-8, max_iter=30000)
    numerical_q = float(2.0 * np.pi * np.trapezoid(numerical * flow.r, flow.r))

    return [
        ComparisonResult(
            name="Circular Poiseuille profile",
            reference="Hagen-Poiseuille analytical solution",
            observed=relative_l2(numerical, analytical),
            expected=0.0,
            relative_error=relative_l2(numerical, analytical),
            tolerance=1e-2,
        ),
        ComparisonResult(
            name="Circular Poiseuille flow rate",
            reference="Q = -pi * dp/dx * R^4 / (8 mu)",
            observed=numerical_q,
            expected=flow.get_flow_rate(),
            relative_error=relative_error(numerical_q, flow.get_flow_rate()),
            tolerance=1e-2,
        ),
    ]


def compare_power_law() -> list[ComparisonResult]:
    results: list[ComparisonResult] = []
    for n in (0.8, 1.4):
        flow = PowerLawPlanePoiseuille(L=0.01, dp_dx=-100.0, K=1e-3, n=n, ny=41)
        analytical = flow.analytical_solution()
        numerical = flow.numerical_solution(tol=1e-7, max_iter=60000, omega=0.1)
        err = relative_l2(numerical, analytical)
        results.append(
            ComparisonResult(
                name=f"Power-law Poiseuille profile n={n}",
                reference="Ostwald-de Waele plane-channel analytical solution",
                observed=err,
                expected=0.0,
                relative_error=err,
                tolerance=5e-3,
            )
        )
    return results


def compare_ghia_lbm(steps: int) -> list[ComparisonResult]:
    sim = LBM(nx=48, ny=48, Re=400.0, u_lid=0.1)
    sim.advance(steps)
    _, ux, uy = sim.get_fields()

    y_norm = np.linspace(0.0, 1.0, sim.ny)
    x_norm = np.linspace(0.0, 1.0, sim.nx)
    u_center = ux[sim.nx // 2, :] / sim.u_lid
    v_center = uy[:, sim.ny // 2] / sim.u_lid

    u_interp = np.interp(GHIA_Y_U, y_norm, u_center)
    v_interp = np.interp(GHIA_X_V, x_norm, v_center)
    u_err = relative_l2(u_interp, GHIA_U_RE400)
    v_err = relative_l2(v_interp, GHIA_V_RE400)

    return [
        ComparisonResult(
            name=f"LBM cavity u-centerline after {steps} steps",
            reference="Ghia et al. (1982), Re=400",
            observed=u_err,
            expected=0.0,
            relative_error=u_err,
            tolerance=0.8,
        ),
        ComparisonResult(
            name=f"LBM cavity v-centerline after {steps} steps",
            reference="Ghia et al. (1982), Re=400",
            observed=v_err,
            expected=0.0,
            relative_error=v_err,
            tolerance=0.8,
        ),
    ]


def collect_results(include_lbm: bool, lbm_steps: int) -> list[ComparisonResult]:
    results = []
    results.extend(compare_plane_poiseuille())
    results.extend(compare_circular_poiseuille())
    results.extend(compare_power_law())
    if include_lbm:
        results.extend(compare_ghia_lbm(lbm_steps))
    return results


def print_results(results: list[ComparisonResult]) -> None:
    header = f"{'status':<6} {'case':<46} {'rel_error':>12} {'tolerance':>12} reference"
    print(header)
    print("-" * len(header))
    for result in results:
        status = "PASS" if result.passed else "FAIL"
        print(
            f"{status:<6} {result.name:<46} "
            f"{result.relative_error:>12.4e} {result.tolerance:>12.4e} "
            f"{result.reference}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--include-lbm",
        action="store_true",
        help="include the slower LBM cavity comparison against Ghia et al. (1982)",
    )
    parser.add_argument(
        "--lbm-steps",
        type=int,
        default=5000,
        help="number of LBM steps for --include-lbm",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results = collect_results(include_lbm=args.include_lbm, lbm_steps=args.lbm_steps)
    print_results(results)
    return 0 if all(result.passed for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())

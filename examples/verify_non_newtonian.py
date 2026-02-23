import logging

import matplotlib.pyplot as plt
import numpy as np

from cfd import PowerLawPlanePoiseuille

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def verify_non_newtonian():
    """
    Verification of the Power-Law fluid solver against analytical solutions.
    """
    logger.info("Starting Detailed Verification of Power-Law Solver")

    L = 0.01
    dp_dx = -500
    K = 0.5

    # Test cases: Shear-thinning, Newtonian, Shear-thickening
    indices = [0.6, 1.0, 1.4]

    plt.figure(figsize=(12, 8))

    for n in indices:
        logger.info(f"--- Case n = {n} ---")
        fluid = PowerLawPlanePoiseuille(L, dp_dx, K, n, ny=101)

        u_ana = fluid.analytical_solution()
        # Using a very small omega for non-linear stability in shear-thinning cases
        current_omega = 0.05 if n < 1 else 0.2
        u_num = fluid.numerical_solution(tol=1e-9, omega=current_omega)

        max_error = np.max(np.abs(u_num - u_ana))
        rel_error = max_error / np.max(u_ana)

        logger.info(f"Max Error: {max_error:.2e} ({rel_error:.2%})")

        plt.plot(u_ana, fluid.y, "-", label=f"Analytical (n={n})", alpha=0.6)
        plt.plot(u_num, fluid.y, "o", markersize=3, label=f"Numerical (n={n})", markevery=4)

    plt.xlabel("Velocity u [m/s]")
    plt.ylabel("y [m]")
    plt.title("Verification: Power-Law Fluid Solver")
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.savefig("non_newtonian_verification.png")
    logger.info("Verification plot saved as 'non_newtonian_verification.png'")


if __name__ == "__main__":
    verify_non_newtonian()

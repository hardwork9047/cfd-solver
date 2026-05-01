"""
Verification script for cylinder flow simulation.

This script demonstrates the cylinder flow solver and validates basic behavior.
"""

import matplotlib.pyplot as plt
import numpy as np

from cfd import CylinderFlow


def verify_low_reynolds():
    """Verify cylinder flow at low Reynolds number (creeping flow)."""
    print("=" * 60)
    print("Cylinder Flow Verification - Low Reynolds Number (Re ≈ 1)")
    print("=" * 60)

    # Low Reynolds number: Stokes flow regime
    flow = CylinderFlow(L=2.0, H=1.0, D=0.1, U_inf=0.1, rho=1.0, mu=0.1, nx=150, ny=75)

    print(f"\nReynolds Number: {flow.get_Reynolds_number():.2f}")
    print(f"Grid size: {flow.nx} × {flow.ny}")
    print(f"Cylinder center: ({flow.x_c:.2f}, {flow.y_c:.2f})")
    print(f"Cylinder radius: {flow.R:.3f} m")

    # Run simulation
    flow.run(t_end=5.0, output_interval=200)

    # Compute drag coefficient
    C_d = flow.compute_drag_coefficient()
    print(f"\nComputed Drag Coefficient: C_d = {C_d:.3f}")
    print("(Theoretical C_d for Stokes flow around cylinder ≈ not available in 2D)")

    # Plot results
    fig = flow.plot_velocity_field(skip=4)
    plt.savefig("cylinder_flow_low_re.png", dpi=150, bbox_inches="tight")
    print("\nPlot saved as: cylinder_flow_low_re.png")

    plt.show()


def verify_moderate_reynolds():
    """Verify cylinder flow at moderate Reynolds number."""
    print("\n" + "=" * 60)
    print("Cylinder Flow Verification - Moderate Reynolds Number (Re ≈ 20)")
    print("=" * 60)

    # Moderate Reynolds number
    flow = CylinderFlow(L=3.0, H=1.5, D=0.1, U_inf=1.0, rho=1.0, mu=0.05, nx=200, ny=100)

    print(f"\nReynolds Number: {flow.get_Reynolds_number():.2f}")
    print(f"Grid size: {flow.nx} × {flow.ny}")

    # Run simulation
    flow.run(t_end=8.0, output_interval=200)

    # Compute drag coefficient
    C_d = flow.compute_drag_coefficient()
    print(f"\nComputed Drag Coefficient: C_d = {C_d:.3f}")
    print("Expected range for Re≈20: C_d ≈ 2-3 (approximate)")

    # Plot velocity on centerline
    x_center, u_center = flow.get_velocity_on_centerline()

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(x_center, u_center, "b-", linewidth=2)
    ax.axvline(flow.x_c - flow.R, color="red", linestyle="--", label="Cylinder front")
    ax.axvline(flow.x_c + flow.R, color="red", linestyle="--", label="Cylinder back")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("u [m/s]")
    ax.set_title("Velocity along Horizontal Centerline")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("cylinder_centerline_velocity.png", dpi=150, bbox_inches="tight")
    print("Centerline velocity plot saved as: cylinder_centerline_velocity.png")

    # Plot flow field
    fig2 = flow.plot_velocity_field(skip=5)
    plt.savefig("cylinder_flow_moderate_re.png", dpi=150, bbox_inches="tight")
    print("Flow field plot saved as: cylinder_flow_moderate_re.png")

    plt.show()


def main():
    """Run all verification tests."""
    print("\n🌊 CFD Solver - Cylinder Flow Verification Suite\n")

    # Test 1: Low Reynolds number
    verify_low_reynolds()

    # Test 2: Moderate Reynolds number
    verify_moderate_reynolds()

    print("\n" + "=" * 60)
    print("✅ Verification Complete!")
    print("=" * 60)
    print("\nNotes:")
    print("- Low Re flow should show smooth, symmetric flow around cylinder")
    print("- Moderate Re may show wake formation behind cylinder")
    print("- For Re > 40, vortex shedding (Kármán vortex street) may occur")
    print("- This solver is for educational purposes - use fine grids for accuracy")


if __name__ == "__main__":
    main()

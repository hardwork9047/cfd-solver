"""
Demo: Flow around a Circular Cylinder

This demo showcases the cylinder flow simulation with visualization.
"""

import matplotlib.pyplot as plt

from cfd import CylinderFlow


def demo_cylinder_flow():
    """Demonstrate cylinder flow simulation with nice visualization."""
    print("🌊 CFD Solver - Cylinder Flow Demo")
    print("=" * 60)
    print("Simulating flow around a circular cylinder...")
    print()

    # Create a cylinder flow simulation
    # Re = rho * U * D / mu = 1.0 * 1.0 * 0.1 / 0.01 = 10
    flow = CylinderFlow(
        L=2.5,  # Domain length [m]
        H=1.2,  # Domain height [m]
        D=0.1,  # Cylinder diameter [m]
        U_inf=1.0,  # Free-stream velocity [m/s]
        rho=1.0,  # Density [kg/m³]
        mu=0.01,  # Viscosity [Pa·s]
        nx=180,  # Grid points in x
        ny=90,  # Grid points in y
    )

    print(f"Reynolds Number: {flow.get_Reynolds_number():.1f}")
    print(f"Domain: {flow.L}m × {flow.H}m")
    print(f"Cylinder: D={flow.D}m at ({flow.x_c:.2f}, {flow.y_c:.2f})")
    print(f"Grid Resolution: {flow.nx} × {flow.ny}")
    print()

    # Run the simulation
    print("Running simulation...")
    flow.run(t_end=6.0, output_interval=150)

    # Compute drag coefficient
    C_d = flow.compute_drag_coefficient()
    print(f"\n✨ Drag Coefficient: C_d = {C_d:.3f}")

    # Create visualization
    print("\nGenerating visualization...")
    fig = flow.plot_velocity_field(skip=4)
    plt.savefig("demo_cylinder_flow.png", dpi=150, bbox_inches="tight")
    print("Plot saved as: demo_cylinder_flow.png")

    # Show the plot
    plt.show()

    print("\n✅ Demo complete!")
    print("\nWhat to observe:")
    print("- Velocity decreases near the cylinder surface (no-slip)")
    print("- Flow accelerates around the cylinder sides")
    print("- Wake region forms behind the cylinder")
    print("- For higher Re, you may see vortex shedding!")


if __name__ == "__main__":
    demo_cylinder_flow()

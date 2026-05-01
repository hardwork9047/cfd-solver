import matplotlib.pyplot as plt

from cfd import CavityFlow
from cfd.result_paths import program_results_dir


def demo_cavity():
    """
    Part 2.1: Demo of Lid-Driven Cavity Flow.
    """
    print("--- Lid-Driven Cavity Flow Demo ---")

    # Re = 100
    cavity = CavityFlow(L=1.0, rho=1.0, mu=0.01, nx=41, ny=41)

    # Run simulation
    cavity.solve_steady_state(max_iterations=200)

    output_dir = program_results_dir(__file__)

    # Plot results
    fig1 = cavity.plot_velocity_field()
    velocity_path = output_dir / "cavity_velocity.png"
    plt.savefig(velocity_path)

    fig2 = cavity.plot_streamlines()
    streamlines_path = output_dir / "cavity_streamlines.png"
    plt.savefig(streamlines_path)

    print(f"\nResults saved as '{velocity_path}' and '{streamlines_path}'")


if __name__ == "__main__":
    demo_cavity()

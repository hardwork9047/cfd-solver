import matplotlib.pyplot as plt

from cfd import CavityFlow


def demo_cavity():
    """
    Part 2.1: Demo of Lid-Driven Cavity Flow.
    """
    print("--- Lid-Driven Cavity Flow Demo ---")

    # Re = 100
    cavity = CavityFlow(L=1.0, rho=1.0, mu=0.01, nx=41, ny=41)

    # Run simulation
    cavity.solve_steady_state(max_iterations=200)

    # Plot results
    fig1 = cavity.plot_velocity_field()
    plt.savefig("cavity_velocity.png")

    fig2 = cavity.plot_streamlines()
    plt.savefig("cavity_streamlines.png")

    print("\nResults saved as 'cavity_velocity.png' and 'cavity_streamlines.png'")


if __name__ == "__main__":
    demo_cavity()

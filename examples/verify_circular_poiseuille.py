import numpy as np
import matplotlib.pyplot as plt
from cfd import CircularPoiseuille

def verify_circular_poiseuille():
    """
    Step 1.2 of the Tutorial: Verification of Circular Poiseuille Flow.
    Compares the numerical solution with the analytical solution (Hagen-Poiseuille).
    """
    print("--- Circular Poiseuille Flow Verification ---")
    
    # Parameters
    R = 0.005     # Pipe radius: 5mm
    dp_dx = -100  # Pressure gradient: -100 Pa/m
    mu = 1e-3     # Dynamic viscosity: 1 mPa·s (water)
    
    flow = CircularPoiseuille(R, dp_dx, mu, nr=51)
    u_ana = flow.analytical_solution()
    u_num = flow.numerical_solution(tol=1e-9, omega=1.7)
    
    # Check max velocity
    u_max_ana = flow.get_max_velocity()
    u_max_num = np.max(u_num)
    
    # Check flow rate
    Q_ana = flow.get_flow_rate()
    # Numerical integration of flow rate: Q = integral(u * 2*pi*r * dr)
    dr = flow.r[1] - flow.r[0]
    Q_num = np.sum(u_num * 2 * np.pi * flow.r * dr)
    
    print(f"\nMax Velocity (Analytical): {u_max_ana:.6f} m/s")
    print(f"Max Velocity (Numerical):  {u_max_num:.6f} m/s")
    print(f"Velocity Error:            {abs(u_max_ana - u_max_num):.2e} m/s")
    
    print(f"\nFlow Rate (Analytical): {Q_ana:.6e} m^3/s")
    print(f"Flow Rate (Numerical):  {Q_num:.6e} m^3/s")
    print(f"Flow Rate Error:        {abs(Q_ana - Q_num):.2e} m^3/s")
    
    # Plotting
    flow.plot()
    plt.savefig('circular_poiseuille_verification.png')
    print("\nVerification plot saved as 'circular_poiseuille_verification.png'")

if __name__ == "__main__":
    verify_circular_poiseuille()

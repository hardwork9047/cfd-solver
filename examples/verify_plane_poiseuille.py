import logging
import numpy as np
import matplotlib.pyplot as plt
from cfd import PlanePoiseuille

# Configure logging to see our beautiful progress
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def verify_plane_poiseuille():
    """
    Step 1.1 of the Tutorial: Verification of Plane Poiseuille Flow.
    Compares the numerical solution with the analytical solution for different grid sizes.
    """
    logger.info("Starting Plane Poiseuille Flow Verification Process")
    
    # Parameters for our simulation
    L = 0.01      # Plate distance: 10mm
    dp_dx = -100  # Pressure gradient: -100 Pa/m
    mu = 1e-3     # Dynamic viscosity: 1 mPa·s (water)
    
    grid_sizes = [11, 21, 51, 101]
    errors = []
    
    plt.figure(figsize=(10, 6))
    
    # Let's get a smooth reference for the analytical solution
    flow_ref = PlanePoiseuille(L, dp_dx, mu, ny=201)
    u_ana_ref = flow_ref.analytical_solution()
    plt.plot(u_ana_ref, flow_ref.y, 'k--', label='Analytical Solution', alpha=0.7)
    
    for ny in grid_sizes:
        try:
            logger.info(f"Analyzing grid size: ny = {ny}")
            flow = PlanePoiseuille(L, dp_dx, mu, ny=ny)
            u_ana = flow.analytical_solution()
            u_num = flow.numerical_solution(tol=1e-8)
            
            # Calculate the maximum deviation
            error = np.max(np.abs(u_num - u_ana))
            errors.append(error)
            logger.info(f"Results for ny={ny}: Max Error = {error:.2e}")
            
            plt.plot(u_num, flow.y, 'o-', markersize=4, label=f'Numerical (ny={ny})')
        except Exception as e:
            logger.error(f"Something went wrong with ny={ny}: {e}")
    
    plt.xlabel('Velocity u [m/s]')
    plt.ylabel('y [m]')
    plt.title('Plane Poiseuille Flow: Grid Convergence')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Save our work as a beautiful image
    output_file = 'plane_poiseuille_verification.png'
    plt.savefig(output_file)
    logger.info(f"Verification plot successfully saved as '{output_file}'")
    
    # Final error summary for you
    logger.info("--- Final Error Analysis Summary ---")
    for i, ny in enumerate(grid_sizes):
        logger.info(f"  Grid {ny:3d}: Error = {errors[i]:.2e}")

if __name__ == "__main__":
    verify_plane_poiseuille()

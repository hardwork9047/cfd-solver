import logging
import matplotlib.pyplot as plt
import numpy as np
from cfd import PowerLawPlanePoiseuille, PlanePoiseuille

# Setting up logging to track our success
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def run_shear_thinning_demo():
    """
    Demonstrates the effect of the flow behavior index (n) on the velocity profile.
    Compare shear-thinning (n=0.5), Newtonian (n=1.0), and shear-thickening (n=1.5).
    """
    logger.info("Starting Shear-Thinning vs Newtonian Comparison")
    
    L = 0.02      # 20mm gap
    dp_dx = -200  # Pressure gradient
    K = 0.1       # Consistency index
    
    # Indices to compare
    indices = [0.5, 1.0, 1.5]
    labels = ["Shear-thinning (n=0.5)", "Newtonian (n=1.0)", "Shear-thickening (n=1.5)"]
    colors = ['r', 'k', 'b']
    
    plt.figure(figsize=(10, 7))
    
    for n, label, color in zip(indices, labels, colors):
        try:
            logger.info(f"Calculating for {label}...")
            # For n=1, PowerLaw fluid is Newtonian with mu=K
            fluid = PowerLawPlanePoiseuille(L, dp_dx, K, n, ny=101)
            
            # Calculate numerical solution
            u_num = fluid.numerical_solution(tol=1e-8)
            
            # Normalize for better visual comparison of the SHAPE
            u_norm = u_num / np.max(u_num)
            
            plt.plot(u_norm, fluid.y, color, label=label, linewidth=2)
            
        except Exception as e:
            logger.error(f"Error during calculation for n={n}: {e}")
            
    plt.xlabel('Normalized Velocity u/u_max')
    plt.ylabel('Position y [m]')
    plt.title('Velocity Profile Comparison: Effect of Flow Behavior Index (n)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    output_image = "shear_thinning_comparison.png"
    plt.savefig(output_image)
    logger.info(f"Comparison plot saved to {output_image}")
    plt.show()

if __name__ == "__main__":
    run_shear_thinning_demo()

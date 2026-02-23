import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from cfd import CavityFlow

# Let's keep track of our progress
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def generate_frames():
    # Create a directory for our frames
    output_dir = "frames"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f"Created directory: {output_dir}")

    # Use parameters that are proven to work
    L = 1.0
    rho = 1.0
    mu = 0.01  # Re = 100
    cavity = CavityFlow(L, rho, mu, nx=65, ny=65)
    
    # We want a fixed scale for the colorbar so it doesn't flicker!
    v_max = 1.0
    
    # We'll run the simulation and save frames along the way
    num_steps = 300
    save_interval = 10
    
    logger.info(f"Starting frame generation: {num_steps} iterations...")
    
    frame_count = 0
    for i in range(num_steps):
        try:
            cavity.step()
        except Exception as e:
            logger.error(f"Failed at step {i}: {e}")
            break
            
        if (i + 1) % save_interval == 0:
            frame_count += 1
            fig, ax = plt.subplots(figsize=(8, 6))
            
            # Speed magnitude
            speed = np.sqrt(cavity.u**2 + cavity.v**2)
            
            # Fixed levels for stability
            levels = np.linspace(0, v_max, 21)
            contour = ax.contourf(cavity.X, cavity.Y, speed, levels=levels, cmap='viridis', extend='both')
            
            # Streamlines
            ax.streamplot(cavity.x, cavity.y, cavity.u, cavity.v, color='white', linewidth=0.5, density=1.0)
            
            ax.set_title(f'Cavity Flow Evolution (Re=100) - Step: {i+1}')
            ax.set_aspect('equal')
            plt.colorbar(contour, ax=ax, label='Velocity Magnitude [m/s]')
            
            # Use the naming convention from our tutorial!
            filename = os.path.join(output_dir, f"frame_{frame_count:04d}.png")
            plt.savefig(filename, dpi=100)
            plt.close(fig)
            
            logger.info(f"Saved frame {frame_count} (Iteration {i+1})")

    logger.info(f"Done! {frame_count} frames saved in '{output_dir}/'")

if __name__ == "__main__":
    generate_frames()

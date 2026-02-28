"""
Analyze particle data from CSV to debug DEM simulation
"""
import csv
import matplotlib.pyplot as plt
import numpy as np

# Load data
data = []
with open('results/particle_data.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        data.append({k: float(v) for k, v in row.items()})

print("Data loaded:", len(data), "rows")
print("Time range:", data[0]['time'], "to", data[-1]['time'])
print()

# Analyze particle 0 (bottom left particle)
print("=== Particle 0 Analysis ===")
print(f"Initial position: ({data[0]['p0_x']:.3f}, {data[0]['p0_y']:.3f})")
print(f"Final position: ({data[-1]['p0_x']:.3f}, {data[-1]['p0_y']:.3f})")
print(f"Position change: ({data[-1]['p0_x'] - data[0]['p0_x']:.6f}, {data[-1]['p0_y'] - data[0]['p0_y']:.6f})")
print()

# Check when particle hits floor (radius = 0.1, so floor contact at y = 0.1)
radius = 0.1
floor_y = radius

first_contact = None
for i, row in enumerate(data):
    if row['p0_y'] <= floor_y:
        first_contact = i
        break

if first_contact is not None:
    print(f"Particle 0 hits floor at step {first_contact}, time {data[first_contact]['time']:.4f}s")
    print(f"Position at contact: ({data[first_contact]['p0_x']:.3f}, {data[first_contact]['p0_y']:.3f})")
    print(f"Velocity at contact: ({data[first_contact]['p0_vx']:.3f}, {data[first_contact]['p0_vy']:.3f})")
    print(f"Force at contact: ({data[first_contact]['p0_fx']:.3f}, {data[first_contact]['p0_fy']:.3f})")
    print()
    
    # Check a few steps after contact
    if first_contact + 100 < len(data):
        step = first_contact + 100
        print(f"After 100 steps (step {step}):")
        print(f"  Position: ({data[step]['p0_x']:.6f}, {data[step]['p0_y']:.6f})")
        print(f"  Velocity: ({data[step]['p0_vx']:.6f}, {data[step]['p0_vy']:.6f})")
        print(f"  Force: ({data[step]['p0_fx']:.3f}, {data[step]['p0_fy']:.3f})")
        
        penetration = radius - data[step]['p0_y']
        print(f"  Floor penetration: {penetration:.6f}m")
        
        # Expected force
        mass = 78.54
        gravity = 9.81
        k_n = 5e4
        damping = 0.8
        
        if penetration > 0:
            f_n_expected = k_n * (penetration ** 1.5)
            v_n = data[step]['p0_vy']
            f_d_expected = -damping * v_n * np.sqrt(k_n * mass)
            total_expected = f_n_expected + f_d_expected
            
            print(f"  Expected normal force: {f_n_expected:.3f}N")
            print(f"  Expected damping force: {f_d_expected:.3f}N")
            print(f"  Expected total contact: {total_expected:.3f}N")
            print(f"  Expected net force: {total_expected - mass * gravity:.3f}N")
            print(f"  Actual net force: {data[step]['p0_fy']:.3f}N")
        print()

else:
    print("Particle 0 never hits floor in simulation")
    print()

# Check final state
print("=== Final State (last 10 steps average) ===")
for i in range(4):
    final_y = np.mean([data[j][f'p{i}_y'] for j in range(len(data)-10, len(data))])
    final_vy = np.mean([data[j][f'p{i}_vy'] for j in range(len(data)-10, len(data))])
    final_fy = np.mean([data[j][f'p{i}_fy'] for j in range(len(data)-10, len(data))])
    penetration = radius - final_y
    print(f"Particle {i}:")
    print(f"  Y position: {final_y:.6f}m, penetration: {penetration:.6f}m")
    print(f"  Y velocity: {final_vy:.6f} m/s")
    print(f"  Y force: {final_fy:.3f} N (gravity: -770.48N)")
    
    if penetration > 0:
        # This particle is penetrating the floor
        # Check if force is increasing with penetration
        print(f"  WARNING: Particle penetrating floor by {penetration*1000:.2f}mm")
        
print()

# Create plots
fig, axes = plt.subplots(3, 1, figsize=(12, 10))

# Extract time series data
times = [row['time'] for row in data]
p0_y = [row['p0_y'] for row in data]
p1_y = [row['p1_y'] for row in data]
p2_y = [row['p2_y'] for row in data]
p3_y = [row['p3_y'] for row in data]

p0_vy = [row['p0_vy'] for row in data]
p1_vy = [row['p1_vy'] for row in data]
p2_vy = [row['p2_vy'] for row in data]
p3_vy = [row['p3_vy'] for row in data]

p0_fy = [row['p0_fy'] for row in data]
p1_fy = [row['p1_fy'] for row in data]
p2_fy = [row['p2_fy'] for row in data]
p3_fy = [row['p3_fy'] for row in data]

# Plot 1: Y position vs time
ax = axes[0]
ax.plot(times, p0_y, label='Particle 0')
ax.plot(times, p1_y, label='Particle 1')
ax.plot(times, p2_y, label='Particle 2')
ax.plot(times, p3_y, label='Particle 3')
ax.axhline(y=radius, color='k', linestyle='--', label='Floor (y=0.1)')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Y Position [m]')
ax.set_title('Particle Y Positions vs Time')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 2: Y velocity vs time
ax = axes[1]
ax.plot(times, p0_vy, label='Particle 0')
ax.plot(times, p1_vy, label='Particle 1')
ax.plot(times, p2_vy, label='Particle 2')
ax.plot(times, p3_vy, label='Particle 3')
ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Y Velocity [m/s]')
ax.set_title('Particle Y Velocities vs Time')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 3: Y force vs time
ax = axes[2]
ax.plot(times, p0_fy, label='Particle 0')
ax.plot(times, p1_fy, label='Particle 1')
ax.plot(times, p2_fy, label='Particle 2')
ax.plot(times, p3_fy, label='Particle 3')
ax.axhline(y=-770.48, color='k', linestyle='--', label='Gravity')
ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Y Force [N]')
ax.set_title('Particle Y Forces vs Time')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('results/images/particle_analysis.png', dpi=150)
print("Analysis plot saved to: results/images/particle_analysis.png")

plt.show()

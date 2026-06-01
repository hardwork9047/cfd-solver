#!/usr/bin/env python3
"""
Graphing utility for 3D LBM-DEM fouling simulation analysis.
Generates three comprehensive figure sets from time_series.csv.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from matplotlib import rcParams

# Font configuration
rcParams['font.sans-serif'] = ['DejaVu Sans']
rcParams['axes.unicode_minus'] = False


def generate_time_series_analysis(result_dir: Path, output_dir: Path = None) -> Path:
    """
    Generate 9-panel time series analysis figure.

    Includes:
    1. Particle count over time
    2. Generation & removal cumulative
    3. Throughput (flux) variation
    4. Maximum flow velocity
    5. Particle velocity profile
    6. Fluid density stability
    7. Capture progress ratio
    8. Velocity decline rate
    9. Permeability decline (normalized flux)
    """
    if output_dir is None:
        output_dir = result_dir

    csv_path = result_dir / "analysis" / "time_series.csv"
    df = pd.read_csv(csv_path)
    df['passed_rate'] = df['removed_particles'].diff() / df['step'].diff()

    fig = plt.figure(figsize=(16, 12))

    # Panel 1: Particle count
    ax1 = plt.subplot(3, 3, 1)
    ax1.plot(df['step'], df['particle_count'], 'o-', color='#2E86AB', linewidth=2, markersize=4)
    ax1.set_xlabel('Step', fontsize=10)
    ax1.set_ylabel('Number of particles', fontsize=10)
    ax1.set_title('Particle Count Over Time', fontsize=11, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    # Panel 2: Generation & removal
    ax2 = plt.subplot(3, 3, 2)
    ax2.plot(df['step'], df['generated_particles'], 'o-', color='#A23B72', linewidth=2, markersize=4, label='Generated')
    ax2.plot(df['step'], df['removed_particles'], 's-', color='#F18F01', linewidth=2, markersize=4, label='Removed')
    ax2.set_xlabel('Step', fontsize=10)
    ax2.set_ylabel('Cumulative count', fontsize=10)
    ax2.set_title('Generation & Removal', fontsize=11, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9)

    # Panel 3: Throughput flux
    ax3 = plt.subplot(3, 3, 3)
    ax3.plot(df['step'][1:], df['passed_rate'][1:], 'o-', color='#C73E1D', linewidth=2, markersize=4)
    mean_rate = df['passed_rate'][1:].mean()
    ax3.axhline(y=mean_rate, color='red', linestyle='--', linewidth=1.5, label=f'Mean: {mean_rate:.4f}')
    ax3.set_xlabel('Step', fontsize=10)
    ax3.set_ylabel('Passed rate (particles/step)', fontsize=10)
    ax3.set_title('Particle Throughput (Flux)', fontsize=11, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=9)

    # Panel 4: Maximum flow velocity
    ax4 = plt.subplot(3, 3, 4)
    ax4.plot(df['step'], df['u_max'], 'o-', color='#06A77D', linewidth=2.5, markersize=5)
    ax4.fill_between(df['step'], df['u_max'], alpha=0.3, color='#06A77D')
    ax4.set_xlabel('Step', fontsize=10)
    ax4.set_ylabel('u_max (LU/step)', fontsize=10)
    ax4.set_title('Maximum Flow Velocity', fontsize=11, fontweight='bold')
    ax4.grid(True, alpha=0.3)

    # Panel 5: Particle velocity profile
    ax5 = plt.subplot(3, 3, 5)
    ax5.plot(df['step'], df['particle_speed_mean'], 'o-', color='#FF006E', linewidth=2, markersize=4, label='Mean')
    ax5.plot(df['step'], df['particle_speed_max'], 's-', color='#8338EC', linewidth=2, markersize=4, label='Max')
    ax5.set_xlabel('Step', fontsize=10)
    ax5.set_ylabel('Particle speed (LU/step)', fontsize=10)
    ax5.set_title('Particle Velocity Profile', fontsize=11, fontweight='bold')
    ax5.grid(True, alpha=0.3)
    ax5.legend(fontsize=9)

    # Panel 6: Density stability
    ax6 = plt.subplot(3, 3, 6)
    ax6.plot(df['step'], df['rho_mean'], 'o-', color='#1F77B4', linewidth=2, markersize=4)
    ax6.axhline(y=1.0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax6.set_xlabel('Step', fontsize=10)
    ax6.set_ylabel('Mean density', fontsize=10)
    ax6.set_title('Fluid Density Stability', fontsize=11, fontweight='bold')
    ax6.grid(True, alpha=0.3)
    ax6.set_ylim([1.0009, 1.0016])

    # Panel 7: Capture progress ratio
    ax7 = plt.subplot(3, 3, 7)
    capture_ratio = 100 * df['particle_count'] / (df['generated_particles'] + 1e-10)
    ax7.plot(df['step'], capture_ratio, 'o-', color='#DC2F02', linewidth=2.5, markersize=5)
    ax7.fill_between(df['step'], capture_ratio, alpha=0.2, color='#DC2F02')
    ax7.set_xlabel('Step', fontsize=10)
    ax7.set_ylabel('Capture ratio (%)', fontsize=10)
    ax7.set_title('Particle Capture Progress', fontsize=11, fontweight='bold')
    ax7.grid(True, alpha=0.3)

    # Panel 8: Velocity decline rate
    ax8 = plt.subplot(3, 3, 8)
    u_max_norm = (df['u_max'].values - df['u_max'].iloc[0]) / df['u_max'].iloc[0] * 100
    ax8.plot(df['step'], u_max_norm, 'o-', color='#FB5607', linewidth=2.5, markersize=5)
    ax8.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax8.fill_between(df['step'], u_max_norm, 0, alpha=0.2, color='#FB5607')
    ax8.set_xlabel('Step', fontsize=10)
    ax8.set_ylabel('Velocity change (%)', fontsize=10)
    ax8.set_title('Flow Velocity Decline Rate', fontsize=11, fontweight='bold')
    ax8.grid(True, alpha=0.3)

    # Panel 9: Permeability decline (normalized throughput)
    ax9 = plt.subplot(3, 3, 9)
    # Calculate normalized permeability (throughput normalized to mean, not first point which is often noisy)
    passed_rate_mean = df['passed_rate'][5:].mean()  # Use stable middle period
    passed_rate_normalized = df['passed_rate'][1:].values / passed_rate_mean * 100
    ax9.plot(df['step'][1:], passed_rate_normalized, 'o-', color='#D62828', linewidth=2.5, markersize=5)
    ax9.axhline(y=100, color='gray', linestyle='--', linewidth=1.5, alpha=0.6, label='Reference')
    ax9.fill_between(df['step'][1:], passed_rate_normalized, 100,
                     where=(passed_rate_normalized < 100), alpha=0.25, color='#D62828', label='Fouling region')
    ax9.set_xlabel('Step', fontsize=10)
    ax9.set_ylabel('Permeability (Relative, %)', fontsize=10)
    ax9.set_title('Membrane Permeability Decline', fontsize=11, fontweight='bold')
    ax9.grid(True, alpha=0.3)
    ax9.legend(fontsize=8, loc='best')
    ax9.set_ylim([85, 115])

    # Add text annotation for final permeability
    final_perm = passed_rate_normalized[-1]
    perm_change = final_perm - 100
    ax9.text(0.98, 0.05, f'Final: {final_perm:.1f}%\nChange: {perm_change:+.1f}%',
            transform=ax9.transAxes, fontsize=9, verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='lightyellow' if perm_change > 0 else 'lightcoral', alpha=0.7))

    plt.tight_layout()
    output_path = output_dir / "time_series_analysis.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def generate_permeability_analysis(result_dir: Path, output_dir: Path = None) -> Path:
    """
    Generate focused permeability decline analysis figure.

    Includes:
    1. Normalized permeability vs time (primary metric)
    2. Particle capture vs permeability loss (dual axis)
    3. Permeability decline rate (dP/dt)
    4. Fouling resistance index
    """
    if output_dir is None:
        output_dir = result_dir

    csv_path = result_dir / "analysis" / "time_series.csv"
    df = pd.read_csv(csv_path)
    df['passed_rate'] = df['removed_particles'].diff() / df['step'].diff()

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel 1: Normalized permeability
    ax1 = axes[0, 0]
    passed_rate_mean = df['passed_rate'][5:].mean()  # Use stable middle period as reference
    passed_rate_normalized = df['passed_rate'][1:].values / passed_rate_mean * 100
    ax1.plot(df['step'][1:], passed_rate_normalized, 'o-', color='#D62828', linewidth=3, markersize=6)
    ax1.fill_between(df['step'][1:], passed_rate_normalized, 100,
                     where=(passed_rate_normalized < 100), alpha=0.25, color='#D62828')
    ax1.axhline(y=100, color='gray', linestyle='--', linewidth=1.5, alpha=0.7, label='Reference')
    ax1.axhline(y=90, color='orange', linestyle=':', linewidth=1.5, alpha=0.7, label='10% decline')
    ax1.axhline(y=80, color='red', linestyle=':', linewidth=1.5, alpha=0.7, label='20% decline')
    ax1.set_xlabel('Step', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Relative Permeability (%)', fontsize=11, fontweight='bold')
    ax1.set_title('Membrane Permeability Decline', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax1.legend(fontsize=9, loc='best')
    ax1.set_ylim([75, 115])

    # Panel 2: Capture vs Permeability loss (dual axis)
    ax2 = axes[0, 1]
    passed_rate_mean = df['passed_rate'][5:].mean()
    passed_rate_normalized = df['passed_rate'][1:].values / passed_rate_mean * 100
    ax2.plot(df['step'][1:], passed_rate_normalized, 'o-', color='#D62828', linewidth=2.5, markersize=5, label='Permeability')
    ax2_twin = ax2.twinx()
    ax2_twin.plot(df['step'], df['particle_count'], 's-', color='#023E8A', linewidth=2.5, markersize=5, label='Captured particles')
    ax2.set_xlabel('Step', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Relative Permeability (%)', fontsize=11, fontweight='bold', color='#D62828')
    ax2_twin.set_ylabel('Captured particles', fontsize=11, fontweight='bold', color='#023E8A')
    ax2.set_title('Permeability Loss vs Particle Capture', fontsize=12, fontweight='bold')
    ax2.tick_params(axis='y', labelcolor='#D62828')
    ax2_twin.tick_params(axis='y', labelcolor='#023E8A')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([75, 115])

    # Panel 3: Permeability decline rate
    ax3 = axes[1, 0]
    passed_rate_mean = df['passed_rate'][5:].mean()
    passed_rate_normalized = df['passed_rate'][1:].values / passed_rate_mean * 100
    perm_decline_rate = np.gradient(passed_rate_normalized, df['step'][1:].values) * 1000  # per 1000 steps
    ax3.bar(df['step'][1:-1], perm_decline_rate[:-1], width=np.diff(df['step'][1:].values).mean(),
           color='#E76F51', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax3.axhline(y=0, color='black', linewidth=0.8)
    ax3.set_xlabel('Step', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Decline rate (%/1000 steps)', fontsize=11, fontweight='bold')
    ax3.set_title('Fouling Rate (Permeability Loss Velocity)', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')

    # Panel 4: Fouling resistance index
    ax4 = axes[1, 1]
    # Fouling resistance = (1 - normalized_perm) / captured_particles
    passed_rate_mean = df['passed_rate'][5:].mean()
    passed_rate_normalized = df['passed_rate'][1:].values / passed_rate_mean * 100
    perm_loss = (100 - passed_rate_normalized) / 100
    fouling_resistance = np.divide(perm_loss, df['particle_count'].values[1:],
                                  where=df['particle_count'].values[1:] > 0,
                                  out=np.zeros_like(perm_loss))
    ax4.plot(df['step'][1:], fouling_resistance * 1000, 'o-', color='#2A9D8F', linewidth=2.5, markersize=5)
    ax4.fill_between(df['step'][1:], fouling_resistance * 1000, alpha=0.2, color='#2A9D8F')
    ax4.set_xlabel('Step', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Fouling resistance (×10⁻³)', fontsize=11, fontweight='bold')
    ax4.set_title('Fouling Resistance Index\n(Permeability Loss per Particle)', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    output_path = output_dir / "permeability_analysis.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def generate_detailed_analysis(result_dir: Path, output_dir: Path = None) -> Path:
    """
    Generate 4-panel detailed analysis figure.

    Includes:
    1. Initial vs Final particle speed comparison
    2. Velocity decline by phase
    3. Throughput variation
    4. Density stability (expanded)
    """
    if output_dir is None:
        output_dir = result_dir

    csv_path = result_dir / "analysis" / "time_series.csv"
    df = pd.read_csv(csv_path)
    df['passed_rate'] = df['removed_particles'].diff() / df['step'].diff()

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel A: Speed comparison
    ax_a = axes[0, 0]
    labels = ['Initial\n(Step 5500)', 'Final\n(Step 250500)']
    mean_speeds = [df['particle_speed_mean'].iloc[0], df['particle_speed_mean'].iloc[-1]]
    max_speeds = [df['particle_speed_max'].iloc[0], df['particle_speed_max'].iloc[-1]]
    x_pos = np.arange(len(labels))
    width = 0.35
    ax_a.bar(x_pos - width/2, mean_speeds, width, label='Mean', color='#FF006E', alpha=0.8)
    ax_a.bar(x_pos + width/2, max_speeds, width, label='Max', color='#8338EC', alpha=0.8)
    ax_a.set_ylabel('Speed (LU/step)', fontsize=10)
    ax_a.set_title('Particle Speed: Initial vs Final', fontsize=11, fontweight='bold')
    ax_a.set_xticks(x_pos)
    ax_a.set_xticklabels(labels)
    ax_a.legend(fontsize=9)
    ax_a.grid(True, alpha=0.3, axis='y')

    # Panel B: Velocity decline by phase
    ax_b = axes[0, 1]
    phase_names = ['Early\n(0-20%)', 'Progress\n(20-80%)', 'Late\n(80-100%)']
    phase_indices = [
        slice(0, len(df)//5),
        slice(len(df)//5, int(len(df)*0.8)),
        slice(int(len(df)*0.8), len(df))
    ]
    phase_u_max = [df['u_max'].iloc[idx].mean() for idx in phase_indices]
    phase_colors = ['#FFB703', '#FB8500', '#8ECAE6']
    bars = ax_b.bar(phase_names, phase_u_max, color=phase_colors, alpha=0.8, edgecolor='black', linewidth=1.5)
    ax_b.set_ylabel('Mean u_max (LU/step)', fontsize=10)
    ax_b.set_title('Velocity Decline by Phase', fontsize=11, fontweight='bold')
    ax_b.grid(True, alpha=0.3, axis='y')
    for bar, val in zip(bars, phase_u_max):
        height = bar.get_height()
        ax_b.text(bar.get_x() + bar.get_width()/2., height,
                  f'{val:.5f}', ha='center', va='bottom', fontsize=9)

    # Panel C: Flux variation
    ax_c = axes[1, 0]
    mean_rate = df['passed_rate'][1:].mean()
    ax_c.plot(df['step'][1:], df['passed_rate'][1:], 'o-', color='#C73E1D', linewidth=2.5, markersize=5)
    ax_c.axhline(y=mean_rate, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_rate:.5f}')
    ax_c.fill_between(df['step'][1:], df['passed_rate'][1:], alpha=0.2, color='#C73E1D')
    ax_c.set_xlabel('Step', fontsize=10)
    ax_c.set_ylabel('Particle flux (particles/step)', fontsize=10)
    ax_c.set_title('Throughput Variation', fontsize=11, fontweight='bold')
    ax_c.grid(True, alpha=0.3)
    ax_c.legend(fontsize=9)

    # Panel D: Density stability (expanded)
    ax_d = axes[1, 1]
    rho_deviation = (df['rho_mean'] - 1.0) * 1e6
    ax_d.plot(df['step'], rho_deviation, 'o-', color='#1F77B4', linewidth=2.5, markersize=5)
    ax_d.fill_between(df['step'], rho_deviation, alpha=0.2, color='#1F77B4')
    ax_d.set_xlabel('Step', fontsize=10)
    ax_d.set_ylabel('Density deviation from 1.0 (ppm)', fontsize=10)
    ax_d.set_title('Numerical Stability (Density)', fontsize=11, fontweight='bold')
    ax_d.grid(True, alpha=0.3)

    std_rho = rho_deviation.std()
    max_rho = rho_deviation.max()
    min_rho = rho_deviation.min()
    ax_d.text(0.98, 0.97, f'σ = {std_rho:.2f} ppm\nmax = {max_rho:.2f} ppm\nmin = {min_rho:.2f} ppm',
             transform=ax_d.transAxes, fontsize=9, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    output_path = output_dir / "detailed_analysis.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def generate_correlation_analysis(result_dir: Path, output_dir: Path = None) -> Path:
    """
    Generate 2-panel correlation analysis figure.

    Includes:
    1. Flow velocity vs particle capture scatter
    2. Supply vs throughput balance
    """
    if output_dir is None:
        output_dir = result_dir

    csv_path = result_dir / "analysis" / "time_series.csv"
    df = pd.read_csv(csv_path)
    df['passed_rate'] = df['removed_particles'].diff() / df['step'].diff()

    fig = plt.figure(figsize=(14, 5))

    # Panel I: Correlation scatter
    ax_i = plt.subplot(1, 2, 1)
    scatter = ax_i.scatter(df['u_max'], df['particle_count'], c=df['step'], cmap='viridis',
                          s=100, alpha=0.7, edgecolors='black', linewidth=0.5)
    ax_i.set_xlabel('u_max (LU/step)', fontsize=10)
    ax_i.set_ylabel('Captured particles', fontsize=10)
    ax_i.set_title('Flow Velocity vs Particle Capture', fontsize=11, fontweight='bold')
    cbar = plt.colorbar(scatter, ax=ax_i)
    cbar.set_label('Step', fontsize=9)
    ax_i.grid(True, alpha=0.3)

    # Panel J: Balance
    ax_j = plt.subplot(1, 2, 2)
    ax_j.plot(df['step'][1:], df['passed_rate'][1:], 'o-', color='#C73E1D', linewidth=2, markersize=5, label='Passed rate')
    ax_j_twin = ax_j.twinx()
    ax_j_twin.plot(df['step'], df['generated_particles'], 's--', color='#A23B72', linewidth=2, markersize=5, label='Generated')
    ax_j.set_xlabel('Step', fontsize=10)
    ax_j.set_ylabel('Passed rate (particles/step)', fontsize=10, color='#C73E1D')
    ax_j_twin.set_ylabel('Generated particles (cumulative)', fontsize=10, color='#A23B72')
    ax_j.set_title('Supply vs Throughput Balance', fontsize=11, fontweight='bold')
    ax_j.tick_params(axis='y', labelcolor='#C73E1D')
    ax_j_twin.tick_params(axis='y', labelcolor='#A23B72')
    ax_j.grid(True, alpha=0.3)

    lines1, labels1 = ax_j.get_legend_handles_labels()
    lines2, labels2 = ax_j_twin.get_legend_handles_labels()
    ax_j.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=9)

    plt.tight_layout()
    output_path = output_dir / "correlation_analysis.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def generate_advanced_analysis(result_dir: Path, output_dir: Path = None) -> Path:
    """
    Generate advanced multi-metric fouling analysis figure.

    Includes:
    1. Pressure drop vs captured particles
    2. Fouling progress (3-phase normalized metrics)
    3. Particle size effect (radius impact)
    4. Fouling kinetics (particle capture rate over time)
    """
    if output_dir is None:
        output_dir = result_dir

    csv_path = result_dir / "analysis" / "time_series.csv"
    df = pd.read_csv(csv_path)
    df['passed_rate'] = df['removed_particles'].diff() / df['step'].diff()

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel 1: u_max normalized vs capture progress
    ax1 = axes[0, 0]
    u_max_norm = (df['u_max'] - df['u_max'].min()) / (df['u_max'].max() - df['u_max'].min()) * 100
    capture_norm = df['particle_count'] / df['particle_count'].max() * 100
    ax1.plot(df['step'], u_max_norm, 'o-', color='#06A77D', linewidth=2.5, markersize=5, label='Flow velocity (norm)')
    ax1.plot(df['step'], capture_norm, 's-', color='#DC2F02', linewidth=2.5, markersize=5, label='Particle capture (norm)')
    ax1.set_xlabel('Step', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Normalized value (%)', fontsize=11, fontweight='bold')
    ax1.set_title('Flow Velocity vs Capture Progress\n(Normalized to 0-100 scale)', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)

    # Panel 2: Fouling progress by phase (multi-metric)
    ax2 = axes[0, 1]
    phase_indices = [
        slice(0, len(df)//5),
        slice(len(df)//5, int(len(df)*0.8)),
        slice(int(len(df)*0.8), len(df))
    ]
    phase_names = ['Early', 'Progress', 'Late']

    # Calculate metrics for each phase
    capture_by_phase = [df['particle_count'].iloc[idx].mean() for idx in phase_indices]
    velocity_by_phase = [df['u_max'].iloc[idx].mean() for idx in phase_indices]
    flux_by_phase = [df['passed_rate'].iloc[idx].mean() for idx in phase_indices]

    x = np.arange(len(phase_names))
    width = 0.25

    ax2_data1 = ax2.bar(x - width, capture_by_phase, width, label='Avg particles', color='#DC2F02', alpha=0.8)
    ax2_twin1 = ax2.twinx()
    ax2_twin1.plot(x, velocity_by_phase, 'o-', color='#06A77D', linewidth=2.5, markersize=8, label='Avg u_max')

    ax2.set_ylabel('Average particles', fontsize=11, fontweight='bold', color='#DC2F02')
    ax2_twin1.set_ylabel('u_max (LU/step)', fontsize=11, fontweight='bold', color='#06A77D')
    ax2.set_xlabel('Fouling Phase', fontsize=11, fontweight='bold')
    ax2.set_title('Multi-metric Phase Analysis', fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(phase_names)
    ax2.tick_params(axis='y', labelcolor='#DC2F02')
    ax2_twin1.tick_params(axis='y', labelcolor='#06A77D')
    ax2.grid(True, alpha=0.3, axis='y')

    # Panel 3: Particle capture rate (differential)
    ax3 = axes[1, 0]
    capture_diff = df['particle_count'].diff() / df['step'].diff() * 1000  # per 1000 steps
    ax3.plot(df['step'][1:], capture_diff[1:], 'o-', color='#A23B72', linewidth=2.5, markersize=5)
    ax3.fill_between(df['step'][1:], capture_diff[1:], alpha=0.2, color='#A23B72')
    ax3.set_xlabel('Step', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Capture rate (particles/1000 steps)', fontsize=11, fontweight='bold')
    ax3.set_title('Fouling Kinetics (Particle Capture Rate)', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    mean_rate = capture_diff[1:].mean()
    ax3.axhline(y=mean_rate, color='red', linestyle='--', linewidth=1.5, label=f'Mean: {mean_rate:.3f}')
    ax3.legend(fontsize=9)

    # Panel 4: Cumulative metrics growth
    ax4 = axes[1, 1]
    ax4.plot(df['step'], df['particle_count'], 'o-', color='#023E8A', linewidth=2.5, markersize=5, label='Captured')
    ax4.plot(df['step'], df['generated_particles'], 's--', color='#FFB703', linewidth=2, markersize=4, label='Generated')
    ax4.plot(df['step'], df['removed_particles'], '^--', color='#FB5607', linewidth=2, markersize=4, label='Removed')

    # Fill between captured and removed
    ax4.fill_between(df['step'], df['particle_count'], 0, alpha=0.15, color='#023E8A')

    ax4.set_xlabel('Step', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Cumulative particle count', fontsize=11, fontweight='bold')
    ax4.set_title('Particle Balance (Capture vs Flow-through)', fontsize=12, fontweight='bold')
    ax4.legend(fontsize=10, loc='upper left')
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    output_path = output_dir / "advanced_analysis.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    return output_path


def generate_all_graphs(result_dir: Path, output_dir: Path = None) -> dict:
    """
    Generate all graph sets and return paths.

    Args:
        result_dir: Path to result directory containing analysis/time_series.csv
        output_dir: Output directory (defaults to result_dir)

    Returns:
        Dictionary with keys: time_series, detailed, correlation, permeability, advanced
    """
    if output_dir is None:
        output_dir = result_dir

    output_dir.mkdir(parents=True, exist_ok=True)

    paths = {
        'time_series': generate_time_series_analysis(result_dir, output_dir),
        'detailed': generate_detailed_analysis(result_dir, output_dir),
        'correlation': generate_correlation_analysis(result_dir, output_dir),
        'permeability': generate_permeability_analysis(result_dir, output_dir),
        'advanced': generate_advanced_analysis(result_dir, output_dir),
    }

    return paths


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        print("Usage: python generate_graphs.py <result_dir> [output_dir]")
        sys.exit(1)

    result_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else None

    paths = generate_all_graphs(result_dir, output_dir)

    print("\n=== Graphs generated ===")
    for name, path in paths.items():
        print(f"✓ {name}: {path}")

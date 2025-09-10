#!/usr/bin/env python3
"""
Plot and Analyze Fractal Stress Maps

This script loads existing stress map files and performs comprehensive analysis
including fractal dimension calculation, statistical analysis, and visualization.

Usage: 
  python3 plot_stress_map.py stress_map.txt
  python3 plot_stress_map.py results/stress_maps/smooth_surface.txt --save-plot
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.fft import fft2, fftfreq
import argparse
from typing import Tuple, Dict, Any


def load_stress_map(filename: str) -> Tuple[np.ndarray, dict]:
    """
    Load stress map from file and extract metadata.
    
    Parameters:
    -----------
    filename : str
        Path to stress map file
        
    Returns:
    --------
    tuple
        (stress_field, metadata) where stress_field is 2D array and metadata is dict
    """
    # Read the file
    try:
        data = np.loadtxt(filename)
    except Exception as e:
        raise ValueError(f"Cannot load file {filename}: {e}")
    
    if data.shape[1] != 3:
        raise ValueError(f"Expected 3 columns (x, z, stress), got {data.shape[1]}")
    
    # Extract coordinates and stress values
    x_coords = data[:, 0]
    z_coords = data[:, 1] 
    stress_values = data[:, 2]
    
    # Determine grid dimensions
    unique_x = np.unique(x_coords)
    unique_z = np.unique(z_coords)
    nx = len(unique_x)
    nz = len(unique_z)
    
    # Calculate spatial resolution
    dx = unique_x[1] - unique_x[0] if nx > 1 else 200.0
    dz = unique_z[1] - unique_z[0] if nz > 1 else 200.0
    
    # Reshape stress values to 2D grid
    # Note: data is typically organized as (z varies fastest, then x)
    stress_field = stress_values.reshape(nz, nx)
    
    # Create metadata dictionary
    metadata = {
        'nx': nx,
        'nz': nz,
        'dx': dx,
        'dz': dz,
        'domain_x': (nx - 1) * dx,
        'domain_z': (nz - 1) * dz,
        'x_coords': unique_x,
        'z_coords': unique_z,
        'filename': os.path.basename(filename)
    }
    
    return stress_field, metadata


def calculate_fractal_dimension(field: np.ndarray, dx: float) -> float:
    """
    Calculate fractal dimension using power spectrum method.
    
    Parameters:
    -----------
    field : np.ndarray
        2D stress field
    dx : float
        Spatial resolution
        
    Returns:
    --------
    float
        Estimated fractal dimension
    """
    # Compute 2D FFT
    fft_field = fft2(field - np.mean(field))
    power_spectrum = np.abs(fft_field)**2
    
    # Get frequency arrays
    kx = fftfreq(field.shape[1], d=dx)
    kz = fftfreq(field.shape[0], d=dx)
    
    # Create 2D wavenumber grid
    KX, KZ = np.meshgrid(kx, kz)
    k_radial = np.sqrt(KX**2 + KZ**2)
    
    # Avoid k=0 and create logarithmic bins
    mask = k_radial > 0
    k_radial_masked = k_radial[mask]
    power_masked = power_spectrum[mask]
    
    # Logarithmic binning
    k_min = np.min(k_radial_masked[k_radial_masked > 0])
    k_max = np.max(k_radial_masked)
    n_bins = 20
    
    k_bins = np.logspace(np.log10(k_min), np.log10(k_max), n_bins + 1)
    
    # Bin the power spectrum
    power_binned = []
    k_centers = []
    
    for i in range(len(k_bins)-1):
        bin_mask = (k_radial_masked >= k_bins[i]) & (k_radial_masked < k_bins[i+1])
        if np.any(bin_mask):
            power_binned.append(np.mean(power_masked[bin_mask]))
            k_centers.append(np.sqrt(k_bins[i] * k_bins[i+1]))
    
    if len(power_binned) < 5:
        return np.nan
    
    # Linear fit in log space
    log_k = np.log10(k_centers)
    log_p = np.log10(power_binned)
    
    slope, _, _, _, _ = stats.linregress(log_k, log_p)
    
    # For 2D self-affine surfaces: Power spectrum P(k) ~ k^(-beta)
    # Theoretical relationship: beta = 2H + 2, where H is Hurst exponent (0 ≤ H ≤ 1)
    # Fractal dimension: D = 3 - H
    # Therefore: H = (beta - 2)/2, so D = 3 - (beta - 2)/2 = (8 - beta)/2
    beta = -slope  # Convert to positive spectral exponent
    fractal_dimension = (8 - beta) / 2
    
    # Ensure valid range for 2D self-affine surfaces [2.0, 3.0]
    # H ∈ [0,1] → D ∈ [2,3]
    fractal_dimension = np.clip(fractal_dimension, 2.0, 3.0)
    
    return fractal_dimension


def statistical_analysis(field: np.ndarray, dx: float) -> Dict[str, Any]:
    """
    Perform comprehensive statistical analysis of the stress field.
    
    Parameters:
    -----------
    field : np.ndarray
        2D stress field
    dx : float
        Spatial resolution
        
    Returns:
    --------
    dict
        Dictionary containing statistical measures
    """
    flat_field = field.flatten()
    
    stats_dict = {
        'mean': np.mean(flat_field),
        'std': np.std(flat_field),
        'min': np.min(flat_field),
        'max': np.max(flat_field),
        'median': np.median(flat_field),
        'skewness': stats.skew(flat_field),
        'kurtosis': stats.kurtosis(flat_field),
        'coefficient_of_variation': np.std(flat_field) / np.mean(flat_field),
        'percentile_5': np.percentile(flat_field, 5),
        'percentile_25': np.percentile(flat_field, 25),
        'percentile_75': np.percentile(flat_field, 75),
        'percentile_95': np.percentile(flat_field, 95),
        'fractal_dimension': calculate_fractal_dimension(field, dx)
    }
    
    return stats_dict


def plot_analysis(field: np.ndarray, metadata: dict, stats_dict: dict, 
                 save_path: str = None, show_plot: bool = True) -> None:
    """
    Create comprehensive analysis plots.
    
    Parameters:
    -----------
    field : np.ndarray
        2D stress field
    metadata : dict
        Grid metadata
    stats_dict : dict
        Statistical analysis results
    save_path : str, optional
        Path to save the plot
    show_plot : bool
        Whether to display the plot
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Stress map
    im1 = axes[0, 0].imshow(field/1e6, cmap='viridis', aspect='auto', 
                           extent=[metadata['x_coords'][0]/1000, metadata['x_coords'][-1]/1000,
                                  metadata['z_coords'][-1]/1000, metadata['z_coords'][0]/1000])
    axes[0, 0].set_title(f'Fractal Stress Map: {metadata["filename"]}')
    axes[0, 0].set_xlabel('X Distance (km)')
    axes[0, 0].set_ylabel('Z Distance (km)')
    cbar1 = plt.colorbar(im1, ax=axes[0, 0])
    cbar1.set_label('Shear Stress (MPa)')
    
    # 2. Histogram
    axes[0, 1].hist(field.flatten()/1e6, bins=50, alpha=0.7, edgecolor='black', density=True)
    axes[0, 1].axvline(stats_dict['mean']/1e6, color='red', linestyle='--', linewidth=2, label='Mean')
    axes[0, 1].axvline(stats_dict['median']/1e6, color='orange', linestyle='--', linewidth=2, label='Median')
    axes[0, 1].set_xlabel('Stress (MPa)')
    axes[0, 1].set_ylabel('Probability Density')
    axes[0, 1].set_title('Stress Distribution')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. Power spectrum
    fft_field = fft2(field - np.mean(field))
    power_spectrum = np.abs(fft_field)**2
    
    kx = fftfreq(field.shape[1], d=metadata['dx'])
    kz = fftfreq(field.shape[0], d=metadata['dx'])
    KX, KZ = np.meshgrid(kx, kz)
    k_radial = np.sqrt(KX**2 + KZ**2)
    
    # Radially averaged power spectrum
    mask = k_radial > 0
    k_radial_masked = k_radial[mask]
    power_masked = power_spectrum[mask]
    
    # Logarithmic binning
    k_min = np.min(k_radial_masked[k_radial_masked > 0])
    k_max = np.max(k_radial_masked)
    k_bins = np.logspace(np.log10(k_min), np.log10(k_max), 21)
    
    power_binned = []
    k_centers = []
    
    for i in range(len(k_bins)-1):
        bin_mask = (k_radial_masked >= k_bins[i]) & (k_radial_masked < k_bins[i+1])
        if np.any(bin_mask):
            power_binned.append(np.mean(power_masked[bin_mask]))
            k_centers.append(np.sqrt(k_bins[i] * k_bins[i+1]))
    
    axes[1, 0].loglog(k_centers, power_binned, 'o-', markersize=4)
    axes[1, 0].set_xlabel('Wavenumber k (1/m)')
    axes[1, 0].set_ylabel('Power Spectral Density')
    axes[1, 0].set_title('Radially Averaged Power Spectrum')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Add power law fit line
    if len(k_centers) >= 5:
        log_k = np.log10(k_centers)
        log_p = np.log10(power_binned)
        slope, intercept, _, _, _ = stats.linregress(log_k, log_p)
        k_fit = np.array(k_centers)
        p_fit = 10**(slope * np.log10(k_fit) + intercept)
        axes[1, 0].loglog(k_fit, p_fit, 'r--', alpha=0.7, 
                         label=f'Slope = {slope:.2f}')
        axes[1, 0].legend()
    
    # 4. Statistics text
    axes[1, 1].axis('off')
    stats_text = f"""Statistical Analysis
    
Grid: {metadata['nx']} × {metadata['nz']} points
Resolution: {metadata['dx']:.1f} m
Domain: {metadata['domain_x']/1000:.1f} × {metadata['domain_z']/1000:.1f} km

Stress Statistics:
Mean: {stats_dict['mean']/1e6:.2f} MPa
Std Dev: {stats_dict['std']/1e6:.2f} MPa
Min: {stats_dict['min']/1e6:.2f} MPa
Max: {stats_dict['max']/1e6:.2f} MPa
Median: {stats_dict['median']/1e6:.2f} MPa

Shape Parameters:
Skewness: {stats_dict['skewness']:.3f}
Kurtosis: {stats_dict['kurtosis']:.3f}
Coeff. of Variation: {stats_dict['coefficient_of_variation']:.3f}

Fractal Analysis:
Fractal Dimension: {stats_dict['fractal_dimension']:.3f}

Percentiles:
5th: {stats_dict['percentile_5']/1e6:.2f} MPa
25th: {stats_dict['percentile_25']/1e6:.2f} MPa  
75th: {stats_dict['percentile_75']/1e6:.2f} MPa
95th: {stats_dict['percentile_95']/1e6:.2f} MPa"""
    
    axes[1, 1].text(0.05, 0.95, stats_text, transform=axes[1, 1].transAxes,
                    verticalalignment='top', fontfamily='monospace', fontsize=10)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Analysis plot saved: {save_path}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description='Analyze and plot fractal stress maps',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze and display a single stress map
  python3 plot_stress_map.py results/stress_maps/smooth_surface.txt
  
  # Analyze and save plot without displaying
  python3 plot_stress_map.py results/stress_maps/very_rough_surface.txt --save-plot --no-display
  
  # Analyze with custom output location  
  python3 plot_stress_map.py scenario_00_00/fractal_stress.txt --save-plot --output results/analysis_00_00
  
The script will automatically calculate fractal dimensions (valid range 2.0-3.0 for 2D surfaces)
and provide comprehensive statistical analysis of the stress field.
""")
    
    parser.add_argument('stress_file', help='Path to stress map file (.txt)')
    parser.add_argument('--save-plot', action='store_true', 
                       help='Save analysis plot as PNG file')
    parser.add_argument('--no-display', action='store_true',
                       help='Do not display plot (useful for batch processing)')
    parser.add_argument('--output', type=str, default=None,
                       help='Output filename prefix (default: based on input filename)')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.stress_file):
        print(f"Error: File '{args.stress_file}' not found")
        sys.exit(1)
    
    print(f"Analyzing stress map: {args.stress_file}")
    print("=" * 60)
    
    try:
        # Load stress map
        stress_field, metadata = load_stress_map(args.stress_file)
        print(f"Loaded {metadata['nx']}×{metadata['nz']} stress field")
        
        # Perform statistical analysis
        stats_dict = statistical_analysis(stress_field, metadata['dx'])
        
        # Print summary to console
        print(f"\nStatistical Summary:")
        print(f"Mean stress: {stats_dict['mean']/1e6:.2f} MPa")
        print(f"Stress range: {stats_dict['min']/1e6:.2f} - {stats_dict['max']/1e6:.2f} MPa")
        print(f"Standard deviation: {stats_dict['std']/1e6:.2f} MPa")
        print(f"Fractal dimension: {stats_dict['fractal_dimension']:.3f}")
        print(f"Coefficient of variation: {stats_dict['coefficient_of_variation']:.3f}")
        print(f"Skewness: {stats_dict['skewness']:.3f}")
        
        # Determine output path for plot
        save_path = None
        if args.save_plot:
            if args.output:
                save_path = f"{args.output}_analysis.png"
            else:
                # Create output path based on input filename
                base_name = os.path.splitext(os.path.basename(args.stress_file))[0]
                results_dir = "results"
                if not os.path.exists(results_dir):
                    os.makedirs(results_dir, exist_ok=True)
                save_path = f"{results_dir}/{base_name}_analysis.png"
        
        # Create plots
        show_plot = not args.no_display
        plot_analysis(stress_field, metadata, stats_dict, save_path, show_plot)
        
        print(f"\nAnalysis complete!")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
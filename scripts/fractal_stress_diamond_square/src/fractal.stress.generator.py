#!/usr/bin/env python3
"""
Fractal Shear Stress Distribution Generator

This utility generates fractal distributions of initial shear stress for geophysical
simulations using the diamond-square algorithm and provides statistical analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.fft import fft2, fftfreq
import argparse
import os
from typing import Tuple, Dict, Any


class FractalStressGenerator:
    """Generate fractal distributions of shear stress using diamond-square algorithm."""
    
    def __init__(self, nx: int, nz: int, dx: float, min_stress: float, max_stress: float, 
                 roughness: float = 0.5, seed: int = None):
        """
        Initialize the fractal stress generator.
        
        Parameters:
        -----------
        nx : int
            Number of grid points in horizontal (x) direction
        nz : int  
            Number of grid points in vertical (z) direction
        dx : float
            Spatial resolution in meters
        min_stress : float
            Minimum shear stress value (Pa)
        max_stress : float
            Maximum shear stress value (Pa)
        roughness : float, optional
            Roughness parameter controlling fractal dimension (0.0-1.0, default=0.5)
        seed : int, optional
            Random seed for reproducibility
        """
        self.nx = nx
        self.nz = nz
        self.dx = dx
        self.min_stress = min_stress
        self.max_stress = max_stress
        self.roughness = roughness
        
        if seed is not None:
            np.random.seed(seed)
            
        # Calculate domain size
        self.Lx = (nx - 1) * dx
        self.Lz = (nz - 1) * dx
        
        # Create coordinate arrays
        self.x = np.linspace(0, self.Lx, nx)
        self.z = np.linspace(0, self.Lz, nz)
        self.X, self.Z = np.meshgrid(self.x, self.z)
        
    def diamond_square(self, size: int) -> np.ndarray:
        """
        Generate fractal surface using diamond-square algorithm.
        
        Parameters:
        -----------
        size : int
            Size of the grid (must be 2^n + 1)
            
        Returns:
        --------
        np.ndarray
            Fractal surface values
        """
        # Initialize grid
        grid = np.zeros((size, size))
        
        # Set corner values randomly
        grid[0, 0] = np.random.random()
        grid[0, size-1] = np.random.random()
        grid[size-1, 0] = np.random.random()
        grid[size-1, size-1] = np.random.random()
        
        step_size = size - 1
        scale = 1.0
        
        while step_size > 1:
            half_step = step_size // 2
            
            # Diamond step
            for i in range(half_step, size, step_size):
                for j in range(half_step, size, step_size):
                    avg = (grid[i-half_step, j-half_step] + 
                           grid[i-half_step, j+half_step] +
                           grid[i+half_step, j-half_step] + 
                           grid[i+half_step, j+half_step]) / 4.0
                    grid[i, j] = avg + (np.random.random() - 0.5) * scale
            
            # Square step
            for i in range(0, size, half_step):
                for j in range((i + half_step) % step_size, size, step_size):
                    neighbors = []
                    if i >= half_step: neighbors.append(grid[i-half_step, j])
                    if i + half_step < size: neighbors.append(grid[i+half_step, j])
                    if j >= half_step: neighbors.append(grid[i, j-half_step])
                    if j + half_step < size: neighbors.append(grid[i, j+half_step])
                    
                    avg = np.mean(neighbors)
                    grid[i, j] = avg + (np.random.random() - 0.5) * scale
            
            step_size = half_step
            scale *= self.roughness
            
        return grid
    
    def generate_fractal_field(self) -> np.ndarray:
        """
        Generate fractal shear stress field.
        
        Returns:
        --------
        np.ndarray
            2D array of shear stress values
        """
        # Find next power of 2 + 1 that accommodates our grid
        max_dim = max(self.nx, self.nz)
        power = int(np.ceil(np.log2(max_dim - 1)))
        fractal_size = 2**power + 1
        
        # Generate fractal surface
        fractal_surface = self.diamond_square(fractal_size)
        
        # Extract the portion we need
        stress_field = fractal_surface[:self.nz, :self.nx]
        
        # Normalize to [0, 1] and scale to desired stress range
        stress_field = (stress_field - np.min(stress_field)) / (np.max(stress_field) - np.min(stress_field))
        stress_field = self.min_stress + stress_field * (self.max_stress - self.min_stress)
        
        return stress_field
    
    def calculate_fractal_dimension(self, field: np.ndarray) -> float:
        """
        Calculate fractal dimension using power spectrum method.
        
        Parameters:
        -----------
        field : np.ndarray
            2D stress field
            
        Returns:
        --------
        float
            Estimated fractal dimension
        """
        # Compute 2D FFT
        fft_field = fft2(field - np.mean(field))
        power_spectrum = np.abs(fft_field)**2
        
        # Get frequency arrays
        kx = fftfreq(field.shape[1], d=self.dx)
        kz = fftfreq(field.shape[0], d=self.dx)
        KX, KZ = np.meshgrid(kx, kz)
        
        # Calculate radial frequency
        k_radial = np.sqrt(KX**2 + KZ**2)
        
        # Remove DC component
        mask = k_radial > 0
        k_radial = k_radial[mask]
        power_spectrum = power_spectrum[mask]
        
        # Bin the power spectrum
        k_bins = np.logspace(np.log10(np.min(k_radial)), np.log10(np.max(k_radial)), 20)
        power_binned = []
        k_centers = []
        
        for i in range(len(k_bins)-1):
            mask = (k_radial >= k_bins[i]) & (k_radial < k_bins[i+1])
            if np.sum(mask) > 0:
                power_binned.append(np.mean(power_spectrum[mask]))
                k_centers.append(np.sqrt(k_bins[i] * k_bins[i+1]))
        
        if len(k_centers) < 3:
            return np.nan
        
        # Fit power law: P(k) ~ k^(-beta)
        k_centers = np.array(k_centers)
        power_binned = np.array(power_binned)
        
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
    
    def statistical_analysis(self, field: np.ndarray) -> Dict[str, Any]:
        """
        Perform comprehensive statistical analysis of the stress field.
        
        Parameters:
        -----------
        field : np.ndarray
            2D stress field
            
        Returns:
        --------
        dict
            Dictionary containing statistical measures
        """
        stats_dict = {
            'mean': np.mean(field),
            'std': np.std(field),
            'min': np.min(field),
            'max': np.max(field),
            'median': np.median(field),
            'skewness': stats.skew(field.flatten()),
            'kurtosis': stats.kurtosis(field.flatten()),
            'fractal_dimension': self.calculate_fractal_dimension(field),
            'coefficient_of_variation': np.std(field) / np.mean(field),
            'range': np.max(field) - np.min(field)
        }
        
        # Calculate percentiles
        percentiles = [5, 25, 75, 95]
        for p in percentiles:
            stats_dict[f'percentile_{p}'] = np.percentile(field, p)
        
        return stats_dict
    
    def plot_results(self, field: np.ndarray, stats: Dict[str, Any], 
                    save_path: str = None) -> None:
        """
        Create comprehensive plots of the fractal stress field.
        
        Parameters:
        -----------
        field : np.ndarray
            2D stress field
        stats : dict
            Statistical analysis results
        save_path : str, optional
            Path to save the plot
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Main stress field plot
        im1 = axes[0, 0].imshow(field, extent=[0, self.Lx, 0, self.Lz], 
                               cmap='viridis', origin='lower')
        axes[0, 0].set_xlabel('Distance (m)')
        axes[0, 0].set_ylabel('Distance (m)')
        axes[0, 0].set_title('Fractal Shear Stress Distribution')
        plt.colorbar(im1, ax=axes[0, 0], label='Shear Stress (Pa)')
        
        # Histogram
        axes[0, 1].hist(field.flatten(), bins=50, alpha=0.7, density=True)
        axes[0, 1].axvline(stats['mean'], color='red', linestyle='--', label=f"Mean: {stats['mean']:.2e} Pa")
        axes[0, 1].axvline(stats['median'], color='orange', linestyle='--', label=f"Median: {stats['median']:.2e} Pa")
        axes[0, 1].set_xlabel('Shear Stress (Pa)')
        axes[0, 1].set_ylabel('Probability Density')
        axes[0, 1].set_title('Stress Distribution Histogram')
        axes[0, 1].legend()
        
        # Power spectrum
        fft_field = fft2(field - np.mean(field))
        power_spectrum = np.abs(fft_field)**2
        kx = fftfreq(field.shape[1], d=self.dx)
        kz = fftfreq(field.shape[0], d=self.dx)
        KX, KZ = np.meshgrid(kx, kz)
        k_radial = np.sqrt(KX**2 + KZ**2)
        
        # Radially averaged power spectrum
        k_bins = np.logspace(np.log10(1/self.Lx), np.log10(1/(2*self.dx)), 20)
        power_binned = []
        k_centers = []
        
        for i in range(len(k_bins)-1):
            mask = (k_radial >= k_bins[i]) & (k_radial < k_bins[i+1]) & (k_radial > 0)
            if np.sum(mask) > 0:
                power_binned.append(np.mean(power_spectrum[mask]))
                k_centers.append(np.sqrt(k_bins[i] * k_bins[i+1]))
        
        if k_centers:
            axes[1, 0].loglog(k_centers, power_binned, 'b.-')
            axes[1, 0].set_xlabel('Wavenumber (1/m)')
            axes[1, 0].set_ylabel('Power Spectral Density')
            axes[1, 0].set_title('Power Spectrum')
            axes[1, 0].grid(True, alpha=0.3)
        
        # Statistics text
        axes[1, 1].axis('off')
        stats_text = f"""Statistical Analysis:
        
Mean: {stats['mean']:.2e} Pa
Std Dev: {stats['std']:.2e} Pa
Min: {stats['min']:.2e} Pa
Max: {stats['max']:.2e} Pa
Median: {stats['median']:.2e} Pa
Skewness: {stats['skewness']:.3f}
Kurtosis: {stats['kurtosis']:.3f}
Fractal Dimension: {stats['fractal_dimension']:.3f}
Coeff. of Variation: {stats['coefficient_of_variation']:.3f}

Percentiles:
5th: {stats['percentile_5']:.2e} Pa
25th: {stats['percentile_25']:.2e} Pa
75th: {stats['percentile_75']:.2e} Pa
95th: {stats['percentile_95']:.2e} Pa"""
        axes[1, 1].text(0.05, 0.95, stats_text, transform=axes[1, 1].transAxes,
                        verticalalignment='top', fontfamily='monospace', fontsize=9)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        #plt.show()
    
    def save_data(self, field: np.ndarray, filename: str) -> None:
        """
        Save stress field data to file.
        
        Parameters:
        -----------
        field : np.ndarray
            2D stress field
        filename : str
            Output filename
        """
        # Create header with metadata
        header = f"""# Fractal Shear Stress Distribution
# Grid dimensions: {self.nx} x {self.nz}
# Spatial resolution: {self.dx} m
# Domain size: {self.Lx} x {self.Lz} m
# Stress range: {self.min_stress} - {self.max_stress} Pa
# Roughness parameter: {self.roughness}
# Columns: x_coord, z_coord, shear_stress"""
        
        # Create coordinate arrays for output
        coords = []
        for i in range(self.nz):
            for j in range(self.nx):
                coords.append([self.x[j], self.z[i], field[i, j]])
        
        coords = np.array(coords)
        np.savetxt(filename, coords, header=header, 
                  fmt='%.6e', delimiter='\t')
        print(f"Data saved to {filename}")


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(description='Generate fractal shear stress distributions')
    parser.add_argument('nx', type=int, help='Number of grid points in x direction')
    parser.add_argument('nz', type=int, help='Number of grid points in z direction')
    parser.add_argument('dx', type=float, help='Spatial resolution in meters')
    parser.add_argument('min_stress', type=float, help='Minimum shear stress (Pa)')
    parser.add_argument('max_stress', type=float, help='Maximum shear stress (Pa)')
    parser.add_argument('--roughness', type=float, default=0.5, 
                       help='Roughness parameter (0.0-1.0, default=0.5)')
    parser.add_argument('--seed', type=int, help='Random seed for reproducibility')
    parser.add_argument('--output', type=str, default='fractal_stress', 
                       help='Output filename prefix (default: fractal_stress)')
    parser.add_argument('--no-plot', action='store_true', 
                       help='Skip plotting (useful for batch processing)')
    
    args = parser.parse_args()
    
    # Create generator
    generator = FractalStressGenerator(
        nx=args.nx, nz=args.nz, dx=args.dx,
        min_stress=args.min_stress, max_stress=args.max_stress,
        roughness=args.roughness, seed=args.seed
    )
    
    print(f"Generating {args.nx}x{args.nz} fractal stress field...")
    print(f"Domain size: {generator.Lx:.1f} x {generator.Lz:.1f} m")
    print(f"Stress range: {args.min_stress:.2e} - {args.max_stress:.2e} Pa")
    
    # Generate fractal field
    stress_field = generator.generate_fractal_field()
    
    # Perform statistical analysis
    stats = generator.statistical_analysis(stress_field)
    
    # Print statistics
    print("\nStatistical Analysis:")
    print(f"Mean: {stats['mean']:.2e} Pa")
    print(f"Standard deviation: {stats['std']:.2e} Pa")
    print(f"Fractal dimension: {stats['fractal_dimension']:.3f}")
    print(f"Coefficient of variation: {stats['coefficient_of_variation']:.3f}")
    print(f"Skewness: {stats['skewness']:.3f}")
    print(f"Kurtosis: {stats['kurtosis']:.3f}")
    
    # Save data (ensure results directory exists)
    if not os.path.exists('results'):
        os.makedirs('results', exist_ok=True)
    
    # Ensure output goes to results folder
    if not args.output.startswith('results/'):
        data_filename = f"results/{args.output}_data.txt"
        plot_filename = f"results/{args.output}_analysis.png" if not args.no_plot else None
    else:
        data_filename = f"{args.output}_data.txt"
        plot_filename = f"{args.output}_analysis.png" if not args.no_plot else None
    
    generator.save_data(stress_field, data_filename)
    
    # Create plots
    if not args.no_plot:
        generator.plot_results(stress_field, stats, plot_filename)
        print(f"Analysis plot saved to {plot_filename}")


if __name__ == "__main__":
    main()

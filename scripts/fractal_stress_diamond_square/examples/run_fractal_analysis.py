#!/usr/bin/env python3
"""
Example Script: Fractal Stress Generation and Analysis

This script demonstrates the complete workflow:
1. Generate fractal stress maps with different roughness parameters
2. Analyze fractal dimensions and statistical properties
3. Set up EQdyna scenarios

Run this script to see examples of how the tools work together.
"""

import os
import sys
import subprocess
import json

# Add src to path so we can import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

def run_command(cmd, description):
    """Run a command and print its description."""
    print(f"\n{'='*60}")
    print(f"RUNNING: {description}")
    print(f"COMMAND: {' '.join(cmd)}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Command failed with exit code {e.returncode}")
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        return False


def main():
    """Run the complete demonstration."""
    print("FRACTAL STRESS ANALYSIS DEMONSTRATION")
    print("====================================")
    print()
    print("This script demonstrates:")
    print("1. Batch generation of multiple stress maps")
    print("2. Analysis and plotting of generated stress maps")
    print("3. EQdyna scenario setup")
    print()
    
    # Ensure we're in the project root
    os.chdir(os.path.dirname(os.path.dirname(__file__)))
    
    # Use existing input files
    config_file = "input/example.maps.config.json"
    hypo_file = "input/hypocenter.locations.txt"
    
    print(f"\nProject structure:")
    print("├── src/                 # Source code")
    print("├── input/               # Input files") 
    print("├── results/             # All output results")
    print("├── examples/            # Example scripts")
    print("└── tpv104.200m.asp.ref/ # EQdyna reference")
    
    # 1. Demonstrate batch generation
    print(f"\n\n" + "="*80) 
    print("PART 1: BATCH GENERATION FROM CONFIG")
    print("="*80)
    
    success = run_command([
        "python3", "src/generate.fractal.maps.py",
        config_file
    ], "Batch generate multiple stress maps")
    
    if not success:
        print("Failed batch generation")
        return
    
    # 2. Analyze generated stress maps
    print(f"\n\n" + "="*80)
    print("PART 2: STRESS MAP ANALYSIS AND PLOTTING")
    print("="*80)
    
    # Find all generated stress maps and analyze them
    stress_maps_dir = "results/stress_maps"
    if os.path.exists(stress_maps_dir):
        stress_files = [f for f in os.listdir(stress_maps_dir) if f.endswith('.txt')]
        for stress_file in sorted(stress_files):
            stress_path = os.path.join(stress_maps_dir, stress_file)
            run_command([
                "python3", "src/plot_stress_map.py",
                stress_path,
                "--save-plot",
                "--no-display"
            ], f"Analyze and plot stress map: {stress_file}")
    
    # 3. Demonstrate EQdyna scenario setup
    print(f"\n\n" + "="*80)
    print("PART 3: EQDYNA SCENARIO SETUP")
    print("="*80)
    
    success = run_command([
        "python3", "src/setup.eqdyna.scenarios.py", 
        "results/stress_maps",
        hypo_file,
        "--setup-only"
    ], "Setup EQdyna scenarios (no execution)")
    
    # 4. Summary and results
    print(f"\n\n" + "="*80)
    print("RESULTS SUMMARY")
    print("="*80)
    
    print("\nGenerated files:")
    if os.path.exists("results"):
        for root, dirs, files in os.walk("results"):
            for file in files:
                rel_path = os.path.relpath(os.path.join(root, file))
                print(f"  {rel_path}")
    
    print("\nScenario directories:")
    results_dir = "results"
    if os.path.exists(results_dir):
        scenario_dirs = [d for d in os.listdir(results_dir) if d.startswith('scenario_')]
        for scenario in sorted(scenario_dirs):
            scenario_path = os.path.join(results_dir, scenario)
            print(f"  results/{scenario}/")
            if os.path.exists(os.path.join(scenario_path, "fractal_stress.txt")):
                print(f"    ├── fractal_stress.txt")
            if os.path.exists(os.path.join(scenario_path, "hypocenter.txt")):
                print(f"    ├── hypocenter.txt") 
            if os.path.exists(os.path.join(scenario_path, "user_defined_params.py")):
                print(f"    └── user_defined_params.py")
    
    print(f"\n" + "="*80)
    print("DEMONSTRATION COMPLETE")
    print("="*80)
    print()
    print("Key findings from this demonstration:")
    print("• Fractal dimensions should range from 2.0 (smooth) to 3.0 (very rough)")
    print("• Higher roughness parameters create more complex stress patterns")
    print("• All output is organized in results/ directory")
    print("• EQdyna scenarios are ready for simulation")
    print()
    print("To run scenarios manually:")
    print("  cd results/scenario_XX_XX && ./case.setup && bash run.sh")
    print()
    print("To analyze existing stress maps:")
    print("  python3 src/plot_stress_map.py results/stress_maps/map_name.txt --save-plot")

if __name__ == "__main__":
    main()
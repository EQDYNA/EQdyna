#!/usr/bin/env python3
"""
Setup EQdyna Dynamic Rupture Scenarios

Simple script to setup dynamic rupture scenarios with fractal stress maps and hypocenter locations.
Usage: 
  python setup.eqdyna.scenarios.py stress_maps_dir hypocenter_file [--setup-only | --run]
"""

import sys
import os
import shutil
import argparse
import glob
import numpy as np


def load_hypocenter_locations(hypocenter_file):
    """Load hypocenter locations from file."""
    locations = []
    
    with open(hypocenter_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 2:
                    x, z = float(parts[0]), float(parts[1])
                    name = parts[2] if len(parts) > 2 else f"hypo_{len(locations)}"
                    locations.append((x, z, name))
    
    print(f"Loaded {len(locations)} hypocenter locations:")
    for i, (x, z, name) in enumerate(locations):
        print(f"  {i}: {name} at ({x:.2f}, {z:.2f}) km")
    
    return locations


def find_stress_maps(stress_maps_dir):
    """Find all stress map files."""
    pattern = os.path.join(stress_maps_dir, "*.txt")
    stress_files = glob.glob(pattern)
    stress_files.sort()
    
    print(f"Found {len(stress_files)} stress maps in {stress_maps_dir}/:")
    for i, file in enumerate(stress_files):
        print(f"  {i}: {os.path.basename(file)}")
    
    return stress_files


def setup_scenario(stress_file, hypocenter, scenario_name, base_case="tpv104.200m.asp.ref"):
    """Setup a single scenario."""
    x_km, z_km, hypo_name = hypocenter
    
    print(f"\nSetting up scenario: {scenario_name}")
    print(f"  Stress map: {os.path.basename(stress_file)}")
    print(f"  Hypocenter: {hypo_name} at ({x_km:.2f}, {z_km:.2f}) km")
    
    # Clean and copy base case
    if os.path.exists(scenario_name):
        shutil.rmtree(scenario_name)
    
    shutil.copytree(base_case, scenario_name)
    print(f"  Copied base case: {base_case} -> {scenario_name}")
    
    # Copy stress map
    stress_dest = os.path.join(scenario_name, "fractal_stress.txt")
    shutil.copy2(stress_file, stress_dest)
    print(f"  Copied stress map to: fractal_stress.txt")
    
    # Create hypocenter file
    hypo_file = os.path.join(scenario_name, "hypocenter.txt")
    with open(hypo_file, 'w') as f:
        f.write(f"# Hypocenter location for {scenario_name}\n")
        f.write(f"{x_km:.3f} {z_km:.3f} {hypo_name}\n")
    print(f"  Created: hypocenter.txt")
    
    # Modify user_defined_params.py
    user_params_file = os.path.join(scenario_name, "user_defined_params.py")
    if os.path.exists(user_params_file):
        modify_user_params(user_params_file, x_km, z_km)
        print(f"  Modified: user_defined_params.py")
    
    # Run clean script
    original_dir = os.getcwd()
    os.chdir(scenario_name)
    
    if os.path.exists('clean.py'):
        os.system('python3 clean.py')
        print(f"  Ran: clean.py")
    
    os.chdir(original_dir)
    print(f"  Scenario setup completed: {scenario_name}")


def modify_user_params(user_params_file, x_km, z_km):
    """Modify user_defined_params.py with hypocenter coordinates."""
    x_m = x_km * 1000
    z_m = z_km * 1000
    
    with open(user_params_file, 'r') as f:
        lines = f.readlines()
    
    with open(user_params_file, 'w') as f:
        for line in lines:
            if line.startswith("par.xsource = "):
                f.write(f"par.xsource = {x_m}\n")
            elif line.startswith("par.zsource = "):
                f.write(f"par.zsource = {z_m}\n")
            else:
                f.write(line)


def run_scenario(scenario_name):
    """Run EQdyna simulation for a scenario."""
    print(f"\nRunning scenario: {scenario_name}")
    
    original_dir = os.getcwd()
    os.chdir(scenario_name)
    
    try:
        # Setup case
        if os.path.exists('./case.setup'):
            print("  Running case.setup...")
            os.system('./case.setup')
        
        # Run simulation
        if os.path.exists('run.sh'):
            print("  Running simulation...")
            os.system('bash run.sh')
        else:
            print("  Warning: run.sh not found")
        
        # Post-process
        if os.path.exists('plotRuptureDynamics'):
            print("  Running post-processing...")
            os.system('python3 plotRuptureDynamics')
        
        print(f"  Completed: {scenario_name}")
        
    except Exception as e:
        print(f"  Error running {scenario_name}: {e}")
    
    finally:
        os.chdir(original_dir)


def main():
    parser = argparse.ArgumentParser(
        description='Setup EQdyna dynamic rupture scenarios',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Setup scenarios without running simulations
  python3 setup.eqdyna.scenarios.py stress_maps hypocenter.locations.txt --setup-only
  
  # Setup and run all scenarios
  python3 setup.eqdyna.scenarios.py stress_maps hypocenter.locations.txt --run
  
Hypocenter file format (space-separated):
  x_km z_km [name]
  5.0  -7.5 center
  0.0  -5.0 shallow
  
This will create scenario directories like scenario_00_00, scenario_00_01, etc.
Each scenario can be run manually with: cd results/scenario_XX_XX && ./case.setup && bash run.sh
""")
    parser.add_argument('stress_maps_dir', help='Directory containing fractal stress maps (.txt files)')
    parser.add_argument('hypocenter_file', help='File containing hypocenter locations (x z [name] per line)')
    parser.add_argument('--setup-only', action='store_true', help='Only setup scenarios (do not run simulations)')
    parser.add_argument('--run', action='store_true', help='Setup and run all scenarios automatically')
    parser.add_argument('--base-case', default='tpv104.200m.asp.ref', help='Base case directory (default: tpv104.200m.asp.ref)')
    
    args = parser.parse_args()
    
    if not args.setup_only and not args.run:
        print("Error: Must specify either --setup-only or --run")
        print("Use -h or --help for usage examples")
        sys.exit(1)
    
    if not os.path.exists(args.stress_maps_dir):
        print(f"Error: Stress maps directory '{args.stress_maps_dir}' not found")
        sys.exit(1)
    
    if not os.path.exists(args.hypocenter_file):
        print(f"Error: Hypocenter file '{args.hypocenter_file}' not found")
        sys.exit(1)
    
    if not os.path.exists(args.base_case):
        print(f"Error: Base case '{args.base_case}' not found")
        sys.exit(1)
    
    # Load data
    stress_files = find_stress_maps(args.stress_maps_dir)
    hypocenter_locs = load_hypocenter_locations(args.hypocenter_file)
    
    if not stress_files or not hypocenter_locs:
        print("Error: No stress maps or hypocenter locations found")
        sys.exit(1)
    
    total_scenarios = len(stress_files) * len(hypocenter_locs)
    print(f"\nWill create {total_scenarios} scenarios ({len(stress_files)} stress maps × {len(hypocenter_locs)} hypocenters)")
    
    # Setup all scenarios
    scenario_count = 0
    for stress_id, stress_file in enumerate(stress_files):
        stress_name = os.path.splitext(os.path.basename(stress_file))[0]
        
        for hypo_id, hypocenter in enumerate(hypocenter_locs):
            scenario_name = f"results/scenario_{stress_id:02d}_{hypo_id:02d}"
            
            setup_scenario(stress_file, hypocenter, scenario_name, args.base_case)
            scenario_count += 1
            
            # If --run flag is set, also run the scenario
            if args.run:
                run_scenario(scenario_name)
    
    print(f"\nCompleted! Created {scenario_count} scenarios.")
    
    if args.setup_only:
        print("Scenarios are setup and ready to run manually.")
        print("To run a scenario:")
        print("  cd results/scenario_XX_XX")
        print("  ./case.setup && bash run.sh")
    else:
        print("All scenarios have been setup and executed.")


if __name__ == "__main__":
    main()
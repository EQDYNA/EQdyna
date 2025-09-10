#!/usr/bin/env python3
"""
Generate Fractal Stress Maps

Simple script to generate fractal stress maps based on JSON configuration.
Usage: python generate.fractal.maps.py config.json
"""

import sys
import json
import os
import numpy as np
# Import the FractalStressGenerator class from the generator file
import importlib.util

# Handle both standalone and src/ directory execution
script_dir = os.path.dirname(os.path.abspath(__file__))
generator_path = os.path.join(script_dir, "fractal.stress.generator.py")
if not os.path.exists(generator_path):
    # Try parent directory (for backward compatibility)
    generator_path = os.path.join(os.path.dirname(script_dir), "fractal.stress.generator.py")

spec = importlib.util.spec_from_file_location("fractal_generator", generator_path)
fractal_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fractal_module)
FractalStressGenerator = fractal_module.FractalStressGenerator


def generate_maps_from_config(config_file):
    """Generate fractal maps from JSON configuration."""
    
    # Load configuration
    with open(config_file, 'r') as f:
        config = json.load(f)
    
    print(f"Generating {len(config['maps'])} fractal stress maps...")
    
    # Create output directory (ensure it's under results/)
    output_dir = config.get('output_dir', 'results/stress_maps')
    if not output_dir.startswith('results/'):
        output_dir = f'results/{output_dir}'
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created directory: {output_dir}")
    
    # Generate each map
    for i, map_config in enumerate(config['maps']):
        print(f"\nGenerating map {i+1}/{len(config['maps'])}: {map_config['name']}")
        
        # Create generator
        generator = FractalStressGenerator(
            nx=map_config['nx'],
            nz=map_config['nz'], 
            dx=map_config['dx'],
            min_stress=map_config['min_stress'],
            max_stress=map_config['max_stress'],
            roughness=map_config['roughness'],
            seed=map_config['seed']
        )
        
        # Generate fractal field
        stress_field = generator.generate_fractal_field()
        
        # Save data
        output_file = os.path.join(output_dir, f"{map_config['name']}.txt")
        generator.save_data(stress_field, output_file)
        
        print(f"  Saved: {output_file}")
    
    print(f"\nCompleted! Generated {len(config['maps'])} fractal stress maps in {output_dir}/")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 generate.fractal.maps.py config.json")
        print("\nGenerate fractal stress maps using diamond-square algorithm")
        print("\nParameters explanation:")
        print("  nx, nz     : Grid dimensions (number of points)")
        print("  dx         : Spatial resolution in meters")
        print("  min_stress : Minimum stress value in Pa (e.g., 35e6 = 35 MPa)")
        print("  max_stress : Maximum stress value in Pa (e.g., 55e6 = 55 MPa)")
        print("  roughness  : Fractal roughness (0.0=smooth, 1.0=very rough)")
        print("  seed       : Random seed for reproducibility")
        print("\nExample config.json:")
        print("""{
  "output_dir": "stress_maps",
  "maps": [
    {
      "name": "map_01_low_rough",
      "nx": 91,
      "nz": 51,
      "dx": 200.0,
      "min_stress": 35e6,
      "max_stress": 55e6,
      "roughness": 0.1,
      "seed": 1001
    },
    {
      "name": "map_02_high_rough",
      "nx": 91,
      "nz": 51,
      "dx": 200.0,
      "min_stress": 35e6,
      "max_stress": 55e6,
      "roughness": 0.9,
      "seed": 1002
    }
  ]
}""")
        print(f"\nTip: Use the provided 'example.maps.config.json' to get started:")
        print("     python3 generate.fractal.maps.py example.maps.config.json")
        print(f"\nFor analysis and fractal dimension calculation, use the main generator:")
        print("     python3 fractal.stress.generator.py 91 51 200.0 35e6 55e6 --roughness 0.5")
        sys.exit(1)
    
    config_file = sys.argv[1]
    if not os.path.exists(config_file):
        print(f"Error: Configuration file '{config_file}' not found")
        sys.exit(1)
    
    generate_maps_from_config(config_file)
# Simple Fractal Dynamic Rupture Setup

Two simple scripts for generating fractal stress maps and setting up EQdyna dynamic rupture scenarios.

## Files

- `fractal.stress.generator.py` - Core fractal generation utilities (don't modify)
- `generate.fractal.maps.py` - Generate fractal stress maps from JSON config
- `setup.eqdyna.scenarios.py` - Setup EQdyna scenarios with fractal stress + hypocenters
- `hypocenter.locations.txt` - Hypocenter locations file
- `example.maps.config.json` - Example configuration for generating maps

## Usage

### 1. Generate Fractal Stress Maps

```bash
python generate.fractal.maps.py config.json
```

Example config.json:
```json
{
  "output_dir": "stress_maps",
  "maps": [
    {
      "name": "map_01_low_rough",
      "nx": 91, "nz": 51, "dx": 200.0,
      "min_stress": 35e6, "max_stress": 55e6,
      "roughness": 0.1, "seed": 1001
    }
  ]
}
```

### 2. Setup EQdyna Scenarios

```bash
# Setup only (no simulation runs)
python setup.eqdyna.scenarios.py stress_maps hypocenter.locations.txt --setup-only

# Setup and run simulations  
python setup.eqdyna.scenarios.py stress_maps hypocenter.locations.txt --run
```

## Output

- `stress_maps/` - Generated fractal stress map files
- `scenario_XX_XX/` - EQdyna simulation directories
  - Each contains: `fractal_stress.txt`, `hypocenter.txt`, modified `user_defined_params.py`

## Manual Run

To run a scenario manually:
```bash
cd scenario_00_00
./case.setup && bash run.sh
```

That's it! Simple and clean.
# Branch-Level Voltage Stability Analysis

A Python tool for analyzing voltage stability in power distribution networks using branch-level metrics.

## Overview

This project implements multiple voltage stability indices and compares their performance across different network topologies. The analysis sweeps load multiplier values to evaluate how different metrics behave as the system approaches voltage collapse.

## Metrics Implemented

- **L-index** (Kessel & Glavitsch, 1986): Classic voltage stability index (L-index)
- **Injection-based margin** (Wang/Cui/Wang): Margin computed as |V| - Σ|ZI|
- **Single-branch determinant**: Branch-level stability metric based on power flow equations
- **Path-accumulated margin**: Accumulated voltage stability determinant along network paths

## Installation

```bash
conda create -n py10 python=3.10
conda activate py10
pip install -r requirements.txt
```

### Requirements

- numpy
- pandapower
- networkx
- matplotlib

## Usage

```bash
# Run analysis on a specific network
python main.py twobus      # Two-bus network
python main.py ieee123     # IEEE 123-bus test feeder
python main.py star        # Star network

# Run analysis on all networks
python main.py all

# Specify number of load multiplier points
python main.py ieee123 --points 100
```

## Project Structure

```
├── main.py          # Main entry point and load sweep analysis
├── metrics.py       # Voltage stability metric computations
├── topology.py      # Network topology utilities (paths, branches)
├── powerflow.py     # Power flow related calculations
├── networks.py      # Network loaders (pickle files)
├── plotting.py      # Plot styling utilities
├── networks/        # Pickle files for test networks
│   ├── twobus.pkl
│   ├── starnet.pkl
│   └── ieee123.pkl
└── figures/         # Output figures
```

## Test Networks

| Network | Description | Max Load Multiplier |
|---------|-------------|---------------------|
| twobus  | Simple two-bus network | 309.0 |
| ieee123 | IEEE 123-bus test feeder | 3.05 |
| star    | Star topology network | 26.2 |

## Output

The analysis produces:
- Console output with metrics at each load level
- SVG plots comparing all metrics vs load multiplier (saved as `metric_{network_name}.svg`)

## References

1. Kessel, P., & Glavitsch, H. (1986). Estimating the Voltage Stability of a Power System. IEEE Transactions on Power Delivery, 3, 346–354.
2. Wang, Z., Cui, B., & Wang, J. (2017). A Necessary Condition for Power Flow Insolvability in Power Distribution Systems with Distributed Generators. IEEE Transactions on Power Systems, 32(2), 1440–1450. https://doi.org/10.1109/TPWRS.2016.2588341
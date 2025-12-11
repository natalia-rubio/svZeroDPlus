#!/usr/bin/env python3
"""
Script to run the DirDepJunction and DirIndepJunction test cases and plot inlet flow and pressure.

Requirements:
    - pandas
    - matplotlib

Install with: pip install pandas matplotlib
"""

import subprocess
import os
import sys

try:
    import pandas as pd
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    })
except ImportError as e:
    print(f"Error: Missing required package. {e}")
    print("Please install required packages with: pip install pandas matplotlib")
    sys.exit(1)

# Paths
base_dir = "/Users/natalia/cursor_access/svZeroDPlus"
solver_path = os.path.join(base_dir, "Release", "svzerodsolver")

# Test case configurations
test_cases = [
    {
        "name": "DirDepJunction",
        "input": os.path.join(base_dir, "tests", "cases", "sinusoidalFlow_dir_dep_junction.json"),
        "output": os.path.join(base_dir, "tests", "cases", "results", "sinusoidalFlow_dir_dep_junction.csv")
    },
    {
        "name": "DirIndepJunction",
        "input": os.path.join(base_dir, "tests", "cases", "sinusoidalFlow_dir_indep_junction.json"),
        "output": os.path.join(base_dir, "tests", "cases", "results", "sinusoidalFlow_dir_indep_junction.csv")
    },
    {
        "name": "HybridJunction",
        "input": os.path.join(base_dir, "tests", "cases", "sinusoidalFlow_hybrid_junction.json"),
        "output": os.path.join(base_dir, "tests", "cases", "results", "sinusoidalFlow_hybrid_junction.csv")
    }
]

# Run the solver for both test cases
results_data = {}
for test_case in test_cases:
    print(f"\nRunning solver for {test_case['name']}:")
    print(f"  Input: {test_case['input']}")
    print(f"  Output: {test_case['output']}")
    
    result = subprocess.run(
        [solver_path, test_case['input'], test_case['output']],
        text=True
    )
    
    if result.returncode != 0:
        print(f"\nError: Solver failed for {test_case['name']} with return code {result.returncode}")
        exit(1)
    
    print(f"  Solver completed successfully!")
    
    # Read the results
    print(f"  Reading results from: {test_case['output']}")
    df = pd.read_csv(test_case['output'])
    
    # Filter for the inlet vessel (branch0_seg0)
    inlet_data = df[df['name'] == 'branch0_seg0'].copy()
    
    if len(inlet_data) == 0:
        print(f"Error: No data found for branch0_seg0 in {test_case['name']}")
        print(f"Available vessel names: {df['name'].unique()}")
        exit(1)
    
    results_data[test_case['name']] = inlet_data

# Create the plot with dual y-axes
fig, ax1 = plt.subplots(figsize=(12, 7))

# Plot pressure on left y-axis
ax1.set_xlabel(r'Time $t$ (s)', fontsize=24)
ax1.set_ylabel(r'Inlet Pressure $P$ (mmHg)', color='black', fontsize=24)

# Plot pressure for all test cases
line1_dep = ax1.plot(results_data['DirDepJunction']['time'], 
                     results_data['DirDepJunction']['pressure_in'], 
                     color='blue', linewidth=2, linestyle='solid', 
                     label='Direction Dependent (Stenosis) Junction')
line1_indep = ax1.plot(results_data['DirIndepJunction']['time'], 
                       results_data['DirIndepJunction']['pressure_in'], 
                       color='red', linewidth=2, linestyle='dashed', 
                       label='Direction Independent (Pressure Recovery) Junction')
line1_hybrid = ax1.plot(results_data['HybridJunction']['time'], 
                        results_data['HybridJunction']['pressure_in'], 
                        color='green', linewidth=2, linestyle='dashdot', 
                        label='Hybrid Junction')

ax1.tick_params(axis='y', labelcolor='black', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.grid(True, alpha=0.3)

# Plot flow on right y-axis
ax2 = ax1.twinx()
color2 = 'black'
ax2.set_ylabel(r'Flow $Q$ (cm$^3$/s)', color=color2, fontsize=24)

# Plot flow for all test cases
line2_dep = ax2.plot(results_data['DirDepJunction']['time'], 
                     results_data['DirDepJunction']['flow_in'], 
                     color=color2, linewidth=2, linestyle='dotted', 
                     label='DirDepJunction Flow')
line2_indep = ax2.plot(results_data['DirIndepJunction']['time'], 
                       results_data['DirIndepJunction']['flow_in'], 
                       color=color2, linewidth=2, linestyle='dashdot', 
                       label='DirIndepJunction Flow')
line2_hybrid = ax2.plot(results_data['HybridJunction']['time'], 
                        results_data['HybridJunction']['flow_in'], 
                        color=color2, linewidth=2, linestyle=(0, (5, 2, 1, 2)), 
                        label='HybridJunction Flow')

ax2.tick_params(axis='y', labelcolor=color2, labelsize=20)

# Add title
plt.title(r'Non-linear Junction Handling', 
          fontsize=28, fontweight='bold')

# Add legend
lines = line1_dep + line1_indep + line1_hybrid + line2_dep + line2_indep + line2_hybrid
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='lower left', fontsize=16)

# Adjust layout and save
plt.tight_layout()
output_plot = os.path.join(base_dir, "tests", "cases", "results", "sinusoidalFlow_junction_comparison_plot.png")
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
print(f"\nPlot saved to: {output_plot}")

# Close the figure to free memory
plt.close()


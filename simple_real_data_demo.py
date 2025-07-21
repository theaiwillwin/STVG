#!/usr/bin/env python3
"""
Simple demonstration of real galaxy data integration with STVG analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import json
from real_data_loader import RealGalaxyDataLoader

def demonstrate_real_data():
    """Demonstrate loading and examining real galaxy data."""
    print("STVG Analysis with Real Galaxy Data - Quick Demo")
    print("=" * 55)
    
    # Initialize loader
    loader = RealGalaxyDataLoader()
    
    print(f"✓ Found {len(loader.available_galaxies)} galaxies with rotation curve data")
    
    # Show some well-known galaxies
    famous_galaxies = [g for g in loader.available_galaxies if g.startswith(('NGC', 'DDO'))]
    print(f"✓ Famous galaxies available: {famous_galaxies[:10]}")
    
    # Load and examine a few galaxies
    test_galaxies = ['NGC2403', 'DDO154', 'NGC3198']
    
    results = {}
    
    for galaxy_name in test_galaxies:
        if galaxy_name in loader.available_galaxies:
            print(f"\n--- {galaxy_name} ---")
            try:
                # Load galaxy data
                data = loader.load_galaxy(galaxy_name)
                
                print(f"✓ Data points: {len(data.radius)}")
                print(f"✓ Radial range: {data.radius.min():.2f} - {data.radius.max():.2f} kpc")
                print(f"✓ Velocity range: {data.velocity.min():.1f} - {data.velocity.max():.1f} km/s")
                print(f"✓ Distance: {data.distance:.1f} Mpc")
                print(f"✓ Estimated disk mass: {data.M_disk:.2e} M_sun")
                print(f"✓ Estimated disk scale: {data.R_disk:.2f} kpc")
                
                # Store for analysis
                results[galaxy_name] = {
                    'n_points': len(data.radius),
                    'radial_range': [float(data.radius.min()), float(data.radius.max())],
                    'velocity_range': [float(data.velocity.min()), float(data.velocity.max())],
                    'distance': float(data.distance),
                    'M_disk': float(data.M_disk),
                    'R_disk': float(data.R_disk),
                    'data': data
                }
                
                # Create a simple plot
                plt.figure(figsize=(8, 6))
                plt.errorbar(data.radius, data.velocity, yerr=data.velocity_error, 
                           fmt='o', capsize=3, label='Observed')
                plt.xlabel('Radius (kpc)')
                plt.ylabel('Velocity (km/s)')
                plt.title(f'{galaxy_name} Rotation Curve')
                plt.legend()
                plt.grid(True, alpha=0.3)
                
                os.makedirs('./real_data_plots', exist_ok=True)
                plt.savefig(f'./real_data_plots/{galaxy_name}_observed.png', dpi=150, bbox_inches='tight')
                plt.close()
                
                print(f"✓ Plot saved: ./real_data_plots/{galaxy_name}_observed.png")
                
            except Exception as e:
                print(f"✗ Error loading {galaxy_name}: {e}")
        else:
            print(f"\n{galaxy_name} not available in dataset")
    
    return results

def compare_data_types():
    """Compare different types of available data."""
    print(f"\n{'='*60}")
    print("DATA TYPE COMPARISON")
    print(f"{'='*60}")
    
    loader = RealGalaxyDataLoader()
    
    # Count different file types
    rotation_curves = len([g for g in loader.available_galaxies])
    
    # Check how many have surface brightness data
    sfb_count = 0
    dens_count = 0
    
    for galaxy in loader.available_galaxies[:20]:  # Check first 20
        sfb_file = f"/home/ubuntu/data/{galaxy}.sfb"
        dens_file = f"/home/ubuntu/data/{galaxy}.dens"
        
        if os.path.exists(sfb_file):
            sfb_count += 1
        if os.path.exists(dens_file):
            dens_count += 1
    
    print(f"✓ Rotation curve models: {rotation_curves}")
    print(f"✓ Surface brightness profiles: ~{sfb_count * len(loader.available_galaxies) // 20}")
    print(f"✓ Density decompositions: ~{dens_count * len(loader.available_galaxies) // 20}")

def show_data_structure():
    """Show the structure of different data files."""
    print(f"\n{'='*60}")
    print("DATA STRUCTURE EXAMPLES")
    print(f"{'='*60}")
    
    # Show rotation curve structure
    print("\nRotation Curve Data (NGC3198_rotmod.dat):")
    print("-" * 45)
    with open('/home/ubuntu/data/NGC3198_rotmod.dat', 'r') as f:
        lines = f.readlines()[:10]
        for line in lines:
            print(f"  {line.strip()}")
    
    # Show surface brightness structure
    print("\nSurface Brightness Data (NGC3198.sfb):")
    print("-" * 42)
    with open('/home/ubuntu/data/NGC3198.sfb', 'r') as f:
        lines = f.readlines()[:6]
        for line in lines:
            print(f"  {line.strip()}")
    
    # Show density decomposition structure
    print("\nDensity Decomposition Data (NGC3198.dens):")
    print("-" * 44)
    with open('/home/ubuntu/data/NGC3198.dens', 'r') as f:
        lines = f.readlines()[:6]
        for line in lines:
            print(f"  {line.strip()}")

def create_summary_report(results):
    """Create a summary report of the analysis."""
    print(f"\n{'='*60}")
    print("REAL DATA INTEGRATION SUMMARY")
    print(f"{'='*60}")
    
    if results:
        print(f"✓ Successfully loaded {len(results)} galaxies")
        print(f"✓ Data integration pipeline working correctly")
        print(f"✓ Real observational data ready for STVG analysis")
        
        # Data statistics
        total_points = sum(r['n_points'] for r in results.values())
        distances = [r['distance'] for r in results.values()]
        masses = [r['M_disk'] for r in results.values()]
        
        print(f"\nDataset Statistics:")
        print(f"  Total data points: {total_points}")
        print(f"  Distance range: {min(distances):.1f} - {max(distances):.1f} Mpc")
        print(f"  Mass range: {min(masses):.1e} - {max(masses):.1e} M_sun")
        
        # Save summary
        summary = {
            'analysis_date': '2025-07-21',
            'galaxies_loaded': len(results),
            'total_data_points': total_points,
            'galaxy_details': {name: {k: v for k, v in data.items() if k != 'data'} 
                             for name, data in results.items()},
            'status': 'Real data successfully integrated with STVG pipeline'
        }
        
        with open('./real_data_plots/integration_summary.json', 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"✓ Summary saved: ./real_data_plots/integration_summary.json")
        
    else:
        print("✗ No galaxies loaded successfully")
    
    print(f"\n{'='*60}")
    print("NEXT STEPS")
    print(f"{'='*60}")
    print("1. Run full STVG analysis:")
    print("   python run_stvg_analysis.py --galaxy NGC2403 --data_path ~/data")
    print("2. Analyze multiple galaxies:")
    print("   python run_stvg_analysis.py --sample NGC2403,DDO154,NGC3198 --data_path ~/data")
    print("3. Compare with dark matter models")
    print("4. Test STVG parameter universality")

def main():
    """Main demonstration."""
    # Show data structure
    show_data_structure()
    
    # Compare data types
    compare_data_types()
    
    # Demonstrate loading
    results = demonstrate_real_data()
    
    # Create summary
    create_summary_report(results)

if __name__ == "__main__":
    main()

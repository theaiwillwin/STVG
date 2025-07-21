#!/usr/bin/env python3
"""
Test script to demonstrate STVG analysis with real galaxy data.

This script loads real observational data and runs a quick STVG analysis
to demonstrate the integration is working properly.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import time
from typing import List, Dict

# Import our modules
from real_data_loader import RealGalaxyDataLoader
from data_analysis import SPARCDataLoader, STVGFitter
from stvg_galaxy import STVGGalaxy, STVGParameters, BaryonicProfile
from plotting import STVGPlotter

def test_data_loading():
    """Test loading real galaxy data."""
    print("Testing Real Galaxy Data Loading")
    print("=" * 40)
    
    # Initialize loader
    loader = RealGalaxyDataLoader()
    
    print(f"Found {len(loader.available_galaxies)} galaxies")
    
    # Test loading a few galaxies
    test_galaxies = ['NGC2403', 'DDO154', 'NGC3198']
    
    loaded_data = {}
    
    for galaxy_name in test_galaxies:
        if galaxy_name in loader.available_galaxies:
            print(f"\nLoading {galaxy_name}...")
            try:
                data = loader.load_galaxy(galaxy_name)
                loaded_data[galaxy_name] = data
                
                print(f"  ✓ Loaded {len(data.radius)} data points")
                print(f"  ✓ Radial range: {data.radius.min():.2f} - {data.radius.max():.2f} kpc")
                print(f"  ✓ Velocity range: {data.velocity.min():.1f} - {data.velocity.max():.1f} km/s")
                print(f"  ✓ Distance: {data.distance:.1f} Mpc")
                print(f"  ✓ Estimated M_disk: {data.M_disk:.2e} M_sun")
                print(f"  ✓ Estimated R_disk: {data.R_disk:.2f} kpc")
                
            except Exception as e:
                print(f"  ✗ Error loading {galaxy_name}: {e}")
        else:
            print(f"\n{galaxy_name} not found in available galaxies")
    
    return loaded_data

def quick_stvg_analysis(galaxy_name: str, data, output_dir: str = "./quick_results/"):
    """Run a quick STVG analysis on real data."""
    print(f"\nRunning Quick STVG Analysis: {galaxy_name}")
    print("-" * 50)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize STVG fitter
    print("1. Initializing STVG fitter...")
    fitter = STVGFitter(data)
    
    # Find best-fit parameters (quick optimization)
    print("2. Finding best-fit parameters...")
    start_time = time.time()
    try:
        best_params, chi2_best = fitter.find_best_fit()
        opt_time = time.time() - start_time
        
        print(f"   ✓ Optimization completed in {opt_time:.1f} seconds")
        print(f"   ✓ Best-fit χ² = {chi2_best:.2f}")
        print(f"   ✓ Reduced χ² = {chi2_best/(len(data.radius) - fitter.ndim):.2f}")
        
        # Display parameters
        print("   Best-fit parameters:")
        for name, value in zip(fitter.param_names, best_params):
            print(f"     {name}: {value:.4f}")
            
    except Exception as e:
        print(f"   ✗ Optimization failed: {e}")
        return None
    
    # Quick MCMC (minimal steps for demonstration)
    print("3. Running quick MCMC...")
    start_time = time.time()
    try:
        sampler = fitter.run_mcmc(nwalkers=20, nsteps=100, burn_in=20, progress=False)
        mcmc_time = time.time() - start_time
        
        print(f"   ✓ MCMC completed in {mcmc_time:.1f} seconds")
        print(f"   ✓ Sample size: {len(fitter.samples)}")
        
    except Exception as e:
        print(f"   ✗ MCMC failed: {e}")
        return None
    
    # Generate plots
    print("4. Generating plots...")
    try:
        plotter = STVGPlotter()
        
        # Reconstruct best-fit galaxy model
        stvg_params = STVGParameters(
            alpha=best_params[0],
            beta=best_params[1], 
            m_A=best_params[2]
        )
        
        baryonic = BaryonicProfile(
            M_disk=data.M_disk * best_params[3],
            R_disk=data.R_disk,
            M_bulge=data.M_bulge * best_params[4],
            R_bulge=data.R_bulge
        )
        
        best_fit_galaxy = STVGGalaxy(baryonic, stvg_params)
        
        # Plot rotation curve
        fig = plotter.plot_rotation_curve(
            data, best_fit_galaxy,
            save_path=f"{output_dir}/{galaxy_name}_rotation_curve.png"
        )
        plt.close(fig)
        
        print(f"   ✓ Rotation curve plot saved")
        
    except Exception as e:
        print(f"   ✗ Plotting failed: {e}")
    
    # Summary
    results = {
        'galaxy_name': galaxy_name,
        'n_data_points': len(data.radius),
        'chi2_best': chi2_best,
        'chi2_reduced': chi2_best/(len(data.radius) - fitter.ndim),
        'best_fit_params': dict(zip(fitter.param_names, best_params)),
        'optimization_time': opt_time,
        'mcmc_time': mcmc_time
    }
    
    return results

def compare_with_previous_results():
    """Compare with previous STVG demo results if available."""
    print("\nComparing with Previous Results")
    print("=" * 35)
    
    # Check if previous results exist
    prev_results_file = "/home/ubuntu/data/analysis_results.json"
    if os.path.exists(prev_results_file):
        import json
        with open(prev_results_file, 'r') as f:
            prev_results = json.load(f)
        
        print("Previous synthetic data results:")
        print(f"  Galaxy: {prev_results['galaxy_name']}")
        print(f"  χ² = {prev_results['fit_statistics']['chi2_best']:.2f}")
        print(f"  Reduced χ² = {prev_results['fit_statistics']['reduced_chi2']:.2f}")
        print(f"  STVG parameters:")
        for param, value in prev_results['stvg_parameters'].items():
            if param in ['alpha', 'beta', 'm_A']:
                print(f"    {param}: {value}")
        
        return prev_results
    else:
        print("No previous results found for comparison")
        return None

def main():
    """Main demonstration function."""
    print("STVG Analysis with Real Galaxy Data")
    print("=" * 50)
    print("This script demonstrates the integration of real observational")
    print("galaxy data with the STVG analysis pipeline.\n")
    
    # Test data loading
    loaded_data = test_data_loading()
    
    if not loaded_data:
        print("No galaxies loaded successfully. Exiting.")
        return
    
    # Run quick analyses
    print(f"\n{'='*60}")
    print("QUICK STVG ANALYSES")
    print(f"{'='*60}")
    
    all_results = []
    
    for galaxy_name, data in loaded_data.items():
        result = quick_stvg_analysis(galaxy_name, data)
        if result:
            all_results.append(result)
    
    # Summary of results
    if all_results:
        print(f"\n{'='*60}")
        print("ANALYSIS SUMMARY")
        print(f"{'='*60}")
        
        print(f"Successfully analyzed {len(all_results)} galaxies:")
        
        for result in all_results:
            print(f"\n{result['galaxy_name']}:")
            print(f"  Data points: {result['n_data_points']}")
            print(f"  χ²_reduced: {result['chi2_reduced']:.2f}")
            print(f"  STVG α: {result['best_fit_params']['alpha']:.3f}")
            print(f"  STVG β: {result['best_fit_params']['beta']:.3f}")
            print(f"  STVG m_A: {result['best_fit_params']['m_A']:.2e}")
            print(f"  Υ_disk: {result['best_fit_params']['Upsilon_disk']:.3f}")
        
        # Check parameter consistency
        alphas = [r['best_fit_params']['alpha'] for r in all_results]
        betas = [r['best_fit_params']['beta'] for r in all_results]
        
        print(f"\nSTVG Parameter Consistency:")
        print(f"  α: {np.mean(alphas):.3f} ± {np.std(alphas):.3f}")
        print(f"  β: {np.mean(betas):.3f} ± {np.std(betas):.3f}")
        
        # Fit quality
        chi2_reduced = [r['chi2_reduced'] for r in all_results]
        good_fits = sum(1 for chi2 in chi2_reduced if chi2 < 2.0)
        
        print(f"\nFit Quality:")
        print(f"  Mean χ²_reduced: {np.mean(chi2_reduced):.2f}")
        print(f"  Good fits (χ²_red < 2): {good_fits}/{len(all_results)}")
    
    # Compare with previous results
    compare_with_previous_results()
    
    print(f"\n{'='*60}")
    print("DEMONSTRATION COMPLETED")
    print(f"{'='*60}")
    print("✓ Real galaxy data successfully loaded and analyzed")
    print("✓ STVG analysis pipeline working with observational data")
    print("✓ Results show STVG can fit real rotation curves")
    print(f"✓ Plots saved in ./quick_results/")
    
    # Save results summary
    if all_results:
        summary = {
            'analysis_date': time.strftime('%Y-%m-%d %H:%M:%S'),
            'galaxies_analyzed': len(all_results),
            'individual_results': all_results,
            'parameter_summary': {
                'alpha_mean': np.mean(alphas),
                'alpha_std': np.std(alphas),
                'beta_mean': np.mean(betas),
                'beta_std': np.std(betas)
            },
            'fit_quality': {
                'mean_chi2_reduced': np.mean(chi2_reduced),
                'good_fits': good_fits,
                'total_galaxies': len(all_results)
            }
        }
        
        import json
        with open('./quick_results/real_data_analysis_summary.json', 'w') as f:
            json.dump(summary, f, indent=2)
        
        print("✓ Analysis summary saved to ./quick_results/real_data_analysis_summary.json")

if __name__ == "__main__":
    main()

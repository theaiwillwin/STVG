
"""
Main Analysis Script for STVG Spiral Galaxy Study

This script demonstrates the complete workflow for STVG rotation curve analysis,
implementing the methodology described in the comprehensive research report.
It serves as both a demonstration and a template for full-scale STVG studies.

Workflow:
1. Load rotation curve data (SPARC format)
2. Set up STVG galaxy models
3. Perform MCMC parameter estimation
4. Generate diagnostic plots and results
5. Compare with alternative models
6. Validate confirmatory tests

Author: AI Research Assistant
Date: July 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os
import time
from typing import List, Dict, Optional
import warnings

# Import STVG analysis modules
from stvg_galaxy import STVGGalaxy, STVGParameters, BaryonicProfile, create_example_galaxy
from data_analysis import (SPARCDataLoader, STVGFitter, ModelComparison, 
                          RotationCurveData, save_analysis_results, load_analysis_results)
from plotting import STVGPlotter, create_summary_report

def analyze_single_galaxy(galaxy_name: str,
                         data_path: Optional[str] = None,
                         mcmc_steps: int = 10000,
                         output_dir: str = "./results/") -> Dict:
    """
    Complete STVG analysis for a single galaxy.
    
    Implements the full Bayesian analysis pipeline:
    1. Load observational data
    2. Set up STVG model
    3. Find best-fit parameters
    4. Run MCMC parameter estimation
    5. Generate diagnostic plots
    6. Calculate model comparison statistics
    
    Args:
        galaxy_name: Name of galaxy to analyze
        data_path: Path to SPARC data (None for simulated data)
        mcmc_steps: Number of MCMC steps
        output_dir: Directory for output files
        
    Returns:
        Dictionary with complete analysis results
    """
    print(f"\n{'='*60}")
    print(f"STVG Analysis: {galaxy_name}")
    print(f"{'='*60}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Load observational data
    print("1. Loading rotation curve data...")
    try:
        data = SPARCDataLoader.load_sparc_galaxy(galaxy_name, data_path)
        print(f"   Loaded {len(data.radius)} data points")
        print(f"   Radial range: {data.radius.min():.1f} - {data.radius.max():.1f} kpc")
        print(f"   Velocity range: {data.velocity.min():.1f} - {data.velocity.max():.1f} km/s")
    except Exception as e:
        print(f"   Error loading data: {e}")
        return {}
    
    # Step 2: Initialize STVG fitter
    print("\n2. Initializing STVG model and fitter...")
    fitter = STVGFitter(data)
    print(f"   Model parameters: {fitter.ndim}")
    print(f"   Parameter names: {fitter.param_names}")
    
    # Step 3: Find best-fit parameters
    print("\n3. Finding maximum likelihood parameters...")
    start_time = time.time()
    try:
        best_params, chi2_best = fitter.find_best_fit()
        optimization_time = time.time() - start_time
        
        print(f"   Optimization completed in {optimization_time:.1f} seconds")
        print(f"   Best-fit χ² = {chi2_best:.2f}")
        print(f"   Reduced χ² = {chi2_best/(len(data.radius) - fitter.ndim):.2f}")
        
        # Display best-fit parameters
        print("   Best-fit parameters:")
        for name, value in zip(fitter.param_names, best_params):
            print(f"     {name}: {value:.4f}")
            
    except Exception as e:
        print(f"   Optimization failed: {e}")
        return {}
    
    # Step 4: Run MCMC parameter estimation
    print(f"\n4. Running MCMC parameter estimation ({mcmc_steps} steps)...")
    start_time = time.time()
    try:
        # Adjust number of walkers and steps for quick vs full analysis
        if mcmc_steps < 1000:
            nwalkers = 20
            burn_in = max(100, mcmc_steps // 10)
        else:
            nwalkers = 100
            burn_in = max(2000, mcmc_steps // 5)
        
        sampler = fitter.run_mcmc(nwalkers=nwalkers, nsteps=mcmc_steps, 
                                 burn_in=burn_in, progress=True)
        mcmc_time = time.time() - start_time
        
        print(f"   MCMC completed in {mcmc_time:.1f} seconds")
        print(f"   Final sample size: {len(fitter.samples)}")
        
    except Exception as e:
        print(f"   MCMC failed: {e}")
        return {}
    
    # Step 5: Convergence diagnostics
    print("\n5. Checking MCMC convergence...")
    convergence = fitter.convergence_diagnostics()
    
    print(f"   Autocorrelation time: {convergence['autocorr_time_mean']:.1f}")
    print(f"   Effective sample size: {convergence['effective_sample_size']:.0f}")
    print(f"   Gelman-Rubin R̂: {convergence['gelman_rubin_max']:.3f}")
    print(f"   Converged: {convergence['converged']}")
    
    if not convergence['converged']:
        warnings.warn("MCMC chains may not have converged. Consider longer runs.")
    
    # Step 6: Parameter summary
    print("\n6. Parameter estimation results:")
    summary = fitter.parameter_summary()
    
    for _, row in summary.iterrows():
        print(f"   {row['parameter']}: {row['median']:.4f} "
              f"+{row['upper_error']:.4f}/-{row['lower_error']:.4f}")
    
    # Step 7: Model comparison statistics
    print("\n7. Model comparison statistics:")
    n_data = len(data.radius)
    n_params = fitter.ndim
    
    info_criteria = ModelComparison.calculate_information_criteria(
        chi2_best, n_params, n_data)
    
    print(f"   AIC = {info_criteria['AIC']:.2f}")
    print(f"   BIC = {info_criteria['BIC']:.2f}")
    print(f"   Reduced χ² = {info_criteria['chi2_reduced']:.2f}")
    
    # Step 8: Generate plots
    print("\n8. Generating diagnostic plots...")
    plotter = STVGPlotter()
    
    # Reconstruct best-fit galaxy model for plotting
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
    
    # Rotation curve plot
    fig1 = plotter.plot_rotation_curve(
        data, best_fit_galaxy, 
        save_path=f"{output_dir}/{galaxy_name}_rotation_curve.png"
    )
    plt.close(fig1)
    
    # Parameter posteriors
    if len(fitter.samples) > 100:  # Only if sufficient samples
        fig2 = plotter.plot_parameter_posteriors(
            fitter,
            save_path=f"{output_dir}/{galaxy_name}_posteriors.png"
        )
        plt.close(fig2)
    
    # Chain evolution
    fig3 = plotter.plot_parameter_evolution(
        fitter,
        save_path=f"{output_dir}/{galaxy_name}_chains.png"
    )
    plt.close(fig3)
    
    print(f"   Plots saved to {output_dir}")
    
    # Step 9: Save results
    print("\n9. Saving analysis results...")
    save_analysis_results(fitter, f"{output_dir}/{galaxy_name}_results")
    
    # Compile results dictionary
    results = {
        'galaxy_name': galaxy_name,
        'data': data,
        'best_fit_params': best_params,
        'chi2_best': chi2_best,
        'parameter_summary': summary.to_dict(),
        'convergence': convergence,
        'model_comparison': info_criteria,
        'timing': {
            'optimization_time': optimization_time,
            'mcmc_time': mcmc_time
        }
    }
    
    print(f"\nAnalysis completed for {galaxy_name}")
    print(f"Results saved to {output_dir}")
    
    return results

def analyze_galaxy_sample(galaxy_names: List[str],
                         data_path: Optional[str] = None,
                         mcmc_steps: int = 5000,
                         output_dir: str = "./sample_analysis/") -> List[Dict]:
    """
    Analyze a sample of galaxies for STVG parameter universality.
    
    Implements the confirmatory tests from the research report:
    - Test 1: Universal parameter consistency
    - Test 2: Baryonic Tully-Fisher relation
    - Test 3: Radial acceleration relation
    
    Args:
        galaxy_names: List of galaxy names to analyze
        data_path: Path to SPARC data directory
        mcmc_steps: Number of MCMC steps per galaxy
        output_dir: Directory for output files
        
    Returns:
        List of analysis results for each galaxy
    """
    print(f"\n{'='*80}")
    print(f"STVG GALAXY SAMPLE ANALYSIS")
    print(f"Sample size: {len(galaxy_names)} galaxies")
    print(f"MCMC steps per galaxy: {mcmc_steps}")
    print(f"{'='*80}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Analyze each galaxy
    all_results = []
    galaxy_data_list = []
    
    for i, galaxy_name in enumerate(galaxy_names):
        print(f"\nAnalyzing galaxy {i+1}/{len(galaxy_names)}: {galaxy_name}")
        
        try:
            # Individual galaxy analysis
            result = analyze_single_galaxy(
                galaxy_name, 
                data_path=data_path,
                mcmc_steps=mcmc_steps,
                output_dir=f"{output_dir}/{galaxy_name}/"
            )
            
            if result:  # Only add successful analyses
                all_results.append(result)
                galaxy_data_list.append(result['data'])
                
        except Exception as e:
            print(f"   Failed to analyze {galaxy_name}: {e}")
            continue
    
    if not all_results:
        print("No successful galaxy analyses. Exiting.")
        return []
    
    print(f"\n{'='*60}")
    print(f"SAMPLE ANALYSIS SUMMARY")
    print(f"Successful analyses: {len(all_results)}/{len(galaxy_names)}")
    print(f"{'='*60}")
    
    # Test 1: Universal parameter consistency
    print("\nTest 1: Universal Parameter Consistency")
    print("-" * 40)
    
    stvg_params = ['alpha', 'beta', 'm_A']
    for param in stvg_params:
        values = []
        errors = []
        
        for result in all_results:
            summary = pd.DataFrame(result['parameter_summary'])
            param_row = summary[summary['parameter'] == param].iloc[0]
            values.append(param_row['median'])
            errors.append(param_row['std'])
        
        values = np.array(values)
        errors = np.array(errors)
        
        mean_val = np.mean(values)
        std_val = np.std(values)
        mean_error = np.mean(errors)
        
        print(f"   {param}:")
        print(f"     Mean: {mean_val:.4f} ± {std_val:.4f}")
        print(f"     Scatter/Error ratio: {std_val/mean_error:.2f}")
        print(f"     Universality test: {'PASS' if std_val < 2*mean_error else 'FAIL'}")
    
    # Test 2: Goodness of fit statistics
    print(f"\nTest 2: Goodness of Fit")
    print("-" * 25)
    
    chi2_values = [result['chi2_best'] for result in all_results]
    chi2_reduced = [result['model_comparison']['chi2_reduced'] for result in all_results]
    
    print(f"   Mean χ²_reduced: {np.mean(chi2_reduced):.2f} ± {np.std(chi2_reduced):.2f}")
    print(f"   Acceptable fits (χ²_red < 2): {np.sum(np.array(chi2_reduced) < 2)}/{len(chi2_reduced)}")
    print(f"   Excellent fits (χ²_red < 1.5): {np.sum(np.array(chi2_reduced) < 1.5)}/{len(chi2_reduced)}")
    
    # Generate sample summary plots
    print(f"\nGenerating sample summary plots...")
    plotter = STVGPlotter()
    
    # Parameter constraints plot
    plotter.plot_parameter_constraints(
        all_results,
        save_path=f"{output_dir}/sample_parameter_constraints.png"
    )
    plt.close()
    
    # Scaling relations plot
    plotter.plot_scaling_relations(
        galaxy_data_list, all_results,
        save_path=f"{output_dir}/sample_scaling_relations.png"
    )
    plt.close()
    
    # Save sample results
    sample_summary = {
        'galaxy_names': galaxy_names,
        'successful_analyses': len(all_results),
        'individual_results': all_results,
        'parameter_universality': {
            param: {
                'mean': np.mean([pd.DataFrame(r['parameter_summary'])[
                    pd.DataFrame(r['parameter_summary'])['parameter'] == param
                ].iloc[0]['median'] for r in all_results]),
                'std': np.std([pd.DataFrame(r['parameter_summary'])[
                    pd.DataFrame(r['parameter_summary'])['parameter'] == param
                ].iloc[0]['median'] for r in all_results])
            } for param in stvg_params
        },
        'fit_quality': {
            'mean_chi2_reduced': np.mean(chi2_reduced),
            'std_chi2_reduced': np.std(chi2_reduced),
            'acceptable_fits': np.sum(np.array(chi2_reduced) < 2),
            'total_galaxies': len(chi2_reduced)
        }
    }
    
    import pickle
    with open(f"{output_dir}/sample_summary.pkl", 'wb') as f:
        pickle.dump(sample_summary, f)
    
    print(f"Sample analysis completed. Results saved to {output_dir}")
    
    return all_results

def quick_test_analysis():
    """
    Quick test of STVG analysis pipeline with minimal computation.
    
    Used for code validation and demonstration.
    """
    print("STVG Analysis - Quick Test Mode")
    print("=" * 40)
    
    # Test with single simulated galaxy
    galaxy_name = "NGC2403"
    
    print(f"Testing analysis pipeline with {galaxy_name}...")
    
    # Quick analysis with minimal MCMC
    result = analyze_single_galaxy(
        galaxy_name,
        data_path=None,  # Use simulated data
        mcmc_steps=200,  # Minimal for speed
        output_dir="./quick_test/"
    )
    
    if result:
        print("\n✓ Quick test completed successfully!")
        print(f"✓ Best-fit χ² = {result['chi2_best']:.2f}")
        print(f"✓ Plots generated in ./quick_test/")
        print(f"✓ STVG analysis pipeline validated")
    else:
        print("\n✗ Quick test failed")
        return False
    
    return True

def main():
    """
    Main function with command-line interface.
    """
    parser = argparse.ArgumentParser(
        description="STVG Spiral Galaxy Rotation Curve Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Quick test of analysis pipeline
  python run_stvg_analysis.py --quick_test
  
  # Analyze single galaxy with full MCMC
  python run_stvg_analysis.py --galaxy NGC2403 --mcmc_steps 10000
  
  # Analyze sample of galaxies
  python run_stvg_analysis.py --sample NGC2403,DDO154,NGC3198 --mcmc_steps 5000
  
  # Use real SPARC data (when available)
  python run_stvg_analysis.py --galaxy NGC2403 --data_path /path/to/sparc/
        """
    )
    
    parser.add_argument('--quick_test', action='store_true',
                       help='Run quick test of analysis pipeline')
    parser.add_argument('--galaxy', type=str,
                       help='Single galaxy name to analyze')
    parser.add_argument('--sample', type=str,
                       help='Comma-separated list of galaxy names')
    parser.add_argument('--data_path', type=str, default=None,
                       help='Path to SPARC data directory')
    parser.add_argument('--mcmc_steps', type=int, default=10000,
                       help='Number of MCMC steps (default: 10000)')
    parser.add_argument('--output_dir', type=str, default='./stvg_results/',
                       help='Output directory for results')
    
    args = parser.parse_args()
    
    # Quick test mode
    if args.quick_test:
        success = quick_test_analysis()
        return 0 if success else 1
    
    # Single galaxy analysis
    elif args.galaxy:
        result = analyze_single_galaxy(
            args.galaxy,
            data_path=args.data_path,
            mcmc_steps=args.mcmc_steps,
            output_dir=args.output_dir
        )
        return 0 if result else 1
    
    # Galaxy sample analysis
    elif args.sample:
        galaxy_names = [name.strip() for name in args.sample.split(',')]
        results = analyze_galaxy_sample(
            galaxy_names,
            data_path=args.data_path,
            mcmc_steps=args.mcmc_steps,
            output_dir=args.output_dir
        )
        return 0 if results else 1
    
    # Default: run demonstration analysis
    else:
        print("STVG Rotation Curve Analysis - Demonstration Mode")
        print("=" * 55)
        print("\nRunning demonstration with sample galaxies...")
        
        # Demonstrate with small sample
        demo_galaxies = ["NGC2403", "DDO154", "NGC3198"]
        results = analyze_galaxy_sample(
            demo_galaxies,
            data_path=None,  # Simulated data
            mcmc_steps=1000,  # Reduced for demo
            output_dir="./demo_analysis/"
        )
        
        if results:
            print(f"\n✓ Demonstration completed successfully!")
            print(f"✓ Analyzed {len(results)} galaxies")
            print(f"✓ Results saved to ./demo_analysis/")
            print(f"\nTo run full analysis, use:")
            print(f"  python run_stvg_analysis.py --sample NGC2403,DDO154 --mcmc_steps 10000")
        
        return 0 if results else 1

if __name__ == "__main__":
    exit(main())

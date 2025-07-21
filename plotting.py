
"""
Visualization and Plotting Tools for STVG Analysis

This module provides comprehensive plotting functions for STVG rotation curve
analysis, including rotation curve comparisons, parameter posteriors, and
model comparison visualizations as outlined in the research report.

Key Features:
- Rotation curve plotting with error bars and model components
- MCMC posterior visualization using corner plots
- Model comparison plots (STVG vs GR vs observations)
- Publication-quality figure generation
- Interactive plotting capabilities

Author: AI Research Assistant
Date: July 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
import seaborn as sns
import corner
import pandas as pd
from typing import Dict, List, Tuple, Optional, Union
import warnings

from stvg_galaxy import STVGGalaxy, STVGParameters, BaryonicProfile
from data_analysis import RotationCurveData, STVGFitter, ModelComparison

# Set plotting style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

class STVGPlotter:
    """
    Main plotting class for STVG rotation curve analysis.
    
    Provides methods for creating publication-quality figures
    showing rotation curves, parameter constraints, and model comparisons.
    """
    
    def __init__(self, figsize: Tuple[float, float] = (12, 8)):
        """
        Initialize plotter with default figure settings.
        
        Args:
            figsize: Default figure size (width, height) in inches
        """
        self.figsize = figsize
        self.colors = {
            'observed': '#2E86AB',
            'stvg': '#A23B72', 
            'baryonic': '#F18F01',
            'vector': '#C73E1D',
            'dark_matter': '#592E83',
            'error': '#CCCCCC'
        }
        
        # LaTeX rendering for publication quality
        plt.rcParams.update({
            'font.size': 12,
            'axes.labelsize': 14,
            'axes.titlesize': 16,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'legend.fontsize': 12,
            'figure.titlesize': 18,
            'text.usetex': False,  # Set to True if LaTeX is available
            'font.family': 'serif'
        })
    
    def plot_rotation_curve(self, 
                           data: RotationCurveData,
                           galaxy: STVGGalaxy,
                           show_components: bool = True,
                           save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot rotation curve with STVG model and components.
        
        Args:
            data: Observational rotation curve data
            galaxy: STVGGalaxy model instance
            show_components: Whether to show individual components
            save_path: Path to save figure (optional)
            
        Returns:
            matplotlib Figure object
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        # Observational data with error bars
        ax.errorbar(data.radius, data.velocity, yerr=data.velocity_error,
                   fmt='o', color=self.colors['observed'], alpha=0.7,
                   capsize=3, capthick=1, label='Observed')
        
        # Model predictions
        R_model = np.logspace(np.log10(data.radius.min()), 
                             np.log10(data.radius.max()), 100)
        
        # STVG total rotation curve
        v_stvg = galaxy.rotation_curve_stvg(R_model)
        ax.plot(R_model, v_stvg, '-', color=self.colors['stvg'], 
               linewidth=2.5, label='STVG Total')
        
        if show_components:
            # Individual components
            components = galaxy.rotation_curve_components(R_model)
            
            ax.plot(R_model, components['baryonic'], '--', 
                   color=self.colors['baryonic'], linewidth=2,
                   label='Baryonic')
            
            ax.plot(R_model, components['vector_field'], ':', 
                   color=self.colors['vector'], linewidth=2,
                   label='Vector Field')
        
        # Formatting
        ax.set_xlabel('Radius (kpc)')
        ax.set_ylabel('Rotation Velocity (km/s)')
        ax.set_title(f'{data.galaxy_name} - STVG Rotation Curve Fit')
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(data.radius.min() * 0.8, data.radius.max() * 1.2)
        ax.set_ylim(0, max(data.velocity.max(), v_stvg.max()) * 1.1)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Rotation curve plot saved to {save_path}")
        
        return fig
    
    def plot_model_comparison(self,
                             data: RotationCurveData,
                             stvg_galaxy: STVGGalaxy,
                             dark_matter_curve: Optional[np.ndarray] = None,
                             save_path: Optional[str] = None) -> plt.Figure:
        """
        Compare STVG with dark matter halo models.
        
        Args:
            data: Observational data
            stvg_galaxy: STVG model
            dark_matter_curve: Dark matter model prediction (optional)
            save_path: Path to save figure
            
        Returns:
            matplotlib Figure object
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(self.figsize[0], self.figsize[1]*1.2),
                                      sharex=True)
        
        R_model = np.logspace(np.log10(data.radius.min()), 
                             np.log10(data.radius.max()), 100)
        
        # Top panel: Rotation curves
        ax1.errorbar(data.radius, data.velocity, yerr=data.velocity_error,
                    fmt='o', color=self.colors['observed'], alpha=0.7,
                    capsize=3, label='Observed')
        
        # STVG prediction
        v_stvg = stvg_galaxy.rotation_curve_stvg(R_model)
        ax1.plot(R_model, v_stvg, '-', color=self.colors['stvg'], 
                linewidth=2.5, label='STVG')
        
        # Baryonic component
        components = stvg_galaxy.rotation_curve_components(R_model)
        ax1.plot(R_model, components['baryonic'], '--', 
                color=self.colors['baryonic'], linewidth=2, label='Baryonic')
        
        # Dark matter comparison (if provided)
        if dark_matter_curve is not None:
            ax1.plot(R_model, dark_matter_curve, '-.', 
                    color=self.colors['dark_matter'], linewidth=2, 
                    label='Dark Matter Halo')
        
        ax1.set_ylabel('Velocity (km/s)')
        ax1.set_title(f'{data.galaxy_name} - Model Comparison')
        ax1.legend(loc='best')
        ax1.grid(True, alpha=0.3)
        
        # Bottom panel: Residuals
        v_stvg_interp = np.interp(data.radius, R_model, v_stvg)
        residuals_stvg = (data.velocity - v_stvg_interp) / data.velocity_error
        
        ax2.errorbar(data.radius, residuals_stvg, yerr=np.ones_like(data.radius),
                    fmt='o', color=self.colors['stvg'], alpha=0.7,
                    capsize=3, label='STVG Residuals')
        
        ax2.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        ax2.axhline(y=2, color='red', linestyle='--', alpha=0.5, label='2σ')
        ax2.axhline(y=-2, color='red', linestyle='--', alpha=0.5)
        
        ax2.set_xlabel('Radius (kpc)')
        ax2.set_ylabel('Residuals (σ)')
        ax2.legend(loc='best')
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(-4, 4)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Model comparison plot saved to {save_path}")
        
        return fig
    
    def plot_parameter_posteriors(self,
                                 fitter: STVGFitter,
                                 save_path: Optional[str] = None) -> plt.Figure:
        """
        Create corner plot of parameter posteriors from MCMC.
        
        Args:
            fitter: STVGFitter with completed MCMC analysis
            save_path: Path to save figure
            
        Returns:
            matplotlib Figure object
        """
        if fitter.samples is None:
            raise ValueError("Must run MCMC analysis first")
        
        # Parameter labels with units
        labels = [
            r'$\alpha$',
            r'$\beta$', 
            r'$m_A$ (kg)',
            r'$\Upsilon_{\rm disk}$',
            r'$\Upsilon_{\rm bulge}$',
            r'$D$ (Mpc)',
            r'$i$ (deg)'
        ]
        
        # Create corner plot
        fig = corner.corner(fitter.samples, 
                           labels=labels,
                           quantiles=[0.16, 0.5, 0.84],
                           show_titles=True,
                           title_kwargs={"fontsize": 12},
                           color=self.colors['stvg'])
        
        # Add title
        fig.suptitle(f'{fitter.data.galaxy_name} - STVG Parameter Posteriors', 
                    fontsize=16, y=0.98)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Parameter posteriors saved to {save_path}")
        
        return fig
    
    def plot_parameter_evolution(self,
                                fitter: STVGFitter,
                                save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot MCMC chain evolution for convergence diagnostics.
        
        Args:
            fitter: STVGFitter with completed MCMC analysis
            save_path: Path to save figure
            
        Returns:
            matplotlib Figure object
        """
        if fitter.sampler is None:
            raise ValueError("Must run MCMC analysis first")
        
        samples = fitter.sampler.get_chain()
        
        fig, axes = plt.subplots(fitter.ndim, figsize=(12, 2*fitter.ndim), 
                                sharex=True)
        
        for i in range(fitter.ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], "k", alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(fitter.param_names[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)
        
        axes[-1].set_xlabel("Step Number")
        fig.suptitle(f'{fitter.data.galaxy_name} - MCMC Chain Evolution', 
                    fontsize=16)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Chain evolution plot saved to {save_path}")
        
        return fig
    
    def plot_scaling_relations(self,
                              galaxy_data: List[RotationCurveData],
                              stvg_results: List[Dict],
                              save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot STVG predictions for scaling relations (Tully-Fisher, RAR).
        
        Args:
            galaxy_data: List of galaxy rotation curve data
            stvg_results: List of STVG fitting results
            save_path: Path to save figure
            
        Returns:
            matplotlib Figure object
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # Extract galaxy properties
        M_bar = []
        V_flat = []
        g_obs = []
        g_bar = []
        
        for data, result in zip(galaxy_data, stvg_results):
            # Baryonic mass
            M_bar.append(data.M_disk + data.M_bulge)
            
            # Flat rotation velocity (average of outer points)
            outer_mask = data.radius > 0.7 * data.radius.max()
            V_flat.append(np.mean(data.velocity[outer_mask]))
            
            # Accelerations for RAR
            # This would require more detailed calculation
            # Placeholder for demonstration
            g_obs.append(V_flat[-1]**2 / (data.radius.max() * 3.086e19))
            g_bar.append(g_obs[-1] * 0.5)  # Simplified
        
        M_bar = np.array(M_bar)
        V_flat = np.array(V_flat)
        g_obs = np.array(g_obs)
        g_bar = np.array(g_bar)
        
        # Baryonic Tully-Fisher Relation
        ax1.scatter(M_bar, V_flat, c=self.colors['stvg'], s=50, alpha=0.7)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel(r'$M_{\rm bar}$ (M$_\odot$)')
        ax1.set_ylabel(r'$V_{\rm flat}$ (km/s)')
        ax1.set_title('Baryonic Tully-Fisher Relation')
        ax1.grid(True, alpha=0.3)
        
        # Theoretical BTFR line
        M_theory = np.logspace(8, 12, 100)
        V_theory = 50 * (M_theory / 1e10)**0.25  # Approximate scaling
        ax1.plot(M_theory, V_theory, '--', color='black', alpha=0.5, 
                label='Theoretical')
        ax1.legend()
        
        # Radial Acceleration Relation
        ax2.scatter(g_bar, g_obs, c=self.colors['stvg'], s=50, alpha=0.7)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_xlabel(r'$g_{\rm bar}$ (m/s$^2$)')
        ax2.set_ylabel(r'$g_{\rm obs}$ (m/s$^2$)')
        ax2.set_title('Radial Acceleration Relation')
        ax2.grid(True, alpha=0.3)
        
        # One-to-one line
        g_range = np.logspace(-12, -9, 100)
        ax2.plot(g_range, g_range, '--', color='black', alpha=0.5, 
                label='1:1 line')
        ax2.legend()
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Scaling relations plot saved to {save_path}")
        
        return fig
    
    def plot_parameter_constraints(self,
                                  results_list: List[Dict],
                                  save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot STVG parameter constraints across galaxy sample.
        
        Args:
            results_list: List of STVG fitting results from multiple galaxies
            save_path: Path to save figure
            
        Returns:
            matplotlib Figure object
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
        
        # Extract parameter values and errors
        params = ['alpha', 'beta', 'm_A']
        param_labels = [r'$\alpha$', r'$\beta$', r'$m_A$ (kg)']
        
        for i, (param, label) in enumerate(zip(params, param_labels)):
            ax = axes[i]
            
            values = []
            errors = []
            galaxy_names = []
            
            for result in results_list:
                summary = pd.DataFrame(result['parameter_summary'])
                param_row = summary[summary['parameter'] == param].iloc[0]
                
                values.append(param_row['median'])
                errors.append(param_row['std'])
                galaxy_names.append(result['galaxy_name'])
            
            values = np.array(values)
            errors = np.array(errors)
            
            # Plot parameter values with error bars
            x_pos = np.arange(len(values))
            ax.errorbar(x_pos, values, yerr=errors, fmt='o', 
                       color=self.colors['stvg'], capsize=5)
            
            # Horizontal line at mean value
            mean_val = np.mean(values)
            ax.axhline(y=mean_val, color='red', linestyle='--', alpha=0.7,
                      label=f'Mean = {mean_val:.3f}')
            
            ax.set_ylabel(label)
            ax.set_title(f'{param.capitalize()} Constraints')
            ax.set_xticks(x_pos)
            ax.set_xticklabels(galaxy_names, rotation=45, ha='right')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        # Fourth panel: Chi-squared distribution
        ax = axes[3]
        chi2_values = [result['chi2_best'] for result in results_list]
        ax.hist(chi2_values, bins=10, alpha=0.7, color=self.colors['stvg'])
        ax.set_xlabel(r'$\chi^2$')
        ax.set_ylabel('Number of Galaxies')
        ax.set_title('Goodness of Fit Distribution')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Parameter constraints plot saved to {save_path}")
        
        return fig

def create_summary_report(galaxy_data: List[RotationCurveData],
                         stvg_results: List[Dict],
                         output_dir: str = "./stvg_analysis_plots/") -> None:
    """
    Generate complete set of analysis plots for STVG study.
    
    Args:
        galaxy_data: List of galaxy rotation curve data
        stvg_results: List of STVG analysis results
        output_dir: Directory to save plots
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    plotter = STVGPlotter()
    
    print("Generating STVG analysis summary plots...")
    
    # Individual galaxy rotation curves
    for data, result in zip(galaxy_data, stvg_results):
        # Reconstruct galaxy model from results
        best_params = result['best_fit_params']
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
        
        galaxy = STVGGalaxy(baryonic, stvg_params)
        
        # Rotation curve plot
        plotter.plot_rotation_curve(
            data, galaxy, 
            save_path=f"{output_dir}/{data.galaxy_name}_rotation_curve.png"
        )
        plt.close()
    
    # Parameter constraints across sample
    plotter.plot_parameter_constraints(
        stvg_results,
        save_path=f"{output_dir}/parameter_constraints.png"
    )
    plt.close()
    
    # Scaling relations
    plotter.plot_scaling_relations(
        galaxy_data, stvg_results,
        save_path=f"{output_dir}/scaling_relations.png"
    )
    plt.close()
    
    print(f"Summary plots saved to {output_dir}")

if __name__ == "__main__":
    # Example usage and testing
    print("STVG Plotting Tools")
    print("=" * 30)
    
    # Create example data and model
    from data_analysis import SPARCDataLoader
    
    data = SPARCDataLoader.load_sparc_galaxy("NGC2403")
    
    # Example STVG model
    baryonic = BaryonicProfile(data.M_disk, data.R_disk, data.M_bulge, data.R_bulge)
    stvg_params = STVGParameters(alpha=0.8, beta=1.1, m_A=3e-31)
    galaxy = STVGGalaxy(baryonic, stvg_params)
    
    # Create plotter
    plotter = STVGPlotter()
    
    # Test rotation curve plot
    print("Creating rotation curve plot...")
    fig = plotter.plot_rotation_curve(data, galaxy, save_path="test_rotation_curve.png")
    plt.close()
    
    # Test model comparison
    print("Creating model comparison plot...")
    fig = plotter.plot_model_comparison(data, galaxy, save_path="test_model_comparison.png")
    plt.close()
    
    print("Plotting tools ready for STVG analysis!")

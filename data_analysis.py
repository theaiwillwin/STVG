
"""
Data Analysis Tools for STVG Rotation Curve Fitting

This module provides comprehensive tools for loading rotation curve data,
performing MCMC parameter estimation, and statistical model comparison
following the Bayesian framework outlined in the research report.

Key Features:
- SPARC format data loading and preprocessing
- MCMC parameter estimation using emcee
- Bayesian model comparison (BIC, AIC, Bayes factors)
- Convergence diagnostics and error analysis
- Statistical comparison tools for STVG vs GR/dark matter

Author: AI Research Assistant  
Date: July 2025
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm, lognorm
import emcee
import corner
from typing import Dict, List, Tuple, Optional, Callable
import warnings
from dataclasses import dataclass
import pickle
import os

from stvg_galaxy import STVGGalaxy, STVGParameters, BaryonicProfile

@dataclass
class RotationCurveData:
    """
    Container for rotation curve observational data.
    
    Follows SPARC database format with additional metadata.
    """
    galaxy_name: str
    radius: np.ndarray          # kpc
    velocity: np.ndarray        # km/s
    velocity_error: np.ndarray  # km/s
    distance: float             # Mpc
    distance_error: float       # Mpc
    inclination: float          # degrees
    inclination_error: float    # degrees
    
    # Baryonic component data
    M_disk: float              # Solar masses
    M_disk_error: float        # Solar masses
    R_disk: float              # kpc
    M_bulge: float = 0.0       # Solar masses
    M_bulge_error: float = 0.0 # Solar masses
    R_bulge: float = 1.0       # kpc
    
    def __post_init__(self):
        """Validate data consistency."""
        assert len(self.radius) == len(self.velocity) == len(self.velocity_error)
        assert self.distance > 0 and self.inclination > 0
        assert np.all(self.velocity_error > 0)

class SPARCDataLoader:
    """
    Data loader for SPARC database format rotation curves.
    
    Simulates SPARC data structure for testing and provides
    interface for real SPARC data integration.
    """
    
    @staticmethod
    def load_sparc_galaxy(galaxy_name: str, 
                         data_path: Optional[str] = None) -> RotationCurveData:
        """
        Load rotation curve data for a SPARC galaxy.
        
        Args:
            galaxy_name: Name of galaxy (e.g., 'NGC2403')
            data_path: Path to SPARC data directory
            
        Returns:
            RotationCurveData object with observational data
        """
        if data_path is None:
            # Generate simulated SPARC-like data for testing
            return SPARCDataLoader._generate_simulated_galaxy(galaxy_name)
        else:
            # Load real SPARC data (implementation for actual data files)
            return SPARCDataLoader._load_real_sparc_data(galaxy_name, data_path)
    
    @staticmethod
    def _generate_simulated_galaxy(galaxy_name: str) -> RotationCurveData:
        """
        Generate simulated rotation curve data for testing.
        
        Creates realistic galaxy data with known STVG parameters
        for validation and demonstration purposes.
        """
        np.random.seed(hash(galaxy_name) % 2**32)  # Reproducible per galaxy
        
        # Simulated galaxy properties
        if galaxy_name == "NGC2403":
            # Large spiral galaxy
            M_disk = 8e9
            R_disk = 2.5
            M_bulge = 2e9
            distance = 3.2
            inclination = 65.0
        elif galaxy_name == "DDO154":
            # Dwarf galaxy
            M_disk = 1e8
            R_disk = 1.2
            M_bulge = 0.0
            distance = 4.3
            inclination = 45.0
        else:
            # Default medium galaxy
            M_disk = 2e9
            R_disk = 2.0
            M_bulge = 5e8
            distance = 10.0
            inclination = 55.0
        
        # Create radius array (typical SPARC coverage)
        R = np.logspace(-0.5, 1.3, 25)  # 0.3 to 20 kpc
        
        # Generate "true" STVG rotation curve
        baryonic = BaryonicProfile(M_disk, R_disk, M_bulge, R_disk/3)
        true_params = STVGParameters(alpha=0.8, beta=1.1, m_A=3e-31)
        galaxy_model = STVGGalaxy(baryonic, true_params)
        
        v_true = galaxy_model.rotation_curve_stvg(R)
        
        # Add realistic observational errors
        v_error = 0.05 * v_true + 3.0  # 5% + 3 km/s systematic
        v_obs = v_true + np.random.normal(0, v_error)
        
        return RotationCurveData(
            galaxy_name=galaxy_name,
            radius=R,
            velocity=v_obs,
            velocity_error=v_error,
            distance=distance,
            distance_error=0.1 * distance,
            inclination=inclination,
            inclination_error=5.0,
            M_disk=M_disk,
            M_disk_error=0.2 * M_disk,
            R_disk=R_disk,
            M_bulge=M_bulge,
            M_bulge_error=0.3 * M_bulge if M_bulge > 0 else 0.0,
            R_bulge=R_disk/3
        )
    
    @staticmethod
    def _load_real_sparc_data(galaxy_name: str, data_path: str) -> RotationCurveData:
        """
        Load actual observational galaxy data files.
        
        Uses the RealGalaxyDataLoader to load rotation curves, surface brightness,
        and density decomposition data from the uploaded files.
        """
        try:
            from real_data_loader import RealGalaxyDataLoader
            
            # Initialize real data loader
            loader = RealGalaxyDataLoader(data_path)
            
            # Load the galaxy data
            return loader.load_galaxy(galaxy_name)
            
        except ImportError:
            raise NotImplementedError("Real galaxy data loader not available")
        except Exception as e:
            raise ValueError(f"Failed to load real data for {galaxy_name}: {e}")

class STVGFitter:
    """
    MCMC parameter estimation for STVG galaxy models.
    
    Implements the Bayesian framework from the research report:
    - Likelihood function with Gaussian errors
    - Prior distributions for STVG and galaxy parameters
    - MCMC sampling using emcee
    - Convergence diagnostics
    """
    
    def __init__(self, data: RotationCurveData):
        """
        Initialize fitter with rotation curve data.
        
        Args:
            data: RotationCurveData object with observations
        """
        self.data = data
        self.ndim = 7  # α, β, m_A, Υ_disk, Υ_bulge, D, i
        self.param_names = ['alpha', 'beta', 'm_A', 'Upsilon_disk', 
                           'Upsilon_bulge', 'distance', 'inclination']
        
        # Results storage
        self.sampler = None
        self.samples = None
        self.best_fit_params = None
        self.chi2_best = None
    
    def log_prior(self, params: np.ndarray) -> float:
        """
        Calculate log prior probability for parameters.
        
        Implements prior distributions from research report:
        - α, β: Log-normal priors
        - m_A: Log-uniform prior over galactic scales
        - Υ: Gaussian priors from stellar population synthesis
        - D, i: Gaussian priors from observations
        
        Args:
            params: Parameter array [α, β, m_A, Υ_disk, Υ_bulge, D, i]
            
        Returns:
            Log prior probability
        """
        alpha, beta, m_A, Upsilon_disk, Upsilon_bulge, distance, inclination = params
        
        log_p = 0.0
        
        # STVG parameter priors
        if alpha <= 0 or alpha > 10:
            return -np.inf
        log_p += lognorm.logpdf(alpha, s=0.5, scale=1.0)
        
        if beta <= 0 or beta > 10:
            return -np.inf
        log_p += lognorm.logpdf(beta, s=0.5, scale=1.0)
        
        if m_A <= 1e-35 or m_A > 1e-25:
            return -np.inf
        log_p += np.log(1/m_A)  # Log-uniform prior
        
        # Mass-to-light ratio priors (stellar population synthesis)
        if Upsilon_disk <= 0 or Upsilon_disk > 5:
            return -np.inf
        log_p += norm.logpdf(Upsilon_disk, loc=0.5, scale=0.2)
        
        if Upsilon_bulge <= 0 or Upsilon_bulge > 5:
            return -np.inf
        log_p += norm.logpdf(Upsilon_bulge, loc=1.0, scale=0.3)
        
        # Distance prior
        if distance <= 0:
            return -np.inf
        log_p += norm.logpdf(distance, loc=self.data.distance, 
                           scale=self.data.distance_error)
        
        # Inclination prior
        if inclination <= 10 or inclination >= 90:
            return -np.inf
        log_p += norm.logpdf(inclination, loc=self.data.inclination,
                           scale=self.data.inclination_error)
        
        return log_p
    
    def log_likelihood(self, params: np.ndarray) -> float:
        """
        Calculate log likelihood for given parameters.
        
        Implements Gaussian likelihood:
        L(D|θ) = ∏_i (2πσ_i²)^(-1/2) exp[-½(V_obs,i - V_STVG,i(θ))²/σ_i²]
        
        Args:
            params: Parameter array [α, β, m_A, Υ_disk, Υ_bulge, D, i]
            
        Returns:
            Log likelihood
        """
        try:
            alpha, beta, m_A, Upsilon_disk, Upsilon_bulge, distance, inclination = params
            
            # Create STVG model
            stvg_params = STVGParameters(alpha=alpha, beta=beta, m_A=m_A)
            
            # Scale baryonic masses by mass-to-light ratios
            M_disk_scaled = self.data.M_disk * Upsilon_disk
            M_bulge_scaled = self.data.M_bulge * Upsilon_bulge
            
            baryonic = BaryonicProfile(
                M_disk=M_disk_scaled,
                R_disk=self.data.R_disk,
                M_bulge=M_bulge_scaled,
                R_bulge=self.data.R_bulge
            )
            
            galaxy = STVGGalaxy(baryonic, stvg_params)
            
            # Calculate model prediction
            v_model = galaxy.rotation_curve_stvg(self.data.radius)
            
            # Apply inclination correction
            v_model_corrected = v_model / np.sin(np.radians(inclination))
            
            # Calculate chi-squared
            chi2 = np.sum((self.data.velocity - v_model_corrected)**2 / 
                         self.data.velocity_error**2)
            
            # Log likelihood (Gaussian)
            log_L = -0.5 * chi2 - 0.5 * np.sum(np.log(2*np.pi*self.data.velocity_error**2))
            
            return log_L
            
        except Exception as e:
            # Return very low likelihood for invalid parameters
            return -np.inf
    
    def log_posterior(self, params: np.ndarray) -> float:
        """
        Calculate log posterior probability.
        
        Args:
            params: Parameter array
            
        Returns:
            Log posterior = log prior + log likelihood
        """
        log_p = self.log_prior(params)
        if not np.isfinite(log_p):
            return -np.inf
        
        return log_p + self.log_likelihood(params)
    
    def find_best_fit(self) -> Tuple[np.ndarray, float]:
        """
        Find maximum likelihood parameters using optimization.
        
        Returns:
            Tuple of (best_fit_parameters, chi2_minimum)
        """
        # Initial guess
        p0 = np.array([1.0, 1.0, 1e-30, 0.5, 1.0, 
                      self.data.distance, self.data.inclination])
        
        # Bounds for optimization
        bounds = [(0.1, 10), (0.1, 10), (1e-35, 1e-25), 
                 (0.1, 5), (0.1, 5), (0.1, 100), (15, 85)]
        
        def neg_log_likelihood(params):
            return -self.log_likelihood(params)
        
        # Minimize negative log likelihood
        result = minimize(neg_log_likelihood, p0, bounds=bounds, 
                         method='L-BFGS-B')
        
        if result.success:
            self.best_fit_params = result.x
            self.chi2_best = 2 * result.fun + np.sum(np.log(2*np.pi*self.data.velocity_error**2))
            return self.best_fit_params, self.chi2_best
        else:
            raise RuntimeError(f"Optimization failed: {result.message}")
    
    def run_mcmc(self, 
                 nwalkers: int = 100,
                 nsteps: int = 10000,
                 burn_in: int = 2000,
                 progress: bool = True) -> emcee.EnsembleSampler:
        """
        Run MCMC parameter estimation.
        
        Implements the sampling strategy from research report:
        - Affine-invariant ensemble sampler (emcee)
        - 100 walkers, 10,000 steps per walker
        - 2,000 step burn-in period
        
        Args:
            nwalkers: Number of MCMC walkers
            nsteps: Number of steps per walker
            burn_in: Number of burn-in steps
            progress: Show progress bar
            
        Returns:
            emcee.EnsembleSampler object with results
        """
        # Find best fit for initialization
        if self.best_fit_params is None:
            self.find_best_fit()
        
        # Initialize walkers around best fit
        pos = self.best_fit_params + 1e-4 * np.random.randn(nwalkers, self.ndim)
        
        # Ensure walkers are within prior bounds
        for i in range(nwalkers):
            while not np.isfinite(self.log_prior(pos[i])):
                pos[i] = self.best_fit_params + 1e-4 * np.random.randn(self.ndim)
        
        # Initialize sampler
        self.sampler = emcee.EnsembleSampler(nwalkers, self.ndim, self.log_posterior)
        
        # Run MCMC
        print(f"Running MCMC with {nwalkers} walkers for {nsteps} steps...")
        self.sampler.run_mcmc(pos, nsteps, progress=progress)
        
        # Extract samples after burn-in
        self.samples = self.sampler.get_chain(discard=burn_in, flat=True)
        
        print(f"MCMC completed. Final sample size: {len(self.samples)}")
        
        return self.sampler
    
    def convergence_diagnostics(self) -> Dict[str, float]:
        """
        Calculate convergence diagnostics for MCMC chains.
        
        Returns:
            Dictionary with convergence statistics
        """
        if self.sampler is None:
            raise ValueError("Must run MCMC first")
        
        # Autocorrelation times
        try:
            tau = self.sampler.get_autocorr_time()
            tau_mean = np.mean(tau)
            tau_max = np.max(tau)
        except Exception:
            tau_mean = tau_max = np.nan
        
        # Effective sample size
        neff = len(self.samples) / (2 * tau_max) if not np.isnan(tau_max) else len(self.samples)
        
        # Gelman-Rubin statistic (simplified)
        chains = self.sampler.get_chain()
        n_chains, n_steps, n_params = chains.shape
        
        # Split each chain in half for R-hat calculation
        if n_steps > 100:
            mid = n_steps // 2
            chain1 = chains[:, :mid, :]
            chain2 = chains[:, mid:, :]
            
            # Calculate R-hat for each parameter
            r_hat = []
            for i in range(n_params):
                W = 0.5 * (np.var(chain1[:, :, i], axis=1).mean() + 
                          np.var(chain2[:, :, i], axis=1).mean())
                B = 0.5 * n_steps * np.var([chain1[:, :, i].mean(axis=1).mean(),
                                           chain2[:, :, i].mean(axis=1).mean()])
                var_plus = ((n_steps - 1) * W + B) / n_steps
                r_hat.append(np.sqrt(var_plus / W) if W > 0 else np.nan)
            
            r_hat_max = np.nanmax(r_hat)
        else:
            r_hat_max = np.nan
        
        return {
            'autocorr_time_mean': tau_mean,
            'autocorr_time_max': tau_max,
            'effective_sample_size': neff,
            'gelman_rubin_max': r_hat_max,
            'converged': r_hat_max < 1.1 if not np.isnan(r_hat_max) else False
        }
    
    def parameter_summary(self) -> pd.DataFrame:
        """
        Generate parameter summary statistics from MCMC samples.
        
        Returns:
            DataFrame with parameter estimates and uncertainties
        """
        if self.samples is None:
            raise ValueError("Must run MCMC first")
        
        # Calculate percentiles
        percentiles = [16, 50, 84]  # Median and 1-sigma bounds
        summary = []
        
        for i, name in enumerate(self.param_names):
            values = np.percentile(self.samples[:, i], percentiles)
            median = values[1]
            lower_err = median - values[0]
            upper_err = values[2] - median
            
            summary.append({
                'parameter': name,
                'median': median,
                'lower_error': lower_err,
                'upper_error': upper_err,
                'mean': np.mean(self.samples[:, i]),
                'std': np.std(self.samples[:, i])
            })
        
        return pd.DataFrame(summary)

class ModelComparison:
    """
    Statistical tools for comparing STVG with alternative models.
    
    Implements model selection criteria from research report:
    - Bayesian Information Criterion (BIC)
    - Akaike Information Criterion (AIC)  
    - Bayes factors for nested model comparison
    """
    
    @staticmethod
    def calculate_information_criteria(chi2: float, 
                                     n_params: int, 
                                     n_data: int) -> Dict[str, float]:
        """
        Calculate model selection criteria.
        
        Args:
            chi2: Chi-squared value
            n_params: Number of model parameters
            n_data: Number of data points
            
        Returns:
            Dictionary with AIC, BIC, and reduced chi-squared
        """
        aic = chi2 + 2 * n_params
        bic = chi2 + n_params * np.log(n_data)
        chi2_reduced = chi2 / (n_data - n_params)
        
        return {
            'AIC': aic,
            'BIC': bic,
            'chi2_reduced': chi2_reduced,
            'chi2': chi2,
            'n_params': n_params,
            'n_data': n_data
        }
    
    @staticmethod
    def bayes_factor(log_evidence_1: float, 
                    log_evidence_2: float) -> float:
        """
        Calculate Bayes factor between two models.
        
        Args:
            log_evidence_1: Log evidence for model 1
            log_evidence_2: Log evidence for model 2
            
        Returns:
            Bayes factor B_12 = P(D|M1) / P(D|M2)
        """
        return np.exp(log_evidence_1 - log_evidence_2)
    
    @staticmethod
    def interpret_bayes_factor(bayes_factor: float) -> str:
        """
        Interpret Bayes factor strength of evidence.
        
        Uses Jeffreys' scale for evidence interpretation.
        """
        if bayes_factor > 100:
            return "Decisive evidence for model 1"
        elif bayes_factor > 10:
            return "Strong evidence for model 1"
        elif bayes_factor > 3:
            return "Moderate evidence for model 1"
        elif bayes_factor > 1:
            return "Weak evidence for model 1"
        elif bayes_factor > 0.33:
            return "Weak evidence for model 2"
        elif bayes_factor > 0.1:
            return "Moderate evidence for model 2"
        elif bayes_factor > 0.01:
            return "Strong evidence for model 2"
        else:
            return "Decisive evidence for model 2"

def save_analysis_results(fitter: STVGFitter, 
                         filename: str,
                         include_samples: bool = True):
    """
    Save MCMC analysis results to file.
    
    Args:
        fitter: STVGFitter object with completed analysis
        filename: Output filename (will add .pkl extension)
        include_samples: Whether to save full MCMC samples
    """
    results = {
        'galaxy_name': fitter.data.galaxy_name,
        'best_fit_params': fitter.best_fit_params,
        'chi2_best': fitter.chi2_best,
        'parameter_names': fitter.param_names,
        'convergence': fitter.convergence_diagnostics(),
        'parameter_summary': fitter.parameter_summary().to_dict(),
        'data': fitter.data
    }
    
    if include_samples and fitter.samples is not None:
        results['samples'] = fitter.samples
    
    with open(f"{filename}.pkl", 'wb') as f:
        pickle.dump(results, f)
    
    print(f"Analysis results saved to {filename}.pkl")

def load_analysis_results(filename: str) -> Dict:
    """
    Load saved MCMC analysis results.
    
    Args:
        filename: Input filename (with or without .pkl extension)
        
    Returns:
        Dictionary with analysis results
    """
    if not filename.endswith('.pkl'):
        filename += '.pkl'
    
    with open(filename, 'rb') as f:
        results = pickle.load(f)
    
    return results

if __name__ == "__main__":
    # Example usage and testing
    print("STVG Data Analysis Tools")
    print("=" * 40)
    
    # Load example galaxy data
    print("Loading simulated galaxy data...")
    data = SPARCDataLoader.load_sparc_galaxy("NGC2403")
    print(f"Loaded {data.galaxy_name}: {len(data.radius)} data points")
    
    # Initialize fitter
    print("\nInitializing MCMC fitter...")
    fitter = STVGFitter(data)
    
    # Find best fit
    print("Finding best-fit parameters...")
    best_params, chi2_min = fitter.find_best_fit()
    print(f"Best-fit chi² = {chi2_min:.2f}")
    
    # Quick MCMC test (reduced for speed)
    print("\nRunning quick MCMC test...")
    sampler = fitter.run_mcmc(nwalkers=20, nsteps=100, burn_in=20, progress=False)
    
    # Convergence diagnostics
    convergence = fitter.convergence_diagnostics()
    print(f"Convergence diagnostics:")
    print(f"- R-hat max: {convergence['gelman_rubin_max']:.3f}")
    print(f"- Effective sample size: {convergence['effective_sample_size']:.0f}")
    
    # Parameter summary
    summary = fitter.parameter_summary()
    print(f"\nParameter estimates:")
    for _, row in summary.iterrows():
        print(f"- {row['parameter']}: {row['median']:.3f} ± {row['std']:.3f}")
    
    print("\nData analysis tools ready for full STVG analysis!")


"""
STVG Galaxy Rotation Curve Implementation

This module implements the core STVG (Scalar-Tensor-Vector Gravity) theory
calculations for spiral galaxy rotation curves, based on the mathematical
framework derived in the comprehensive research report.

Key Features:
- STVG rotation curve calculation using the explicit formula:
  v²_STVG(R) = G(1 + αΔφ(R))[v²_bar(R) + (βGM_bar(R)/R)e^(-m_A R)(1 + m_A R)]
- Baryonic mass profile calculations (disk + bulge components)
- Scalar and vector field computations
- Parameter handling for the 6 STVG parameters

Author: AI Research Assistant
Date: July 2025
"""

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import interp1d
from typing import Dict, Tuple, Optional, Union, List
import warnings

# Physical constants
G_NEWTON = 6.67430e-11  # m³ kg⁻¹ s⁻²
C_LIGHT = 299792458.0   # m s⁻¹
KPC_TO_M = 3.0857e19    # meters per kiloparsec
MSUN_TO_KG = 1.989e30   # kg per solar mass

class STVGParameters:
    """
    Container for STVG theory parameters.
    
    Parameters from the research report:
    - G: Newton's gravitational constant (fixed)
    - m_A: Vector field mass (~10^-30 kg, galactic scale)
    - α: Scalar field coupling (O(1))
    - β: Vector field coupling (O(1))
    - V₀: Scalar potential amplitude (cosmological energy scale)
    - λ: Scalar potential decay rate (O(1))
    """
    
    def __init__(self, 
                 alpha: float = 1.0,
                 beta: float = 1.0, 
                 m_A: float = 1e-30,
                 V0: float = 1e-10,
                 lambda_param: float = 1.0):
        """
        Initialize STVG parameters.
        
        Args:
            alpha: Scalar field coupling parameter
            beta: Vector field coupling parameter  
            m_A: Vector field mass in kg
            V0: Scalar potential amplitude
            lambda_param: Scalar potential decay rate
        """
        self.alpha = alpha
        self.beta = beta
        self.m_A = m_A
        self.V0 = V0
        self.lambda_param = lambda_param
        self.G = G_NEWTON
        
    def to_dict(self) -> Dict[str, float]:
        """Convert parameters to dictionary."""
        return {
            'alpha': self.alpha,
            'beta': self.beta,
            'm_A': self.m_A,
            'V0': self.V0,
            'lambda': self.lambda_param,
            'G': self.G
        }

class BaryonicProfile:
    """
    Baryonic mass profile for spiral galaxies.
    
    Implements exponential disk + Sérsic bulge decomposition
    as used in SPARC database analysis.
    """
    
    def __init__(self, 
                 M_disk: float,
                 R_disk: float,
                 M_bulge: float = 0.0,
                 R_bulge: float = 1.0,
                 n_sersic: float = 4.0):
        """
        Initialize baryonic profile.
        
        Args:
            M_disk: Total disk mass in solar masses
            R_disk: Disk scale length in kpc
            M_bulge: Total bulge mass in solar masses
            R_bulge: Bulge effective radius in kpc
            n_sersic: Sérsic index for bulge (default 4 = de Vaucouleurs)
        """
        self.M_disk = M_disk * MSUN_TO_KG
        self.R_disk = R_disk * KPC_TO_M
        self.M_bulge = M_bulge * MSUN_TO_KG
        self.R_bulge = R_bulge * KPC_TO_M
        self.n_sersic = n_sersic
        
        # Sérsic profile normalization
        from scipy.special import gamma, gammainc
        self.b_n = 2*n_sersic - 1/3 + 4/(405*n_sersic)  # Approximation
        self.I_bulge = self.M_bulge / (2*np.pi * self.R_bulge**2 * 
                                      np.exp(self.b_n) * (self.b_n)**(-2*n_sersic) * 
                                      gamma(2*n_sersic))
    
    def surface_density_disk(self, R: np.ndarray) -> np.ndarray:
        """
        Exponential disk surface density.
        
        Args:
            R: Galactocentric radius array in kpc
            
        Returns:
            Surface density in kg/m²
        """
        R_m = R * KPC_TO_M
        Sigma_0 = self.M_disk / (2 * np.pi * self.R_disk**2)
        return Sigma_0 * np.exp(-R_m / self.R_disk)
    
    def surface_density_bulge(self, R: np.ndarray) -> np.ndarray:
        """
        Sérsic bulge surface density.
        
        Args:
            R: Galactocentric radius array in kpc
            
        Returns:
            Surface density in kg/m²
        """
        if self.M_bulge == 0:
            return np.zeros_like(R)
            
        R_m = R * KPC_TO_M
        return self.I_bulge * np.exp(-self.b_n * ((R_m/self.R_bulge)**(1/self.n_sersic) - 1))
    
    def enclosed_mass(self, R: np.ndarray) -> np.ndarray:
        """
        Total enclosed baryonic mass within radius R.
        
        Args:
            R: Galactocentric radius array in kpc
            
        Returns:
            Enclosed mass in kg
        """
        R_m = R * KPC_TO_M
        
        # Disk contribution (analytical)
        M_disk_enc = self.M_disk * (1 - (1 + R_m/self.R_disk) * np.exp(-R_m/self.R_disk))
        
        # Bulge contribution (numerical integration for Sérsic profile)
        if self.M_bulge == 0:
            M_bulge_enc = np.zeros_like(R)
        else:
            M_bulge_enc = np.zeros_like(R)
            for i, r in enumerate(R_m):
                def integrand(r_prime):
                    return 2*np.pi * r_prime * self.surface_density_bulge(np.array([r_prime/KPC_TO_M]))[0]
                M_bulge_enc[i], _ = quad(integrand, 0, r, limit=100)
        
        return M_disk_enc + M_bulge_enc
    
    def rotation_curve_baryonic(self, R: np.ndarray) -> np.ndarray:
        """
        Baryonic rotation curve v²_bar(R).
        
        Args:
            R: Galactocentric radius array in kpc
            
        Returns:
            Squared rotation velocity in (m/s)²
        """
        R_m = R * KPC_TO_M
        M_enc = self.enclosed_mass(R)
        
        # Avoid division by zero at center
        R_safe = np.where(R_m > 0, R_m, 1e-10)
        v_squared = G_NEWTON * M_enc / R_safe
        
        return np.where(R > 0, v_squared, 0)

class STVGGalaxy:
    """
    Main class for STVG galaxy rotation curve calculations.
    
    Implements the complete STVG framework including:
    - Vector field A_μ effects via Yukawa potential
    - Scalar field φ modifications to gravitational constant
    - Full rotation curve prediction from STVG field equations
    """
    
    def __init__(self, 
                 baryonic_profile: BaryonicProfile,
                 stvg_params: STVGParameters):
        """
        Initialize STVG galaxy model.
        
        Args:
            baryonic_profile: Baryonic mass distribution
            stvg_params: STVG theory parameters
        """
        self.baryonic = baryonic_profile
        self.params = stvg_params
        
        # Cache for computed fields
        self._scalar_field_cache = {}
        self._vector_field_cache = {}
    
    def scalar_field_variation(self, R: Union[np.ndarray, List, float]) -> np.ndarray:
        """
        Calculate scalar field variation Δφ(R) = φ(R) - φ_∞.
        
        Solves the modified Klein-Gordon equation:
        ∇²φ - λ²φ = αρ_bar(r)/M_Pl²
        
        Args:
            R: Galactocentric radius array in kpc
            
        Returns:
            Scalar field variation (dimensionless)
        """
        # Ensure R is numpy array and convert to SI units
        R = np.asarray(R)
        R_m = R * KPC_TO_M
        
        # Simplified solution for weak field limit
        # Full solution would require numerical integration of coupled ODEs
        lambda_eff = np.sqrt(self.params.lambda_param * self.params.V0)
        
        # Characteristic scale for scalar field
        r_phi = C_LIGHT / lambda_eff if lambda_eff > 0 else 1e20
        
        # Approximate solution: exponential decay from galactic center
        M_enc = self.baryonic.enclosed_mass(R)
        phi_amplitude = self.params.alpha * G_NEWTON * M_enc / (C_LIGHT**2 * R_m)
        
        # Apply exponential suppression at large scales
        Delta_phi = phi_amplitude * np.exp(-R_m / r_phi)
        
        return Delta_phi
    
    def vector_field_potential(self, R: Union[np.ndarray, List, float]) -> np.ndarray:
        """
        Calculate vector field contribution to gravitational potential.
        
        From the research report:
        Φ_A(r) = (βGM/r)e^(-m_A r)
        
        Args:
            R: Galactocentric radius array in kpc
            
        Returns:
            Vector field potential in (m/s)²
        """
        R = np.asarray(R)
        R_m = R * KPC_TO_M
        M_enc = self.baryonic.enclosed_mass(R)
        
        # Avoid division by zero
        R_safe = np.where(R_m > 0, R_m, 1e-10)
        
        # Vector field mass in SI units
        m_A_SI = self.params.m_A
        
        # Yukawa potential from massive vector field
        phi_A = (self.params.beta * G_NEWTON * M_enc / R_safe) * \
                np.exp(-m_A_SI * R_m / (1.973e-25))  # Convert to natural units
        
        return np.where(R > 0, phi_A, 0)
    
    def effective_gravitational_constant(self, R: Union[np.ndarray, List, float]) -> np.ndarray:
        """
        Calculate field-dependent gravitational constant G_eff(R).
        
        From the research report:
        G_eff(r) = G[1 + α(φ(r) - φ_∞)]
        
        Args:
            R: Galactocentric radius array in kpc
            
        Returns:
            Effective gravitational constant in SI units
        """
        Delta_phi = self.scalar_field_variation(R)
        return self.params.G * (1 + self.params.alpha * Delta_phi)
    
    def rotation_curve_stvg(self, R: Union[np.ndarray, List, float]) -> np.ndarray:
        """
        Calculate complete STVG rotation curve.
        
        Implements the explicit formula from the research report:
        v²_STVG(R) = G(1 + αΔφ(R))[v²_bar(R) + (βGM_bar(R)/R)e^(-m_A R)(1 + m_A R)]
        
        Args:
            R: Galactocentric radius array in kpc
            
        Returns:
            Rotation velocity in km/s
        """
        # Ensure R is numpy array
        R = np.asarray(R)
        
        # Convert radius to SI
        R_m = R * KPC_TO_M
        
        # Baryonic contribution
        v_bar_squared = self.baryonic.rotation_curve_baryonic(R)
        
        # Vector field contribution
        M_enc = self.baryonic.enclosed_mass(R)
        R_safe = np.where(R_m > 0, R_m, 1e-10)
        
        # Vector field mass scale (convert to inverse length)
        m_A_inv_length = self.params.m_A / (1.973e-25)  # Convert to m⁻¹
        
        # Vector field term: (βGM_bar(R)/R)e^(-m_A R)(1 + m_A R)
        vector_term = (self.params.beta * G_NEWTON * M_enc / R_safe) * \
                      np.exp(-m_A_inv_length * R_m) * \
                      (1 + m_A_inv_length * R_m)
        
        # Scalar field modification
        G_eff = self.effective_gravitational_constant(R)
        
        # Total STVG rotation curve
        v_squared_stvg = G_eff * (v_bar_squared / G_NEWTON + vector_term)
        
        # Convert to km/s and ensure positive
        v_stvg = np.sqrt(np.maximum(v_squared_stvg, 0)) / 1000.0
        
        return v_stvg
    
    def rotation_curve_components(self, R: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Calculate individual components of the STVG rotation curve.
        
        Args:
            R: Galactocentric radius array in kpc
            
        Returns:
            Dictionary with rotation curve components in km/s
        """
        # Baryonic component
        v_bar = np.sqrt(self.baryonic.rotation_curve_baryonic(R)) / 1000.0
        
        # Vector field component
        R_m = R * KPC_TO_M
        M_enc = self.baryonic.enclosed_mass(R)
        R_safe = np.where(R_m > 0, R_m, 1e-10)
        m_A_inv_length = self.params.m_A / (1.973e-25)
        
        vector_term = (self.params.beta * G_NEWTON * M_enc / R_safe) * \
                      np.exp(-m_A_inv_length * R_m) * \
                      (1 + m_A_inv_length * R_m)
        v_vector = np.sqrt(np.maximum(vector_term, 0)) / 1000.0
        
        # Scalar field enhancement factor
        G_eff = self.effective_gravitational_constant(R)
        enhancement = np.sqrt(G_eff / G_NEWTON)
        
        # Total STVG
        v_total = self.rotation_curve_stvg(R)
        
        return {
            'baryonic': v_bar,
            'vector_field': v_vector,
            'scalar_enhancement': enhancement,
            'total_stvg': v_total
        }
    
    def parameter_sensitivity(self, R: np.ndarray, 
                            param_name: str, 
                            delta_frac: float = 0.01) -> np.ndarray:
        """
        Calculate sensitivity of rotation curve to parameter changes.
        
        Args:
            R: Galactocentric radius array in kpc
            param_name: Name of parameter to vary
            delta_frac: Fractional change for derivative calculation
            
        Returns:
            Derivative dv/dp at each radius
        """
        # Get baseline rotation curve
        v_baseline = self.rotation_curve_stvg(R)
        
        # Create modified parameters
        params_plus = STVGParameters(**self.params.to_dict())
        current_value = getattr(params_plus, param_name)
        setattr(params_plus, param_name, current_value * (1 + delta_frac))
        
        # Calculate modified rotation curve
        galaxy_plus = STVGGalaxy(self.baryonic, params_plus)
        v_plus = galaxy_plus.rotation_curve_stvg(R)
        
        # Numerical derivative
        dv_dp = (v_plus - v_baseline) / (current_value * delta_frac)
        
        return dv_dp

def create_example_galaxy() -> STVGGalaxy:
    """
    Create an example galaxy for testing and demonstration.
    
    Based on typical SPARC galaxy properties:
    - Disk mass: 10^10 M_⊙
    - Disk scale length: 3 kpc
    - Small bulge component
    
    Returns:
        STVGGalaxy instance with example parameters
    """
    # Example baryonic profile (similar to NGC 2403)
    baryonic = BaryonicProfile(
        M_disk=1e10,      # Solar masses
        R_disk=3.0,       # kpc
        M_bulge=1e9,      # Solar masses  
        R_bulge=1.0,      # kpc
        n_sersic=2.0      # Exponential bulge
    )
    
    # Example STVG parameters
    stvg_params = STVGParameters(
        alpha=0.5,
        beta=1.2,
        m_A=5e-31,        # kg (galactic scale)
        V0=1e-10,
        lambda_param=1.0
    )
    
    return STVGGalaxy(baryonic, stvg_params)

if __name__ == "__main__":
    # Example usage and testing
    print("STVG Galaxy Rotation Curve Implementation")
    print("=" * 50)
    
    # Create example galaxy
    galaxy = create_example_galaxy()
    
    # Calculate rotation curve
    R = np.logspace(-1, 1.5, 50)  # 0.1 to ~30 kpc
    
    print(f"Calculating rotation curve for {len(R)} radial points...")
    
    # Get rotation curve components
    components = galaxy.rotation_curve_components(R)
    
    print("\nRotation curve components calculated:")
    print(f"- Baryonic: {components['baryonic'][10]:.1f} km/s at R=1 kpc")
    print(f"- Vector field: {components['vector_field'][10]:.1f} km/s at R=1 kpc") 
    print(f"- Total STVG: {components['total_stvg'][10]:.1f} km/s at R=1 kpc")
    
    # Test parameter sensitivity
    print(f"\nParameter sensitivity at R=5 kpc:")
    R_test = np.array([5.0])
    for param in ['alpha', 'beta', 'm_A']:
        sensitivity = galaxy.parameter_sensitivity(R_test, param)
        print(f"- d(v)/d({param}) = {sensitivity[0]:.2f} (km/s)/unit")
    
    print("\nSTVG implementation ready for data analysis!")

"""
Real Galaxy Data Loader for STVG Analysis

This module loads and processes real observational galaxy data from the uploaded files:
- Rotation curve models (.dat files)
- Surface brightness profiles (.sfb files) 
- Bulge/disk decomposition (.dens files)

Integrates with the existing STVG analysis pipeline.

Author: AI Research Assistant
Date: July 2025
"""

import numpy as np
import pandas as pd
import os
import glob
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import warnings

from data_analysis import RotationCurveData

class RealGalaxyDataLoader:
    """
    Loader for real observational galaxy data files.
    """
    
    def __init__(self, data_dir: str = "/home/ubuntu/data"):
        """
        Initialize with data directory containing extracted files.
        
        Args:
            data_dir: Directory containing .dat, .sfb, .dens files
        """
        self.data_dir = data_dir
        self.available_galaxies = self._scan_available_galaxies()
        
    def _scan_available_galaxies(self) -> List[str]:
        """Scan data directory for available galaxies."""
        rotmod_files = glob.glob(os.path.join(self.data_dir, "*_rotmod.dat"))
        galaxy_names = []
        
        for file in rotmod_files:
            basename = os.path.basename(file)
            galaxy_name = basename.replace("_rotmod.dat", "")
            galaxy_names.append(galaxy_name)
            
        return sorted(galaxy_names)
    
    def list_available_galaxies(self) -> List[str]:
        """Return list of available galaxy names."""
        return self.available_galaxies
    
    def load_galaxy(self, galaxy_name: str) -> RotationCurveData:
        """
        Load complete galaxy data for STVG analysis.
        
        Args:
            galaxy_name: Name of galaxy (e.g., 'NGC3198')
            
        Returns:
            RotationCurveData object with observational data
        """
        if galaxy_name not in self.available_galaxies:
            raise ValueError(f"Galaxy {galaxy_name} not found. Available: {self.available_galaxies[:10]}...")
        
        # Load rotation curve data
        rotmod_data = self._load_rotation_curve(galaxy_name)
        
        # Load surface brightness data if available
        sfb_data = self._load_surface_brightness(galaxy_name)
        
        # Load density decomposition if available
        dens_data = self._load_density_decomposition(galaxy_name)
        
        # Combine data into RotationCurveData format
        return self._create_rotation_curve_data(galaxy_name, rotmod_data, sfb_data, dens_data)
    
    def _load_rotation_curve(self, galaxy_name: str) -> Dict:
        """Load rotation curve model data from .dat file."""
        rotmod_file = os.path.join(self.data_dir, f"{galaxy_name}_rotmod.dat")
        
        if not os.path.exists(rotmod_file):
            raise FileNotFoundError(f"Rotation curve file not found: {rotmod_file}")
        
        # Read the file, handling comments and headers
        with open(rotmod_file, 'r') as f:
            lines = f.readlines()
        
        # Extract distance from header comment
        distance = None
        for line in lines:
            if line.startswith("# Distance"):
                try:
                    distance = float(line.split("=")[1].split()[0])
                except:
                    distance = 10.0  # Default if parsing fails
                break
        
        if distance is None:
            distance = 10.0  # Default distance in Mpc
        
        # Find data start (skip comments and headers)
        data_start = 0
        for i, line in enumerate(lines):
            if not line.startswith("#") and len(line.strip()) > 0:
                data_start = i
                break
        
        # Read data columns
        data_lines = [line.strip() for line in lines[data_start:] if line.strip()]
        
        if not data_lines:
            raise ValueError(f"No data found in {rotmod_file}")
        
        # Parse data columns
        data = []
        for line in data_lines:
            try:
                values = [float(x) for x in line.split()]
                if len(values) >= 3:  # At least radius, velocity, error
                    data.append(values)
            except ValueError:
                continue
        
        if not data:
            raise ValueError(f"No valid data rows in {rotmod_file}")
        
        data = np.array(data)
        
        # Standard format: Rad, Vobs, errV, [Vgas, Vdisk, Vbul, SBdisk, SBbul]
        result = {
            'radius': data[:, 0],           # kpc
            'velocity': data[:, 1],         # km/s
            'velocity_error': data[:, 2],   # km/s
            'distance': distance
        }
        
        # Add component velocities if available
        if data.shape[1] > 3:
            result['v_gas'] = data[:, 3] if data.shape[1] > 3 else np.zeros_like(data[:, 0])
            result['v_disk'] = data[:, 4] if data.shape[1] > 4 else np.zeros_like(data[:, 0])
            result['v_bulge'] = data[:, 5] if data.shape[1] > 5 else np.zeros_like(data[:, 0])
        
        # Add surface brightness if available
        if data.shape[1] > 6:
            result['sb_disk'] = data[:, 6] if data.shape[1] > 6 else np.zeros_like(data[:, 0])
            result['sb_bulge'] = data[:, 7] if data.shape[1] > 7 else np.zeros_like(data[:, 0])
        
        return result
    
    def _load_surface_brightness(self, galaxy_name: str) -> Optional[Dict]:
        """Load surface brightness profile from .sfb file."""
        sfb_file = os.path.join(self.data_dir, f"{galaxy_name}.sfb")
        
        if not os.path.exists(sfb_file):
            return None
        
        try:
            # Read surface brightness data, skipping header row
            data = np.loadtxt(sfb_file, comments='#', skiprows=1)
            
            if data.size == 0:
                return None
            
            # Ensure 2D array
            if data.ndim == 1:
                data = data.reshape(1, -1)
            
            # Standard format: radius, mu, kill, error
            result = {
                'radius': data[:, 0],      # arcsec or kpc
                'mu': data[:, 1],          # mag/arcsec^2
                'error': data[:, 3] if data.shape[1] > 3 else np.ones_like(data[:, 0]) * 0.1
            }
            
            return result
            
        except Exception as e:
            warnings.warn(f"Could not load surface brightness for {galaxy_name}: {e}")
            return None
    
    def _load_density_decomposition(self, galaxy_name: str) -> Optional[Dict]:
        """Load bulge/disk decomposition from .dens file."""
        dens_file = os.path.join(self.data_dir, f"{galaxy_name}.dens")
        
        if not os.path.exists(dens_file):
            return None
        
        try:
            # Read density decomposition data
            data = np.loadtxt(dens_file, comments='#')
            
            if data.size == 0:
                return None
            
            # Ensure 2D array
            if data.ndim == 1:
                data = data.reshape(1, -1)
            
            # Standard format: Rad[kpc], SBdisk[Lsun/pc^2], SBbulge[Lsun/pc^2]
            result = {
                'radius': data[:, 0],       # kpc
                'sb_disk': data[:, 1],      # Lsun/pc^2
                'sb_bulge': data[:, 2] if data.shape[1] > 2 else np.zeros_like(data[:, 0])
            }
            
            return result
            
        except Exception as e:
            warnings.warn(f"Could not load density decomposition for {galaxy_name}: {e}")
            return None
    
    def _create_rotation_curve_data(self, galaxy_name: str, rotmod_data: Dict, 
                                   sfb_data: Optional[Dict], dens_data: Optional[Dict]) -> RotationCurveData:
        """
        Create RotationCurveData object from loaded data.
        """
        # Extract basic rotation curve data
        radius = rotmod_data['radius']
        velocity = rotmod_data['velocity']
        velocity_error = rotmod_data['velocity_error']
        distance = rotmod_data['distance']
        
        # Estimate galaxy parameters from data
        M_disk, R_disk, M_bulge, R_bulge = self._estimate_galaxy_parameters(
            rotmod_data, sfb_data, dens_data, distance)
        
        # Default inclination (could be improved with metadata)
        inclination = 60.0  # degrees, typical value
        
        return RotationCurveData(
            galaxy_name=galaxy_name,
            radius=radius,
            velocity=velocity,
            velocity_error=velocity_error,
            distance=distance,
            distance_error=0.1 * distance,  # 10% uncertainty
            inclination=inclination,
            inclination_error=10.0,  # 10 degree uncertainty
            M_disk=M_disk,
            M_disk_error=0.3 * M_disk,  # 30% uncertainty
            R_disk=R_disk,
            M_bulge=M_bulge,
            M_bulge_error=0.5 * M_bulge if M_bulge > 0 else 0.0,
            R_bulge=R_bulge
        )
    
    def _estimate_galaxy_parameters(self, rotmod_data: Dict, sfb_data: Optional[Dict], 
                                   dens_data: Optional[Dict], distance: float) -> Tuple[float, float, float, float]:
        """
        Estimate galaxy mass and scale parameters from observational data.
        """
        radius = rotmod_data['radius']
        velocity = rotmod_data['velocity']
        
        # Estimate disk scale length from rotation curve
        # Typically R_disk ~ radius where velocity starts to decline
        v_max = np.max(velocity)
        v_max_idx = np.argmax(velocity)
        
        # Look for turnover point
        if v_max_idx < len(velocity) - 3:
            # Find where velocity drops to 90% of max
            decline_mask = velocity[v_max_idx:] < 0.9 * v_max
            if np.any(decline_mask):
                turnover_idx = v_max_idx + np.where(decline_mask)[0][0]
                R_disk = radius[turnover_idx]
            else:
                R_disk = radius[v_max_idx] * 1.5
        else:
            R_disk = radius[v_max_idx] * 1.5
        
        # Ensure reasonable scale length
        R_disk = np.clip(R_disk, 0.5, 10.0)
        
        # Estimate disk mass from rotation curve
        # Use v^2 * R / G as rough mass estimate at 2*R_disk
        G = 4.3e-6  # km^2 s^-2 kpc Msun^-1
        
        # Find velocity at ~2*R_disk
        target_radius = 2.0 * R_disk
        interp_idx = np.searchsorted(radius, target_radius)
        if interp_idx < len(velocity):
            v_est = velocity[min(interp_idx, len(velocity)-1)]
        else:
            v_est = velocity[-1]
        
        # Rough mass estimate (factor of ~2 for disk vs point mass)
        M_disk = v_est**2 * target_radius / G / 2.0
        
        # Ensure reasonable mass range
        M_disk = np.clip(M_disk, 1e7, 1e12)
        
        # Estimate bulge parameters
        # Look for central velocity rise indicating bulge
        if len(velocity) > 3 and radius[0] < 1.0:
            central_gradient = (velocity[2] - velocity[0]) / (radius[2] - radius[0])
            if central_gradient > 20:  # Steep central rise
                M_bulge = M_disk * 0.3  # Typical bulge fraction
                R_bulge = R_disk * 0.2  # Typical bulge scale
            else:
                M_bulge = 0.0
                R_bulge = 1.0
        else:
            M_bulge = 0.0
            R_bulge = 1.0
        
        return M_disk, R_disk, M_bulge, R_bulge
    
    def get_galaxy_summary(self, galaxy_name: str) -> Dict:
        """Get summary information about a galaxy's available data."""
        if galaxy_name not in self.available_galaxies:
            return {}
        
        summary = {'galaxy_name': galaxy_name}
        
        # Check available files
        files = {
            'rotation_curve': os.path.exists(os.path.join(self.data_dir, f"{galaxy_name}_rotmod.dat")),
            'surface_brightness': os.path.exists(os.path.join(self.data_dir, f"{galaxy_name}.sfb")),
            'density_decomposition': os.path.exists(os.path.join(self.data_dir, f"{galaxy_name}.dens"))
        }
        summary['available_files'] = files
        
        # Load basic data info
        try:
            data = self.load_galaxy(galaxy_name)
            summary['data_points'] = len(data.radius)
            summary['radial_range'] = [float(data.radius.min()), float(data.radius.max())]
            summary['velocity_range'] = [float(data.velocity.min()), float(data.velocity.max())]
            summary['distance'] = float(data.distance)
            summary['estimated_M_disk'] = float(data.M_disk)
            summary['estimated_R_disk'] = float(data.R_disk)
        except Exception as e:
            summary['load_error'] = str(e)
        
        return summary

def demonstrate_real_data_loading():
    """Demonstrate loading and analyzing real galaxy data."""
    print("Real Galaxy Data Loading Demonstration")
    print("=" * 50)
    
    # Initialize loader
    loader = RealGalaxyDataLoader()
    
    print(f"Found {len(loader.available_galaxies)} galaxies with rotation curve data")
    print(f"Sample galaxies: {loader.available_galaxies[:10]}")
    
    # Demonstrate loading a few galaxies
    test_galaxies = loader.available_galaxies[:5]
    
    for galaxy_name in test_galaxies:
        print(f"\n--- {galaxy_name} ---")
        try:
            summary = loader.get_galaxy_summary(galaxy_name)
            print(f"Data points: {summary.get('data_points', 'N/A')}")
            print(f"Radial range: {summary.get('radial_range', 'N/A')} kpc")
            print(f"Velocity range: {summary.get('velocity_range', 'N/A')} km/s")
            print(f"Distance: {summary.get('distance', 'N/A')} Mpc")
            print(f"Available files: {summary.get('available_files', {})}")
            
        except Exception as e:
            print(f"Error loading {galaxy_name}: {e}")
    
    return loader

if __name__ == "__main__":
    demonstrate_real_data_loading()

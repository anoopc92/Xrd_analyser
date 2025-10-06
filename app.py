import streamlit as st
import json
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass, asdict, field
from datetime import datetime
from enum import Enum

import numpy as np
from scipy import signal
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, savgol_filter
from lmfit import Model, Parameters
from lmfit.models import VoigtModel, LinearModel

from pymatgen.core import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from mp_api.client import MPRester

import plotly.graph_objects as go
from plotly.subplots import make_subplots

warnings.filterwarnings('ignore')

# ============================================================================
# Configuration and Data Models
# ============================================================================

class FileFormat(Enum):
    """Supported XRD file formats"""
    CSV = "csv"
    TXT = "txt"
    XY = "xy"
    RAW = "raw"
    UXD = "uxd"


@dataclass
class XRDMetadata:
    """Metadata for XRD measurement"""
    wavelength: float = 1.5406  # Cu K-alpha default
    instrument: str = "Unknown"
    sample_id: str = "Unknown"
    scan_speed: Optional[float] = None
    step_size: Optional[float] = None
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())


@dataclass
class PeakParameters:
    """Detected peak parameters"""
    position: float  # 2θ position
    intensity: float
    fwhm: float
    area: float
    voigt_sigma: Optional[float] = None
    voigt_gamma: Optional[float] = None
    fit_quality: Optional[float] = None  # R²


@dataclass
class CrystalliteAnalysis:
    """Crystallite size and strain analysis results"""
    scherrer_size: float  # nm
    scherrer_std: Optional[float] = None
    wh_size: float = None  # nm
    wh_strain: float = None  # dimensionless
    wh_r_squared: float = None


@dataclass
class PhaseMatch:
    """Matched reference phase"""
    material_id: str
    formula: str
    space_group: str
    matched_peaks: int
    total_peaks: int
    match_score: float
    tolerance: float


@dataclass
class XRDAnalysisState:
    """Complete analysis state with provenance"""
    run_id: str
    metadata: XRDMetadata
    raw_data: Dict[str, np.ndarray]
    preprocessed_data: Dict[str, np.ndarray]
    detected_peaks: List[PeakParameters]
    crystallite_analysis: Optional[CrystalliteAnalysis]
    phase_matches: List[PhaseMatch]
    hyperparameters: Dict[str, Any]
    provenance: List[Dict[str, Any]]
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())

    def to_dict(self) -> Dict:
        """Convert to JSON-serializable dict"""
        return {
            'run_id': self.run_id,
            'timestamp': self.timestamp,
            'metadata': asdict(self.metadata),
            'detected_peaks': [asdict(p) for p in self.detected_peaks],
            'crystallite_analysis': asdict(self.crystallite_analysis) if self.crystallite_analysis else None,
            'phase_matches': [asdict(m) for m in self.phase_matches],
            'hyperparameters': self.hyperparameters,
            'provenance': self.provenance
        }


# ============================================================================
# File I/O and Format Handlers
# ============================================================================

class XRDFileLoader:
    """Heterogeneous XRD file format loader"""
    
    @staticmethod
    def load(filepath: Union[str, Path]) -> Tuple[np.ndarray, np.ndarray, XRDMetadata]:
        """
        Load XRD data from various file formats
        
        Returns:
            two_theta, intensity, metadata
        """
        filepath = Path(filepath)
        
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        ext = filepath.suffix.lower()
        
        # Dispatch to appropriate loader
        if ext in ['.csv']:
            return XRDFileLoader._load_csv(filepath)
        elif ext in ['.txt', '.xy']:
            return XRDFileLoader._load_txt(filepath)
        elif ext in ['.uxd']:
            return XRDFileLoader._load_uxd(filepath)
        elif ext in ['.raw']:
            return XRDFileLoader._load_raw(filepath)
        else:
            # Try generic loader
            return XRDFileLoader._load_generic(filepath)
    
    @staticmethod
    def _load_csv(filepath: Path) -> Tuple[np.ndarray, np.ndarray, XRDMetadata]:
        """Load CSV format"""
        try:
            data = np.loadtxt(filepath, delimiter=',', skiprows=1)
        except:
            data = np.loadtxt(filepath, delimiter=',')
        
        two_theta = data[:, 0]
        intensity = data[:, 1]
        metadata = XRDMetadata(sample_id=filepath.stem)
        
        return two_theta, intensity, metadata
    
    @staticmethod
    def _load_txt(filepath: Path) -> Tuple[np.ndarray, np.ndarray, XRDMetadata]:
        """Load TXT/XY format with header parsing"""
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Parse metadata from header
        metadata = XRDMetadata(sample_id=filepath.stem)
        data_start = 0
        
        for i, line in enumerate(lines):
            if 'wavelength' in line.lower():
                try:
                    metadata.wavelength = float(line.split()[-1])
                except:
                    pass
            if not line.strip().startswith('#') and line.strip():
                try:
                    float(line.split()[0])
                    data_start = i
                    break
                except:
                    continue
        
        # Load numerical data
        data = np.loadtxt(lines[data_start:])
        two_theta = data[:, 0]
        intensity = data[:, 1]
        
        return two_theta, intensity, metadata
    
    @staticmethod
    def _load_uxd(filepath: Path) -> Tuple[np.ndarray, np.ndarray, XRDMetadata]:
        """Load Bruker UXD format"""
        with open(filepath, 'r', encoding='latin-1') as f:
            content = f.read()
        
        metadata = XRDMetadata(sample_id=filepath.stem, instrument="Bruker")
        
        # Parse data section
        data_lines = []
        in_data = False
        for line in content.split('\n'):
            if '_2THETACOUNTS' in line or '_COUNTS' in line:
                in_data = True
                continue
            if in_data and line.strip():
                data_lines.append(line.strip())
        
        # Parse data
        data = np.array([list(map(float, line.split())) for line in data_lines])
        two_theta = data[:, 0]
        intensity = data[:, 1]
        
        return two_theta, intensity, metadata
    
    @staticmethod
    def _load_raw(filepath: Path) -> Tuple[np.ndarray, np.ndarray, XRDMetadata]:
        """Load Rigaku RAW format"""
        # Simplified RAW loader - extend for full binary format
        return XRDFileLoader._load_txt(filepath)
    
    @staticmethod
    def _load_generic(filepath: Path) -> Tuple[np.ndarray, np.ndarray, XRDMetadata]:
        """Generic loader for unknown formats"""
        data = np.loadtxt(filepath)
        two_theta = data[:, 0]
        intensity = data[:, 1]
        metadata = XRDMetadata(sample_id=filepath.stem)
        
        return two_theta, intensity, metadata


# ============================================================================
# Data Preprocessing Pipeline
# ============================================================================

class XRDPreprocessor:
    """Robust preprocessing for noisy XRD signals"""
    
    def __init__(self, 
                 smoothing_window: int = 5,
                 smoothing_order: int = 2,
                 baseline_lambda: float = 1e5,
                 noise_threshold: float = 3.0):
        self.smoothing_window = smoothing_window
        self.smoothing_order = smoothing_order
        self.baseline_lambda = baseline_lambda
        self.noise_threshold = noise_threshold
    
    def preprocess(self, two_theta: np.ndarray, intensity: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Complete preprocessing pipeline
        
        Returns:
            Dict with 'two_theta', 'intensity', 'baseline', 'corrected'
        """
        # Remove NaN and sort
        mask = ~(np.isnan(two_theta) | np.isnan(intensity))
        two_theta = two_theta[mask]
        intensity = intensity[mask]
        
        sort_idx = np.argsort(two_theta)
        two_theta = two_theta[sort_idx]
        intensity = intensity[sort_idx]
        
        # Smooth
        if len(intensity) >= self.smoothing_window:
            intensity_smooth = savgol_filter(intensity, 
                                            self.smoothing_window, 
                                            self.smoothing_order)
        else:
            intensity_smooth = intensity.copy()
        
        # Baseline correction (asymmetric least squares)
        baseline = self._baseline_als(intensity_smooth, 
                                      lam=self.baseline_lambda, 
                                      p=0.01)
        
        # Background subtraction
        intensity_corrected = intensity_smooth - baseline
        intensity_corrected = np.maximum(intensity_corrected, 0)
        
        # Normalize
        if intensity_corrected.max() > 0:
            intensity_corrected = intensity_corrected / intensity_corrected.max() * 100
        
        return {
            'two_theta': two_theta,
            'intensity': intensity_corrected,
            'baseline': baseline,
            'raw_intensity': intensity
        }
    
    @staticmethod
    def _baseline_als(y: np.ndarray, lam: float = 1e5, p: float = 0.01, 
                     niter: int = 10) -> np.ndarray:
        """Asymmetric Least Squares baseline correction"""
        from scipy import sparse
        from scipy.sparse.linalg import spsolve
        
        L = len(y)
        D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L-2))
        D = lam * D.dot(D.transpose())
        w = np.ones(L)
        W = sparse.spdiags(w, 0, L, L)
        
        for _ in range(niter):
            W.setdiag(w)
            Z = W + D
            z = spsolve(Z, w * y)
            w = p * (y > z) + (1 - p) * (y < z)
        
        return z


# ============================================================================
# Peak Detection and Voigt Fitting
# ============================================================================

class PeakDetector:
    """Robust peak detection with Voigt profile fitting"""
    
    def __init__(self,
                 prominence: float = 5.0,
                 min_distance: int = 5,
                 width_range: Tuple[float, float] = (0.1, 2.0)):
        self.prominence = prominence
        self.min_distance = min_distance
        self.width_range = width_range
    
    def detect_peaks(self, two_theta: np.ndarray, intensity: np.ndarray) -> List[PeakParameters]:
        """
        Detect and fit peaks with Voigt profiles
        
        Returns:
            List of PeakParameters with fitted values
        """
        # Initial peak detection
        peaks, properties = find_peaks(
            intensity,
            prominence=self.prominence,
            distance=self.min_distance,
            width=1
        )
        
        if len(peaks) == 0:
            return []
        
        # Fit each peak with Voigt profile
        fitted_peaks = []
        for peak_idx in peaks:
            try:
                params = self._fit_voigt_peak(two_theta, intensity, peak_idx)
                if params is not None:
                    fitted_peaks.append(params)
            except Exception as e:
                # Fallback to simple parameters
                fitted_peaks.append(PeakParameters(
                    position=two_theta[peak_idx],
                    intensity=intensity[peak_idx],
                    fwhm=0.2,
                    area=intensity[peak_idx] * 0.2
                ))
        
        return fitted_peaks
    
    def _fit_voigt_peak(self, two_theta: np.ndarray, intensity: np.ndarray, 
                       peak_idx: int) -> Optional[PeakParameters]:
        """Fit single peak with Voigt profile using lmfit"""
        # Define fitting window
        window = 30
        start = max(0, peak_idx - window)
        end = min(len(two_theta), peak_idx + window)
        
        x_fit = two_theta[start:end]
        y_fit = intensity[start:end]
        
        if len(x_fit) < 5:
            return None
        
        # Setup Voigt model
        model = VoigtModel() + LinearModel()
        
        params = model.make_params()
        params['center'].set(value=two_theta[peak_idx], min=x_fit.min(), max=x_fit.max())
        params['amplitude'].set(value=intensity[peak_idx], min=0)
        params['sigma'].set(value=0.1, min=0.01, max=1.0)
        params['gamma'].set(value=0.1, min=0.01, max=1.0)
        params['slope'].set(value=0)
        params['intercept'].set(value=0)
        
        try:
            result = model.fit(y_fit, params, x=x_fit)
            
            # Extract parameters
            center = result.params['center'].value
            amplitude = result.params['amplitude'].value
            sigma = result.params['sigma'].value
            gamma = result.params['gamma'].value
            
            # Calculate FWHM (approximation for Voigt)
            fwhm = 0.5346 * 2 * gamma + np.sqrt(0.2166 * (2 * gamma)**2 + (2.355 * sigma)**2)
            
            # Calculate area
            area = amplitude * np.pi * gamma * (1 + (sigma/gamma))
            
            return PeakParameters(
                position=center,
                intensity=amplitude,
                fwhm=fwhm,
                area=area,
                voigt_sigma=sigma,
                voigt_gamma=gamma,
                fit_quality=result.rsquared if hasattr(result, 'rsquared') else None
            )
        except:
            return None


# ============================================================================
# Scherrer and Williamson-Hall Analysis
# ============================================================================

class CrystalliteSizeAnalyzer:
    """Crystallite size and strain analysis"""
    
    def __init__(self, wavelength: float = 1.5406, shape_factor: float = 0.9):
        self.wavelength = wavelength  # Angstroms
        self.shape_factor = shape_factor  # Scherrer constant
    
    def analyze(self, peaks: List[PeakParameters]) -> CrystalliteAnalysis:
        """
        Perform Scherrer and Williamson-Hall analysis
        
        Returns:
            CrystalliteAnalysis with size and strain
        """
        if not peaks:
            return CrystalliteAnalysis(scherrer_size=0)
        
        # Scherrer analysis
        scherrer_sizes = []
        for peak in peaks:
            if peak.fwhm > 0:
                size = self._scherrer_equation(peak.position, peak.fwhm)
                if 0.5 < size < 1000:  # Reasonable range (nm)
                    scherrer_sizes.append(size)
        
        if not scherrer_sizes:
            return CrystalliteAnalysis(scherrer_size=0)
        
        scherrer_size = np.mean(scherrer_sizes)
        scherrer_std = np.std(scherrer_sizes) if len(scherrer_sizes) > 1 else None
        
        # Williamson-Hall analysis
        wh_size, wh_strain, wh_r2 = self._williamson_hall(peaks)
        
        return CrystalliteAnalysis(
            scherrer_size=scherrer_size,
            scherrer_std=scherrer_std,
            wh_size=wh_size,
            wh_strain=wh_strain,
            wh_r_squared=wh_r2
        )
    
    def _scherrer_equation(self, two_theta: float, fwhm: float) -> float:
        """
        Calculate crystallite size using Scherrer equation
        
        D = K*λ / (β*cos(θ))
        """
        theta_rad = np.radians(two_theta / 2)
        beta_rad = np.radians(fwhm)
        
        size_angstrom = (self.shape_factor * self.wavelength) / (beta_rad * np.cos(theta_rad))
        size_nm = size_angstrom / 10.0
        
        return size_nm
    
    def _williamson_hall(self, peaks: List[PeakParameters]) -> Tuple[float, float, float]:
        """
        Williamson-Hall analysis: β*cos(θ) = K*λ/D + 4*ε*sin(θ)
        
        Returns:
            size (nm), strain, R²
        """
        if len(peaks) < 3:
            return None, None, None
        
        x_data = []  # 4*sin(θ)
        y_data = []  # β*cos(θ)
        
        for peak in peaks:
            theta_rad = np.radians(peak.position / 2)
            beta_rad = np.radians(peak.fwhm)
            
            x = 4 * np.sin(theta_rad)
            y = beta_rad * np.cos(theta_rad)
            
            if 0 < x < 2 and 0 < y < 0.1:  # Reasonable ranges
                x_data.append(x)
                y_data.append(y)
        
        if len(x_data) < 3:
            return None, None, None
        
        x_data = np.array(x_data)
        y_data = np.array(y_data)
        
        # Linear fit
        coeffs = np.polyfit(x_data, y_data, 1)
        slope = coeffs[0]  # strain
        intercept = coeffs[1]  # K*λ/D
        
        # Calculate size
        size_angstrom = (self.shape_factor * self.wavelength) / intercept
        size_nm = size_angstrom / 10.0
        
        # Calculate R²
        y_pred = np.polyval(coeffs, x_data)
        ss_res = np.sum((y_data - y_pred) ** 2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        return size_nm, slope, r_squared


# ============================================================================
# Materials Project Integration
# ============================================================================

class MaterialsProjectMatcher:
    """Phase identification using Materials Project"""
    
    def __init__(self, api_key: str, wavelength: float = 1.5406):
        self.api_key = api_key
        self.wavelength = wavelength
        self.calculator = XRDCalculator(wavelength=wavelength)
    
    def match_phases(self, 
                    experimental_peaks: List[PeakParameters],
                    material_ids: Optional[List[str]] = None,
                    formulas: Optional[List[str]] = None,
                    base_tolerance: float = 0.3,
                    adaptive: bool = True) -> List[PhaseMatch]:
        """
        Match experimental peaks against Materials Project references
        
        Args:
            experimental_peaks: Detected peaks from experiment
            material_ids: Specific MP IDs to check (e.g., ['mp-149'])
            formulas: Chemical formulas to search (e.g., ['Si', 'SiO2'])
            base_tolerance: Base tolerance for peak matching (degrees 2θ)
            adaptive: Use adaptive tolerance based on peak width
        
        Returns:
            List of PhaseMatch objects sorted by match score
        """
        matches = []
        
        with MPRester(self.api_key) as mpr:
            # Get structures
            structures = []
            
            if material_ids:
                for mp_id in material_ids:
                    try:
                        struct = mpr.get_structure_by_material_id(mp_id)
                        structures.append((mp_id, struct))
                    except:
                        print(f"Warning: Could not retrieve {mp_id}")
            
            elif formulas:
                for formula in formulas:
                    try:
                        results = mpr.materials.summary.search(
                            formula=formula,
                            fields=["material_id", "structure"]
                        )
                        for doc in results[:5]:  # Limit to top 5
                            structures.append((doc.material_id, doc.structure))
                    except:
                        print(f"Warning: Could not search formula {formula}")
            
            # Match each structure
            for mp_id, structure in structures:
                match = self._match_structure(
                    experimental_peaks,
                    structure,
                    mp_id,
                    base_tolerance,
                    adaptive
                )
                if match:
                    matches.append(match)
        
        # Sort by match score
        matches.sort(key=lambda m: m.match_score, reverse=True)
        return matches
    
    def _match_structure(self,
                        experimental_peaks: List[PeakParameters],
                        structure: Structure,
                        mp_id: str,
                        base_tolerance: float,
                        adaptive: bool) -> Optional[PhaseMatch]:
        """Match experimental peaks against a single structure"""
        
        # Generate simulated pattern
        pattern = self.calculator.get_pattern(structure, two_theta_range=(10, 90))
        
        # Extract simulated peaks (above threshold)
        threshold = pattern.y.max() * 0.05
        sim_peaks = [(pattern.x[i], pattern.y[i]) 
                    for i in range(len(pattern.x)) 
                    if pattern.y[i] > threshold]
        
        if not sim_peaks:
            return None
        
        # Match peaks
        matched_count = 0
        exp_positions = [p.position for p in experimental_peaks]
        
        for sim_pos, sim_int in sim_peaks:
            # Adaptive tolerance based on peak width
            if adaptive and experimental_peaks:
                avg_fwhm = np.mean([p.fwhm for p in experimental_peaks])
                tolerance = max(base_tolerance, avg_fwhm / 2)
            else:
                tolerance = base_tolerance
            
            # Check if any experimental peak matches
            for exp_pos in exp_positions:
                if abs(exp_pos - sim_pos) <= tolerance:
                    matched_count += 1
                    break
        
        # Calculate match score
        match_score = matched_count / len(sim_peaks) if sim_peaks else 0
        
        # Only return if reasonable match
        if match_score < 0.3:
            return None
        
        return PhaseMatch(
            material_id=mp_id,
            formula=structure.composition.reduced_formula,
            space_group=structure.get_space_group_info()[0],
            matched_peaks=matched_count,
            total_peaks=len(sim_peaks),
            match_score=match_score,
            tolerance=tolerance if adaptive else base_tolerance
        )


# ============================================================================
# Visualization and Reporting
# ============================================================================

class XRDVisualizer:
    """Generate interactive and static plots"""
    
    @staticmethod
    def create_analysis_plot(state: XRDAnalysisState, 
                           output_path: Optional[Path] = None) -> go.Figure:
        """Create comprehensive analysis plot"""
        
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('XRD Pattern with Detected Peaks', 
                          'Peak Fitting Quality',
                          'Williamson-Hall Plot', 
                          'Phase Matching'),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"type": "table"}]]
        )
        
        # Plot 1: XRD pattern
        two_theta = state.preprocessed_data['two_theta']
        intensity = state.preprocessed_data['intensity']
        
        fig.add_trace(
            go.Scatter(x=two_theta, y=intensity, 
                      mode='lines', name='XRD Pattern',
                      line=dict(color='blue')),
            row=1, col=1
        )
        
        # Add detected peaks
        if state.detected_peaks:
            peak_positions = [p.position for p in state.detected_peaks]
            peak_intensities = [p.intensity for p in state.detected_peaks]
            
            fig.add_trace(
                go.Scatter(x=peak_positions, y=peak_intensities,
                          mode='markers', name='Detected Peaks',
                          marker=dict(color='red', size=10, symbol='x')),
                row=1, col=1
            )
        
        # Plot 2: Peak fit quality
        if state.detected_peaks:
            positions = [p.position for p in state.detected_peaks]
            fwhm = [p.fwhm for p in state.detected_peaks]
            
            fig.add_trace(
                go.Scatter(x=positions, y=fwhm,
                          mode='markers+lines', name='FWHM',
                          marker=dict(size=8)),
                row=1, col=2
            )
        
        # Plot 3: Williamson-Hall
        if state.crystallite_analysis and len(state.detected_peaks) >= 3:
            x_wh = []
            y_wh = []
            
            for peak in state.detected_peaks:
                theta_rad = np.radians(peak.position / 2)
                beta_rad = np.radians(peak.fwhm)
                x_wh.append(4 * np.sin(theta_rad))
                y_wh.append(beta_rad * np.cos(theta_rad))
            
            fig.add_trace(
                go.Scatter(x=x_wh, y=y_wh,
                          mode='markers', name='W-H Data',
                          marker=dict(size=10)),
                row=2, col=1
            )
            
            # Add fit line
            if state.crystallite_analysis.wh_strain:
                x_fit = np.linspace(min(x_wh), max(x_wh), 100)
                # Reconstruct fit from size and strain
                intercept = 0.9 * 1.5406 / (state.crystallite_analysis.wh_size * 10)
                y_fit = state.crystallite_analysis.wh_strain * x_fit + intercept
                
                fig.add_trace(
                    go.Scatter(x=x_fit, y=y_fit,
                              mode='lines', name='W-H Fit',
                              line=dict(dash='dash')),
                    row=2, col=1
                )
        
        # Table: Phase matches
        if state.phase_matches:
            headers = ['Material ID', 'Formula', 'Match Score', 'Peaks Matched']
            values = [
                [m.material_id for m in state.phase_matches[:5]],
                [m.formula for m in state.phase_matches[:5]],
                [f"{m.match_score:.2f}" for m in state.phase_matches[:5]],
                [f"{m.matched_peaks}/{m.total_peaks}" for m in state.phase_matches[:5]]
            ]
            
            fig.add_trace(
                go.Table(header=dict(values=headers),
                        cells=dict(values=values)),
                row=2, col=2
            )
        
        # Update layout
        fig.update_xaxes(title_text="2θ (degrees)", row=1, col=1)
        fig.update_yaxes(title_text="Intensity (a.u.)", row=1, col=1)
        fig.update_xaxes(title_text="2θ (degrees)", row=1, col=2)
        fig.update_yaxes(title_text="FWHM (degrees)", row=1, col=2)
        fig.update_xaxes(title_text="4sin(θ)", row=2, col=1)
        fig.update_yaxes(title_text="βcos(θ)", row=2, col=1)
        
        fig.update_layout(
            height=800,
            showlegend=True,
            title_text=f"XRD Analysis Report - {state.metadata.sample_id}"
        )
        
        if output_path:
            fig.write_html(str(output_path))
        
        return fig


# ============================================================================
# Main Pipeline Orchestrator
# ============================================================================

class XRDAnalysisPipeline:
    """Complete end-to-end XRD analysis pipeline"""
    
    def __init__(self, mp_api_key: Optional[str] = None):
        self.mp_api_key = mp_api_key
        self.file_loader = XRDFileLoader()
        
    def run(self,
            filepath: Union[str, Path],
            material_ids: Optional[List[str]] = None,
            formulas: Optional[List[str]] = None,
            hyperparameters: Optional[Dict[str, Any]] = None,
            output_dir: Optional[Path] = None) -> XRDAnalysisState:
        """
        Execute complete analysis pipeline
        
        Args:
            filepath: Path to XRD data file
            material_ids: List of Materials Project IDs to match
            formulas: List of chemical formulas to search
            hyperparameters: Optional dict of hyperparameters for tuning
            output_dir: Directory for output files
        
        Returns:
            XRDAnalysisState with complete analysis results
        """
        import uuid
        
        # Initialize state
        run_id = str(uuid.uuid4())[:8]
        provenance = []
        
        # Default hyperparameters
        if hyperparameters is None:
            hyperparameters = {
                'smoothing_window': 5,
                'smoothing_order': 2,
                'baseline_lambda': 1e5,
                'peak_prominence': 5.0,
                'peak_min_distance': 5,
                'match_tolerance': 0.3,
                'adaptive_tolerance': True
            }
        
        # Step 1: Load data
        provenance.append({
            'step': 'load_data',
            'timestamp': datetime.now().isoformat(),
            'filepath': str(filepath)
        })
        
        two_theta, intensity, metadata = self.file_loader.load(filepath)
        raw_data = {'two_theta': two_theta, 'intensity': intensity}
        
        # Step 2: Preprocess
        provenance.append({
            'step': 'preprocess',
            'timestamp': datetime.now().isoformat(),
            'params': {
                'smoothing_window': hyperparameters['smoothing_window'],
                'smoothing_order': hyperparameters['smoothing_order'],
                'baseline_lambda': hyperparameters['baseline_lambda']
            }
        })
        
        preprocessor = XRDPreprocessor(
            smoothing_window=hyperparameters['smoothing_window'],
            smoothing_order=hyperparameters['smoothing_order'],
            baseline_lambda=hyperparameters['baseline_lambda']
        )
        
        preprocessed_data = preprocessor.preprocess(two_theta, intensity)
        
        # Step 3: Detect peaks
        provenance.append({
            'step': 'detect_peaks',
            'timestamp': datetime.now().isoformat(),
            'params': {
                'prominence': hyperparameters['peak_prominence'],
                'min_distance': hyperparameters['peak_min_distance']
            }
        })
        
        peak_detector = PeakDetector(
            prominence=hyperparameters['peak_prominence'],
            min_distance=hyperparameters['peak_min_distance']
        )
        
        detected_peaks = peak_detector.detect_peaks(
            preprocessed_data['two_theta'],
            preprocessed_data['intensity']
        )
        
        # Step 4: Crystallite size/strain analysis
        provenance.append({
            'step': 'crystallite_analysis',
            'timestamp': datetime.now().isoformat(),
            'wavelength': metadata.wavelength
        })
        
        size_analyzer = CrystalliteSizeAnalyzer(wavelength=metadata.wavelength)
        crystallite_analysis = size_analyzer.analyze(detected_peaks)
        
        # Step 5: Phase matching (if MP API key provided)
        phase_matches = []
        if self.mp_api_key and (material_ids or formulas):
            provenance.append({
                'step': 'phase_matching',
                'timestamp': datetime.now().isoformat(),
                'material_ids': material_ids,
                'formulas': formulas,
                'tolerance': hyperparameters['match_tolerance']
            })
            
            matcher = MaterialsProjectMatcher(
                api_key=self.mp_api_key,
                wavelength=metadata.wavelength
            )
            
            phase_matches = matcher.match_phases(
                detected_peaks,
                material_ids=material_ids,
                formulas=formulas,
                base_tolerance=hyperparameters['match_tolerance'],
                adaptive=hyperparameters['adaptive_tolerance']
            )
        
        # Create analysis state
        state = XRDAnalysisState(
            run_id=run_id,
            metadata=metadata,
            raw_data=raw_data,
            preprocessed_data=preprocessed_data,
            detected_peaks=detected_peaks,
            crystallite_analysis=crystallite_analysis,
            phase_matches=phase_matches,
            hyperparameters=hyperparameters,
            provenance=provenance
        )
        
        # Step 6: Generate outputs
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Save JSON report
            json_path = output_dir / f"analysis_{run_id}.json"
            with open(json_path, 'w') as f:
                json.dump(state.to_dict(), f, indent=2, default=str)
            
            provenance.append({
                'step': 'save_json',
                'timestamp': datetime.now().isoformat(),
                'path': str(json_path)
            })
            
            # Generate plots
            html_path = output_dir / f"analysis_{run_id}.html"
            visualizer = XRDVisualizer()
            visualizer.create_analysis_plot(state, html_path)
            
            provenance.append({
                'step': 'save_plot',
                'timestamp': datetime.now().isoformat(),
                'path': str(html_path)
            })
        
        return state
    
    def optimize_hyperparameters(self,
                                filepath: Union[str, Path],
                                parameter_grid: Dict[str, List],
                                metric: str = 'num_peaks',
                                max_iterations: int = 10) -> Dict[str, Any]:
        """
        Iterative hyperparameter optimization
        
        Args:
            filepath: Path to XRD data
            parameter_grid: Dict of parameter names to lists of values
            metric: Optimization metric ('num_peaks', 'peak_quality', 'wh_r2')
            max_iterations: Maximum optimization iterations
        
        Returns:
            Best hyperparameters found
        """
        from itertools import product
        
        # Generate parameter combinations
        keys = parameter_grid.keys()
        values = parameter_grid.values()
        combinations = [dict(zip(keys, v)) for v in product(*values)]
        
        best_score = -np.inf
        best_params = None
        
        for i, params in enumerate(combinations[:max_iterations]):
            if i >= max_iterations:
                break
            
            try:
                # Run analysis with these parameters
                state = self.run(filepath, hyperparameters=params)
                
                # Calculate score based on metric
                if metric == 'num_peaks':
                    score = len(state.detected_peaks)
                elif metric == 'peak_quality':
                    quality_scores = [p.fit_quality for p in state.detected_peaks 
                                    if p.fit_quality is not None]
                    score = np.mean(quality_scores) if quality_scores else 0
                elif metric == 'wh_r2':
                    score = (state.crystallite_analysis.wh_r_squared 
                           if state.crystallite_analysis.wh_r_squared else 0)
                else:
                    score = len(state.detected_peaks)
                
                if score > best_score:
                    best_score = score
                    best_params = params.copy()
                    
            except Exception as e:
                continue
        
        return best_params if best_params else parameter_grid


# ============================================================================
# Utility Functions
# ============================================================================

def batch_analyze(file_paths: List[Union[str, Path]],
                 mp_api_key: Optional[str] = None,
                 output_dir: Optional[Path] = None,
                 **kwargs) -> List[XRDAnalysisState]:
    """
    Batch process multiple XRD files
    
    Args:
        file_paths: List of file paths to analyze
        mp_api_key: Materials Project API key
        output_dir: Output directory for results
        **kwargs: Additional arguments passed to pipeline.run()
    
    Returns:
        List of XRDAnalysisState objects
    """
    pipeline = XRDAnalysisPipeline(mp_api_key=mp_api_key)
    results = []
    
    for filepath in file_paths:
        try:
            state = pipeline.run(filepath, output_dir=output_dir, **kwargs)
            results.append(state)
            print(f"✓ Analyzed: {filepath}")
        except Exception as e:
            print(f"✗ Failed: {filepath} - {e}")
    
    return results


def compare_samples(states: List[XRDAnalysisState]) -> go.Figure:
    """
    Create comparative visualization of multiple samples
    
    Args:
        states: List of XRDAnalysisState objects to compare
    
    Returns:
        Plotly figure with overlaid patterns
    """
    fig = go.Figure()
    
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown']
    
    for i, state in enumerate(states):
        two_theta = state.preprocessed_data['two_theta']
        intensity = state.preprocessed_data['intensity']
        
        # Normalize for comparison
        intensity_norm = intensity / intensity.max() * 100
        
        color = colors[i % len(colors)]
        
        fig.add_trace(go.Scatter(
            x=two_theta,
            y=intensity_norm + i * 20,  # Offset for clarity
            mode='lines',
            name=state.metadata.sample_id,
            line=dict(color=color)
        ))
    
    fig.update_layout(
        title="XRD Pattern Comparison",
        xaxis_title="2θ (degrees)",
        yaxis_title="Normalized Intensity (a.u.)",
        height=600
    )
    
    return fig

# Add this Streamlit wrapper at the very end:
st.title("XRD Analysis Pipeline App")

# User inputs
uploaded_file = st.file_uploader("Upload XRD file (CSV, TXT, XY, UXD, RAW)", type=["csv", "txt", "xy", "uxd", "raw"])
mp_api_key = st.text_input("Materials Project API Key (optional for phase matching)", type="password")
formulas = st.text_input("Chemical Formulas (comma-separated, e.g., Si,SiO2)", "")
material_ids = st.text_input("Material IDs (comma-separated, e.g., mp-149)", "")
output_dir = Path("outputs")  # Local dir for temp outputs; Streamlit Cloud uses ephemeral storage

# Hyperparameters (optional, with defaults)
st.sidebar.header("Hyperparameters")
smoothing_window = st.sidebar.slider("Smoothing Window", 3, 10, 5)
peak_prominence = st.sidebar.slider("Peak Prominence", 1.0, 10.0, 5.0)
match_tolerance = st.sidebar.slider("Match Tolerance", 0.1, 1.0, 0.3)

hyperparameters = {
    'smoothing_window': smoothing_window,
    'peak_prominence': peak_prominence,
    'match_tolerance': match_tolerance,
    'smoothing_order': 2,  # Default
    'baseline_lambda': 1e5,  # Default
    'peak_min_distance': 5,  # Default
    'adaptive_tolerance': True  # Default
}

if uploaded_file and st.button("Run Analysis"):
    # Save uploaded file temporarily
    file_path = Path(uploaded_file.name)
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getvalue())
    
    # Parse inputs
    formulas_list = [f.strip() for f in formulas.split(",") if f.strip()] or None
    material_ids_list = [m.strip() for m in material_ids.split(",") if m.strip()] or None
    
    # Run pipeline
    with st.spinner("Running XRD analysis..."):
        pipeline = XRDAnalysisPipeline(mp_api_key=mp_api_key if mp_api_key else None)
        state = pipeline.run(
            filepath=file_path,
            formulas=formulas_list,
            material_ids=material_ids_list,
            hyperparameters=hyperparameters,
            output_dir=output_dir
        )
    
    # Display results
    st.subheader("Analysis Summary")
    st.write(f"Run ID: {state.run_id}")
    st.write(f"Peaks Detected: {len(state.detected_peaks)}")
    if state.crystallite_analysis:
        st.write(f"Scherrer Size: {state.crystallite_analysis.scherrer_size:.2f} nm")
        if state.crystallite_analysis.wh_size:
            st.write(f"W-H Size: {state.crystallite_analysis.wh_size:.2f} nm")
            st.write(f"W-H Strain: {state.crystallite_analysis.wh_strain:.4f}")
    
    if state.phase_matches:
        st.subheader("Phase Matches")
        match_data = [[m.material_id, m.formula, f"{m.match_score:.2f}", f"{m.matched_peaks}/{m.total_peaks}"] for m in state.phase_matches]
        st.table({"Material ID": [row[0] for row in match_data], "Formula": [row[1] for row in match_data], "Match Score": [row[2] for row in match_data], "Peaks Matched": [row[3] for row in match_data]})
    
    # Display plot
    st.subheader("Interactive Plot")
    visualizer = XRDVisualizer()
    fig = visualizer.create_analysis_plot(state)
    st.plotly_chart(fig, use_container_width=True)
    
    # Download JSON
    json_data = json.dumps(state.to_dict(), indent=2, default=str)
    st.download_button("Download JSON Report", json_data, file_name=f"analysis_{state.run_id}.json")
    
    # Clean up temp file
    file_path.unlink()

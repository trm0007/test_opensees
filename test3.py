# import numpy as np
# import matplotlib.pyplot as plt
# import pandas as pd

# class SeismicDesignTool:
#     def __init__(self, S, TB, TC, TD, xi, Z, I, R):
#         self.S = S
#         self.TB = TB
#         self.TC = TC
#         self.TD = TD
#         self.xi = xi
#         self.Z = Z
#         self.I = I
#         self.R = R
        
#         # Generate response spectrum
#         self.T_values, self.Sa_values = self._generate_response_spectrum()
    
#     def _generate_response_spectrum(self):
#         """Generate the response spectrum curve"""
#         def calculate_Cs(S, T, TB, TC, TD, xi):
#             mu = max((10 / (5 + xi)) ** 0.5, 0.55)
            
#             if 0 <= T <= TB:
#                 Cs = S * (1 + (T / TB) * (2.5 * mu - 1))
#             elif TB < T <= TC:
#                 Cs = 2.5 * S * mu
#             elif TC < T <= TD:
#                 Cs = 2.5 * S * mu * (TC / T)
#             elif TD < T <= 4:
#                 Cs = 2.5 * S * mu * (TC * TD / T ** 2)
#             else:
#                 Cs = 0
            
#             return Cs
        
#         T_values = np.arange(0, 4.005, 0.01)
#         Cs_values = [calculate_Cs(self.S, T, self.TB, self.TC, self.TD, self.xi) for T in T_values]
#         Sa_values = (2/3) * (self.Z * self.I / self.R) * np.array(Cs_values)
        
#         return T_values, Sa_values
    
#     def get_spectral_acceleration(self, period):
#         """Get spectral acceleration for a given period using interpolation"""
#         return np.interp(period, self.T_values, self.Sa_values)
    
#     def calculate_base_shear(self, weight, period, structure_type="building"):
#         """Calculate base shear force"""
#         Sa = self.get_spectral_acceleration(period)
        
#         # Base shear formula: V = Sa * W
#         base_shear = Sa * weight
        
#         return {
#             'period': period,
#             'spectral_acceleration': Sa,
#             'weight': weight,
#             'base_shear': base_shear,
#             'structure_type': structure_type
#         }
    
#     def design_building_elements(self, building_data):
#         """Design structural elements for a building"""
#         results = {}
        
#         # Calculate base shear
#         base_shear_data = self.calculate_base_shear(
#             building_data['total_weight'], 
#             building_data['natural_period'],
#             "building"
#         )
        
#         base_shear = base_shear_data['base_shear']
#         num_stories = building_data['num_stories']
#         story_height = building_data['story_height']
        
#         # Distribute base shear to stories (using code-based distribution)
#         story_weights = [building_data['total_weight'] / num_stories] * num_stories
#         story_heights = [(i + 1) * story_height for i in range(num_stories)]
        
#         # Calculate story forces using ASCE 7 distribution
#         total_wh = sum(w * h for w, h in zip(story_weights, story_heights))
#         story_forces = [(w * h * base_shear) / total_wh for w, h in zip(story_weights, story_heights)]
        
#         # Calculate shear and moment at each level
#         cumulative_shear = []
#         overturning_moment = []
        
#         for i in range(num_stories):
#             # Shear at level i (sum of forces above)
#             shear = sum(story_forces[i:])
#             cumulative_shear.append(shear)
            
#             # Overturning moment at level i
#             moment = sum(f * (j - i) * story_height for j, f in enumerate(story_forces) if j >= i)
#             overturning_moment.append(moment)
        
#         # Column design (simplified)
#         max_axial_force = building_data['total_weight'] / building_data['num_columns']
#         max_shear_force = max(cumulative_shear) / building_data['num_columns']
        
#         # Beam design (simplified)
#         max_beam_moment = max(overturning_moment) / building_data['num_bays']
        
#         # Foundation design
#         foundation_load = building_data['total_weight'] + base_shear * 0.1  # Including overturning
#         foundation_pressure = foundation_load / building_data['foundation_area']
        
#         results = {
#             'base_shear_analysis': base_shear_data,
#             'story_analysis': {
#                 'story_forces': story_forces,
#                 'story_shears': cumulative_shear,
#                 'overturning_moments': overturning_moment
#             },
#             'element_design': {
#                 'column_design': {
#                     'max_axial_force': max_axial_force,
#                     'max_shear_force': max_shear_force,
#                     'required_column_size': self._estimate_column_size(max_axial_force, max_shear_force)
#                 },
#                 'beam_design': {
#                     'max_moment': max_beam_moment,
#                     'required_beam_size': self._estimate_beam_size(max_beam_moment)
#                 },
#                 'foundation_design': {
#                     'total_load': foundation_load,
#                     'bearing_pressure': foundation_pressure,
#                     'required_foundation_area': foundation_load / building_data.get('soil_bearing_capacity', 200)
#                 }
#             }
#         }
        
#         return results
    
#     def design_bridge_elements(self, bridge_data):
#         """Design structural elements for a bridge"""
#         results = {}
        
#         # Calculate base shear for bridge
#         base_shear_data = self.calculate_base_shear(
#             bridge_data['total_weight'],
#             bridge_data['natural_period'],
#             "bridge"
#         )
        
#         base_shear = base_shear_data['base_shear']
        
#         # Bridge-specific calculations
#         span_length = bridge_data['span_length']
#         num_spans = bridge_data['num_spans']
#         deck_width = bridge_data['deck_width']
        
#         # Pier design
#         num_piers = num_spans + 1
#         pier_shear = base_shear / num_piers
#         pier_moment = pier_shear * bridge_data['pier_height']
        
#         # Deck design
#         deck_moment = (bridge_data['total_weight'] * span_length**2) / (8 * num_spans)
#         seismic_deck_moment = deck_moment * 1.2  # Amplification for seismic
        
#         # Foundation design for piers
#         pier_axial_load = bridge_data['total_weight'] / num_piers
#         foundation_load = pier_axial_load + pier_moment / bridge_data['pier_spacing']
        
#         results = {
#             'base_shear_analysis': base_shear_data,
#             'pier_design': {
#                 'pier_shear': pier_shear,
#                 'pier_moment': pier_moment,
#                 'required_pier_size': self._estimate_pier_size(pier_shear, pier_moment)
#             },
#             'deck_design': {
#                 'deck_moment': seismic_deck_moment,
#                 'required_deck_thickness': self._estimate_deck_thickness(seismic_deck_moment, span_length)
#             },
#             'foundation_design': {
#                 'pier_foundation_load': foundation_load,
#                 'required_foundation_size': foundation_load / bridge_data.get('soil_bearing_capacity', 300)
#             }
#         }
        
#         return results
    
#     def _estimate_column_size(self, axial_force, shear_force):
#         """Simplified column size estimation"""
#         # Simplified approach - real design requires detailed analysis
#         required_area = axial_force / 15000  # Assuming 15 MPa allowable stress
#         column_size = np.sqrt(required_area)
#         return max(300, column_size)  # Minimum 300mm
    
#     def _estimate_beam_size(self, moment):
#         """Simplified beam size estimation"""
#         # Simplified approach
#         required_section_modulus = moment / 150  # Assuming 150 MPa allowable stress
#         beam_depth = (required_section_modulus * 6) ** (1/3)
#         return max(300, beam_depth)  # Minimum 300mm depth
    
#     def _estimate_pier_size(self, shear, moment):
#         """Simplified pier size estimation"""
#         required_area = moment / (200 * 1000)  # Simplified calculation
#         pier_diameter = np.sqrt(4 * required_area / np.pi)
#         return max(800, pier_diameter)  # Minimum 800mm diameter
    
#     def _estimate_deck_thickness(self, moment, span):
#         """Simplified deck thickness estimation"""
#         thickness = span / 25  # Rule of thumb: span/25
#         return max(200, thickness)  # Minimum 200mm
    
#     def generate_design_report(self, structure_data, structure_type="building"):
#         """Generate comprehensive design report"""
#         if structure_type == "building":
#             results = self.design_building_elements(structure_data)
#         else:
#             results = self.design_bridge_elements(structure_data)
        
#         return results
    
#     def plot_design_results(self, building_results=None, bridge_results=None):
#         """Plot design results"""
#         fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        
#         # Plot 1: Response Spectrum
#         axes[0,0].plot(self.T_values, self.Sa_values, 'b-', linewidth=2)
#         axes[0,0].set_xlabel('Period T (sec)')
#         axes[0,0].set_ylabel('Spectral Acceleration Sa (g)')
#         axes[0,0].set_title('Response Spectrum')
#         axes[0,0].grid(True, alpha=0.3)
        
#         if building_results:
#             period = building_results['base_shear_analysis']['period']
#             sa = building_results['base_shear_analysis']['spectral_acceleration']
#             axes[0,0].plot(period, sa, 'ro', markersize=8, label=f'Building T={period:.2f}s')
#             axes[0,0].legend()
        
#         # Plot 2: Story Forces (if building)
#         if building_results:
#             story_forces = building_results['story_analysis']['story_forces']
#             stories = list(range(1, len(story_forces) + 1))
#             axes[0,1].barh(stories, story_forces, color='red', alpha=0.7)
#             axes[0,1].set_xlabel('Story Force (kN)')
#             axes[0,1].set_ylabel('Story Level')
#             axes[0,1].set_title('Story Forces Distribution')
#             axes[0,1].grid(True, alpha=0.3)
        
#         # Plot 3: Shear Diagram
#         if building_results:
#             story_shears = building_results['story_analysis']['story_shears']
#             axes[1,0].plot(story_shears, stories, 'g-o', linewidth=2, markersize=6)
#             axes[1,0].set_xlabel('Story Shear (kN)')
#             axes[1,0].set_ylabel('Story Level')
#             axes[1,0].set_title('Story Shear Diagram')
#             axes[1,0].grid(True, alpha=0.3)
        
#         # Plot 4: Moment Diagram
#         if building_results:
#             moments = building_results['story_analysis']['overturning_moments']
#             axes[1,1].plot(moments, stories, 'purple', linewidth=2, marker='s', markersize=6)
#             axes[1,1].set_xlabel('Overturning Moment (kN-m)')
#             axes[1,1].set_ylabel('Story Level')
#             axes[1,1].set_title('Overturning Moment Diagram')
#             axes[1,1].grid(True, alpha=0.3)
        
#         plt.tight_layout()
#         plt.show()

# # Example Usage
# def main():
#     # Initialize seismic design tool with your parameters
#     design_tool = SeismicDesignTool(
#         S=1.5, TB=0.5, TC=1.5, TD=2.0, xi=5, Z=0.2, I=1.5, R=5.0
#     )
    
#     # Example Building Data
#     building_data = {
#         'total_weight': 50000,  # kN
#         'natural_period': 1.2,  # seconds
#         'num_stories': 10,
#         'story_height': 3.0,  # meters
#         'num_columns': 20,
#         'num_bays': 5,
#         'foundation_area': 1000,  # m²
#         'soil_bearing_capacity': 200  # kN/m²
#     }
    
#     # Example Bridge Data
#     bridge_data = {
#         'total_weight': 30000,  # kN
#         'natural_period': 2.5,  # seconds
#         'span_length': 40,  # meters
#         'num_spans': 3,
#         'deck_width': 12,  # meters
#         'pier_height': 15,  # meters
#         'pier_spacing': 40,  # meters
#         'soil_bearing_capacity': 300  # kN/m²
#     }
    
#     # Perform building design
#     print("=== BUILDING DESIGN ANALYSIS ===")
#     building_results = design_tool.generate_design_report(building_data, "building")
    
#     print(f"Base Shear: {building_results['base_shear_analysis']['base_shear']:.1f} kN")
#     print(f"Spectral Acceleration: {building_results['base_shear_analysis']['spectral_acceleration']:.3f} g")
#     print(f"Required Column Size: {building_results['element_design']['column_design']['required_column_size']:.0f} mm")
#     print(f"Required Beam Depth: {building_results['element_design']['beam_design']['required_beam_size']:.0f} mm")
#     print(f"Foundation Pressure: {building_results['element_design']['foundation_design']['bearing_pressure']:.1f} kN/m²")
    
#     # Perform bridge design
#     print("\n=== BRIDGE DESIGN ANALYSIS ===")
#     bridge_results = design_tool.generate_design_report(bridge_data, "bridge")
    
#     print(f"Base Shear: {bridge_results['base_shear_analysis']['base_shear']:.1f} kN")
#     print(f"Pier Shear: {bridge_results['pier_design']['pier_shear']:.1f} kN")
#     print(f"Pier Moment: {bridge_results['pier_design']['pier_moment']:.1f} kN-m")
#     print(f"Required Pier Diameter: {bridge_results['pier_design']['required_pier_size']:.0f} mm")
#     print(f"Required Deck Thickness: {bridge_results['deck_design']['required_deck_thickness']:.0f} mm")
    
#     # Plot results
#     design_tool.plot_design_results(building_results=building_results)
    
#     return design_tool, building_results, bridge_results

# # Run the analysis
# if __name__ == "__main__":
#     design_tool, building_results, bridge_results = main()


# Three dimensional Frame: Eigenvalue Analysis & Effective Modal Mass Participation Ratios
# REFERENCES:
# Used in verification by SAP2000 and SeismoStruct:
# SAP2000 Integrated Finite Element Analysis and Design of Structures, Verification Manual,
# Computers and Structures, 2009. Example 1-024.
# SeismoStruct, Verification Report, 2020. Example 12.




import openseespy.opensees as op
import numpy as np
from typing import Dict, Tuple, List

# Constants and Unit Conversions
class Units:
    # Basic Units
    m = 1.0
    kN = 1.0
    sec = 1.0
    
    # Length
    mm = m/1000.0
    cm = m/100.0
    inch = 25.4*mm
    ft = 12.0*inch
    
    # Area
    m2 = m**2
    cm2 = cm**2
    mm2 = mm**2
    
    # Second Moment of Area
    m4 = m**4
    cm4 = cm**4
    
    # Force
    N = kN/1000.0
    kip = kN*4.448221615
    
    # Stress
    Pa = N/m2
    kPa = Pa*1.0e3
    MPa = Pa*1.0e6
    ksi = 6.8947573*MPa

def create_3d_frame_model():
    """Create the 3D frame model for modal analysis."""
    op.wipe()
    op.model('basic', '-ndm', 3, '-ndf', 6)
    
    # Frame grid
    Xs = [0, 35*Units.ft, 70*Units.ft]
    Ys = [0, 25*Units.ft, 50*Units.ft]
    Zs = [0, 13*Units.ft, 26*Units.ft]
    
    # Center of mass at each floor
    Xcm = 38*Units.ft
    Ycm = 27*Units.ft
    
    # Lumped floor masses
    massX = 6.2112*Units.kip*Units.sec**2/Units.ft
    massY = 6.2112*Units.kip*Units.sec**2/Units.ft
    
    # Material properties
    v = 0.2  # Poisson's ratio
    
    # Beam properties
    Eb = 500000*Units.kip/Units.ft**2
    Gb = Eb/(2*(1+v))
    Ab = 5*Units.ft**2
    Iyb = 2.61*Units.ft**4
    Izb = 1.67*Units.ft**4
    Jb = 0
    
    # Column properties
    Ec = 350000*Units.kip/Units.ft**2
    Gc = Ec/(2*(1+v))
    Ac = 4*Units.ft**2
    Izc = 1.25*Units.ft**4
    Iyc = 1.25*Units.ft**4
    Jc = 0
    
    # Define transformations
    ColTransf = 1
    BeamXTransf = 2
    BeamYTransf = 3
    op.geomTransf('Linear', ColTransf, 0, 1, 0)
    op.geomTransf('Linear', BeamXTransf, 0, 0, 1)
    op.geomTransf('Linear', BeamYTransf, 0, 0, 1)
    
    # Create nodes
    create_nodes(Xs, Ys, Zs, Xcm, Ycm, massX, massY)
    
    # Create elements
    create_columns(Xs, Ys, Zs, Ac, Ec, Gc, Jc, Iyc, Izc, ColTransf)
    beamEles = create_beams(Xs, Ys, Zs, Ab, Eb, Gb, Jb, Iyb, Izb, BeamXTransf, BeamYTransf)

def create_nodes(Xs: List[float], Ys: List[float], Zs: List[float], 
                Xcm: float, Ycm: float, massX: float, massY: float):
    """Create nodes for the 3D frame model."""
    storey = 0
    for k in range(len(Zs)):
        no = 1
        constrained = []
        
        for i in range(len(Xs)):
            for j in range(len(Ys)):
                nodeID = int(f"{no}00{storey}")
                op.node(nodeID, Xs[i], Ys[j], Zs[k])
                if k == 0:
                    op.fix(nodeID, 1, 1, 1, 1, 1, 1)     
                else:
                    constrained.append(nodeID)
                no += 1
                
        if k != 0:  # Center of mass
            nodeID = int(f"{no}00{storey}")
            op.node(nodeID, Xcm, Ycm, Zs[k])
            op.mass(nodeID, massX, massY, 0, 0, 0, 0)
            op.fix(nodeID, 0, 0, 1, 1, 1, 0)        
            op.rigidDiaphragm(3, nodeID, *constrained)
        storey += 1

def calculate_correlation(Ti: float, Tj: float, damping: float = 0.05) -> float:
    """
    Calculate correlation coefficient ρ_ij for CQC method.
    
    Parameters:
    -----------
    Ti, Tj : float
        Periods of modes i and j
    damping : float
        Damping ratio (default 5%)
    
    Returns:
    --------
    float
        Correlation coefficient between modes i and j
    """
    beta = Ti/Tj
    if beta < 1:
        beta = 1/beta  # Ensure beta >= 1
    
    rho = (8 * damping**2 * (1 + beta) * beta**(3/2) / 
          ((1 - beta**2)**2 + 4 * damping**2 * beta * (1 + beta)**2))
    return rho

def modal_combination(modal_responses: np.ndarray, 
                     periods: np.ndarray = None, 
                     method: str = 'SRSS',
                     damping: float = 0.05) -> float:
    """
    Combine modal responses using SRSS or CQC method.
    
    Parameters:
    -----------
    modal_responses : np.ndarray
        Array of modal responses to combine
    periods : np.ndarray, optional
        Required for CQC method
    method : str ('SRSS' or 'CQC')
        Combination method
    damping : float
        Damping ratio for CQC (default 5%)
    
    Returns:
    --------
    float
        Combined response
    """
    if method.upper() == 'SRSS':
        return np.sqrt(np.sum([r**2 for r in modal_responses]))
    
    elif method.upper() == 'CQC':
        if periods is None:
            raise ValueError("Periods array required for CQC method")
        
        total = 0.0
        nmodes = len(modal_responses)
        for i in range(nmodes):
            for j in range(nmodes):
                rho = calculate_correlation(periods[i], periods[j], damping)
                total += modal_responses[i] * modal_responses[j] * rho
        return np.sqrt(total)
    
    else:
        raise ValueError("Invalid method. Use 'SRSS' or 'CQC'")

def create_columns(Xs: List[float], Ys: List[float], Zs: List[float], 
                  Ac: float, Ec: float, Gc: float, Jc: float, 
                  Iyc: float, Izc: float, ColTransf: int):
    """Create column elements for the 3D frame model."""
    colTag = '00'
    for no in range(1, (len(Xs))*(len(Ys))+1):
        for storey in range(1, len(Zs)):
            nodeI = int(f"{no}00{storey-1}")
            nodeJ = int(f"{no}00{storey}")     
            eleTag = int(f"{no}{colTag}{storey}")
            op.element('elasticBeamColumn', eleTag, nodeI, nodeJ, 
                      Ac, Ec, Gc, Jc, Iyc, Izc, ColTransf)

def create_beams(Xs: List[float], Ys: List[float], Zs: List[float], 
                Ab: float, Eb: float, Gb: float, Jb: float, 
                Iyb: float, Izb: float, BeamXTransf: int, BeamYTransf: int) -> List[int]:
    """Create beam elements for the 3D frame model."""
    beamEles = []
    
    # Beams in X direction
    beamXtag = '01'    
    for storey in range(1, len(Zs)):
        no = 1
        for i in range(len(Xs)-1):
            for j in range(1, len(Ys)+1):
                nodeI = int(f"{i*len(Ys)+j}00{storey}")
                nodeJ = int(f"{(i+1)*len(Ys)+j}00{storey}")
                eleTag = int(f"{no}{beamXtag}{storey}")
                beamEles.append(eleTag)
                op.element('elasticBeamColumn', eleTag, nodeI, nodeJ, 
                          Ab, Eb, Gb, Jb, Iyb, Izb, BeamXTransf)
                no += 1
    
    # Beams in Y direction
    beamYtag = '02'    
    for storey in range(1, len(Zs)):
        no = 1
        for i in range(len(Xs)):
            for j in range(1, len(Ys)):
                nodeI = int(f"{i*len(Ys)+j}00{storey}")
                nodeJ = int(f"{i*len(Ys)+(j+1)}00{storey}")
                eleTag = int(f"{no}{beamYtag}{storey}")
                beamEles.append(eleTag)
                op.element('elasticBeamColumn', eleTag, nodeI, nodeJ, 
                          Ab, Eb, Gb, Jb, Iyb, Izb, BeamYTransf)
                no += 1
    return beamEles

def test_ModalAnalysis_3DFrame():
    """Test function for 3D frame modal analysis."""
    print("=====================================================================")
    print("ModalFrame3d: ModalAnalysis of 3D Frame Structure, using nodal masses")    
    
    create_3d_frame_model()
    
    # Number of eigenvalues to calculate
    numEigen = 4
    
    # Perform modal analysis
    T, Mratios, Mfactors, Mtots = ModalAnalysis(numEigen, outname='OpenSeespy', pflag=1)
    
    # Validate results against reference solutions
    validate_results(T, Mratios)
    Sa = np.array([0.5, 0.45, 0.3, 0.25])  # Example values for each mode
    
    # Calculate modal displacements (simplified example)
    modal_displacements = Mfactors[1] * Sa  # Using X-direction factors
    # Combine using different methods
    print("\nModal Combination Results:")
    print(f"SRSS combination: {modal_combination(modal_displacements):.4f}")
    print(f"CQC combination: {modal_combination(modal_displacements, T, 'CQC'):.4f}")

def validate_results(T: np.ndarray, Mratios: Dict[int, np.ndarray]):
    """Validate modal analysis results against reference solutions."""
    ok = 0
    
    # Reference results from SAP2000 and SeismoStruct
    comparisonResults = [
        [0.227062, 0.215633, 0.073345, 0.072005],  # SAP2000
        [0.22706191, 0.21563345, 0.07334548, 0.07200536]  # SeismoStruct
    ]
    
    print("\n\nComparisons of Periods [sec]:")
    print('{:>10}{:>15}{:>15}{:>15}'.format('Mode', 'OpenSees', 'SAP2000', 'SeismoStruct'))
    
    for i in range(len(T)):
        print('{:>10}{:>15.5f}{:>15.4f}{:>15.4f}'.format(
            i + 1, T[i], comparisonResults[0][i], comparisonResults[1][i]))
        if abs(T[i] - comparisonResults[0][i]) > 1e-5:
            ok -= 1
    
    # Validate mass participation ratios
    comparisonResults = [
        [90.258, 0.19, 9.388, 0.164],  # SAP2000 U1
        [90.257941, 0.189952, 9.38773, 0.164377]  # SeismoStruct U1
    ]
    print("\n\nComparisons for Effective Modal Mass Participating in U₁ [%]:")
    print('{:>10}{:>15}{:>15}{:>15}'.format('Mode', 'OpenSees', 'SAP2000', 'SeismoStruct'))
    
    for i in range(len(T)):
        print('{:>10}{:>15.4f}{:>15.4f}{:>15.4f}'.format(
            i + 1, Mratios[1][i], comparisonResults[0][i], comparisonResults[1][i]))
        if abs(Mratios[1][i] - comparisonResults[0][i]) > 1e-3:
            ok -= 1
            
    comparisonResults = [
        [0.192, 91.046, 0.151, 8.612],  # SAP2000 U2
        [0.191706, 91.046194, 0.150526, 8.611574]  # SeismoStruct U2
    ]
    print("\n\nComparisons for Effective Modal Mass Participating in U₂ [%]:")
    print('{:>10}{:>15}{:>15}{:>15}'.format('Mode', 'OpenSees', 'SAP2000', 'SeismoStruct'))
    
    for i in range(len(T)):
        print('{:>10}{:>15.4f}{:>15.4f}{:>15.4f}'.format(
            i + 1, Mratios[2][i], comparisonResults[0][i], comparisonResults[1][i]))
        if abs(Mratios[2][i] - comparisonResults[0][i]) > 1e-3:
            ok -= 1

    assert ok == 0

def ModalAnalysis(numEigen: int, pflag: int = 1, outname: str = None) -> Tuple[
    np.ndarray, Dict[int, np.ndarray], Dict[int, np.ndarray], Dict[int, float]]:
    """
    Perform modal analysis of an OpenSees model.
    
    Parameters:
    -----------
    numEigen : int
        Number of eigenvalues to calculate
    pflag : int (1 or 0)
        Flag to print output information on screen
    outname : str, optional
        If not None and pFlag==1, writes modal properties to outname.csv
    
    Returns:
    --------
    T : np.ndarray
        Period array for the first numEigen modes
    Mratios : dict
        Effective modal mass participation ratios
    Mfactors : dict
        Modal participation factors
    Mtots : dict
        Total activated masses
    """
    print("=== MODAL ANALYSIS ===")

    # Set up the system for eigenvalue analysis
    # Análisis modal
    op.wipeAnalysis()

    # Set up the system for eigenvalue analysis
    op.system('BandGeneral')  # or 'FullGeneral' for smaller models
    op.numberer('RCM')
    op.constraints('Transformation')
    op.integrator('LoadControl', 1.0)
    op.algorithm('Linear')
    op.analysis('Static')


    N = op.systemSize()  # Number of equations
    Mmatrix = np.array(op.printA('-ret')).reshape((N, N))
    
    print('\n************************************************************')
    print('Extracting the mass matrix, ignore the warnings...')
        
    # Determine maximum number of DOFs/node
    NDF = max(len(op.nodeDOFs(node)) for node in op.getNodeTags())

    # Initialize dictionaries for results
    DOFs = []  # Indices of unrestrained DOFs
    used = {}  # Nodes and associated unrestrained DOFs
    ldict = {i: np.zeros((N, 1)) for i in range(1, NDF+1)}
    Mratios = {i: np.zeros(numEigen) for i in range(1, NDF+1)}
    Mfactors = {i: np.zeros(numEigen) for i in range(1, NDF+1)}
    
    # Create influence vectors and get unrestrained DOFs
    idx = 0
    for node in op.getNodeTags():
        used[node] = []
        ndof = len(op.nodeDOFs(node))
        for j in range(ndof):
            temp = op.nodeDOFs(node)[j]
            if temp not in DOFs and temp >= 0:
                DOFs.append(temp)
                used[node].append(j+1)
                ldict[j+1][idx, 0] = 1
                idx += 1

    Mmatrix = Mmatrix[DOFs, :][:, DOFs]  # Reorganize mass matrix

    # Calculate total masses assigned to unrestrained DOFs
    Mtots = {i: (ldict[i].T @ Mmatrix @ ldict[i])[0, 0] for i in range(1, NDF+1)}

    # Perform eigenvalue analysis
    op.wipeAnalysis()
    eigenValues = perform_eigenvalue_analysis(numEigen)
    
    # Calculate modal properties
    Lambda = np.asarray(eigenValues)
    Omega = Lambda**0.5
    T = 2*np.pi/Omega
    frq = 1/T

    # Obtain modal properties
    for mode in range(1, numEigen+1):
        idx = 0
        phi = np.zeros((N, 1))
        for node in used:
            for dof in used[node]:
                phi[idx, 0] = op.nodeEigenvector(node, mode, dof)
                idx += 1
                
        phi = phi/(phi.T @ Mmatrix @ phi)**0.5  # Normalize eigenvector
        Mn = phi.T @ Mmatrix @ phi              # Modal mass (should be 1)

        for j in range(1, NDF+1):
            if Mtots[j] != 0:
                Ln = phi.T @ Mmatrix @ ldict[j]                # Modal excitation factor
                Mnstar = (Ln**2/Mn)[0, 0]                      # Effective modal mass
                Mfactors[j][mode-1] = float(Ln/Mn)                   # Modal participation factor
                Mratios[j][mode-1] = (Mnstar/Mtots[j] * 100)   # Effective modal mass ratio [%]

    # Remove rotational DOFs (not correctly implemented yet)
    for j in [4, 5, 6]:
        Mratios.pop(j, None)
        Mfactors.pop(j, None)

    # Calculate cumulative modal mass participation ratio
    sM1 = np.cumsum(Mratios[1])
    sM2 = np.cumsum(Mratios[2])
    sM3 = np.cumsum(Mratios[3])
    
    # Print results if requested
    if pflag == 1:
        print_results(T, frq, Omega, Lambda, Mtots, Mfactors, Mratios, sM1, sM2, sM3, outname)
        
        # Additional validation
        validate_modal_participation(Mratios)
        analyze_structural_behavior(T, Mfactors)

    return T, Mratios, Mfactors, Mtots

# Try this more robust approach for eigenvalue analysis
def perform_eigenvalue_analysis(numEigen):
    try:
        # First try with faster solver
        eigenValues = op.eigen('-genBandArpack', numEigen)
        if all(val > 0 for val in eigenValues):
            return eigenValues
    except:
        pass
    
    # Fallback to more reliable solver
    print("Warning: Using slower fullGenLapack solver")
    try:
        eigenValues = op.eigen('-fullGenLapack', numEigen)
        if all(val > 0 for val in eigenValues):
            return eigenValues
    except:
        print("Error: All eigenvalue solvers failed")
        raise

def print_results(T, frq, Omega, Lambda, Mtots, Mfactors, Mratios, sM1, sM2, sM3, outname):
    """Print modal analysis results."""
    arguments = []
    arguments.append('Modal Periods and Frequencies')
    arguments.append('%4s|%8s|%10s|%12s|%12s' %
          ('Mode', 'T [sec]', 'f [Hz]', 'ω [rad/sec]', 'λ [rad²/sec²]'))
    
    for mode in range(len(T)):      
        arguments.append('%4s|%8s|%10s|%12s|%12s' %
              (f"{mode+1:.0f}", f"{T[mode]:.4f}", f"{frq[mode]:.3f}",
               f"{Omega[mode]:.2f}", f"{Lambda[mode]:.2f}"))
    
    arguments.append('Total Activated Masses')
    arguments.append('%8s|%8s|%8s' % ('M₁', 'M₂', 'M₃'))
    arguments.append('%8s|%8s|%8s' %
          (f"{Mtots[1]:.2f}", f"{Mtots[2]:.2f}", f"{Mtots[3]:.2f}"))
    
    arguments.append('Modal Mass Participation Factors') 
    arguments.append('%4s|%7s|%7s|%7s' %
        ('Mode', 'Γ₁', 'Γ₂', 'Γ₃'))             
    
    for mode in range(len(T)):
        arguments.append('%4s|%7s|%7s|%7s' % (f"{mode+1:.0f}",
            f"{Mfactors[1][mode]:.3f}", f"{Mfactors[2][mode]:.3f}", f"{Mfactors[3][mode]:.3f}"))  
    
    arguments.append('Effective Modal Mass Participation Ratios [%]') 
    arguments.append('%4s|%7s|%7s|%7s' %
        ('Mode', 'U₁', 'U₂', 'U₃'))              
    
    for mode in range(len(T)):
        arguments.append('%4s|%7s|%7s|%7s' % (f"{mode+1:.0f}",
            f"{Mratios[1][mode]:.3f}", f"{Mratios[2][mode]:.3f}", f"{Mratios[3][mode]:.3f}"))  
    
    arguments.append('Cumulative Effective Modal Mass Participation Ratios [%]') 
    arguments.append('%4s|%7s|%7s|%7s' %
        ('Mode', '∑U₁', '∑U₂', '∑U₃'))              
    
    for mode in range(len(T)):
        arguments.append('%4s|%7s|%7s|%7s' % (f"{mode+1:.0f}",
            f"{sM1[mode]:.3f}", f"{sM2[mode]:.3f}", f"{sM3[mode]:.3f}"))  

    # Print to screen
    print('\n'.join(arguments))
 
    # Write to file if requested
    if outname is not None:
        with open(f"{outname}.csv", 'w', encoding='utf-8') as f:
            f.write('\n'.join(arguments))

def validate_modal_participation(Mratios: Dict[int, np.ndarray], threshold: float = 0.90):
    """
    Validate if enough modes were captured based on effective modal mass participation.
    
    Parameters:
    -----------
    Mratios : dict
        Dictionary from ModalAnalysis containing modal mass ratios
    threshold : float (default=0.90)
        Minimum required mass participation (e.g., 0.90 for 90%)
    """
    directions = {1: "U₁ (X)", 2: "U₂ (Y)", 3: "U₃ (Z)"}
    
    print("\nModal Participation Validation:")
    for dir_key, dir_name in directions.items():
        cumulative_mass = np.cumsum(Mratios[dir_key])
        last_mode_participation = cumulative_mass[-1]
        
        if last_mode_participation >= threshold * 100:
            print(f"✓ {dir_name}: Captured {last_mode_participation:.1f}% of mass (≥{threshold*100:.0f}%)")
        else:
            print(f"⚠ {dir_name}: Only {last_mode_participation:.1f}% of mass captured (<{threshold*100:.0f}%)")

def analyze_structural_behavior(periods: np.ndarray, Mfactors: Dict[int, np.ndarray]):
    """
    Analyze structural behavior based on modal properties.
    
    Parameters:
    -----------
    periods : np.ndarray
        Array of modal periods
    Mfactors : dict
        Dictionary of modal participation factors
    """
    fundamental_period = periods[0]
    fundamental_participation = Mfactors[1][0]
    
    print("\nStructural Behavior Analysis:")
    print(f"Fundamental period: {fundamental_period:.3f} seconds")
    
    # Structural stiffness indication
    if fundamental_period > 1.0:
        print("→ Flexible structure (low frequency)")
    elif fundamental_period < 0.3:
        print("→ Stiff structure (high frequency)")
    else:
        print("→ Moderately flexible structure")
    
    # Participation factor significance
    if abs(fundamental_participation) > 0.5:
        print("→ First mode dominates response")
    else:
        print("→ Higher modes are important")

# Run the test
if __name__ == "__main__":
    test_ModalAnalysis_3DFrame()
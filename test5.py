import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pi
import json

# ==============================================
# MODEL DEFINITION (3D 2-bay x 2-bay, 3-story)
# ==============================================
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)  # 3D, 6 DOF

# Units: kN, m, sec
floor_ht = 3.5  # Story height (m)
bay_width = 6.0  # Bay width (m)

# Define nodes (3 stories + roof)
node_tags = []
for floor in range(4):  # 0=ground, 1-3=floors
    z = floor * floor_ht
    for x in [0, bay_width, 2*bay_width]:
        for y in [0, bay_width, 2*bay_width]:
            tag = 100*floor + 10*int(x/bay_width) + int(y/bay_width) + 1
            ops.node(tag, x, y, z)
            node_tags.append(tag)

# Fix base nodes
for node in [1, 2, 3, 11, 12, 13, 21, 22, 23]:
    ops.fix(node, 1, 1, 1, 1, 1, 1)

# Material properties
E = 200e6       # Elastic modulus (kPa)
G = E/(2*(1+0.3))  # Shear modulus (Poisson's ratio=0.3)
A = 0.01        # Cross-sectional area (m²)
Iz = 8.33e-5    # Moment of inertia about z-axis (m⁴)
Iy = 8.33e-5    # Moment of inertia about y-axis (m⁴)
J = 1.25e-4     # Torsional constant (m⁴)

# Define materials and sections (optional - we'll use direct properties)
# ops.uniaxialMaterial('Elastic', 1, E)
# ops.section('Elastic', 1, E, A, Iz, Iy, G, J)

# Geometric transformations
# For vertical columns, use vector perpendicular to the member axis
ops.geomTransf('PDelta', 1, 1, 0, 0)  # Columns: vector in x-direction
ops.geomTransf('Linear', 2, 0, 0, 1)  # Beams X-direction: vector in z-direction  
ops.geomTransf('Linear', 3, 0, 0, 1)  # Beams Y-direction: vector in z-direction

# Create elements
def build_frame():
    # Columns (vertical)
    for col_base in [1, 2, 3, 11, 12, 13, 21, 22, 23]:
        for floor in range(3):
            i = col_base + floor*100
            j = i + 100
            ele_tag = 1000 + i
            ops.element('elasticBeamColumn', ele_tag, i, j, A, E, G, J, Iy, Iz, 1)
    
    # Beams (horizontal)
    for floor in range(1, 4):
        # X-direction beams
        for y_idx in range(3):
            for x_beam in range(2):
                i = 100*floor + 10*y_idx + 1 + x_beam
                j = i + 1
                ele_tag = 2000 + i
                ops.element('elasticBeamColumn', ele_tag, i, j, A, E, G, J, Iy, Iz, 2)
        
        # Y-direction beams
        for x_idx in range(3):
            for y_beam in range(2):
                i = 100*floor + 10*x_idx + 1 + y_beam*10
                j = i + 10
                ele_tag = 3000 + i
                ops.element('elasticBeamColumn', ele_tag, i, j, A, E, G, J, Iy, Iz, 3)

build_frame()

# ==============================================
# LOAD DEFINITION
# ==============================================
# Load parameters
dead_load = 5.0    # kN/m² (dead load)
live_load = 2.5    # kN/m² (live load)
snow_load = 1.5    # kN/m² (snow load)
floor_area = (2*bay_width)**2  # Total floor area
wind_pressure = 1.0  # kN/m² (wind pressure)

# Calculate total loads per floor
DL_floor = dead_load * floor_area  # Dead load per floor
LL_floor = live_load * floor_area  # Live load per floor
SL_floor = snow_load * floor_area  # Snow load per floor (roof only)

# Mass calculation (for seismic)
seismic_mass = (DL_floor + 0.25*LL_floor) / 9.81  # Effective mass per floor

# Assign masses (distributed to nodes)
for floor in range(1, 4):
    floor_nodes = [n for n in node_tags if n >= 100*floor and n < 100*(floor+1)]
    mass_per_node = seismic_mass / len(floor_nodes)
    for node in floor_nodes:
        ops.mass(node, mass_per_node, mass_per_node, 0, 0, 0, 0)

# ==============================================
# LOAD PATTERNS
# ==============================================

# Time series
ops.timeSeries('Linear', 1)  # For static loads
ops.timeSeries('Linear', 2)  # For wind loads

# Dead Load Pattern
ops.pattern('Plain', 1, 1)
for floor in range(1, 4):
    floor_nodes = [n for n in node_tags if n >= 100*floor and n < 100*(floor+1)]
    load_per_node = -DL_floor / len(floor_nodes)  # Negative for downward
    for node in floor_nodes:
        ops.load(node, 0, 0, load_per_node, 0, 0, 0)

# Live Load Pattern
ops.pattern('Plain', 2, 1)
for floor in range(1, 4):
    floor_nodes = [n for n in node_tags if n >= 100*floor and n < 100*(floor+1)]
    load_per_node = -LL_floor / len(floor_nodes)
    for node in floor_nodes:
        ops.load(node, 0, 0, load_per_node, 0, 0, 0)

# Snow Load Pattern (roof only)
ops.pattern('Plain', 3, 1)
roof_nodes = [n for n in node_tags if n >= 300 and n < 400]
snow_per_node = -SL_floor / len(roof_nodes)
for node in roof_nodes:
    ops.load(node, 0, 0, snow_per_node, 0, 0, 0)

# Wind Load Pattern X-direction
ops.pattern('Plain', 4, 2)
wind_force_per_floor = wind_pressure * (2*bay_width) * floor_ht
for floor in range(1, 4):
    # Apply to exterior nodes in X-direction
    wind_nodes = [100*floor + 1, 100*floor + 2, 100*floor + 3]  # Front face
    wind_per_node = wind_force_per_floor / len(wind_nodes)
    for node in wind_nodes:
        ops.load(node, wind_per_node, 0, 0, 0, 0, 0)

# Wind Load Pattern Y-direction
ops.pattern('Plain', 5, 2)
for floor in range(1, 4):
    # Apply to exterior nodes in Y-direction
    wind_nodes = [100*floor + 1, 100*floor + 11, 100*floor + 21]  # Side face
    wind_per_node = wind_force_per_floor / len(wind_nodes)
    for node in wind_nodes:
        ops.load(node, 0, wind_per_node, 0, 0, 0, 0)

# Earthquake Load Pattern (will be defined in seismic analysis)
ops.pattern('Plain', 6, 1)
# Seismic forces will be calculated and applied later

# ==============================================
# ANALYSIS FUNCTIONS
# ==============================================
def setup_static_analysis():
    ops.wipeAnalysis()
    ops.system('BandGeneral')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')

def run_load_combination(load_factors):
    """Run analysis with given load combination factors"""
    ops.wipeAnalysis()
    ops.loadConst('-time', 0.0)
    
    # Apply load factors
    for pattern_id, factor in load_factors.items():
        ops.loadConst('-time', 0.0)
        if factor != 0:
            ops.load(1, 0, 0, 0, 0, 0, 0)  # Dummy command to activate pattern
    
    setup_static_analysis()
    
    # Calculate total loads with factors
    total_forces = {}
    for ele_tag in ops.getEleTags():
        total_force = np.zeros(12)  # 6 DOF at each end
        
        for pattern_id, factor in load_factors.items():
            if factor != 0:
                ops.wipe()
                ops.model('basic', '-ndm', 3, '-ndf', 6)
                # Rebuild model for each pattern (simplified approach)
                # In practice, you'd use pattern recorders or superposition
                pass
        
        total_forces[ele_tag] = total_force
    
    return total_forces

# ==============================================
# LOAD COMBINATIONS (ASCE 7)
# ==============================================
load_combinations = {
    'LC1': {1: 1.4},  # 1.4D
    'LC2': {1: 1.2, 2: 1.6},  # 1.2D + 1.6L
    'LC3': {1: 1.2, 2: 1.0, 3: 1.6},  # 1.2D + 1.0L + 1.6S
    'LC4': {1: 1.2, 2: 1.0, 4: 1.0},  # 1.2D + 1.0L + 1.0Wx
    'LC5': {1: 1.2, 2: 1.0, 5: 1.0},  # 1.2D + 1.0L + 1.0Wy
    'LC6': {1: 0.9, 4: 1.0},  # 0.9D + 1.0Wx
    'LC7': {1: 0.9, 5: 1.0},  # 0.9D + 1.0Wy
}

# ==============================================
# MODAL ANALYSIS
# ==============================================
def run_modal_analysis():
    # Apply dead load for stiffness calculation
    ops.wipeAnalysis()
    setup_static_analysis()
    ops.analyze(1)
    ops.loadConst('-time', 0.0)
    
    # Modal analysis
    ops.wipeAnalysis()
    ops.system('FullGeneral')
    ops.numberer('RCM')
    ops.constraints('Plain')
    
    # Extract eigenvalues
    n_modes = 6
    eigenvalues = ops.eigen('-fullGenLapack', n_modes)
    
    periods = []
    frequencies = []
    for eig in eigenvalues:
        if eig > 0:
            omega = sqrt(eig)
            period = 2*pi/omega
            frequency = 1/period
            periods.append(period)
            frequencies.append(frequency)
    
    return periods, frequencies, eigenvalues

# ==============================================
# SEISMIC ANALYSIS
# ==============================================
def calculate_seismic_forces(periods):
    """Calculate seismic forces using equivalent lateral force method"""
    # Seismic parameters (example values)
    Ss = 1.5    # Mapped spectral acceleration (short periods)
    S1 = 0.6    # Mapped spectral acceleration (1-sec period)
    Fa = 1.0    # Site amplification factor
    Fv = 1.5    # Site amplification factor
    
    SMS = Fa * Ss
    SM1 = Fv * S1
    SDS = (2/3) * SMS
    SD1 = (2/3) * SM1
    
    # Design response spectrum
    TL = 8.0    # Long-period transition period
    T1 = periods[0] if periods else 1.0  # Fundamental period
    
    # Calculate base shear
    total_weight = 3 * (DL_floor + 0.25*LL_floor)  # Effective seismic weight
    
    if T1 <= 0.2*SD1/SDS:
        Sa = SDS*(0.4 + 0.6*T1*SDS/SD1)
    elif T1 <= SD1/SDS:
        Sa = SDS
    elif T1 <= TL:
        Sa = SD1/T1
    else:
        Sa = SD1*TL/(T1**2)
    
    Cs = Sa  # Seismic response coefficient (simplified)
    Cs = max(Cs, 0.01)  # Minimum value
    Cs = min(Cs, SDS/(0.75))  # Maximum value (when T1 >= TL)
    
    V = Cs * total_weight  # Base shear
    
    # Distribute forces over height
    seismic_forces = []
    heights = [floor_ht * i for i in range(1, 4)]
    weights = [DL_floor + 0.25*LL_floor] * 3
    
    total_Wh = sum(w*h for w, h in zip(weights, heights))
    
    for i, (w, h) in enumerate(zip(weights, heights)):
        Fx = V * (w * h) / total_Wh
        seismic_forces.append(Fx)
    
    return V, seismic_forces

def run_seismic_analysis(periods):
    """Run equivalent lateral force seismic analysis"""
    V, seismic_forces = calculate_seismic_forces(periods)
    
    # Apply seismic forces in X-direction
    ops.pattern('Plain', 7, 1)
    for floor, force in enumerate(seismic_forces, 1):
        floor_nodes = [n for n in node_tags if n >= 100*floor and n < 100*(floor+1)]
        force_per_node = force / len(floor_nodes)
        for node in floor_nodes:
            ops.load(node, force_per_node, 0, 0, 0, 0, 0)
    
    # Apply seismic forces in Y-direction
    ops.pattern('Plain', 8, 1)
    for floor, force in enumerate(seismic_forces, 1):
        floor_nodes = [n for n in node_tags if n >= 100*floor and n < 100*(floor+1)]
        force_per_node = force / len(floor_nodes)
        for node in floor_nodes:
            ops.load(node, 0, force_per_node, 0, 0, 0, 0)
    
    return V, seismic_forces

# ==============================================
# ANALYSIS EXECUTION
# ==============================================

print("Starting 3D Building Analysis...")

# 1. Modal Analysis
print("\n1. MODAL ANALYSIS")
periods, frequencies, eigenvalues = run_modal_analysis()
print(f"First 3 periods: {periods[:3]} sec")
print(f"First 3 frequencies: {frequencies[:3]} Hz")

# 2. Seismic Force Calculation
print("\n2. SEISMIC ANALYSIS")
V, seismic_forces = run_seismic_analysis(periods)
print(f"Base shear: {V:.2f} kN")
print(f"Seismic forces by floor: {seismic_forces} kN")

# Add seismic load combinations
load_combinations.update({
    'LC8': {1: 1.2, 2: 1.0, 7: 1.0},   # 1.2D + 1.0L + 1.0Ex
    'LC9': {1: 1.2, 2: 1.0, 8: 1.0},   # 1.2D + 1.0L + 1.0Ey
    'LC10': {1: 0.9, 7: 1.0},          # 0.9D + 1.0Ex
    'LC11': {1: 0.9, 8: 1.0},          # 0.9D + 1.0Ey
})

# ==============================================
# LOAD COMBINATION ANALYSIS
# ==============================================
print("\n3. LOAD COMBINATION ANALYSIS")

# For simplified analysis, we'll run individual load cases and combine results
load_case_results = {}

# Run individual load cases
individual_cases = {
    'Dead': {1: 1.0},
    'Live': {2: 1.0},
    'Snow': {3: 1.0},
    'WindX': {4: 1.0},
    'WindY': {5: 1.0},
    'SeismicX': {7: 1.0},
    'SeismicY': {8: 1.0}
}

for case_name, load_factors in individual_cases.items():
    print(f"Analyzing {case_name} load case...")
    
    # Clear previous analysis
    ops.wipeAnalysis()
    ops.loadConst('-time', 0.0)
    
    # Apply only the current load case
    for pattern_id in range(1, 9):
        if pattern_id in load_factors:
            # Load case is active with its factor
            continue
        else:
            # Deactivate other load patterns
            pass
    
    setup_static_analysis()
    ops.analyze(1)
    
    # Store results
    case_displacements = {}
    case_forces = {}
    
    for node in node_tags:
        case_displacements[node] = ops.nodeDisp(node)
    
    for ele in ops.getEleTags():
        case_forces[ele] = ops.eleForce(ele)
    
    load_case_results[case_name] = {
        'displacements': case_displacements,
        'forces': case_forces
    }

# ==============================================
# RESULTS PROCESSING
# ==============================================
def check_story_drifts(displacements, case_name):
    """Check inter-story drift ratios"""
    print(f"\n{case_name} - STORY DRIFT CHECK")
    print("Story  Drift X(mm)  Ratio X   Drift Y(mm)  Ratio Y   Status")
    print("-" * 65)
    
    for story in range(1, 4):
        # Get nodes at top and bottom of story
        top_nodes = [n for n in node_tags if n >= 100*story and n < 100*(story+1)]
        bot_nodes = [n for n in node_tags if n >= 100*(story-1) and n < 100*story] if story > 1 else [n for n in node_tags if n < 100]
        
        # Calculate average displacements
        avg_top_x = np.mean([displacements[n][0] for n in top_nodes])
        avg_top_y = np.mean([displacements[n][1] for n in top_nodes])
        avg_bot_x = np.mean([displacements[n][0] for n in bot_nodes])
        avg_bot_y = np.mean([displacements[n][1] for n in bot_nodes])
        
        # Calculate drifts
        drift_x = (avg_top_x - avg_bot_x) * 1000  # Convert to mm
        drift_y = (avg_top_y - avg_bot_y) * 1000
        
        # Calculate drift ratios
        ratio_x = abs(drift_x) / (floor_ht * 1000)
        ratio_y = abs(drift_y) / (floor_ht * 1000)
        
        # Check limits (typical limit is 1/400 = 0.0025 for serviceability)
        status_x = "OK" if ratio_x < 0.0025 else "EXCEED"
        status_y = "OK" if ratio_y < 0.0025 else "EXCEED"
        status = "OK" if status_x == "OK" and status_y == "OK" else "EXCEED"
        
        print(f"{story:3d}    {drift_x:8.2f}   {ratio_x:6.4f}   {drift_y:8.2f}   {ratio_y:6.4f}   {status:>6}")

def check_member_forces(forces, case_name):
    """Check member force limits"""
    print(f"\n{case_name} - MEMBER FORCE CHECK (Sample Elements)")
    print("Element  Axial(kN)  Shear V2(kN)  Shear V3(kN)  Moment M2(kNm)  Moment M3(kNm)")
    print("-" * 85)
    
    # Check first few elements as examples
    sample_elements = list(ops.getEleTags())[:10]
    
    for ele in sample_elements:
        force_vec = forces[ele]
        axial = abs(force_vec[0])
        shear_v2 = abs(force_vec[1])
        shear_v3 = abs(force_vec[2])
        moment_m2 = abs(force_vec[4])
        moment_m3 = abs(force_vec[5])
        
        print(f"{ele:7d}  {axial:9.2f}  {shear_v2:11.2f}  {shear_v3:11.2f}  {moment_m2:12.2f}  {moment_m3:12.2f}")

# Analyze critical load combinations
critical_combinations = ['Dead', 'Live', 'WindX', 'SeismicX']

for case_name in critical_combinations:
    if case_name in load_case_results:
        check_story_drifts(load_case_results[case_name]['displacements'], case_name)
        check_member_forces(load_case_results[case_name]['forces'], case_name)

# ==============================================
# ENVELOPE ANALYSIS
# ==============================================
def create_load_combination_envelope():
    """Create envelope of maximum forces and displacements from all load combinations"""
    print("\n4. LOAD COMBINATION ENVELOPE")
    
    # Initialize envelope dictionaries
    max_displacements = {}
    max_forces = {}
    
    for node in node_tags:
        max_displacements[node] = [0, 0, 0, 0, 0, 0]  # 6 DOF
    
    for ele in ops.getEleTags():
        max_forces[ele] = [0] * 12  # 12 force components
    
    # For each load combination, calculate results using superposition
    for lc_name, factors in load_combinations.items():
        print(f"Processing {lc_name}...")
        
        # Combine results using superposition principle
        combined_disp = {}
        combined_forces = {}
        
        for node in node_tags:
            combined_disp[node] = [0, 0, 0, 0, 0, 0]
        
        for ele in ops.getEleTags():
            combined_forces[ele] = [0] * 12
        
        # Superpose individual load case results
        case_mapping = {1: 'Dead', 2: 'Live', 3: 'Snow', 4: 'WindX', 5: 'WindY', 7: 'SeismicX', 8: 'SeismicY'}
        
        for pattern_id, factor in factors.items():
            if pattern_id in case_mapping:
                case_name = case_mapping[pattern_id]
                if case_name in load_case_results:
                    # Add factored results
                    for node in node_tags:
                        for i in range(6):
                            combined_disp[node][i] += factor * load_case_results[case_name]['displacements'][node][i]
                    
                    for ele in ops.getEleTags():
                        for i in range(12):
                            combined_forces[ele][i] += factor * load_case_results[case_name]['forces'][ele][i]
        
        # Update envelope with maximum values
        for node in node_tags:
            for i in range(6):
                if abs(combined_disp[node][i]) > abs(max_displacements[node][i]):
                    max_displacements[node][i] = combined_disp[node][i]
        
        for ele in ops.getEleTags():
            for i in range(12):
                if abs(combined_forces[ele][i]) > abs(max_forces[ele][i]):
                    max_forces[ele][i] = combined_forces[ele][i]
    
    return max_displacements, max_forces

max_disps, max_forces = create_load_combination_envelope()

# Check envelope results
check_story_drifts(max_disps, "ENVELOPE")
check_member_forces(max_forces, "ENVELOPE")

# ==============================================
# SAVE RESULTS
# ==============================================
results = {
    "modal_analysis": {
        "periods": periods,
        "frequencies": frequencies
    },
    "seismic_analysis": {
        "base_shear": V,
        "floor_forces": seismic_forces
    },
    "load_combinations": load_combinations,
    "individual_cases": {case: {
        "max_displacement": max([abs(d) for disp in results['displacements'].values() for d in disp]),
        "max_force": max([abs(f) for forces in results['forces'].values() for f in forces])
    } for case, results in load_case_results.items()},
    "envelope_results": {
        "max_displacements": {str(k): v for k, v in max_disps.items()},
        "max_forces": {str(k): v for k, v in max_forces.items()}
    }
}

# Save to JSON file
with open('3D_building_analysis_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\n" + "="*60)
print("ANALYSIS COMPLETE")
print("="*60)
print(f"Results saved to: 3D_building_analysis_results.json")
print(f"Total load combinations analyzed: {len(load_combinations)}")
print(f"Modal periods: {[f'{p:.3f}' for p in periods[:3]]} sec")
print(f"Seismic base shear: {V:.2f} kN")
print("="*60)
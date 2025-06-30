from import_ import *     
from input import *
from load_combinations import *
from Grid_and_structure_creation import *
from wall_meshing import *
from sections_function import *
from units import *
import matplotlib.pyplot as plt
import opsvis as opsv
from ipywidgets import widgets
from opsvis.model import get_Ew_data_from_ops_domain_3d
from opsvis.secforces import section_force_distribution_3d
from typing import Dict, List, Tuple
'''
do not change the filepaths. make separate the json files creating and plotting creating.
1. create json filepath for beam forces, stresses, strains by using ops.eleResponse(ele_tag, 'forces'), ops.eleResponse(ele_tag, 'stresses'), ops.eleResponse(ele_tag, 'strains') . 2. create json filepath for shell elements forces, stresses, strains by using ops.eleResponse(ele_tag, 'forces'), ops.eleResponse(ele_tag, 'stresses'), ops.eleResponse(ele_tag, 'strains'). remove the repeted codings. and try to make these simple. 3. also calculate the deflection and relative deflection for beams beams and shell elements.
4. make separate the json files creating and plotting creating.
'''



# Extract Element Results
def structural_model_plot(OUTPUT_FOLDER="postprocessing_folder", load_combination="combo"):
    # Define the main output directory and subdirectories
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    
    # Create subdirectories
    JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
    JSON_FOLDER = os.path.join(JSON_FOLDER, load_combination)
    IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
    IMAGE_FOLDER = os.path.join(IMAGE_FOLDER, load_combination)
    STRUCTURAL_MODEL_FOLDER = os.path.join(IMAGE_FOLDER, "structural_model")
    
    os.makedirs(JSON_FOLDER, exist_ok=True)
    os.makedirs(IMAGE_FOLDER, exist_ok=True)
    os.makedirs(STRUCTURAL_MODEL_FOLDER, exist_ok=True)
    
    # Model plot
    fig= plt.figure(figsize=(10, 8))
    # Get all shell elements (in this case, just element 1)
    ax = fig.add_subplot(111, projection='3d')

    # Get all shell elements (in this case, just element 1)
    shell_elements = ops.getEleTags()  # Returns [1]

    for ele_tag in shell_elements:
        # Get node coordinates of the element
        ele_nodes = ops.eleNodes(ele_tag)
        node_coords = np.array([ops.nodeCoord(node) for node in ele_nodes])
        
        # Create a filled polygon
        poly = Poly3DCollection([node_coords], alpha=0.5, linewidth=1, edgecolor='k')
        
        # Assign color (modify logic as needed)
        poly.set_facecolor('yellow')  # Single color for all elements
        
        ax.add_collection3d(poly)

    # Overlay the original model edges (optional)
    opsv.plot_model(element_labels=0, node_labels=0, ax=ax, fmt_model={'color': 'k', 'linewidth': 1})
    # opsv.plot_model()
    plt.title(f"Model - {load_combination}")
    filepath = os.path.join(STRUCTURAL_MODEL_FOLDER, f"model_{load_combination}.png")
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    # print(f"Saved: {filepath}")
    plt.close()
    
    # Deformation plot
    plt.figure(figsize=(10, 8))
    opsv.plot_model()
    plt.title(f"Deformation - {load_combination}")
    filepath = os.path.join(STRUCTURAL_MODEL_FOLDER, f"deformation_{load_combination}.png")
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    # print(f"Saved: {filepath}")
    plt.close()
    
    # Load plot
    plt.figure(figsize=(10, 8))
    opsv.plot_load()
    plt.title(f"Load - {load_combination}")
    filepath = os.path.join(STRUCTURAL_MODEL_FOLDER, f"load_{load_combination}.png")
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    # print(f"Saved: {filepath}")
    plt.close()
    
    # Extruded shapes plot

    element_data_path = os.path.join("output_folder", "json_files", "element_data.json")

    with open(element_data_path, 'r') as f:
        data = json.load(f)
    ele_shapes = {}

    for ele in data["elements"]:
        tag = ele["eleTag"]
        sec_tag = ele["secTag"]

        for section_type in section_definitions:
            sections = section_definitions[section_type]
            for section_name, section_data in sections.items():
                if section_data["section_tag"] != sec_tag:
                    continue

                if section_data["type"] == "rectangular":
                    B = section_data["B"]
                    H = section_data["H"]
                    ele_shapes[tag] = ["rect", [B, H]]

                elif section_data["type"] == "circular":
                    D = section_data["D_Sec"]
                    ele_shapes[tag] = ["circ", [D]]

                elif section_data["type"] == "L":
                    H = section_data["H"]
                    B = section_data["B"]
                    t = section_data["t"]
                    ele_shapes[tag] = ["L", [H, B, t]]

                elif section_data["type"] == "I":
                    B = section_data["B"]
                    H = section_data["H"]
                    tf = section_data.get("tf", B / 10.)
                    tw = section_data.get("tw", H / 6.)
                    ele_shapes[tag] = ["I", [B, H, tf, tw]]
                    
    plt.figure(figsize=(10, 8))
    opsv.plot_extruded_shapes_3d(ele_shapes)
    plt.title(f"Extruded Shapes - {load_combination}")
    filepath = os.path.join(STRUCTURAL_MODEL_FOLDER, f"extruded_shapes_{load_combination}.png")
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    # print(f"Saved: {filepath}")
    plt.close()
    
    print(f"All plots saved to {STRUCTURAL_MODEL_FOLDER}")


def apply_3d_beam_loads(eleTags, Wy=0.0, Wz=0.0, Wx=0.0, Py=0.0, Pz=0.0, Px=0.0, xL=0.5, pattern_tag=1001, time_series_tag=1):
    ops.pattern("Plain", pattern_tag, time_series_tag)
    for eleTag in eleTags:
        if Wy != 0.0 or Wz != 0.0 or Wx != 0.0:
            ops.eleLoad("-ele", eleTag, "-type", "-beamUniform", Wy, Wz, Wx)
        if Py != 0.0 or Pz != 0.0 or Px != 0.0:
            ops.eleLoad("-ele", eleTag, "-type", "-beamPoint", Py, Pz, xL, Px)


def get_node_results(OUTPUT_FOLDER = "postprocessing_folder", load_combination="combo"):
    node_tags = ops.getNodeTags()
    node_data = []
    print(f"Total nodes: {len(node_tags)}")
    for node_tag in node_tags:
        # Get node coordinates
        coords = ops.nodeCoord(node_tag)
        
        # Get displacements
        disp = ops.nodeDisp(node_tag)
        # print(f"Node {node_tag}: disp = {disp}")       
        # Get reactions (if available)
        reaction = ops.nodeReaction(node_tag)
        # print(f"Node {node_tag}: Reaction = {reaction}")
        node_data.append({
            'id': node_tag,
            'coordinates': coords,
            'displacements': disp,
            'reactions': reaction,
            'load_combination': load_combination
        })
    
    OUTPUT_FOLDER = "postprocessing_folder"  # Main output directory
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist

    # Subdirectories
    # Change to:
    JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
    JSON_FOLDER = os.path.join(JSON_FOLDER, load_combination)
    output_file = os.path.join(JSON_FOLDER, f"nodes_results_{load_combination}.json")
    os.makedirs(JSON_FOLDER, exist_ok=True)
    # Save the data to JSON file
    with open(output_file, 'w') as f:
        json.dump(node_data, f, indent=4)
    
    print(f"Node results saved to: {output_file}")
    return node_data

def deflection(disp_i, disp_j, direction, x, L):
    """
    Calculate deflection at position x along beam using linear interpolation.
    
    Parameters:
    -----------
    disp_i : list/array
        Displacement at node i [ux, uy, uz, rx, ry, rz]
    disp_j : list/array
        Displacement at node j [ux, uy, uz, rx, ry, rz]
    direction : str
        Direction of deflection ('dx', 'dy', 'dz')
    x : float
        Position along beam from node i
    L : float
        Total length of beam
    """
    direction_map = {'dx': 0, 'dy': 1, 'dz': 2}
    idx = direction_map[direction]
    
    # Linear interpolation between nodes
    return disp_i[idx] + (disp_j[idx] - disp_i[idx]) * (x / L)


def rel_deflection(disp_i, disp_j, direction, x, L):
    """
    Calculate relative deflection (deflection relative to chord between end points).
    
    Parameters:
    -----------
    disp_i : list/array
        Displacement at node i [ux, uy, uz, rx, ry, rz]
    disp_j : list/array
        Displacement at node j [ux, uy, uz, rx, ry, rz]
    direction : str
        Direction of deflection ('dy', 'dz')
    x : float
        Position along beam from node i
    L : float
        Total length of beam
    """
    direction_map = {'dy': 1, 'dz': 2}
    idx = direction_map[direction]
    
    # Deflection at position x
    deflection_x = disp_i[idx] + (disp_j[idx] - disp_i[idx]) * (x / L)
    
    # Chord deflection (linear between end points)
    chord_deflection = disp_i[idx] + (disp_j[idx] - disp_i[idx]) * (x / L)
    
    # For relative deflection, we need to consider the beam curvature
    # This is a simplified approach - in reality, you'd need the actual beam equation
    # For now, we'll calculate relative to the straight line between end points
    return deflection_x - chord_deflection


def extract_element_results(JSON_FOLDER, nep=5):
    """
    Extract forces, stresses, and strains for all elements (beams and shells)
    
    Args:
        nep: Number of evaluation points along beam elements
        
    Returns:
        Dictionary containing all element results organized by element type
    """
    # Initialize output dictionary
    results = {
        'beam_elements': {},
        'shell_elements': {}
    }
    
    # Get all element tags
    ele_tags = ops.getEleTags()
    
    # Process each element
    for ele_tag in ele_tags:
        ele_type = ops.eleType(ele_tag)
        
        if 'ForceBeamColumn' in ele_type or 'beam' in ele_type.lower():
            # Process beam element
            results['beam_elements'][ele_tag] = extract_beam_results(ele_tag, nep, JSON_FOLDER)
            
        elif 'Shell' in ele_type or 'shell' in ele_type.lower():
            # Process shell element
            results['shell_elements'][ele_tag] = extract_shell_results(ele_tag)
            
        else:
            print(f"Warning: Unsupported element type {ele_type} for element {ele_tag}")
    
    return results


# def extract_beam_results(ele_tag, nep, JSON_FOLDER):
#     """
#     Extract forces, stresses, and strains for a beam element
    
#     Args:
#         ele_tag: Element tag
#         nep: Number of evaluation points
#         element_data: Dictionary containing pre-loaded element properties
        
#     Returns:
#         Dictionary containing beam results
        
#     Raises:
#         ValueError: If required section properties are missing
#     """
#     json_file_path = os.path.join(JSON_FOLDER, "element_data.json")
#     with open(json_file_path, 'r') as f:
#         element_data = json.load(f)
#     # Find the element in the pre-loaded data
#     element = None
#     for el in element_data["elements"]:
#         if el["eleTag"] == ele_tag:
#             element = el
#             break
    
#     if element is None:
#         raise ValueError(f"Element with tag {ele_tag} not found in element data")
    
#     # Get required section properties (no defaults)
#     try:
#         A = element["area"]
#         Iy = element["Iy"]
#         Iz = element["Iz"]
#         B = element["B"]  # Width
#         H = element["H"]  # Depth
#     except KeyError as e:
#         raise ValueError(f"Missing required section property: {str(e)}")
    
#     # Get element geometry
#     node_tags = ops.eleNodes(ele_tag)
#     ecrd = np.array([ops.nodeCoord(node_tags[0]), ops.nodeCoord(node_tags[1])])
#     L = np.linalg.norm(ecrd[1] - ecrd[0])
    
#     # Get element forces
#     forces = ops.eleResponse(ele_tag, 'localForces')
    
#     # Process distributed loads
#     Ew = get_Ew_data_from_ops_domain_3d()
#     eload_data = Ew.get(ele_tag, [['-beamUniform', 0., 0., 0.]])
    
#     # Get force distribution
#     s, xl, _ = section_force_distribution_3d(ecrd, forces, nep, eload_data)
    
#     # Organize forces
#     force_results = []
#     for i in range(len(xl)):
#         force_results.append({
#             'position': float(xl[i]),
#             'N': float(s[i,0]),
#             'Vy': float(s[i,1]),
#             'Vz': float(s[i,2]),
#             'T': float(s[i,3]),
#             'My': float(s[i,4]),
#             'Mz': float(s[i,5])
#         })
    
#     # Calculate torsional constant for rectangular section
#     if B > 0 and H > 0:
#         if B == H:  # Square section
#             J = 0.141 * B * H**3
#         else:  # Rectangular section
#             J = (B * H**3) * (1/3 - 0.21*(H/B)*(1 - (H**4)/(12*B**4)))
#     else:
#         raise ValueError("Section dimensions B and H must be positive")
    
#     # Material properties (steel defaults)
#     E = 2.1e11  # Pa
#     G = 8.1e10  # Pa
    
#     # Calculate stresses and strains at each point
#     stress_results = []
#     strain_results = []
    
#     for i in range(len(xl)):
#         # Current forces
#         N = s[i,0]
#         Vy = s[i,1]
#         Vz = s[i,2]
#         T = s[i,3]
#         My = s[i,4]
#         Mz = s[i,5]
        
#         # Calculate stresses
#         σ_axial = N/A
#         σ_bending_y = My * (H/2) / Iy
#         σ_bending_z = Mz * (B/2) / Iz
#         τ_torsion = T * (H/2) / J
        
#         # Shear stresses
#         Qy = A*H/8  # First moment of area
#         Qz = A*B/8
#         τ_shear_y = Vy*Qz/(Iz*B)
#         τ_shear_z = Vz*Qy/(Iy*H)
        
#         # Von Mises stress
#         σ_vm = np.sqrt((σ_axial + σ_bending_y + σ_bending_z)**2 + 
#                      3*(max(τ_shear_y, τ_shear_z, τ_torsion))**2)
        
#         # Calculate strains
#         ε_axial = σ_axial/E
#         ε_bending_y = σ_bending_y/E
#         ε_bending_z = σ_bending_z/E
#         γ_shear_y = τ_shear_y/G
#         γ_shear_z = τ_shear_z/G
#         γ_torsion = τ_torsion/G
        
#         # Store results
#         stress_results.append({
#             'position': float(xl[i]),
#             'σ_axial': float(σ_axial),
#             'σ_bending_y': float(σ_bending_y),
#             'σ_bending_z': float(σ_bending_z),
#             'τ_shear_y': float(τ_shear_y),
#             'τ_shear_z': float(τ_shear_z),
#             'τ_torsion': float(τ_torsion),
#             'von_mises': float(σ_vm)
#         })
        
#         strain_results.append({
#             'position': float(xl[i]),
#             'ε_axial': float(ε_axial),
#             'ε_bending_y': float(ε_bending_y),
#             'ε_bending_z': float(ε_bending_z),
#             'γ_shear_y': float(γ_shear_y),
#             'γ_shear_z': float(γ_shear_z),
#             'γ_torsion': float(γ_torsion)
#         })
    
#     return {
#         'length': float(L),
#         'forces': force_results,
#         'stresses': stress_results,
#         'strains': strain_results
#     }

def extract_beam_results(ele_tag, nep, JSON_FOLDER):
    """
    Extract forces, stresses, and strains for a beam element
    
    Args:
        ele_tag: Element tag
        nep: Number of evaluation points
        JSON_FOLDER: Path to folder containing element_data.json
        
    Returns:
        Dictionary containing beam results
        
    Raises:
        ValueError: If required section properties are missing
    """
    json_file_path = os.path.join(JSON_FOLDER, "element_data.json")
    with open(json_file_path, 'r') as f:
        element_data = json.load(f)
    
    # Find the element in the pre-loaded data
    element = None
    for el in element_data["elements"]:
        if el["eleTag"] == ele_tag:
            element = el
            break
    
    if element is None:
        raise ValueError(f"Element with tag {ele_tag} not found in element data")
    
    # Get required section properties (no defaults)
    try:
        A = element["area"]
        Iy = element["Iy"]
        Iz = element["Iz"]
        section_type = element["type"]
    except KeyError as e:
        raise ValueError(f"Missing required section property: {str(e)}")
    
    # Handle section-specific properties
    if section_type == "rectangular":
        B = element["B"]
        H = element["H"]
    elif section_type == "circular":
        B = element["B"]
        H = element["H"]
        D = B
    elif section_type == "L":
        B = element.get("B")
        H = element.get("H")
        t = element.get("t")
        # print(f"Debug L: section_type={section_type}, B={B}, H={H}, t={t}, element={element}")

    else:
        raise ValueError(f"Unsupported section type: {section_type}")
    
    # Validate dimensions
    if B is None or H is None or B <= 0 or H <= 0:
        raise ValueError(f"Invalid section dimensions for element {ele_tag}: B={B}, H={H}")
    
    # Get element geometry
    node_tags = ops.eleNodes(ele_tag)
    ecrd = np.array([ops.nodeCoord(node_tags[0]), ops.nodeCoord(node_tags[1])])
    L = np.linalg.norm(ecrd[1] - ecrd[0])
    
    # Get element forces
    forces = ops.eleResponse(ele_tag, 'localForces')
    
    # Process distributed loads
    Ew = get_Ew_data_from_ops_domain_3d()
    eload_data = Ew.get(ele_tag, [['-beamUniform', 0., 0., 0.]])
    
    # Get force distribution
    s, xl, _ = section_force_distribution_3d(ecrd, forces, nep, eload_data)
    
    # Organize forces
    force_results = []
    for i in range(len(xl)):
        force_results.append({
            'position': float(xl[i]),
            'N': float(s[i,0]),
            'Vy': float(s[i,1]),
            'Vz': float(s[i,2]),
            'T': float(s[i,3]),
            'My': float(s[i,4]),
            'Mz': float(s[i,5])
        })
    
    # Calculate torsional constant based on section type
    if section_type == "rectangular":
        if B == H:  # Square section
            J = 0.141 * B * H**3
        else:  # Rectangular section
            J = (B * H**3) * (1/3 - 0.21 * (H/B) * (1 - (H**4) / (12 * B**4)))

    elif section_type == "circular":
        J = 2 * Iz  # For circular sections, J ≈ 2*I
    elif section_type == "L":
        # Approximate torsion constant for L-section (Roark's formula)
        a = max(B, H)
        b = min(B, H)
        t_val = t
        J = (1/3) * (a * t_val**3 + b * t_val**3)
    
    # Material properties (steel defaults)
    E = 2.1e11  # Pa (Young's modulus)
    G = 8.1e10  # Pa (Shear modulus)
    
    # Calculate stresses and strains at each point
    stress_results = []
    strain_results = []
    
    for i in range(len(xl)):
        # Current forces
        N = s[i,0]
        Vy = s[i,1]
        Vz = s[i,2]
        T = s[i,3]
        My = s[i,4]
        Mz = s[i,5]
        
        # Calculate stresses
        σ_axial = N/A
        
        # Bending stresses depend on section type
        if section_type == "rectangular":
            σ_bending_y = My * (H/2) / Iy
            σ_bending_z = Mz * (B/2) / Iz
        elif section_type == "circular":
            σ_bending_y = My * (D/2) / Iy
            σ_bending_z = Mz * (D/2) / Iz
        elif section_type == "L":
            # For L-section, use extreme fiber distances
            σ_bending_y = My * (H/2) / Iy
            σ_bending_z = Mz * (B/2) / Iz
        
        # Torsional stress
        if section_type == "rectangular":
            τ_torsion = T * (H/2) / J
        elif section_type == "circular":
            τ_torsion = T * (D/2) / J
        elif section_type == "L":
            τ_torsion = T * t / J  # Approximate for thin-walled
        
        # Shear stresses
        if section_type == "rectangular":
            Qy = A*H/8  # First moment of area
            Qz = A*B/8
            τ_shear_y = Vy*Qz/(Iz*B)
            τ_shear_z = Vz*Qy/(Iy*H)
        elif section_type == "circular":
            τ_shear_y = 4/3 * Vy/A  # Average shear stress for circular
            τ_shear_z = 4/3 * Vz/A
        elif section_type == "L":
            # Simplified shear stress for L-sections
            τ_shear_y = Vy/A
            τ_shear_z = Vz/A
        
        # Von Mises stress
        σ_vm = np.sqrt((σ_axial + σ_bending_y + σ_bending_z)**2 + 
                      3*(max(τ_shear_y, τ_shear_z, τ_torsion))**2)
        
        # Calculate strains
        ε_axial = σ_axial/E
        ε_bending_y = σ_bending_y/E
        ε_bending_z = σ_bending_z/E
        γ_shear_y = τ_shear_y/G
        γ_shear_z = τ_shear_z/G
        γ_torsion = τ_torsion/G
        
        # Store results
        stress_results.append({
            'position': float(xl[i]),
            'σ_axial': float(σ_axial),
            'σ_bending_y': float(σ_bending_y),
            'σ_bending_z': float(σ_bending_z),
            'τ_shear_y': float(τ_shear_y),
            'τ_shear_z': float(τ_shear_z),
            'τ_torsion': float(τ_torsion),
            'von_mises': float(σ_vm)
        })
        
        strain_results.append({
            'position': float(xl[i]),
            'ε_axial': float(ε_axial),
            'ε_bending_y': float(ε_bending_y),
            'ε_bending_z': float(ε_bending_z),
            'γ_shear_y': float(γ_shear_y),
            'γ_shear_z': float(γ_shear_z),
            'γ_torsion': float(γ_torsion)
        })
    
    return {
        'length': float(L),
        'section_type': section_type,
        'forces': force_results,
        'stresses': stress_results,
        'strains': strain_results
    }

def extract_shell_results(ele_tag):
    """
    Extract forces, stresses, and strains for a shell element
    
    Args:
        ele_tag: Element tag
        
    Returns:
        Dictionary containing shell results
    """
    # Initialize results
    results = {
        'forces': None,
        'stresses': [],
        'strains': []
    }
    
    try:
        # Get element responses
        forces = ops.eleResponse(ele_tag, 'forces')       # Membrane and bending forces
        stresses = ops.eleResponse(ele_tag, 'stresses')   # Stresses at integration points
        strains = ops.eleResponse(ele_tag, 'strains')     # Strains at integration points
        
        # Organize forces (stress resultants)
        if forces:
            # Different shell elements might return different numbers of forces
            force_results = {}
            force_components = ['Nxx', 'Nyy', 'Nxy', 'Mxx', 'Myy', 'Mxy', 'Qxz', 'Qyz']
            
            for i, component in enumerate(force_components):
                if i < len(forces):
                    force_results[component] = float(forces[i])
                else:
                    force_results[component] = 0.0  # Default value if not available
            
            results['forces'] = force_results
        
        # Organize stresses (typically at integration points through thickness)
        if stresses:
            # Stresses are usually reported for each integration point
            # Format can vary - we'll handle different cases
            num_stress_components = len(stresses)
            num_integration_points = num_stress_components // 5  # Most common case
            
            for ip in range(num_integration_points):
                start_idx = ip * 5
                end_idx = start_idx + 5
                
                if end_idx <= num_stress_components:
                    stress_data = {
                        'integration_point': ip + 1,
                        'σ_xx': float(stresses[start_idx]),
                        'σ_yy': float(stresses[start_idx + 1]),
                        'τ_xy': float(stresses[start_idx + 2]),
                        'τ_xz': float(stresses[start_idx + 3]) if (start_idx + 3) < num_stress_components else 0.0,
                        'τ_yz': float(stresses[start_idx + 4]) if (start_idx + 4) < num_stress_components else 0.0
                    }
                    results['stresses'].append(stress_data)
        
        # Organize strains
        if strains:
            num_strain_components = len(strains)
            num_integration_points = num_strain_components // 5  # Most common case
            
            for ip in range(num_integration_points):
                start_idx = ip * 5
                end_idx = start_idx + 5
                
                if end_idx <= num_strain_components:
                    strain_data = {
                        'integration_point': ip + 1,
                        'ε_xx': float(strains[start_idx]),
                        'ε_yy': float(strains[start_idx + 1]),
                        'γ_xy': float(strains[start_idx + 2]),
                        'γ_xz': float(strains[start_idx + 3]) if (start_idx + 3) < num_strain_components else 0.0,
                        'γ_yz': float(strains[start_idx + 4]) if (start_idx + 4) < num_strain_components else 0.0
                    }
                    results['strains'].append(strain_data)
                    
    except Exception as e:
        print(f"Error processing shell element {ele_tag}: {str(e)}")
        # Return partial results if available
        pass
    
    return results


def extract_all_element_data(JSON_FOLDER, nep=5, output_folder="postprocessing_folder", load_combination="combo2"):
    """
    Extract all beam and shell element results from extract_element_results function
    and save them into JSON files (one file per result type).

    Args:
        nep: Number of evaluation points for beam elements.
        output_folder: Folder to store output JSON files.

    Returns:
        Dict with beam and shell results split into forces, stresses, and strains.
    """
    full_results = extract_element_results(JSON_FOLDER,nep)

    beam_forces = {}
    beam_stresses = {}
    beam_strains = {}
    shell_forces = {}
    shell_stresses = {}
    shell_strains = {}

    for ele_tag, result in full_results['beam_elements'].items():
        beam_forces[ele_tag] = result['forces']
        beam_stresses[ele_tag] = result['stresses']
        beam_strains[ele_tag] = result['strains']

    for ele_tag, result in full_results['shell_elements'].items():
        shell_forces[ele_tag] = result['forces']
        shell_stresses[ele_tag] = result['stresses']
        shell_strains[ele_tag] = result['strains']

    # Define output paths
    # Change to:
    JSON_FOLDER = os.path.join(output_folder, "json_files")
    JSON_FOLDER = os.path.join(JSON_FOLDER, load_combination)
    os.makedirs(JSON_FOLDER, exist_ok=True)

    beam_forces_file = os.path.join(JSON_FOLDER, "beam_forces.json")
    beam_stresses_file = os.path.join(JSON_FOLDER, "beam_stress.json")
    beam_strains_file = os.path.join(JSON_FOLDER, "beam_strains.json")

    shell_forces_file = os.path.join(JSON_FOLDER, "shell_forces.json")
    shell_stresses_file = os.path.join(JSON_FOLDER, "shell_stress.json")
    shell_strains_file = os.path.join(JSON_FOLDER, "shell_strains.json")

    # Save all data into corresponding files
    with open(beam_forces_file, 'w') as f:
        json.dump(beam_forces, f, indent=2)
    with open(beam_stresses_file, 'w') as f:
        json.dump(beam_stresses, f, indent=2)
    with open(beam_strains_file, 'w') as f:
        json.dump(beam_strains, f, indent=2)

    with open(shell_forces_file, 'w') as f:
        json.dump(shell_forces, f, indent=2)
    with open(shell_stresses_file, 'w') as f:
        json.dump(shell_stresses, f, indent=2)
    with open(shell_strains_file, 'w') as f:
        json.dump(shell_strains, f, indent=2)

    return {
        'beam_forces': beam_forces,
        'beam_stresses': beam_stresses,
        'beam_strains': beam_strains,
        'shell_forces': shell_forces,
        'shell_stresses': shell_stresses,
        'shell_strains': shell_strains
    }


def calculate_beam_stresses_strains(JSON_FOLDER, nep=10):
    """
    Calculate beam stresses and strains from element data file
    
    Args:
        JSON_FOLDER: Path to folder containing element_data.json
        nep: Number of evaluation points along each beam
        
    Returns:
        Dictionary of results for all beam elements in format:
        {
            element_tag: {
                'length': element length,
                'forces': force distribution,
                'stresses': calculated stresses,
                'strains': calculated strains
            },
            ...
        }
    """
    # 1. Load element data from JSON file
    json_path = os.path.join(JSON_FOLDER, "element_data.json")
    
    try:
        with open(json_path) as f:
            data = json.load(f)
            element_data = data["elements"]
    except FileNotFoundError:
        raise FileNotFoundError(f"Element data file not found at: {json_path}")
    except KeyError:
        raise KeyError("'elements' key not found in the JSON file")
    except Exception as e:
        raise Exception(f"Error loading element data: {str(e)}")
    
    # 2. Initialize results dictionary
    results = {}
    
    # 3. Process all beam elements
    for element in element_data:
        if element["type"].lower() != "beam":
            continue
            
        ele_tag = element["eleTag"]
        
        try:
            # 4. Get element properties
            A = element.get("area")
            Iy = element.get("Iy")
            Iz = element.get("Iz")
            J = element.get("J")
            B = element.get("B")
            H = element.get("H")
            Qy = A*H/8 if A > 0 and H > 0 else 0
            Qz = A*B/8 if A > 0 and B > 0 else 0
            thickness = min(B, H)/10 if min(B, H) > 0 else 0.01
            
            # 5. Get element geometry and forces
            node_tags = ops.eleNodes(ele_tag)
            ecrd = np.array([ops.nodeCoord(node_tags[0]), ops.nodeCoord(node_tags[1])])
            L = np.linalg.norm(ecrd[1] - ecrd[0])
            forces = ops.eleResponse(ele_tag, 'forces')
            
            # 6. Process loads including self-weight
            Ew = get_Ew_data_from_ops_domain_3d()
            eload_data = Ew.get(ele_tag, [['-beamUniform', 0., 0., 0.]])
            if element.get("unit_weight", 0) > 0:
                eload_data.append(['-beamUniform', 0, -element["unit_weight"], 0])
            
            # 7. Get force distribution
            s, xl, _ = section_force_distribution_3d(ecrd, forces, nep, eload_data)
            
            # 8. Get material properties
            try:
                secTag = element.get("secTag")
                E = ops.sectionProperty(secTag, 'E') if secTag else 2.1e11
                G = ops.sectionProperty(secTag, 'G') if secTag else E/2.6
            except:
                E, G = 2.1e11, 8.1e10  # Default steel properties
            
            # 9. Initialize results for this element
            force_results = []
            stress_results = []
            strain_results = []
            
            # 10. Calculate at each evaluation point
            for i in range(len(xl)):
                # Current forces
                N, Vy, Vz, T, My, Mz = [float(s[i,j]) for j in range(6)]
                
                # Store forces
                force_results.append({
                    'position': float(xl[i]),
                    'N': N, 'Vy': Vy, 'Vz': Vz, 'T': T, 'My': My, 'Mz': Mz
                })
                
                # Calculate stresses
                σ_axial = N/A if A != 0 else 0
                σ_bending_z = (Mz*H/2)/Iz if Iz != 0 else 0
                σ_bending_y = (My*B/2)/Iy if Iy != 0 else 0
                τ_shear_y = Vy*Qz/(Iz*thickness) if Iz != 0 else 0
                τ_shear_z = Vz*Qy/(Iy*thickness) if Iy != 0 else 0
                τ_torsion = T*H/(2*J) if J != 0 else 0
                
                # Store stresses
                stress_results.append({
                    'position': float(xl[i]),
                    'σ_axial': σ_axial,
                    'σ_bending_y': σ_bending_y,
                    'σ_bending_z': σ_bending_z,
                    'τ_shear_y': τ_shear_y,
                    'τ_shear_z': τ_shear_z,
                    'τ_torsion': τ_torsion,
                    'σ_max': σ_axial + σ_bending_y + σ_bending_z,
                    'σ_min': σ_axial - σ_bending_y - σ_bending_z,
                    'von_mises': np.sqrt((σ_axial + σ_bending_y + σ_bending_z)**2 + 
                                      3*max(τ_shear_y, τ_shear_z, τ_torsion)**2)
                })
                
                # Calculate and store strains
                strain_results.append({
                    'position': float(xl[i]),
                    'ε_axial': σ_axial/E,
                    'ε_bending_y': σ_bending_y/E,
                    'ε_bending_z': σ_bending_z/E,
                    'γ_shear_y': τ_shear_y/G,
                    'γ_shear_z': τ_shear_z/G,
                    'γ_torsion': τ_torsion/G
                })
            
            # 11. Store results for this element
            results[ele_tag] = {
                'length': float(L),
                'forces': force_results,
                'stresses': stress_results,
                'strains': strain_results
            }
            
        except Exception as e:
            print(f"Warning: Could not process element {ele_tag}: {str(e)}")
            continue
    
    return results





def calculate_slab_reinforcement_from_shell_forces(JSON_FOLDER, output_folder="postprocessing_folder", nep=5, load_combination="combo2"):
    """Calculate slab reinforcement from shell forces according to ACI 318 (FPS units)."""
    # Create output directories
    json_folder = "postprocessing_folder/json_files/slab_reinforcement"
    json_folder = os.path.join("postprocessing_folder", "json_files", load_combination, "slab_reinforcement")
    os.makedirs(json_folder, exist_ok=True)

    # Extract shell forces (placeholder - replace with your actual data extraction)
    results = extract_all_element_data(JSON_FOLDER, nep=nep, output_folder=output_folder, load_combination=load_combination)
    shell_forces = results.get("shell_forces", {})

    # Material properties (FPS units - lb, in, psi)
    phi = 0.9  # Strength reduction factor
    fy = 60000  # psi (yield strength of reinforcement)
    fc = 4000  # psi (concrete compressive strength)
    slab_thickness = 6  # in (total slab thickness)
    cover = 0.75  # in (clear cover to reinforcement)

    # Calculate effective depth (assuming #4 bars - 0.5 in diameter)
    d = slab_thickness - cover - 0.25  # in (0.25 = half of #4 bar diameter)

    # Initialize storage
    reinforcement_results = []
    coords = []
    As_x_bottom = []
    As_y_bottom = []
    As_x_top = []
    As_y_top = []

    for ele_tag, force in shell_forces.items():
        if not force:
            continue

        # Convert moments from lb-ft/ft to lb-in/in (typical FPS units)
        Mx = force.get("Mxx", 0) * 12  # lb-ft/ft to lb-in/ft
        My = force.get("Myy", 0) * 12  # lb-ft/ft to lb-in/ft
        Mxy = force.get("Mxy", 0) * 12  # lb-ft/ft to lb-in/ft

        # Calculate reinforcement for both top and bottom layers
        # Bottom reinforcement (positive moments)
        As_x_b, spacing_x_b, bar_size_x_b = calculate_reinforcement_with_spacing(Mx, d, fc, fy)
        As_y_b, spacing_y_b, bar_size_y_b = calculate_reinforcement_with_spacing(My, d, fc, fy)
        
        # Top reinforcement (negative moments)
        As_x_t, spacing_x_t, bar_size_x_t = calculate_reinforcement_with_spacing(-Mx, d, fc, fy)
        As_y_t, spacing_y_t, bar_size_y_t = calculate_reinforcement_with_spacing(-My, d, fc, fy)

        reinforcement_results.append({
            "element_id": ele_tag,
            "moments": {"Mxx": Mx/12, "Myy": My/12, "Mxy": Mxy/12},  # Return as lb-ft/ft
            "reinforcement": {
                "bottom_x": {
                    "As": As_x_b,
                    "spacing": spacing_x_b,
                    "bar_size": bar_size_x_b
                },
                "bottom_y": {
                    "As": As_y_b,
                    "spacing": spacing_y_b,
                    "bar_size": bar_size_y_b
                },
                "top_x": {
                    "As": As_x_t,
                    "spacing": spacing_x_t,
                    "bar_size": bar_size_x_t
                },
                "top_y": {
                    "As": As_y_t,
                    "spacing": spacing_y_t,
                    "bar_size": bar_size_y_t
                }
            }
        })

        coords.append([0.0, 0.0])  # Placeholder for centroid
        As_x_bottom.append(As_x_b)
        As_y_bottom.append(As_y_b)
        As_x_top.append(As_x_t)
        As_y_top.append(As_y_t)

    # Convert to arrays
    coords = np.array(coords)
    As_x_bottom = np.array(As_x_bottom)
    As_y_bottom = np.array(As_y_bottom)
    As_x_top = np.array(As_x_top)
    As_y_top = np.array(As_y_top)

    # Save results
    with open(os.path.join(json_folder, "results.json"), 'w') as f:
        json.dump(reinforcement_results, f, indent=2)

    print("Slab reinforcement results saved.")
    return reinforcement_results

def calculate_reinforcement_with_spacing(Mu: float, d: float, fc: float, fy: float) -> Tuple[float, float, str]:
    """
    Calculate required reinforcement and spacing per ACI 318 (FPS units).
    
    Args:
        Mu: Factored moment per unit width (lb-in/ft)
        d: Effective depth (in)
        fc: Concrete compressive strength (psi)
        fy: Steel yield strength (psi)
        
    Returns:
        Tuple containing:
        - Required reinforcement area (in²/ft)
        - Center-to-center bar spacing (in)
        - Selected bar size
    """
    # Constants
    phi = 0.9  # Strength reduction factor for flexure
    beta1 = 0.85 if fc <= 4000 else max(0.65, 0.85 - 0.05*(fc-4000)/1000)
    b = 12  # Unit width for calculation (12 in = 1 ft)
    
    # Calculate required Rn (moment parameter)
    Rn = Mu / (phi * b * d**2)
    
    # Calculate required steel ratio
    rho = (0.85 * fc / fy) * (1 - math.sqrt(1 - (2 * Rn) / (0.85 * fc)))
    
    # Check minimum reinforcement (ACI 7.6.1.1)
    rho_min = max(0.0018, 3 * math.sqrt(fc) / fy)
    rho = max(rho, rho_min)
    
    # Calculate required steel area (in²/ft)
    As_req = rho * b * d
    
    # Check maximum reinforcement (ACI 21.2.2)
    epsilon_t = 0.005  # Net tensile strain
    rho_max = (0.85 * beta1 * fc / fy) * (epsilon_t / (epsilon_t + 0.002))
    if rho > rho_max:
        raise ValueError("Error: Reinforcement exceeds maximum allowed by ACI")
    
    # Available bar sizes (#3 to #6)
    bar_sizes = {
        '#3': 0.11,  # in²
        # '#4': 0.20,
        # '#5': 0.31,
        # '#6': 0.44
    }
    
    # Find the most economical bar size and spacing
    selected_bar = None
    for bar_size, bar_area in sorted(bar_sizes.items(), key=lambda x: x[1], reverse=True):
        spacing = (bar_area * b) / As_req
        
        # Check minimum spacing requirements (ACI 25.2.1)
        min_spacing = max(1,  # 1 inch minimum
                         {'#3': 0.375, '#4': 0.5, '#5': 0.625, '#6': 0.75}[bar_size],  # Bar diameter
                         1.33 * 0.75)  # 1.33 * max aggregate size (assume 3/4")
        
        if spacing >= min_spacing:
            selected_bar = (bar_size, spacing)
            break
    
    if not selected_bar:
        # If no bar satisfies spacing requirements, use smallest bar at minimum spacing
        bar_size = '#3'
        bar_area = bar_sizes[bar_size]
        min_spacing = max(1, 0.375, 1.33 * 0.75)
        As_provided = (bar_area * b) / min_spacing
        return (As_provided, min_spacing, bar_size)
    
    # Round spacing to nearest 1/2 inch for practicality
    practical_spacing = round(selected_bar[1] * 2) / 2
    
    return (As_req, practical_spacing, selected_bar[0])

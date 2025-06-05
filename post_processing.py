from import_ import *     
from input import *
from load_combinations import merge_structures
from Grid_and_structure_creation import load_structure
from wall_meshing import create_combined_structure_json
from sections_function import *
from units import *
import matplotlib.pyplot as plt
import opsvis as opsv
from ipywidgets import widgets

    
def plot_beam_forces(force_type='N', sfac=1.0, num_points=10, view_angle=None, figsize=(12, 10), load_combination="combo"):
    """
    Plot beam forces with detailed visualization using multiple points along each element
    and save the data to JSON files.
    
    Parameters:
    -----------
    force_type : str
        Type of force to plot ('N' for axial, 'Vy'/'Vz' for shear, 'T' for torsion, 'My'/'Mz' for bending moments)
    sfac : float
        Scale factor for force visualization
    num_points : int
        Number of points to sample along each element (minimum 10)
    view_angle : tuple
        (elevation, azimuth) for 3D plot view angle. If None, default view is used.
    figsize : tuple
        Figure size as (width, height) in inches
    load_combination : str
        Name of the load combination being analyzed (used for file naming)
    """

    OUTPUT_FOLDER = "postprocessing_folder"  # Main output directory
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist

    # Subdirectories
    JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
    IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
    os.makedirs(JSON_FOLDER, exist_ok=True)
    os.makedirs(IMAGE_FOLDER, exist_ok=True)
    # Ensure minimum number of points
    num_points = max(10, num_points)
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')
    
    # Set specific view angle if provided
    if view_angle:
        ax.view_init(elev=view_angle[0], azim=view_angle[1])
    
    beam_elements = [e for e in ops.getEleTags() if 'ForceBeamColumn3d' in ops.eleType(e)]
    
    # Dictionary to map force types to indices and labels
    force_dict = {
        'N': (0, 'Axial Force'),
        'Vy': (1, 'Shear Force (Y)'),
        'Vz': (2, 'Shear Force (Z)'),
        'T': (3, 'Torsion'),
        'My': (4, 'Bending Moment (Y)'),
        'Mz': (5, 'Bending Moment (Z)')
    }
    
    if force_type not in force_dict:
        print(f"Unknown force type: {force_type}")
        return
    
    force_index, force_title = force_dict[force_type]
    
    max_force = 0  # Track maximum force for color scaling
    all_force_data = {}  # Dictionary to store all force data for JSON output
    
    # First pass to get max force value for consistent scaling and collect data
    for ele_tag in beam_elements:
        # Initialize data structure for this element
        element_data = {
            "element_tag": ele_tag,
            "nodes": ops.eleNodes(ele_tag),
            "node_coords": [ops.nodeCoord(n) for n in ops.eleNodes(ele_tag)],
            "force_type": force_type,
            "force_values": [],
            "locations": [],
            "load_combination": load_combination
        }
        
        # Sample forces along the element
        for loc in np.linspace(0, 1, num_points):
            # Get force at specific location (fraction along element)
            section_forces = ops.eleResponse(ele_tag, 'section', [loc*num_points, 'force'])
            if section_forces:
                force_value = section_forces[force_index]
                max_force = max(max_force, abs(force_value))
                
                # Store data for JSON output
                element_data["force_values"].append(float(force_value))
                element_data["locations"].append(float(loc))
            else:
                # If section forces can't be retrieved, use element end forces and interpolate
                end_forces = ops.eleForce(ele_tag)
                if force_type in ['N', 'T']:  # Constant along element
                    force_value = end_forces[force_index]
                elif force_type in ['Vy', 'Vz']:  # Shear is often constant for simple beams
                    force_value = end_forces[force_index]
                elif force_type in ['My', 'Mz']:  # Linear moment variation for uniform load
                    start_moment = end_forces[force_index]
                    end_moment = end_forces[force_index + 6] if len(end_forces) > 6 else -start_moment
                    force_value = start_moment * (1 - loc) + end_moment * loc
                
                max_force = max(max_force, abs(force_value))
                element_data["force_values"].append(float(force_value))
                element_data["locations"].append(float(loc))
        
        # Add element data to the main dictionary
        all_force_data[f"element_{ele_tag}"] = element_data
    
    if max_force == 0:
        max_force = 1.0  # Avoid division by zero
    
    # Save data to JSON file
    # Ensure the directory exists (if JSON_FOLDER is defined)
    os.makedirs(os.path.join(JSON_FOLDER, "member"), exist_ok=True)

    OUTPUT_FOLDER = "postprocessing_folder"  # Main output directory
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist

    # Subdirectories
    

    # Define the JSON file path with load combination
    json_filename = os.path.join(JSON_FOLDER, f"member/beam_forces_{force_type}_{load_combination}.json")
    with open(json_filename, 'w') as json_file:
        json.dump({
            "force_type": force_type,
            "force_title": force_title,
            "scale_factor": sfac,
            "max_force": float(max_force),
            "load_combination": load_combination,
            "elements": all_force_data
        }, json_file, indent=4)
    
    # Color map for force intensity
    cmap = plt.cm.jet
    
    # Now plot everything
    for ele_tag in beam_elements:
        element_data = all_force_data[f"element_{ele_tag}"]
        coords = element_data["node_coords"]
        forces = element_data["force_values"]
        
        start_point = np.array(coords[0])
        end_point = np.array(coords[1])
        
        # Element vector and length
        element_vector = end_point - start_point
        element_length = np.linalg.norm(element_vector)
        
        # Calculate unit vectors for the local coordinate system
        x_local = element_vector / element_length  # Along the beam axis
        
        # Find perpendicular vectors for y_local and z_local
        if np.abs(x_local[2]) < 0.9:  # If not vertical
            y_local = np.cross(x_local, [0, 0, 1])
            y_local = y_local / np.linalg.norm(y_local)
        else:  # If vertical
            y_local = np.cross(x_local, [0, 1, 0])
            y_local = y_local / np.linalg.norm(y_local)
            
        z_local = np.cross(x_local, y_local)
        z_local = z_local / np.linalg.norm(z_local)
        
        # Choose offset direction based on force type
        if force_type in ['N']:  # Axial force - offset in y_local direction
            offset_dir = y_local
        elif force_type in ['Vy']:  # Shear in y - offset in y_local direction
            offset_dir = y_local
        elif force_type in ['Vz']:  # Shear in z - offset in z_local direction
            offset_dir = z_local
        elif force_type in ['T']:  # Torsion - offset in radial direction
            offset_dir = y_local  # Simplified, ideally would be tangential
        elif force_type in ['My']:  # Moment about y - offset in z_local
            offset_dir = z_local
        elif force_type in ['Mz']:  # Moment about z - offset in y_local
            offset_dir = y_local
        
        # Plot the beam element itself
        ax.plot([coords[0][0], coords[1][0]], 
               [coords[0][1], coords[1][1]], 
               [coords[0][2], coords[1][2]], 'k-', linewidth=2, alpha=0.6)
        
        # Sample points along the element
        points = []
        norm_forces = [f * sfac / max_force for f in forces]
        
        for i, loc in enumerate(np.linspace(0, 1, num_points)):
            # Calculate point position along the beam
            point = start_point + loc * element_vector
            points.append(point)
        
        # Create offset points for the force diagram
        diagram_points = []
        for i, point in enumerate(points):
            # Create offset point based on force magnitude and direction
            offset_point = point + offset_dir * norm_forces[i]
            diagram_points.append(offset_point)
        
        # Plot the force diagram line with color gradient based on force intensity
        for i in range(len(diagram_points) - 1):
            # Normalize force for color mapping (0 to 1)
            norm_force = abs(forces[i]) / max_force
            color = cmap(norm_force)
            
            # Plot segment of the force diagram
            ax.plot([diagram_points[i][0], diagram_points[i+1][0]],
                   [diagram_points[i][1], diagram_points[i+1][1]],
                   [diagram_points[i][2], diagram_points[i+1][2]], 
                   c=color, linewidth=2)
            
            # Plot lines connecting the beam to the force diagram
            ax.plot([points[i][0], diagram_points[i][0]],
                   [points[i][1], diagram_points[i][1]],
                   [points[i][2], diagram_points[i][2]], 
                   'k--', alpha=0.3)
        
        # Connect the last point
        i = len(diagram_points) - 1
        ax.plot([points[i][0], diagram_points[i][0]],
               [points[i][1], diagram_points[i][1]],
               [points[i][2], diagram_points[i][2]], 
               'k--', alpha=0.3)
        
        # Add value labels at a few key points (start, middle, end)
        label_indices = [0, num_points//2, num_points-1]
        for idx in label_indices:
            ax.text(diagram_points[idx][0], diagram_points[idx][1], diagram_points[idx][2],
                   f'{forces[idx]:.1f}', color='red', fontsize=8)
    
    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, max_force))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.1)
    cbar.set_label(f'{force_title} Magnitude')
    
    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    # Move title to the top of the figure rather than the top of the axes
    plt.suptitle(f'{force_title} Diagram (Scale Factor: {sfac}) - {load_combination}', 
                fontsize=14, fontweight='bold', y=0.98)
    
    # Get actual model dimensions instead of forcing equal scaling
    all_nodes = ops.getNodeTags()
    x_coords = []
    y_coords = []
    z_coords = []
    
    for node in all_nodes:
        coords = ops.nodeCoord(node)
        x_coords.append(coords[0])
        y_coords.append(coords[1])
        z_coords.append(coords[2])
    
    # Calculate model dimensions
    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)
    z_min, z_max = min(z_coords), max(z_coords)
    
    # Add some padding (10%)
    x_pad = 0.1 * (x_max - x_min)
    y_pad = 0.1 * (y_max - y_min)
    z_pad = 0.1 * (z_max - z_min)
    
    # Only add force diagram offset to limits if needed
    force_pad = max_force * sfac * 1.2  # Extra 20% for labels
    
    ax.set_xlim(x_min - x_pad - force_pad, x_max + x_pad + force_pad)
    ax.set_ylim(y_min - y_pad - force_pad, y_max + y_pad + force_pad)
    ax.set_zlim(z_min - z_pad - force_pad, z_max + z_pad + force_pad)
    
    # Set aspect ratio based on model dimensions to avoid distortion
    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min
    
    # Avoid zero division
    x_range = max(x_range, 0.1)
    y_range = max(y_range, 0.1)
    z_range = max(z_range, 0.1)
    
    ax.set_box_aspect([x_range/min(x_range, y_range, z_range), 
                      y_range/min(x_range, y_range, z_range), 
                      z_range/min(x_range, y_range, z_range)])
    
    plt.tight_layout()
    
    # Save the figure
    # Ensure the directory exists
    os.makedirs(os.path.join(IMAGE_FOLDER, "member"), exist_ok=True)
    # Construct the path and save with load combination
    image_filename = os.path.join(IMAGE_FOLDER, f"member/beam_forces_{force_type}_{load_combination}.png")
    plt.savefig(image_filename, dpi=300, bbox_inches='tight')
    plt.close()

# Extract Node Results
def get_node_results(load_combination="combo"):
    node_tags = ops.getNodeTags()
    node_data = []
    
    for node_tag in node_tags:
        # Get node coordinates
        coords = ops.nodeCoord(node_tag)
        
        # Get displacements
        disp = ops.nodeDisp(node_tag)
        
        # Get reactions (if available)
        reaction = ops.nodeReaction(node_tag)
        
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
    JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
    IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
    output_file = os.path.join(JSON_FOLDER, f"nodes_results_{load_combination}.json")
    os.makedirs(JSON_FOLDER, exist_ok=True)
    os.makedirs(IMAGE_FOLDER, exist_ok=True)
    # Save the data to JSON file
    with open(output_file, 'w') as f:
        json.dump(node_data, f, indent=4)
    
    print(f"Node results saved to: {output_file}")
    return node_data

# Extract Element Results
def structural_model_plot(ele_shapes, OUTPUT_FOLDER="postprocessing_folder", load_combination="combo"):
    # Define the main output directory and subdirectories
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    
    # Create subdirectories
    JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
    IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
    STRUCTURAL_MODEL_FOLDER = os.path.join(IMAGE_FOLDER, "structural_model")
    
    os.makedirs(JSON_FOLDER, exist_ok=True)
    os.makedirs(IMAGE_FOLDER, exist_ok=True)
    os.makedirs(STRUCTURAL_MODEL_FOLDER, exist_ok=True)
    
    # Model plot
    plt.figure(figsize=(10, 8))
    opsv.plot_model()
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
    plt.figure(figsize=(10, 8))
    opsv.plot_extruded_shapes_3d(ele_shapes)
    plt.title(f"Extruded Shapes - {load_combination}")
    filepath = os.path.join(STRUCTURAL_MODEL_FOLDER, f"extruded_shapes_{load_combination}.png")
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    # print(f"Saved: {filepath}")
    plt.close()
    
    print(f"All plots saved to {STRUCTURAL_MODEL_FOLDER}")


# Modified stress plotting function with save capability
def plot_shell_stress(jstr, element_list):

    os.makedirs(r"postprocessing_folder\images\shell_element", exist_ok=True)
    os.makedirs(r"postprocessing_folder\json_files\shell_element", exist_ok=True)
    plt.figure()
    stresses = []
    coords = []
    element_data = []
    
    # Collect element data
    for ele_tag in element_list:
        try:
            # Get element response (forces in N/m)
            response = ops.eleResponse(ele_tag, 'forces')
            membrane_stress = response[:3]  # Nxx, Nyy, Nxy
            
            # Calculate centroid coordinates
            nodes = ops.eleNodes(ele_tag)
            node_coords = np.array([ops.nodeCoord(n)[:2] for n in nodes])
            centroid = np.mean(node_coords, axis=0)
            
            stresses.append(membrane_stress)
            coords.append(centroid)
            
            # Store element data for JSON
            element_data.append({
                "element_id": ele_tag,
                "nodes": [int(n) for n in nodes],
                "centroid": centroid.tolist(),
                "stresses": {
                    "sxx": float(membrane_stress[0]),
                    "syy": float(membrane_stress[1]),
                    "sxy": float(membrane_stress[2])
                }
            })
        except Exception as e:
            print(f"Warning: Could not process element {ele_tag}. Error: {str(e)}")
            continue
    
    # Convert to numpy arrays
    coords = np.array(coords)
    stresses = np.array(stresses)
    
    # Calculate stress component values
    if jstr == 'sxx':
        values = stresses[:, 0]
    elif jstr == 'syy':
        values = stresses[:, 1]
    elif jstr == 'sxy':
        values = stresses[:, 2]
    elif jstr == 'vmis':
        sxx = stresses[:, 0]
        syy = stresses[:, 1]
        sxy = stresses[:, 2]
        values = np.sqrt(sxx**2 + syy**2 - sxx*syy + 3*sxy**2)
    else:
        raise ValueError(f"Unknown stress component: {jstr}")
    
    # Create triangulation
    triangulation = mtri.Triangulation(coords[:, 0], coords[:, 1])
    
    # Create contour plot
    plt.tricontourf(triangulation, values, levels=20, cmap='jet')
    plt.colorbar()
    plt.title(f'{jstr} Stress Distribution')
    plt.xlabel('X [m]')
    plt.ylabel('Y [m]')
    plt.axis('equal')
    
    # Save image
    image_path = fr"postprocessing_folder\images\shell_element\{jstr}_stress.png"
    plt.savefig(image_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save JSON data
    json_path = fr"postprocessing_folder\json_files\shell_element\{jstr}_stress.json"
    with open(json_path, 'w') as f:
        json.dump({
            "stress_type": jstr,
            "elements": element_data,
            "coordinates": coords.tolist(),
            "stress_values": values.tolist()
        }, f, indent=2)

# Modified moment plotting function with save capability
def plot_shell_moment(jstr, element_list):
    os.makedirs(r"postprocessing_folder\images\shell_element", exist_ok=True)
    os.makedirs(r"postprocessing_folder\json_files\shell_element", exist_ok=True)
    plt.figure()
    moments = []
    coords = []
    element_data = []
    
    # Collect element data
    for ele_tag in element_list:
        try:
            # Get element forces response
            response = ops.eleResponse(ele_tag, 'forces')
            element_moments = response[3:6]  # Mxx, Myy, Mxy
            
            # Calculate element centroid
            nodes = ops.eleNodes(ele_tag)
            node_coords = np.array([ops.nodeCoord(n)[:2] for n in nodes])
            centroid = np.mean(node_coords, axis=0)
            
            moments.append(element_moments)
            coords.append(centroid)
            
            # Store element data for JSON
            element_data.append({
                "element_id": ele_tag,
                "nodes": [int(n) for n in nodes],
                "centroid": centroid.tolist(),
                "moments": {
                    "Mxx": float(element_moments[0]),
                    "Myy": float(element_moments[1]),
                    "Mxy": float(element_moments[2])
                }
            })
        except Exception as e:
            print(f"Warning: Could not process element {ele_tag}. Error: {str(e)}")
            continue
    
    # Convert to numpy arrays
    coords = np.array(coords)
    moments = np.array(moments)
    
    # Select moment component
    if jstr == 'Mxx':
        values = moments[:, 0]
    elif jstr == 'Myy':
        values = moments[:, 1]
    elif jstr == 'Mxy':
        values = moments[:, 2]
    else:
        raise ValueError(f"Invalid moment component: {jstr}")
    
    # Create triangulation
    triangulation = mtri.Triangulation(coords[:, 0], coords[:, 1])
    
    # Create contour plot
    plt.tricontourf(triangulation, values, levels=20, cmap='jet')
    plt.colorbar()
    plt.title(f'{jstr} Moment Distribution [N·m/m]')
    plt.xlabel('X [m]')
    plt.ylabel('Y [m]')
    plt.axis('equal')
    
    # Save image
    image_path = fr"postprocessing_folder\images\shell_element\{jstr}_moment.png"
    plt.savefig(image_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save JSON data
    json_path = fr"postprocessing_folder\json_files\shell_element\{jstr}_moment.json"
    with open(json_path, 'w') as f:
        json.dump({
            "moment_type": jstr,
            "elements": element_data,
            "coordinates": coords.tolist(),
            "moment_values": values.tolist()
        }, f, indent=2)

# Main execution
def shell_analysis(JSON_FOLDER):
    # Load structure data
    # Call the first function to load structure data
    nodes_dict, members_dict, structure = load_structure(JSON_FOLDER)

    # Call the second function to create combined structure JSON
    combined_data = create_combined_structure_json(JSON_FOLDER)
    # 2. Load the combined structure data to get all node coordinates
    structure_data = merge_structures(structure, combined_data)
    
    # Get all shell elements
    all_elements = [shell['id'] for shell in structure_data['shell_elements']]
    
    # Plot and save stress components
    stress_components = ['sxx', 'syy', 'sxy', 'vmis']
    for comp in stress_components:
        plot_shell_stress(comp, all_elements)
    
    # Plot and save moment components
    moment_components = ['Mxx', 'Myy', 'Mxy']
    for comp in moment_components:
        plot_shell_moment(comp, all_elements)
    
    print("Postprocessing complete. Results saved in postprocessing_folder")

def calculate_beam_stress(load_combination="combo", num_points=10, view_angle=None, figsize=(12, 10)):
    """
    Calculate normal stresses in beam elements considering axial and bending effects.
    
    Parameters:
    -----------
    load_combination : str
        Name of the load combination being analyzed (used for file naming)
    num_points : int
        Number of points to sample along each element (minimum 10)
    view_angle : tuple
        (elevation, azimuth) for 3D plot view angle. If None, default view is used.
    figsize : tuple
        Figure size as (width, height) in inches
    """
    
    OUTPUT_FOLDER = "postprocessing_folder"
    JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
    IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
    os.makedirs(JSON_FOLDER, exist_ok=True)
    os.makedirs(IMAGE_FOLDER, exist_ok=True)
    
    # Load element and section data
    element_data_path = os.path.join("output_folder", "json_files", "element_data.json")
    try:
        with open(element_data_path, 'r') as f:
            model_data = json.load(f)
    except Exception as e:
        print(f"Error loading element data: {e}")
        return
    
    elements = model_data["elements"]
    sections = model_data["sections"]
    
    # Prepare stress data storage
    all_stress_data = {}
    max_stress_val = 0
    min_stress_val = 0
    
    # Process each element
    for element in elements:
        ele_tag = element["eleTag"]
        section_name = element.get("section_name", "")
        
        # Find section properties
        section = None
        for category in sections.values():
            if section_name in category:
                section = category[section_name]
                break
        
        if not section:
            print(f"Section {section_name} not found for element {ele_tag}")
            continue
            
        # Get geometric properties
        A = element["area"]
        Iy = element["Iy"]
        Iz = element["Iz"]
        section_type = section["type"]
        
        # Get section dimensions
        if section_type == "rectangular":
            B = section["B"]
            H = section["H"]
            y_ext = B/2  # Extreme y distance
            z_ext = H/2  # Extreme z distance
        elif section_type == "circular":
            D = section["D_Sec"]
            radius = D/2
        else:
            print(f"Unsupported section type {section_type} for element {ele_tag}")
            continue
        
        # Initialize element data
        element_stress = {
            "element_tag": ele_tag,
            "section_type": section_type,
            "stress_values": [],
            "max_stress": 0,
            "min_stress": 0
        }
        
        # Sample points along element
        for loc in np.linspace(0, 1, num_points):
            # Get section forces
            forces = ops.eleResponse(ele_tag, 'section', [loc*num_points, 'force'])
            if not forces:
                end_forces = ops.eleForce(ele_tag)
                N = end_forces[0]
                My = end_forces[4]
                Mz = end_forces[5]
            else:
                N = forces[0]
                My = forces[4]
                Mz = forces[5]
            
            # Calculate axial stress component
            sigma_axial = N / A if A != 0 else 0
            
            # Calculate bending stress components
            if section_type == "rectangular":
                sigma_my = (My * z_ext) / Iy if Iy != 0 else 0
                sigma_mz = (Mz * y_ext) / Iz if Iz != 0 else 0
                
                # Calculate corner stresses
                stresses = [
                    sigma_axial + sigma_my + sigma_mz,  # Top-right
                    sigma_axial + sigma_my - sigma_mz,  # Top-left
                    sigma_axial - sigma_my + sigma_mz,  # Bottom-right
                    sigma_axial - sigma_my - sigma_mz   # Bottom-left
                ]
                
            elif section_type == "circular":
                # Resultant moment and stress
                Mr = np.sqrt(My**2 + Mz**2)
                sigma_bending = (Mr * radius) / Iy if Iy != 0 else 0
                stresses = [sigma_axial + sigma_bending, sigma_axial - sigma_bending]
            
            current_max = max(stresses)
            current_min = min(stresses)
            
            # Update global max/min
            max_stress_val = max(max_stress_val, current_max)
            min_stress_val = min(min_stress_val, current_min)
            
            # Store stress data
            element_stress["stress_values"].append({
                "location": float(loc),
                "max_stress": float(current_max),
                "min_stress": float(current_min),
                "N": float(N),
                "My": float(My),
                "Mz": float(Mz)
            })
        
        # Add element data to main dict
        all_stress_data[f"element_{ele_tag}"] = element_stress
    
    OUTPUT_FOLDER = "postprocessing_folder"
    JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
    IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")

    # Create the main directories
    os.makedirs(JSON_FOLDER, exist_ok=True)
    os.makedirs(IMAGE_FOLDER, exist_ok=True)

    # Create the member subdirectory and store its path
    member_folder = os.path.join(JSON_FOLDER, "member")
    os.makedirs(member_folder, exist_ok=True)

    # Save stress data to JSON
    stress_json_path = os.path.join(member_folder, f"beam_stress_{load_combination}.json")

    with open(stress_json_path, 'w') as f:
        json.dump({
            "load_combination": load_combination,
            "max_stress": float(max_stress_val),
            "min_stress": float(min_stress_val),
            "elements": all_stress_data
        }, f, indent=4)
    
    # Visualization (similar to plot_beam_forces but for stress)
    # ... (Add visualization code here using max_stress_val/min_stress_val for color scaling)
    
    print(f"Stress calculation complete for {load_combination}")




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


def plot_beam_deflections(deflection_data, load_combination, image_folder, view_angle=None, figsize=(12, 10)):
    """
    Create visualization plots for beam deflections.
    """
    fig = plt.figure(figsize=figsize)
    
    # Create subplots for different deflection components
    ax1 = plt.subplot(2, 2, 1)
    ax2 = plt.subplot(2, 2, 2)
    ax3 = plt.subplot(2, 2, 3)
    ax4 = plt.subplot(2, 2, 4, projection='3d')
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(deflection_data)))
    
    for i, (element_name, element_data) in enumerate(deflection_data.items()):
        positions = [d['position'] for d in element_data['deflections']]
        dx_vals = [d['dx'] for d in element_data['deflections']]
        dy_vals = [d['dy'] for d in element_data['deflections']]
        dz_vals = [d['dz'] for d in element_data['deflections']]
        
        color = colors[i]
        label = f"Element {element_data['element_tag']}"
        
        # Plot deflections in each direction
        ax1.plot(positions, dx_vals, color=color, label=label, linewidth=2)
        ax2.plot(positions, dy_vals, color=color, label=label, linewidth=2)
        ax3.plot(positions, dz_vals, color=color, label=label, linewidth=2)
        
        # 3D plot
        ax4.plot(positions, dx_vals, dy_vals, color=color, label=label, linewidth=2)
    
    # Formatting
    ax1.set_xlabel('Position along beam')
    ax1.set_ylabel('X Deflection')
    ax1.set_title('X-Direction Deflections')
    ax1.grid(True)
    ax1.legend()
    
    ax2.set_xlabel('Position along beam')
    ax2.set_ylabel('Y Deflection')
    ax2.set_title('Y-Direction Deflections')
    ax2.grid(True)
    ax2.legend()
    
    ax3.set_xlabel('Position along beam')
    ax3.set_ylabel('Z Deflection')
    ax3.set_title('Z-Direction Deflections')
    ax3.grid(True)
    ax3.legend()
    
    ax4.set_xlabel('Position')
    ax4.set_ylabel('X Deflection')
    ax4.set_zlabel('Y Deflection')
    ax4.set_title('3D Deflection Plot')
    if view_angle:
        ax4.view_init(elev=view_angle[0], azim=view_angle[1])
    
    plt.tight_layout()
    
    # Save plot
    plot_path = os.path.join(image_folder, f"beam_deflections_{load_combination}.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Deflection plot saved to: {plot_path}")


def calculate_beam_deflection(load_combination="combo2", num_points=20, view_angle=None, figsize=(12, 10)):
    """
    Calculate deflections and relative deflections for beam elements.
    
    Parameters:
    -----------
    load_combination : str
        Name of the load combination being analyzed (used for file naming)
    num_points : int
        Number of points to sample along each element (minimum 10)
    view_angle : tuple
        (elevation, azimuth) for 3D plot view angle. If None, default view is used.
    figsize : tuple
        Figure size as (width, height) in inches
    """
    
    OUTPUT_FOLDER = "postprocessing_folder"
    JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
    IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
    os.makedirs(JSON_FOLDER, exist_ok=True)
    os.makedirs(IMAGE_FOLDER, exist_ok=True)
    
    # Load element and node data
    element_data_path = os.path.join("output_folder", "json_files", "element_data.json")
    node_data_path = os.path.join("output_folder", "json_files", "combined_structure_data.json")
    INPUT_FOLDER = os.path.join("output_folder", "json_files")
    try:
        with open(element_data_path, 'r') as f:
            element_data = json.load(f)
        # with open(node_data_path, 'r') as f:
        #     node_data = json.load(f)

        # Call the first function to load structure data
        nodes_dict, members_dict, structure = load_structure(INPUT_FOLDER)

        # Call the second function to create combined structure JSON
        combined_data = create_combined_structure_json(INPUT_FOLDER)
        # 2. Load the combined structure data to get all node coordinates
        node_data = merge_structures(structure, combined_data)
        # print(f'data={node_data}')

    except Exception as e:
        print(f"Error loading data: {e}")
        return None
    
    elements = element_data["elements"]
    nodes = node_data["nodes"]
    
    # Prepare deflection data storage
    all_deflection_data = {}
    max_deflection = {'dx': 0, 'dy': 0, 'dz': 0}
    min_deflection = {'dx': 0, 'dy': 0, 'dz': 0}
    max_rel_deflection = {'dy': 0, 'dz': 0}
    min_rel_deflection = {'dy': 0, 'dz': 0}
    
    # Process each element
    for element in elements:
        ele_tag = element["eleTag"]
        i_node = element["node_i_id"]
        j_node = element["node_j_id"]
        L = element["length"]
        
        # Get node coordinates - Fixed the key reference
        node_i = next((n for n in nodes if n["id"] == i_node), None)
        node_j = next((n for n in nodes if n["id"] == j_node), None)
        
        if not node_i or not node_j:
            print(f"Nodes not found for element {ele_tag}: looking for nodes {i_node}, {j_node}")
            continue
        
        try:
            # Get node displacements
            disp_i = ops.nodeDisp(i_node)
            disp_j = ops.nodeDisp(j_node)
            
            # Ensure we have 6 DOF (pad with zeros if necessary)
            if len(disp_i) < 6:
                disp_i = list(disp_i) + [0.0] * (6 - len(disp_i))
            if len(disp_j) < 6:
                disp_j = list(disp_j) + [0.0] * (6 - len(disp_j))
                
        except Exception as e:
            print(f"Error getting displacements for element {ele_tag}: {e}")
            continue
            
        # Initialize element data
        element_deflection = {
            "element_tag": ele_tag,
            "node_i": i_node,
            "node_j": j_node,
            "length": L,
            "deflections": [],
            "max_deflection": {'dx': 0, 'dy': 0, 'dz': 0},
            "min_deflection": {'dx': 0, 'dy': 0, 'dz': 0},
            "max_rel_deflection": {'dy': 0, 'dz': 0},
            "min_rel_deflection": {'dy': 0, 'dz': 0}
        }
        
        # Sample points along element
        for x in np.linspace(0, L, num_points):
            # Calculate deflections
            dx = deflection(disp_i, disp_j, 'dx', x, L)
            dy = deflection(disp_i, disp_j, 'dy', x, L)
            dz = deflection(disp_i, disp_j, 'dz', x, L)
            
            # Calculate relative deflections (only for bending directions)
            rel_dy = rel_deflection(disp_i, disp_j, 'dy', x, L)
            rel_dz = rel_deflection(disp_i, disp_j, 'dz', x, L)
            
            # Update element max/min deflections
            element_deflection["max_deflection"]['dx'] = max(element_deflection["max_deflection"]['dx'], abs(dx))
            element_deflection["max_deflection"]['dy'] = max(element_deflection["max_deflection"]['dy'], abs(dy))
            element_deflection["max_deflection"]['dz'] = max(element_deflection["max_deflection"]['dz'], abs(dz))
            
            element_deflection["min_deflection"]['dx'] = min(element_deflection["min_deflection"]['dx'], dx)
            element_deflection["min_deflection"]['dy'] = min(element_deflection["min_deflection"]['dy'], dy)
            element_deflection["min_deflection"]['dz'] = min(element_deflection["min_deflection"]['dz'], dz)
            
            element_deflection["max_rel_deflection"]['dy'] = max(element_deflection["max_rel_deflection"]['dy'], abs(rel_dy))
            element_deflection["max_rel_deflection"]['dz'] = max(element_deflection["max_rel_deflection"]['dz'], abs(rel_dz))
            
            element_deflection["min_rel_deflection"]['dy'] = min(element_deflection["min_rel_deflection"]['dy'], rel_dy)
            element_deflection["min_rel_deflection"]['dz'] = min(element_deflection["min_rel_deflection"]['dz'], rel_dz)
            
            # Update global max/min deflections
            max_deflection['dx'] = max(max_deflection['dx'], abs(dx))
            max_deflection['dy'] = max(max_deflection['dy'], abs(dy))
            max_deflection['dz'] = max(max_deflection['dz'], abs(dz))
            
            min_deflection['dx'] = min(min_deflection['dx'], dx)
            min_deflection['dy'] = min(min_deflection['dy'], dy)
            min_deflection['dz'] = min(min_deflection['dz'], dz)
            
            max_rel_deflection['dy'] = max(max_rel_deflection['dy'], abs(rel_dy))
            max_rel_deflection['dz'] = max(max_rel_deflection['dz'], abs(rel_dz))
            
            min_rel_deflection['dy'] = min(min_rel_deflection['dy'], rel_dy)
            min_rel_deflection['dz'] = min(min_rel_deflection['dz'], rel_dz)
            
            # Store deflection data
            element_deflection["deflections"].append({
                "position": float(x),
                "dx": float(dx),
                "dy": float(dy),
                "dz": float(dz),
                "rel_dy": float(rel_dy),
                "rel_dz": float(rel_dz)
            })
        
        # Add element data to main dict
        all_deflection_data[f"element_{ele_tag}"] = element_deflection
    
    # Create member subdirectory
    member_folder = os.path.join(JSON_FOLDER, "member")
    os.makedirs(member_folder, exist_ok=True)
    
    # Save deflection data to JSON
    deflection_json_path = os.path.join(member_folder, f"beam_deflection_{load_combination}.json")
    
    with open(deflection_json_path, 'w') as f:
        json.dump({
            "load_combination": load_combination,
            "max_deflection": {k: float(v) for k, v in max_deflection.items()},
            "min_deflection": {k: float(v) for k, v in min_deflection.items()},
            "max_rel_deflection": {k: float(v) for k, v in max_rel_deflection.items()},
            "min_rel_deflection": {k: float(v) for k, v in min_rel_deflection.items()},
            "elements": all_deflection_data
        }, f, indent=4)
    
    # Visualization
    plot_beam_deflections(all_deflection_data, load_combination, IMAGE_FOLDER, view_angle, figsize)
    
    print(f"Deflection calculation complete for {load_combination}")
    print(f"Results saved to: {deflection_json_path}")
    
    return all_deflection_data








# --------------------------
# Geometry Calculations
# --------------------------
def calculate_shell_area(coords):
    """Calculate area for both triangular and quadrilateral shell elements"""
    if len(coords) == 3:  # Triangle
        v1, v2, v3 = [np.array(c) for c in coords]
        return 0.5 * np.linalg.norm(np.cross(v2-v1, v3-v1))
    elif len(coords) == 4:  # Quad
        p1, p2, p3, p4 = [np.array(c) for c in coords]
        return 0.5 * np.linalg.norm(np.cross(p3-p1, p4-p2))
    return 0

# --------------------------
# Interpolation Functions
# --------------------------
def triangular_interpolation(corner_values, xi, eta):
    """Shape functions for triangular elements"""
    return (1 - xi - eta)*corner_values[0] + xi*corner_values[1] + eta*corner_values[2]

def bilinear_interpolation(corner_values, xi, eta):
    """Shape functions for quadrilateral elements"""
    return 0.25*(
        (1-xi)*(1-eta)*corner_values[0] + 
        (1+xi)*(1-eta)*corner_values[1] + 
        (1+xi)*(1+eta)*corner_values[2] + 
        (1-xi)*(1+eta)*corner_values[3]
    )



def calculate_shell_deflections(shell_element, nodes, num_points=5):
    """Calculate deflections across shell element surface"""
    # print(f"DEBUG: Starting calculation for shell element {shell_element.get('id')}")
    element_id = shell_element["id"]
    node_names = shell_element["node_names"]
    node_data = []
    
    # print(f"DEBUG: Processing {len(node_names)} nodes for element {element_id}")
    
    # Get node data
    for node_name in node_names:
        # print(f"DEBUG: Looking for node {node_name}")
        node = next((n for n in nodes if n["name"] == node_name), None)
        if not node:
            print(f"ERROR: Node {node_name} not found for element {element_id}")
            print(f"DEBUG: Available nodes: {[n['name'] for n in nodes]}")
            return None
        
        try:
            # print(f"DEBUG: Getting displacements for node {node['id']}")
            disp = ops.nodeDisp(node["id"])
            # print(f"DEBUG: Node {node['id']} displacements: {disp}")
            padded_disp = list(disp) + [0]*(6-len(disp))  # Pad to 6 DOF
            
            # Create coordinates array from x, y, z keys
            coords = [node["x"], node["y"], node["z"]]
            # print(f"DEBUG: Node coordinates: {coords}")
            
            node_data.append({
                "id": node["id"],
                "coords": coords,
                "disp": padded_disp
            })
            
        except KeyError as ke:
            print(f"ERROR: Missing key in node data: {str(ke)}")
            print(f"DEBUG: Node keys: {node.keys()}")
            return None
        except Exception as e:
            print(f"ERROR: Could not process node {node['id']}")
            print(f"DEBUG: Exception: {str(e)}")
            return None
    
    n_nodes = len(node_data)
    # print(f"DEBUG: Found {n_nodes} valid nodes")
    if n_nodes not in [3, 4]:
        print(f"ERROR: Element {element_id} has {n_nodes} nodes (only 3 or 4 supported)")
        return None
    
    # Generate sampling points
    # print(f"DEBUG: Generating sampling points for {'triangle' if n_nodes == 3 else 'quad'}")
    try:
        if n_nodes == 3:  # Triangular element
            xi, eta = np.meshgrid(np.linspace(0, 1, num_points), np.linspace(0, 1, num_points))
            mask = xi + eta <= 1  # Only points inside triangle
            xi, eta = xi[mask], eta[mask]
        else:  # Quadrilateral element
            xi, eta = np.meshgrid(np.linspace(-1, 1, num_points), np.linspace(-1, 1, num_points))
            xi, eta = xi.flatten(), eta.flatten()
        # print(f"DEBUG: Generated {len(xi)} sampling points")
    except Exception as e:
        print(f"ERROR: Failed to generate sampling points")
        print(f"DEBUG: Exception: {str(e)}")
        return None
    
    # Interpolate at each point
    deflection_points = []
    # print(f"DEBUG: Starting interpolation at sampling points")
    for i, (x, e) in enumerate(zip(xi, eta)):
        try:
            # print(f"DEBUG: Point {i+1}/{len(xi)} - natural coords: ({x:.2f}, {e:.2f})")
            interp_fn = triangular_interpolation if n_nodes == 3 else bilinear_interpolation
            
            # Interpolate coordinates and displacements
            x_coords = interp_fn([n["coords"][0] for n in node_data], x, e)
            y_coords = interp_fn([n["coords"][1] for n in node_data], x, e)
            z_coords = interp_fn([n["coords"][2] for n in node_data], x, e)
            dx = interp_fn([n["disp"][0] for n in node_data], x, e)
            dy = interp_fn([n["disp"][1] for n in node_data], x, e)
            dz = interp_fn([n["disp"][2] for n in node_data], x, e)
            
            deflection_points.append({
                "natural_coords": [float(x), float(e)],
                "physical_coords": [x_coords, y_coords, z_coords],
                "deflections": {"dx": dx, "dy": dy, "dz": dz}
            })
        except Exception as ex:
            print(f"ERROR: Failed interpolation at point {i} (natural coords: {x}, {e})")
            print(f"DEBUG: Exception: {str(ex)}")
            return None
    
    # print(f"DEBUG: Successfully calculated deflections for element {element_id}")
    return {
        "element_id": element_id,
        "element_name": shell_element.get("name", f"Shell_{element_id}"),
        "node_data": node_data,
        "area": calculate_shell_area([n["coords"] for n in node_data]),
        "deflection_points": deflection_points,
        "element_type": "triangular" if n_nodes == 3 else "quadrilateral"
    }


# --------------------------
# Plotting Functions (Combined Plots for All Elements)
# --------------------------
def plot_combined_deformed_shape(all_element_data, output_dir, figsize=(12, 10)):
    """Create and save combined deformed shape plot for all elements"""
    if not all_element_data:
        return
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')
    
    all_mag = []
    for element_data in all_element_data.values():
        points = element_data["deflection_points"]
        x = [p["physical_coords"][0] for p in points]
        y = [p["physical_coords"][1] for p in points]
        z = [p["physical_coords"][2] for p in points]
        dx = [p["deflections"]["dx"] for p in points]
        dy = [p["deflections"]["dy"] for p in points]
        dz = [p["deflections"]["dz"] for p in points]
        mag = np.sqrt(np.array(dx)**2 + np.array(dy)**2 + np.array(dz)**2)
        all_mag.extend(mag)
        
        # Plot deformed shape for this element
        ax.scatter(np.array(x)+np.array(dx), 
                  np.array(y)+np.array(dy), 
                  np.array(z)+np.array(dz), 
                  c=mag, cmap='viridis', s=30, alpha=0.8)
    
    # Create colorbar with combined range
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=min(all_mag), vmax=max(all_mag)))
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label='Deflection Magnitude')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title("Combined Deformed Shape - All Shell Elements")
    
    # Save figure
    plot_file = os.path.join(output_dir, "combined_deformed_shape.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    return plot_file

def plot_combined_deflection_component(all_element_data, component, output_dir, figsize=(12, 10)):
    """Create and save combined deflection component plot for all elements"""
    if not all_element_data:
        return
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    
    all_deflections = []
    all_x = []
    all_y = []
    
    for element_data in all_element_data.values():
        points = element_data["deflection_points"]
        x = [p["physical_coords"][0] for p in points]
        y = [p["physical_coords"][1] for p in points]
        deflections = [p["deflections"][component] for p in points]
        
        all_x.extend(x)
        all_y.extend(y)
        all_deflections.extend(deflections)
    
    # Create triangulation for all points
    tri = Triangulation(all_x, all_y)
    cs = ax.tricontourf(tri, all_deflections, levels=20, cmap='RdBu_r')
    
    plt.colorbar(cs, label=f'{component.upper()} Deflection')
    ax.set_title(f"Combined {component.upper()} Deflections - All Shell Elements")
    ax.set_aspect('equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    
    # Save figure
    plot_file = os.path.join(output_dir, f"combined_{component}_deflection.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    return plot_file

# --------------------------
# Main Analysis Function
# --------------------------
def analyze_shell_deflections(load_case="combo1", samples_per_side=5):
    """Main analysis workflow"""
    # Setup directories
    output_root = "postprocessing_folder"
    os.makedirs(output_root, exist_ok=True)
    
    json_dir = os.path.join(output_root, "json_files", "shell_elements")
    img_dir = os.path.join(output_root, "images", "shell_deflections")
    os.makedirs(json_dir, exist_ok=True)
    os.makedirs(img_dir, exist_ok=True)
    INPUT_FOLDER = os.path.join("OUTPUT_FOLDER", "json_files")
    # Load structure data
    # Call the first function to load structure data
    nodes_dict, members_dict, structure = load_structure(INPUT_FOLDER)

    # Call the second function to create combined structure JSON
    combined_data = create_combined_structure_json(INPUT_FOLDER)
    # 2. Load the combined structure data to get all node coordinates
    model = merge_structures(structure, combined_data)
    # print(f'model["nodes"]={model["nodes"]}')
    
    # Process each shell element
    results = {}
    for element in model.get("shell_elements", []):
        element_id = element["id"]
        # print(f"Processing element {element_id}...")
        
        # Calculate deflections
        deflection_data = calculate_shell_deflections(element, model["nodes"], samples_per_side)
        # print(f'deflection_data={deflection_data}')
        if not deflection_data:
            continue
        
        results[f"shell_{element_id}"] = deflection_data
    
    # Generate combined plots for all elements
    plot_files = {
        "deformed": plot_combined_deformed_shape(results, img_dir),
        "dx": plot_combined_deflection_component(results, "dx", img_dir),
        "dy": plot_combined_deflection_component(results, "dy", img_dir),
        "dz": plot_combined_deflection_component(results, "dz", img_dir)
    }
    print(f'plot_files={plot_files}')
    
    # Save combined results
    output_file = os.path.join(json_dir, f"shell_deflections_{load_case}.json")
    with open(output_file, 'w') as f:
        json.dump({
            "load_case": load_case,
            "elements": results,
            "plot_files": {k: os.path.relpath(v, output_root) for k,v in plot_files.items()},
            "metadata": {
                "sampling_points": samples_per_side,
                "timestamp": str(np.datetime64('now'))
            }
        }, f, indent=2)
    
    print(f"\nAnalysis complete. Results saved to:")
    print(f"- JSON data: {output_file}")
    print(f"- Combined plots directory: {img_dir}")
    print(f"- Generated plots:")
    for plot_type, plot_file in plot_files.items():
        print(f"  - {plot_type}: {os.path.basename(plot_file)}")


def calculate_slab_reinforcement_from_json():
    # Create output directories
    os.makedirs(r"postprocessing_folder/images/slab_reinforcement", exist_ok=True)
    os.makedirs(r"postprocessing_folder/json_files/slab_reinforcement", exist_ok=True)
    
    # Material properties (SI units)
    phi = 0.9      # Strength reduction factor
    fy = 420e6      # Yield strength of steel (Pa)
    d = 0.2         # Effective depth (m)
    fc = 28e6       # Concrete compressive strength (Pa)
    b = 1.0         # Unit width (m)
    
    # Load all moment components
    moment_types = ['Mxx', 'Myy', 'Mxy']
    moment_data = {}
    
    for mt in moment_types:
        json_path = f"postprocessing_folder/json_files/shell_element/{mt}_moment.json"
        try:
            with open(json_path, 'r') as f:
                moment_data[mt] = json.load(f)
        except FileNotFoundError:
            print(f"Error: Moment file not found at {json_path}")
            return
    
    # Initialize results storage
    reinforcement_results = []
    coords = []
    As_x_bottom = []
    As_y_bottom = []
    As_x_top = []
    As_y_top = []
    
    # Process each element
    for i, element in enumerate(moment_data['Mxx']['elements']):
        ele_tag = element['element_id']
        centroid = element['centroid']
        
        # Get moments for this element from all components
        Mx = element['moments']['Mxx']
        My = moment_data['Myy']['elements'][i]['moments']['Myy']
        Mxy = moment_data['Mxy']['elements'][i]['moments']['Mxy']
        
        # Calculate reinforcement
        As_x_b, As_y_b, As_x_t, As_y_t = calculate_slab_moment(
            Mx, My, Mxy, phi, fy, d, fc, b)
        
        # Store results
        reinforcement_results.append({
            "element_id": ele_tag,
            "centroid": centroid,
            "moments": {"Mxx": Mx, "Myy": My, "Mxy": Mxy},
            "reinforcement": {
                "As_x_bottom": As_x_b,
                "As_y_bottom": As_y_b,
                "As_x_top": As_x_t,
                "As_y_top": As_y_t
            }
        })
        
        coords.append(centroid)
        As_x_bottom.append(As_x_b)
        As_y_bottom.append(As_y_b)
        As_x_top.append(As_x_t)
        As_y_top.append(As_y_t)
    
    # Convert to numpy arrays
    coords = np.array(coords)
    As_x_bottom = np.array(As_x_bottom)
    As_y_bottom = np.array(As_y_bottom)
    As_x_top = np.array(As_x_top)
    As_y_top = np.array(As_y_top)
    
    # Plot results
    plot_reinforcement_contours(coords, As_x_bottom, As_y_bottom, As_x_top, As_y_top)
    
    # Save results to JSON
    output_path = "postprocessing_folder/json_files/slab_reinforcement/results.json"
    with open(output_path, 'w') as f:
        json.dump(reinforcement_results, f, indent=2)
    
    print(f"Reinforcement calculation complete. Results saved to {output_path}")

def calculate_slab_moment(Mx, My, Mxy, phi, fy, d, fc, b):
    """Calculate required reinforcement for a slab element."""
    # Bottom reinforcement
    Mx_bottom = Mx + abs(Mxy)
    My_bottom = My + abs(Mxy)
    
    # Handle negative moments
    if My_bottom < 0:
        My_bottom = 0
        Mx_bottom = Mx + abs(Mxy**2 / My) if My != 0 else Mx
    if Mx_bottom < 0:
        Mx_bottom = 0
        My_bottom = My + abs(Mxy**2 / Mx) if Mx != 0 else My
    
    # Top reinforcement
    Mx_top = Mx - abs(Mxy)
    My_top = My - abs(Mxy)
    
    if My_top < 0:
        My_top = 0
        Mx_top = Mx - abs(Mxy**2 / My) if My != 0 else Mx
    if Mx_top < 0:
        Mx_top = 0
        My_top = My - abs(Mxy**2 / Mx) if Mx != 0 else My
    
    # Calculate required steel area
    As_x_bottom = solve_As(Mx_bottom, phi, fy, d, fc, b)
    As_y_bottom = solve_As(My_bottom, phi, fy, d, fc, b)
    As_x_top = solve_As(Mx_top, phi, fy, d, fc, b)
    As_y_top = solve_As(My_top, phi, fy, d, fc, b)
    
    return As_x_bottom, As_y_bottom, As_x_top, As_y_top

def solve_As(Mu, phi, fy, d, fc, b, tolerance=1e-6, max_iterations=100):
    """Iteratively solve for required steel area."""
    if Mu <= 0:
        return 0.0
    
    As = Mu / (phi * fy * d)  # Initial estimate
    for _ in range(max_iterations):
        a = As * fy / (0.85 * fc * b)
        As_new = Mu / (phi * fy * (d - a/2))
        
        if abs(As_new - As) < tolerance:
            return max(As_new, 0)  # Ensure non-negative
        As = As_new
    
    print(f"Warning: Max iterations reached for Mu={Mu}")
    return max(As, 0)

def plot_reinforcement_contours(coords, As_x_bottom, As_y_bottom, As_x_top, As_y_top):
    """Create contour plots of reinforcement."""
    # Ensure all values are non-negative and replace zeros with small positive numbers
    def clean_values(values):
        values = np.array(values)
        values[values < 0] = 0  # Set negative values to zero
        values[values == 0] = 1e-10  # Replace zeros with small positive value for log scale
        return values
    
    As_x_bottom = clean_values(As_x_bottom)
    As_y_bottom = clean_values(As_y_bottom)
    As_x_top = clean_values(As_x_top)
    As_y_top = clean_values(As_y_top)
    
    # Create triangulation
    triangulation = mtri.Triangulation(coords[:, 0], coords[:, 1])
    
    # Common plot settings
    plot_settings = {
        'levels': 20,
        'cmap': 'jet',
        'norm': colors.LogNorm(vmin=1e-6, vmax=1e-3),  # Adjusted for typical reinforcement
        'extend': 'both'
    }
    
    # Plot each component
    components = {
        'As_x_bottom': As_x_bottom,
        'As_y_bottom': As_y_bottom,
        'As_x_top': As_x_top,
        'As_y_top': As_y_top
    }
    
    for name, values in components.items():
        plt.figure(figsize=(10, 8))
        contour = plt.tricontourf(triangulation, values, **plot_settings)
        
        # Add contour lines
        plt.tricontour(triangulation, values, levels=plot_settings['levels'],
                      colors='k', linewidths=0.5, alpha=0.5)
        
        plt.colorbar(contour, label='Reinforcement Area (m²/m)')
        plt.title(f'{name.replace("_", " ").title()}')
        plt.xlabel('X coordinate (m)')
        plt.ylabel('Y coordinate (m)')
        plt.axis('equal')
        
        # Save figure
        output_path = f"postprocessing_folder/images/slab_reinforcement/{name}.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
   
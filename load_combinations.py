# from build_model import build_model
from import_ import *
from Grid_and_structure_creation import load_structure
# from test import create_combined_structure_json
from wall_meshing import create_combined_structure_json
from units import *


def merge_structures(structure, combined_data):
    """
    Merges the structure data (nodes and members) with combined data (nodes, members, and shell elements).
    
    Args:
        structure: Dictionary from load_structure() with keys "nodes" and "members"
        combined_data: Dictionary from create_combined_structure_json() with keys "nodes", "members", "shell_elements"
        
    Returns:
        A single merged dictionary with all nodes, members, and shell elements
    """
    merged = {
        "nodes": [],
        "members": [],
        "shell_elements": combined_data.get("shell_elements", [])
    }
    
    # Create sets to track existing IDs and avoid duplicates
    existing_node_ids = set()
    existing_member_ids = set()
    
    # First add nodes from structure (original nodes)
    for node in structure.get("nodes", []):
        merged["nodes"].append(node)
        existing_node_ids.add(node["id"])
    
    # Then add nodes from combined_data (may include shell nodes)
    for node in combined_data.get("nodes", []):
        # Convert combined_data format to match structure format if needed
        if node["id"] not in existing_node_ids:
            merged["nodes"].append({
                "id": node["id"],
                "name": node["name"],
                "x": node.get("x", node.get("coordinates", [0,0,0])[0]),
                "y": node.get("y", node.get("coordinates", [0,0,0])[1]),
                "z": node.get("z", node.get("coordinates", [0,0,0])[2])
            })
            existing_node_ids.add(node["id"])
    
    # Add members from structure (original members)
    for member in structure.get("members", []):
        merged["members"].append(member)
        existing_member_ids.add(member["id"])
    
    # Add members from combined_data (if any)
    for member in combined_data.get("members", []):
        if member["id"] not in existing_member_ids:
            merged["members"].append(member)
    
    return merged

    
def surface_element_calculator(node_coords, pressure=None, rho_h=0.0):
    """
    Unified function to calculate surface area, surface loads, and mass matrix 
    for quadrilateral or triangular elements
    
    Parameters:
    -----------
    node_coords : array-like
        List of coordinates for the nodes (4 for quad, 3 for triangle)
    pressure : float, optional
        Pressure value to apply to the surface (pressure per unit area)
    rho_h : float, optional
        Density per unit area for mass calculation, default is 0.0
        
    Returns:
    --------
    results : dict
        Dictionary containing:
        - 'area': surface area
        - 'forces': array of forces (if pressure is provided)
        - 'mass_matrix': mass matrix (if rho_h > 0)
    """
    # Convert to numpy array if not already
    node_coords = np.array(node_coords)
    n_nodes = len(node_coords)
    results = {'area': 0.0}
    
    # Determine element type based on number of nodes
    if n_nodes == 4:
        element_type = 'quad'
        n_dofs = 12  # 4 nodes × 3 DOFs
    elif n_nodes == 3:
        element_type = 'tri'
        n_dofs = 9   # 3 nodes × 3 DOFs
    else:
        raise ValueError("Element must have 3 (triangle) or 4 (quadrilateral) nodes")
    
    # Initialize mass matrix if needed
    if rho_h > 0:
        mass_matrix = np.zeros((n_dofs, n_dofs))
    
    # Quadrilateral element calculations
    if element_type == 'quad':
        # Gauss points for 2×2 integration
        one_over_root3 = 1.0 / np.sqrt(3.0)
        gauss_points = [
            [-one_over_root3, -one_over_root3],
            [one_over_root3, -one_over_root3],
            [one_over_root3, one_over_root3],
            [-one_over_root3, one_over_root3]
        ]
        gauss_weights = [1.0, 1.0, 1.0, 1.0]  # Equal weights for 2×2 integration
        
        # Initialize variables
        area_contributions = []
        normals = []
        shape_functions = []
        
        # First pass: calculate all geometry information
        for gp_idx, (xi, eta) in enumerate(gauss_points):
            # Calculate tangent vectors
            g1 = 0.25 * ((1 - eta) * (node_coords[1] - node_coords[0]) + 
                         (1 + eta) * (node_coords[2] - node_coords[3]))
            g2 = 0.25 * ((1 + xi) * (node_coords[2] - node_coords[1]) + 
                         (1 - xi) * (node_coords[3] - node_coords[0]))
            
            # Calculate normal vector as cross product of g1 and g2
            normal = np.cross(g1, g2)
            normals.append(normal)
            
            # Calculate area contribution
            area_contribution = np.linalg.norm(normal) * gauss_weights[gp_idx]
            area_contributions.append(area_contribution)
            results['area'] += area_contribution
            
            # Store shape functions if pressure or mass is needed
            if pressure is not None or rho_h > 0:
                N = np.zeros(4)
                N[0] = 0.25 * (1 - xi) * (1 - eta)
                N[1] = 0.25 * (1 + xi) * (1 - eta)
                N[2] = 0.25 * (1 + xi) * (1 + eta)
                N[3] = 0.25 * (1 - xi) * (1 + eta)
                shape_functions.append((N, gauss_weights[gp_idx]))
        
        # Calculate forces if pressure is provided
        if pressure is not None:
            forces = np.zeros(n_dofs)
            
            # Calculate total pressure (pressure × area)
            total_pressure = pressure
            # Loop through integration points
            for gp_idx in range(4):
                # Get stored values for this integration point
                normal = normals[gp_idx]
                N, weight = shape_functions[gp_idx]
                
                # Apply weighted pressure to each node
                for i in range(4):
                    for j in range(3):
                        forces[i*3 + j] -= total_pressure * area_contributions[gp_idx] * normal[j]/np.linalg.norm(normal) * N[i]
            
            results['forces'] = forces
    return results


def get_node_name_to_id_mapping(JSON_FOLDER):
        """Helper function to get node name to ID mapping"""
        # Call the first function to load structure data
        nodes_dict, members_dict, structure = load_structure(JSON_FOLDER)

        # Call the second function to create combined structure JSON
        combined_data = create_combined_structure_json(JSON_FOLDER)
        # 2. Load the combined structure data to get all node coordinates
        structure_data = merge_structures(structure, combined_data)
        
        node_name_to_id = {}
        for node in structure_data.get("nodes", []):
            node_name_to_id[node["name"]] = node["id"]
        
        print(f"Loaded {len(node_name_to_id)} nodes from combined structure file")
        return node_name_to_id
    
def get_element_name_to_id_mapping(JSON_FOLDER):
    """Helper function to get element name to ID mapping"""
    # Call the first function to load structure data
    nodes_dict, members_dict, structure = load_structure(JSON_FOLDER)

    # Call the second function to create combined structure JSON
    combined_data = create_combined_structure_json(JSON_FOLDER)
    # 2. Load the combined structure data to get all node coordinates
    structure_data = merge_structures(structure, combined_data)
    
    element_name_to_id = {}
    for member in structure_data.get("members", []):
        element_name_to_id[member["name"]] = member["id"]
    
    print(f"Loaded {len(element_name_to_id)} elements from combined structure file")
    return element_name_to_id


def apply_shell_load(JSON_FOLDER, load_case_names, pressures, numbering=1):
    """
    Apply shell load on each element of mesh_data_with_predefined_points{numbering}.json for multiple load cases,
    convert the surface_load to nodal load, and save them in separate files for each load case.
    Also saves element areas with corresponding pressures.
    
    Args:
        JSON_FOLDER (str): Path to folder containing the mesh data
        load_case_names (list): List of load case names (default is ['DL', 'LL', 'WX'])
        pressures (list): List of pressure values to apply (positive is outward normal)
        numbering (int): Number identifier for the mesh file
    """
    if len(load_case_names) != len(pressures):
        raise ValueError("load_case_names and pressures must have the same length")

    # 1. Load the mesh data file
    mesh_file = os.path.join(JSON_FOLDER, 'wall', f'mesh_data_with_predefined_points{numbering}.json')
    if not os.path.exists(mesh_file):
        raise FileNotFoundError(f"Mesh file not found at: {mesh_file}")
    
    with open(mesh_file, 'r') as f:
        mesh_data = json.load(f)
    
    # Call the first function to load structure data
    nodes_dict, members_dict, structure = load_structure(JSON_FOLDER)

    # Call the second function to create combined structure JSON
    combined_data = create_combined_structure_json(JSON_FOLDER)
    # 2. Load the combined structure data to get all node coordinates
    combined_data = merge_structures(structure, combined_data)
    
    # Create a dictionary of all nodes with their coordinates
    all_nodes = {}
    for node in combined_data["nodes"]:
        all_nodes[node["name"]] = {
            "x": node["x"],
            "y": node["y"],
            "z": node["z"]
        }
    
    # Initialize output data structure that will contain all load cases
    output_data_all = {}
    element_areas = {}  # Store area and pressure info for each element
    
    # First pass to calculate element areas (same for all load cases)
    for elem_name, elem_data in mesh_data["elements"].items():
        # Get the node names for this element
        node_names = elem_data["nodes"]
        
        # Skip if not quad or triangle
        if len(node_names) not in [3, 4]:
            continue
        
        # Get the node coordinates in proper order
        node_coords = [
                np.array([all_nodes[name]['x'], all_nodes[name]['y'], all_nodes[name]['z']])
                for name in node_names  # Maintain the element's node ordering
            ]

        # Calculate element area
        result = surface_element_calculator(node_coords, pressure=1.0)  # Use 1.0 to get area
            
        element_areas[elem_name] = {
            "area": result['area'],
            "load_cases": {}
        }
    
    # Process each load case
    for load_case_name, pressure in zip(load_case_names, pressures):
        # Track loads for this specific load case
        nodal_loads = {}
        total_area = 0.0
        total_force = [0.0, 0.0, 0.0]
        
        for elem_name, elem_data in mesh_data["elements"].items():
            # Get the node names for this element
            node_names = elem_data["nodes"]
            
            # Skip if not quad or triangle
            if len(node_names) not in [3, 4]:
                continue
            
            # Get the node coordinates in proper order
            node_coords = [
                np.array([all_nodes[name]['x'], all_nodes[name]['y'], all_nodes[name]['z']])
                for name in node_names  # Maintain the element's node ordering
            ]

            # Calculate element forces with actual pressure
            result = surface_element_calculator(node_coords, pressure=pressure)
            
            # Store pressure information for this load case
            element_areas[elem_name]["load_cases"][load_case_name] = {
                "pressure": pressure,
                "total_force": -pressure * element_areas[elem_name]["area"]  # Negative for outward normal
            }
            
            # Distribute loads to nodes
            for i, node_name in enumerate(node_names):
                # Get the force components for this node
                fx = result['forces'][i*3]
                fy = result['forces'][i*3 + 1]
                fz = result['forces'][i*3 + 2]
                
                # Track the load in our nodal_loads dictionary
                if node_name not in nodal_loads:
                    nodal_loads[node_name] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                
                nodal_loads[node_name][0] += fx
                nodal_loads[node_name][1] += fy
                nodal_loads[node_name][2] += fz
                total_force[0] += fx
                total_force[1] += fy
                total_force[2] += fz
        
        # Calculate theoretical total force (pressure * total area)
        total_area = sum(elem["area"] for elem in element_areas.values())
        theoretical_force = -pressure * total_area  # Negative because pressure is positive outward
        
        # Print verification information
        print("\n=== Force Balance Verification ===")
        print(f"Total Surface Area: {total_area:.6f}")
        print(f"Theoretical Total Force (pressure × area): {theoretical_force:.6f} (should be in normal direction)")
        print(f"\nLoad Case: {load_case_name}")
        print(f"Sum of Nodal Forces: {total_force}")
        print(f"Difference: {total_force[2] - theoretical_force:.6f} (z-component)")
        
        # Create output data for this specific load case
        output_data = {}
        for node_name, forces in nodal_loads.items():
            output_key = f"{node_name}"  # Unique key for node + load case
            output_data[output_key] = [
                node_name,      # node_id
                load_case_name, # load_case_name
                forces[0],      # Fx
                forces[1],      # Fy
                forces[2],      # Fz
                forces[3],      # Mx
                forces[4],      # My
                forces[5]       # Mz
            ]
        
        # Save this load case to a separate JSON file
        output_file = f"nodal_loads_{load_case_name}.{numbering}.json"
        output_path = os.path.join(JSON_FOLDER, "load_data", output_file)
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(output_data, f, indent=4)
        
        print(f"Processed and saved load case: {load_case_name} to {output_path}")
        
        # Also store in combined output (optional)
        output_data_all.update(output_data)
    

    # Save element areas with pressure information
    areas_file = f"element_areas_{numbering}.json"
    areas_path = os.path.join(JSON_FOLDER, "load_data", areas_file)
    with open(areas_path, 'w') as f:
        json.dump(element_areas, f, indent=4)
    
    print(f"\nElement areas with pressure information saved to: {areas_path}")
    print(f"Shell loads have been applied to all load cases")
    return True


def apply_self_weight(JSON_FOLDER, load_case_names="self_weight"):
    """
    Calculate the self-weight of beam/column members and convert them into nodal loads.
    Save the results in a new JSON file in the load_data folder.
    
    Args:
        JSON_FOLDER (str): Path to folder containing the structural data
        load_case_names (str or list): Load case name(s) for self-weight (default is "self_weight")
    """
    # print("11 is working")
    # Convert single load case name to list for consistent processing
    if isinstance(load_case_names, str):
        load_case_names = [load_case_names]
    
    # 1. Load the element data file
    element_data_path = os.path.join(JSON_FOLDER, "element_data.json")
    if not os.path.exists(element_data_path):
        raise FileNotFoundError(f"Element data file not found at: {element_data_path}")
    
    with open(element_data_path, 'r') as f:
        element_data = json.load(f)
    
    # 2. Load the combined structure data to get node information
    # Call the first function to load structure data
    nodes_dict, members_dict, structure = load_structure(JSON_FOLDER)

    # Call the second function to create combined structure JSON
    combined_data = create_combined_structure_json(JSON_FOLDER)
    # 2. Load the combined structure data to get all node coordinates
    combined_data = merge_structures(structure, combined_data)
    
    # Create a mapping from node ID to node name (if needed)
    node_id_to_name = {node["id"]: node["name"] for node in combined_data["nodes"]}
    
    # 3. Initialize data structures for nodal loads
    nodal_loads = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # Fx, Fy, Fz, Mx, My, Mz
    
    # 4. Process each element to calculate self-weight and distribute to nodes
    total_weight = 0.0
    for element in element_data["elements"]:
        # Skip if element doesn't have unit_weight (like shell elements)
        if "unit_weight" not in element or element["unit_weight"] is None:
            continue
        
        # Calculate total weight of the element
        weight = element["unit_weight"] * element["area"] * element["length"]
        total_weight += weight
        
        # Distribute half the weight to each node (uniform distribution)
        half_weight = weight / 2
        
        # Add to start node (negative z-direction for gravity)
        node_i_id = element["node_i_id"]
        nodal_loads[node_i_id][2] -= half_weight  # Fz component
        
        # Add to end node (negative z-direction for gravity)
        node_j_id = element["node_j_id"]
        nodal_loads[node_j_id][2] -= half_weight  # Fz component
    
    print(f"Total structural self-weight calculated: {total_weight:.2f} N")
    
    # 5. Format the output data for all load cases
    output_data = {}
    for load_case in load_case_names:
        for node_id, forces in nodal_loads.items():
            # Get node name from ID if needed
            node_name = node_id_to_name.get(node_id, str(node_id))
            
            output_key = f"{node_name}"
            output_data[output_key] = [
                node_name,      # node_name (or node_id if name not available)
                load_case,      # load_case_name
                forces[0],      # Fx
                forces[1],      # Fy
                forces[2],      # Fz
                forces[3],     # Mx
                forces[4],     # My
                forces[5]      # Mz
            ]
    
    # 6. Save the output to a JSON file with load case names in filename
    load_case_str = "_".join(load_case_names)
    output_file = f"nodal_loads_{load_case_str}.json"
    output_path = os.path.join(JSON_FOLDER, "load_data", output_file)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    # print(f"output_path={output_path}")
    # print(f"output_data={output_data}")
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=4)
    
    print(f"Successfully saved self-weight nodal loads to: {output_path}")
    return True

def save_nodal_load_cases(JSON_FOLDER, load_entries):
    """
    Save load cases to separate files in the load_data subfolder.
    
    Args:
        JSON_FOLDER (str): Path to the main JSON folder
        load_entries (list): List of all load entries with format:
                            [node_id, load_case, Fx, Fy, Fz, Mx, My, Mz]
    """
    # Create the load_data directory path
    load_data_dir = os.path.join(JSON_FOLDER, "load_data")
    os.makedirs(load_data_dir, exist_ok=True)
    
    # Reorganize the data by load case
    load_cases = defaultdict(list)
    
    for entry in load_entries:
        load_case = entry[1]  # The load case is at index 1
        load_cases[load_case].append(entry)
    
    # Save each load case to separate files
    for load_case, entries in load_cases.items():
        filename = f"nodal_loads_{load_case}.json"
        filepath = os.path.join(load_data_dir, filename)
        
        # Convert to dictionary format with node_id as keys
        output_data = {entry[0]: entry for entry in entries}
        
        with open(filepath, 'w') as f:
            json.dump(output_data, f, indent=4)
        
        print(f"Saved {len(entries)} node loads to {filepath}")
    
    print("\nAll load cases processed successfully!")
    return True


def create_load_combinations(JSON_FOLDER, load_combinations):
    """
    Create load combinations by multiplying nodal loads with appropriate load factors
    and summing them up for each combination. Only handles these filename patterns:
    1. nodal_loads_loadcasename.number.json (e.g., nodal_loads_DL.1.json)
    2. nodal_loads_loadcasename.json (e.g., nodal_loads_LL.json)

    Args:
        JSON_FOLDER (str): Path to folder containing the load data
        load_combinations (dict): Dictionary of load combinations with their factors
                                 Example: {"Comb1": [("DL", 1.4)], 
                                          "Comb2": [("DL", 1.2), ("LL", 1.6)],
                                          "Mass": [("DL", 1.0)]}
    """
    # Constants
    g = 32.17  # ft/s² (gravitational acceleration)
    
    # Get node name to ID mapping
    node_name_to_id = get_node_name_to_id_mapping(JSON_FOLDER)
    
    # 1. Collect all unique load cases from all combinations
    all_load_cases = set()
    for combo_load_cases in load_combinations.values():
        for load_case, _ in combo_load_cases:
            all_load_cases.add(load_case)
    
    # 2. Load all the required load case files (only matching the two patterns)
    load_case_data = {}
    load_data_dir = os.path.join(JSON_FOLDER, "load_data")
    
    # Debug print: show the load_data_dir path we're trying to access
    print(f"Attempting to access load data directory: {load_data_dir}")
    print(f"Directory exists: {os.path.exists(load_data_dir)}")
    if os.path.exists(load_data_dir):
        print(f"Directory contents: {os.listdir(load_data_dir)}")
    
    # Find all nodal load files in the directory with the two allowed patterns
    nodal_load_files = []
    for f in os.listdir(load_data_dir):
        if f.startswith("nodal_loads_") and f.endswith(".json"):
            # Check for pattern 1: nodal_loads_loadcasename.number.json
            if re.match(r"nodal_loads_(.+)\.\d+\.json", f):
                nodal_load_files.append(f)
            # Check for pattern 2: nodal_loads_loadcasename.json (without number)
            elif re.match(r"nodal_loads_(.+)\.json", f) and not '.' in f[len("nodal_loads_"):-len(".json")]:
                nodal_load_files.append(f)
    
    print(f'Filtered nodal_load_files found (only matching allowed patterns): {nodal_load_files}')
    
    # Map each load case to its file(s)
    load_case_to_files = {}
    for filename in nodal_load_files:
        # Check for pattern 1: nodal_loads_loadcasename.number.json
        pattern1_match = re.match(r"nodal_loads_(.+)\.(\d+)\.json", filename)
        if pattern1_match:
            load_case = pattern1_match.group(1)
            if load_case in all_load_cases:
                if load_case not in load_case_to_files:
                    load_case_to_files[load_case] = []
                load_case_to_files[load_case].append(filename)
                print(f"Mapped load case '{load_case}' to file '{filename}' (pattern 1)")
        
        # Check for pattern 2: nodal_loads_loadcasename.json (without number)
        else:
            load_case = filename[len("nodal_loads_"):-len(".json")]
            if load_case in all_load_cases:
                if load_case not in load_case_to_files:
                    load_case_to_files[load_case] = []
                load_case_to_files[load_case].append(filename)
                print(f"Mapped load case '{load_case}' to file '{filename}' (pattern 2)")
    
    # Verify we have files for all required load cases
    missing_load_cases = all_load_cases - set(load_case_to_files.keys())
    print(f"Required load cases: {all_load_cases}")
    print(f"Found load cases: {set(load_case_to_files.keys())}")
    print(f"Missing load cases: {missing_load_cases}")
    
    if missing_load_cases:
        raise FileNotFoundError(f"Missing load case files for: {', '.join(missing_load_cases)}")
    
    # Load the data from each file
    for load_case, filenames in load_case_to_files.items():
        for filename in filenames:
            filepath = os.path.join(load_data_dir, filename)
            print(f"Loading data for load case '{load_case}' from file: {filepath}")
            
            if not os.path.exists(filepath):
                print(f"ERROR: File not found at path: {filepath}")
                continue
                
            with open(filepath, 'r') as f:
                data = json.load(f)
                
                # For files with pattern 1: nodal_loads_loadcasename.number.json
                pattern1_match = re.match(r"nodal_loads_(.+)\.(\d+)\.json", filename)
                if pattern1_match:
                    print(f"File {filename} contains single load case data (pattern 1)")
                    if load_case not in load_case_data:
                        load_case_data[load_case] = {}
                    load_case_data[load_case].update(data)
                
                # For files with pattern 2: nodal_loads_loadcasename.json (without number)
                else:
                    print(f"File {filename} contains single load case data (pattern 2)")
                    if load_case not in load_case_data:
                        load_case_data[load_case] = {}
                    load_case_data[load_case].update(data)
    
    # 3. Process each load combination (rest remains the same)
    for combo_name, combo_load_cases in load_combinations.items():
        combined_loads = {}
        
        # Initialize all nodes with zero loads
        all_nodes_in_combo = set()
        for load_case, _ in combo_load_cases:
            if load_case not in load_case_data:
                raise ValueError(f"Load case '{load_case}' not found in loaded data")
            for node_key in load_case_data[load_case].keys():
                node_name = load_case_data[load_case][node_key][0]  # First element is node name
                all_nodes_in_combo.add(node_name)
        
        # Initialize with zeros
        for node_name in all_nodes_in_combo:
            combined_loads[node_name] = {
                "Fx": 0.0,
                "Fy": 0.0,
                "Fz": 0.0,
                "Mx": 0.0,
                "My": 0.0,
                "Mz": 0.0
            }
        
        # Apply each load case with its factor
        for load_case, factor in combo_load_cases:
            if load_case not in load_case_data:
                raise ValueError(f"Load case '{load_case}' not found in loaded data")
            
            for node_key, load_values in load_case_data[load_case].items():
                node_name = load_values[0]  # First element is node name
                
                if node_name not in combined_loads:
                    continue
                
                # Fixed indexing: load_values contains [node_name, load_case, Fx, Fy, Fz, Mx, My, Mz]
                combined_loads[node_name]["Fx"] += load_values[2] * factor  # Fx (index 2)
                combined_loads[node_name]["Fy"] += load_values[3] * factor  # Fy (index 3)
                combined_loads[node_name]["Fz"] += load_values[4] * factor  # Fz (index 4)
                combined_loads[node_name]["Mx"] += load_values[5] * factor  # Mx (index 5)
                combined_loads[node_name]["My"] += load_values[6] * factor  # My (index 6)
                combined_loads[node_name]["Mz"] += load_values[7] * factor  # Mz (index 7)
        
        # 4. Format the output data differently for mass combinations
        output_data = {}
        
        if combo_name.lower() == "mass":
            # Convert forces to masses (m = F/g)
            for node_name, loads in combined_loads.items():
                output_data[node_name] = [
                    node_name_to_id.get(node_name),  # node_id_tag
                    combo_name,                      # load_combination_name
                    loads["Fz"]/g,                  # mass in x-direction (kip·s²/ft)
                    loads["Fz"]/g,                  # mass in y-direction
                    0.0,                            # mass in z-direction
                    0.0,                            # rotational mass (typically zero)
                    0.0,
                    0.0
                ]
        else:
            # Regular load combination format
            for node_name, loads in combined_loads.items():
                output_data[node_name] = [
                    node_name_to_id.get(node_name),  # node_id_tag
                    combo_name,         # load_combination_name
                    loads["Fx"],       # Fx
                    loads["Fy"],       # Fy
                    loads["Fz"],       # Fz
                    loads["Mx"],      # Mx
                    loads["My"],      # My
                    loads["Mz"],      # Mz
                ]
        
        # Create the output filename using both patterns
        # Pattern 1: nodal_loads_loadcasename.number.json
        
        # Pattern 2: nodal_loads_loadcasename.json
        output_file = f"nodal_loads_{combo_name}.json"
        output_path = os.path.join(load_data_dir, output_file)
        
        # Ensure directory exists
        os.makedirs(load_data_dir, exist_ok=True)
        
        # Save using both patterns
        for output_path in [output_path]:
            print(f"Saving output to: {output_path}")
            with open(output_path, 'w') as f:
                json.dump(output_data, f, indent=4)
            print(f"Created load combination file: {output_path}")
    
    return True


def apply_nodal_loads(JSON_FOLDER, load_combination="Comb2", load_pattern_tag=1, time_series_tag=None):
    """
    Simplified function to apply nodal loads from pre-combined load files.
    Handles mass loads differently by applying them as nodal masses.
    
    Args:
        JSON_FOLDER (str): Path to folder containing load_data
        load_combination (str): Name of load combination to apply (default: "Comb2")
        load_pattern_tag (int): Tag for the load pattern (default: 1)
        time_series_tag (int): Optional time series tag for dynamic analysis
    """
    # 1. Open the load combination file
    nodal_loads_file = os.path.join(JSON_FOLDER, "load_data", f"nodal_loads_{load_combination}.json")
    
    if not os.path.exists(nodal_loads_file):
        raise FileNotFoundError(f"Load combination file not found: {nodal_loads_file}")
    
    with open(nodal_loads_file, 'r') as f:
        nodal_loads = json.load(f)
    
    # Check if this is a mass load case
    is_mass_load = load_combination.lower() == "mass"
    
    if not is_mass_load:
        # 2. Create load pattern for force loads
        pattern_type = 'Plain'
        if time_series_tag is not None:
            ops.pattern(pattern_type, load_pattern_tag, time_series_tag)
        else:
            ops.pattern(pattern_type, load_pattern_tag, 'Linear')
    
    # 3. Apply all nodal loads/masses
    for load_values in nodal_loads.values():
        node_tag = int(load_values[0])  # First value is node name/ID
        
        if is_mass_load:
            # For mass loads: [node_id, "mass", m_x, m_y, m_z, 0, 0, 0]
            mass_components = load_values[2:5]  # Only take m_x, m_y, m_z
            try:
                ops.mass(node_tag, *mass_components)
            except Exception as e:
                print(f"Error applying mass to node {node_tag}: {str(e)}")
        else:
            # For regular loads: [node_id, combo_name, Fx, Fy, Fz, Mx, My, Mz]
            load_components = load_values[2:]  # Skip node name and load case name
            try:
                ops.load(node_tag, *load_components)
            except Exception as e:
                print(f"Error applying load to node {node_tag}: {str(e)}")
    
    if is_mass_load:
        print(f"Applied {len(nodal_loads)} nodal masses from {load_combination}")
    else:
        print(f"Applied {len(nodal_loads)} nodal loads from {load_combination}")
    
    return True


def create_combined_loads(load_combinations, *element_loads, JSON_FOLDER):
    """Combine loads from multiple element_loads dictionaries"""
    combined_loads = {}
    
    # First combine all element loads into a single structure
    all_loads = {}
    for load_dict in element_loads:
        for elem_ids, load_list in load_dict.items():
            if elem_ids not in all_loads:
                all_loads[elem_ids] = []
            all_loads[elem_ids].extend(load_list)
    
    # Process each load combination
    for comb_name, case_factors in load_combinations.items():
        combined_loads[comb_name] = {}
        
        # Initialize structure for each element group
        for elem_ids in all_loads.keys():
            combined_loads[comb_name][elem_ids] = [{
                "LoadCase": [comb_name],
                "uniform": {"x": 0.0, "y": 0.0, "z": 0.0},
                "point": {"x": 0.0, "y": 0.0, "z": 0.0, "location": 0.0},
                "temperature_points": []
            }]
        
        # Apply load factors - skip if load case doesn't exist
        for case_name, factor in case_factors:
            found_case = False
            # Check if this load case exists in any element loads
            for elem_ids, load_list in all_loads.items():
                for load_config in load_list:
                    if case_name in load_config["LoadCase"]:
                        found_case = True
                        break
                if found_case:
                    break
            
            if not found_case:
                print(f"Warning: Load case '{case_name}' not found - skipping")
                continue
                
            print(f"Processing {case_name} with factor {factor} for combination {comb_name}")
            
            # Find matching load cases in all_loads
            for elem_ids, load_list in all_loads.items():
                for load_config in load_list:
                    if case_name in load_config["LoadCase"]:
                        target_config = combined_loads[comb_name][elem_ids][0]
                        
                        # Add uniform loads
                        if "uniform" in load_config:
                            for dim in ["x", "y", "z"]:
                                target_config["uniform"][dim] += load_config["uniform"].get(dim, 0.0) * factor
                        
                        # Add point loads (take location from first non-zero point load)
                        if "point" in load_config:
                            for dim in ["x", "y", "z"]:
                                target_config["point"][dim] += load_config["point"].get(dim, 0.0) * factor
                            # Set location if not already set and current load has a location
                            if (target_config["point"]["location"] == 0.0 and 
                                load_config["point"].get("location", 0.0) != 0.0):
                                target_config["point"]["location"] = load_config["point"]["location"]
                        
                        # Add temperature loads
                        if "temperature_points" in load_config:
                            if not target_config["temperature_points"]:
                                # Initialize temperature points structure
                                target_config["temperature_points"] = [
                                    {"temp": 0.0, "y": tp["y"]} 
                                    for tp in load_config["temperature_points"]
                                ]
                            
                            # Add temperature values
                            for i, tp in enumerate(load_config["temperature_points"]):
                                if i < len(target_config["temperature_points"]):
                                    target_config["temperature_points"][i]["temp"] += tp["temp"] * factor
    
    # Directory to save load data
    load_data_dir = JSON_FOLDER
    os.makedirs(load_data_dir, exist_ok=True)
    load_data_dir = os.path.join(JSON_FOLDER, "load_data")
    os.makedirs(load_data_dir, exist_ok=True)
    # Save each load combination to a separate JSON file
    for combo_name, load_data in combined_loads.items():
        # Prepare output data
        output_data = {
            "load_case": combo_name,
            "load_data": load_data
        }
        
        # Create the output filename
        output_file = f"member_loads_{combo_name}.json"
        output_path = os.path.join(load_data_dir, output_file)
        
        # Save the file
        with open(output_path, 'w') as f:
            json.dump(output_data, f, indent=4)
        print(f"Saved load combination {combo_name} to {output_path}")
    
    return combined_loads

def convert_member_loads_to_nodal(JSON_FOLDER):
    """
    Converts member loads (uniform and point) into equivalent nodal loads.
    Properly sums contributions from multiple members meeting at common nodes.
    Creates separate output files for each load combination.
    
    Args:
        JSON_FOLDER: Path to the directory containing load data
        
    Returns:
        Dictionary of nodal loads for each load combination
    """

    # Load structure data
    nodes_dict, members_dict, structure = load_structure(JSON_FOLDER)
    combined_data = create_combined_structure_json(JSON_FOLDER)
    combined_data = merge_structures(structure, combined_data)
    # Create node ID to coordinates mapping
    node_coords = {}
    for node in combined_data["nodes"]:
        node_coords[node["id"]] = {
            "x": node.get("x", node.get("coordinates", [0,0,0])[0]),
            "y": node.get("y", node.get("coordinates", [0,0,0])[1]),
            "z": node.get("z", node.get("coordinates", [0,0,0])[2])
        }
    
    # Create member ID to node IDs mapping
    member_nodes = {}
    for member in combined_data["members"]:
        member_nodes[member["id"]] = {
            "start_node": member["start_node_id"],
            "end_node": member["end_node_id"]
        }
    
    # Create node to members mapping to track all members connected to each node
    node_members = defaultdict(list)
    for member_id, nodes in member_nodes.items():
        node_members[nodes["start_node"]].append(member_id)
        node_members[nodes["end_node"]].append(member_id)
    
    # Directory to load member load data
    load_data_dir = os.path.join(JSON_FOLDER, "load_data")
    if not os.path.exists(load_data_dir):
        print(f"No load data directory found at {load_data_dir}")
        return {}
    
    # Create output directory for nodal loads
    nodal_loads_dir = os.path.join(JSON_FOLDER, "load_data")
    os.makedirs(nodal_loads_dir, exist_ok=True)
    
    # Initialize dictionary to store all nodal loads
    all_nodal_loads = {}
    
    # Process each load combination file
    for filename in sorted(os.listdir(load_data_dir)):
        if filename.startswith("member_loads_") and filename.endswith(".json") and not "converted_into_nodal_loads" in filename:
            combo_name = filename[len("member_loads_"):-len(".json")]
            
            # Load the member load data
            with open(os.path.join(load_data_dir, filename), 'r') as f:
                try:
                    load_data = json.load(f)
                except json.JSONDecodeError as e:
                    print(f"Error loading {filename}: {e}")
                    continue
            
            # Initialize nodal loads for this combination using defaultdict
            nodal_loads = defaultdict(lambda: {"x": 0.0, "y": 0.0, "z": 0.0})
            
            # Process each member group in the load data
            for member_ids, load_configs in load_data.get("load_data", {}).items():
                if not load_configs:
                    continue
                    
                # Parse member IDs (handling both single ID and list of IDs)
                try:
                    if isinstance(member_ids, str):
                        if member_ids.startswith('[') and member_ids.endswith(']'):
                            member_id_list = [int(id_str.strip()) for id_str in member_ids[1:-1].split(',')]
                        else:
                            member_id_list = [int(member_ids)]
                    else:
                        member_id_list = member_ids
                except (ValueError, TypeError) as e:
                    print(f"Invalid member ID format in {filename}: {member_ids} - {e}")
                    continue
                
                for member_id in member_id_list:
                    if member_id not in member_nodes:
                        print(f"Warning: Member ID {member_id} not found in structure - skipping")
                        continue
                        
                    start_node = member_nodes[member_id]["start_node"]
                    end_node = member_nodes[member_id]["end_node"]
                    
                    # Get node coordinates
                    start_coords = node_coords.get(start_node, {"x": 0, "y": 0, "z": 0})
                    end_coords = node_coords.get(end_node, {"x": 0, "y": 0, "z": 0})
                    
                    # Calculate member length
                    dx = end_coords["x"] - start_coords["x"]
                    dy = end_coords["y"] - start_coords["y"]
                    dz = end_coords["z"] - start_coords["z"]
                    length = math.sqrt(dx**2 + dy**2 + dz**2)
                    
                    if length == 0:
                        print(f"Warning: Member {member_id} has zero length - skipping")
                        continue
                    
                    # Process each load configuration for this member group
                    for load_config in load_configs:
                        # Process uniform loads
                        if "uniform" in load_config:
                            ux = load_config["uniform"].get("x", 0.0)
                            uy = load_config["uniform"].get("y", 0.0)
                            uz = load_config["uniform"].get("z", 0.0)
                            
                            # For uniform load, each node gets half the total load
                            total_x = ux * length
                            total_y = uy * length
                            total_z = uz * length
                            
                            nodal_loads[start_node]["x"] += total_x / 2
                            nodal_loads[start_node]["y"] += total_y / 2
                            nodal_loads[start_node]["z"] += total_z / 2
                            
                            nodal_loads[end_node]["x"] += total_x / 2
                            nodal_loads[end_node]["y"] += total_y / 2
                            nodal_loads[end_node]["z"] += total_z / 2
                        
                        # Process point loads
                        if "point" in load_config:
                            px = load_config["point"].get("x", 0.0)
                            py = load_config["point"].get("y", 0.0)
                            pz = load_config["point"].get("z", 0.0)
                            location = load_config["point"].get("location", 0.0)
                            
                            if location < 0 or location > length:
                                print(f"Warning: Point load location {location} out of bounds for member {member_id} - skipping")
                                continue
                            
                            # Calculate distribution factors
                            a = location
                            b = length - location
                            
                            # Distribute point load to nodes
                            if length > 0:
                                nodal_loads[start_node]["x"] += px * b / length
                                nodal_loads[end_node]["x"] += px * a / length
                                nodal_loads[start_node]["y"] += py * b / length
                                nodal_loads[end_node]["y"] += py * a / length
                                nodal_loads[start_node]["z"] += pz * b / length
                                nodal_loads[end_node]["z"] += pz * a / length
            
            # Convert defaultdict to regular dict for JSON serialization
            nodal_loads_dict = dict(nodal_loads)
            
            # Store nodal loads for this combination
            all_nodal_loads[combo_name] = nodal_loads_dict
            
            # Save this combination's nodal loads to a separate file
            output_file = f"member_loads_converted_into_nodal_loads_{combo_name}.json"
            output_path = os.path.join(nodal_loads_dir, output_file)
            # Create node ID to name mapping from combined_data
            node_id_to_name = {}
            for node in combined_data["nodes"]:
                node_id_to_name[node["id"]] = node["name"]

            output_data = {}
            for node_id, loads in nodal_loads_dict.items():
                # Get the actual node name from the mapping
                node_name = node_id_to_name.get(node_id, f"N{node_id}")  # fallback to N{id} if name not found
                
                output_data[node_name] = [
                    node_id,
                    combo_name,
                    loads["x"],
                    loads["y"], 
                    loads["z"],
                    0.0,  # Mx (moment about x-axis)
                    0.0,  # My (moment about y-axis)
                    0.0   # Mz (moment about z-axis)
                ]
            
            with open(output_path, 'w') as f:
                json.dump(output_data, f, indent=4)
            print(f"Saved nodal loads for {combo_name} to {output_path}")
    

    # print(f"Saved all nodal loads to {combined_output_path}")
    
    return all_nodal_loads

def nodal_mass_combo(JSON_FOLDER, load_combinations):

    # Directory to save load data
    load_data_dir = JSON_FOLDER
    os.makedirs(load_data_dir, exist_ok=True)
    load_data_dir = os.path.join(JSON_FOLDER, "load_data")
    os.makedirs(load_data_dir, exist_ok=True)

    for combo_name, combo_load_cases in load_combinations.items():
        # Define file paths
        member_loads_file = f"nodal_loads_{combo_name}.json"
        member_loads_path = os.path.join(load_data_dir, member_loads_file)
        
        converted_loads_file = f"member_loads_converted_into_nodal_loads_{combo_name}.json"
        converted_loads_path = os.path.join(load_data_dir, converted_loads_file)
        
        # Initialize combined nodal loads dictionary
        combined_nodal_loads = {}
        
        # Read member loads file (if it exists and contains nodal loads)
        if os.path.exists(member_loads_path):
            try:
                with open(member_loads_path, 'r') as f:
                    member_data = json.load(f)
                    # The data is already in nodal loads format, so just add all entries
                    for node_name, loads in member_data.items():
                        combined_nodal_loads[node_name] = loads.copy()
            except (json.JSONDecodeError, KeyError) as e:
                print(f"Error reading {member_loads_file}: {e}")
        
        # Read converted nodal loads file
        if os.path.exists(converted_loads_path):
            try:
                with open(converted_loads_path, 'r') as f:
                    converted_data = json.load(f)
                    print(f'converted_data12={converted_data}')
                    # Handle different possible data structures
                    nodal_loads_data = converted_data
                    print(f'nodal_loads_data={nodal_loads_data}')
                    
                    # Sum up the loads
                    for node_name, loads in nodal_loads_data.items():
                        if node_name not in combined_nodal_loads:
                            # If it's array format [id, combo, fx, fy, fz, mx, my, mz]
                            if isinstance(loads, list) and len(loads) >= 8:
                                combined_nodal_loads[node_name] = [
                                    loads[0],  # node_id
                                    loads[1],  # combo_name
                                    loads[2],  # fx
                                    loads[3],  # fy
                                    loads[4],  # fz
                                    loads[5],  # mx
                                    loads[6],  # my
                                    loads[7]   # mz
                                ]
                            # If it's dict format {"x": fx, "y": fy, "z": fz}
                            elif isinstance(loads, dict):
                                # Get node_id from node_name (extract number)
                                node_id = int(node_name.replace('N', '').replace('n', ''))
                                combined_nodal_loads[node_name] = [
                                    node_id,
                                    combo_name,
                                    loads.get("x", 0.0),
                                    loads.get("y", 0.0),
                                    loads.get("z", 0.0),
                                    0.0,  # mx
                                    0.0,  # my
                                    0.0   # mz
                                ]
                        else:
                            # Sum existing loads
                            existing = combined_nodal_loads[node_name]
                            if isinstance(loads, list) and len(loads) >= 8:
                                # Sum forces and moments
                                existing[2] += loads[2]  # fx
                                existing[3] += loads[3]  # fy
                                existing[4] += loads[4]  # fz
                                existing[5] += loads[5]  # mx
                                existing[6] += loads[6]  # my
                                existing[7] += loads[7]  # mz
                            elif isinstance(loads, dict):
                                # Sum dict format loads
                                existing[2] += loads.get("x", 0.0)  # fx
                                existing[3] += loads.get("y", 0.0)  # fy
                                existing[4] += loads.get("z", 0.0)  # fz
                            
            except (json.JSONDecodeError, KeyError) as e:
                print(f"Error reading {converted_loads_file}: {e}")
        
        # Save combined nodal loads
        if combined_nodal_loads:
            output_file = f"nodal_mass_loads_{combo_name}.json"
            output_path = os.path.join(load_data_dir, output_file)
            
            with open(output_path, 'w') as f:
                json.dump(combined_nodal_loads, f, indent=4)
            
            print(f"Saved combined nodal loads for {combo_name} to {output_path}")
        else:
            print(f"No nodal loads found for combination {combo_name}")


# def process_loads(JSON_FOLDER, surface_configurations=None, nodal_load_entries=None, 
#                  load_combinations=None, all_element_loads=None, run_parts=None):
#     """
#     Consolidated function to handle all load processing:
#     - Applies shell loads and self-weight
#     - Saves nodal load cases
#     - Creates load combinations
#     - Converts member loads to nodal loads
#     - Handles nodal mass combinations
    
#     Args:
#         JSON_FOLDER: Path to folder containing structural data
#         surface_configurations: Dictionary of surface load configurations
#         nodal_load_entries: List of nodal load entries
#         load_combinations: Dictionary of load combinations
#         all_element_loads: List of element load dictionaries
#         run_parts: List of parts being run (for logging)
#     """
    
#     # 1. Process surface configurations (shell loads)
#     if surface_configurations is not None and surface_configurations:
#         for config_name, config_data in surface_configurations.items():           
#             numbering = list(surface_configurations.keys()).index(config_name) + 1
            
#             # Load mesh data
#             mesh_file = os.path.join(JSON_FOLDER, 'wall', f'mesh_data_with_predefined_points{numbering}.json')
#             if not os.path.exists(mesh_file):
#                 raise FileNotFoundError(f"Mesh file not found at: {mesh_file}")
            
#             with open(mesh_file, 'r') as f:
#                 mesh_data = json.load(f)
            
#             # Load structure data
#             nodes_dict, members_dict, structure = load_structure(JSON_FOLDER)
#             combined_data = create_combined_structure_json(JSON_FOLDER)
#             combined_data = merge_structures(structure, combined_data)
            
#             # Create node coordinate mapping
#             all_nodes = {node["name"]: {"x": node["x"], "y": node["y"], "z": node["z"]} 
#                         for node in combined_data["nodes"]}
            
#             # Initialize output structures
#             output_data_all = {}
#             element_areas = {}
            
#             # First pass to calculate element areas
#             for elem_name, elem_data in mesh_data["elements"].items():
#                 node_names = elem_data["nodes"]
#                 if len(node_names) not in [3, 4]: continue
                
#                 node_coords = [np.array([all_nodes[name]['x'], all_nodes[name]['y'], all_nodes[name]['z']])
#                              for name in node_names]
#                 result = surface_element_calculator(node_coords, pressure=1.0)
#                 element_areas[elem_name] = {"area": result['area'], "load_cases": {}}
            
#             # Process each load case
#             for load_case_name, pressure in zip(config_data["load_case_names"], config_data["pressures"]):
#                 nodal_loads = {}
#                 total_force = [0.0, 0.0, 0.0]
                
#                 for elem_name, elem_data in mesh_data["elements"].items():
#                     node_names = elem_data["nodes"]
#                     if len(node_names) not in [3, 4]: continue
                    
#                     node_coords = [np.array([all_nodes[name]['x'], all_nodes[name]['y'], all_nodes[name]['z']])
#                                  for name in node_names]
#                     result = surface_element_calculator(node_coords, pressure=pressure)
                    
#                     # Store pressure info
#                     element_areas[elem_name]["load_cases"][load_case_name] = {
#                         "pressure": pressure,
#                         "total_force": -pressure * element_areas[elem_name]["area"]
#                     }
                    
#                     # Distribute to nodes
#                     for i, node_name in enumerate(node_names):
#                         if node_name not in nodal_loads:
#                             nodal_loads[node_name] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        
#                         nodal_loads[node_name][0] += result['forces'][i*3]
#                         nodal_loads[node_name][1] += result['forces'][i*3 + 1]
#                         nodal_loads[node_name][2] += result['forces'][i*3 + 2]
#                         total_force[0] += result['forces'][i*3]
#                         total_force[1] += result['forces'][i*3 + 1]
#                         total_force[2] += result['forces'][i*3 + 2]
                
#                 # Save load case
#                 output_data = {}
#                 for node_name, forces in nodal_loads.items():
#                     output_data[node_name] = [
#                         node_name, load_case_name, 
#                         forces[0], forces[1], forces[2], 
#                         forces[3], forces[4], forces[5]
#                     ]
                
#                 output_file = f"nodal_loads_{load_case_name}.{numbering}.json"
#                 output_path = os.path.join(JSON_FOLDER, "load_data", output_file)
#                 os.makedirs(os.path.dirname(output_path), exist_ok=True)
#                 with open(output_path, 'w') as f:
#                     json.dump(output_data, f, indent=4)
                
#                 output_data_all.update(output_data)
            
#             # Save element areas
#             areas_file = f"element_areas_{numbering}.json"
#             areas_path = os.path.join(JSON_FOLDER, "load_data", areas_file)
#             with open(areas_path, 'w') as f:
#                 json.dump(element_areas, f, indent=4)
    
#     # 2. Process self-weight
#     load_case_names = "self_weight"
#     if isinstance(load_case_names, str):
#         load_case_names = [load_case_names]
    
#     # Load element data
#     element_data_path = os.path.join(JSON_FOLDER, "element_data.json")
#     if not os.path.exists(element_data_path):
#         raise FileNotFoundError(f"Element data file not found at: {element_data_path}")
    
#     with open(element_data_path, 'r') as f:
#         element_data = json.load(f)
    
#     # Load structure data
#     nodes_dict, members_dict, structure = load_structure(JSON_FOLDER)
#     combined_data = create_combined_structure_json(JSON_FOLDER)
#     combined_data = merge_structures(structure, combined_data)
    
#     # Create node mapping
#     node_id_to_name = {node["id"]: node["name"] for node in combined_data["nodes"]}
    
#     # Calculate self-weight
#     nodal_loads = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#     total_weight = 0.0
    
#     for element in element_data["elements"]:
#         if "unit_weight" not in element or element["unit_weight"] is None:
#             continue
        
#         weight = element["unit_weight"] * element["area"] * element["length"]
#         total_weight += weight
#         half_weight = weight / 2
        
#         nodal_loads[element["node_i_id"]][2] -= half_weight
#         nodal_loads[element["node_j_id"]][2] -= half_weight
    
#     # Save self-weight loads
#     output_data = {}
#     for load_case in load_case_names:
#         for node_id, forces in nodal_loads.items():
#             node_name = node_id_to_name.get(node_id, str(node_id))
#             output_data[f"{node_name}"] = [
#                 node_name, load_case, 
#                 forces[0], forces[1], forces[2], 
#                 forces[3], forces[4], forces[5]
#             ]
    
#     load_case_str = "_".join(load_case_names)
#     output_file = f"nodal_loads_{load_case_str}.json"
#     output_path = os.path.join(JSON_FOLDER, "load_data", output_file)
#     os.makedirs(os.path.dirname(output_path), exist_ok=True)
#     with open(output_path, 'w') as f:
#         json.dump(output_data, f, indent=4)
    
#     # 3. Save nodal load cases
#     if nodal_load_entries:
#         load_data_dir = os.path.join(JSON_FOLDER, "load_data")
#         os.makedirs(load_data_dir, exist_ok=True)
        
#         load_cases = defaultdict(list)
#         for entry in nodal_load_entries:
#             load_cases[entry[1]].append(entry)
        
#         for load_case, entries in load_cases.items():
#             filename = f"nodal_loads_{load_case}.json"
#             filepath = os.path.join(load_data_dir, filename)
#             output_data = {entry[0]: entry for entry in entries}
#             with open(filepath, 'w') as f:
#                 json.dump(output_data, f, indent=4)
    
#     # 4. Create load combinations
#     if load_combinations:
#         node_name_to_id = get_node_name_to_id_mapping(JSON_FOLDER)
#         load_data_dir = os.path.join(JSON_FOLDER, "load_data")
        
#         # Find all nodal load files
#         nodal_load_files = []
#         for f in os.listdir(load_data_dir):
#             if f.startswith("nodal_loads_") and f.endswith(".json"):
#                 if re.match(r"nodal_loads_(.+)\.\d+\.json", f) or \
#                    (re.match(r"nodal_loads_(.+)\.json", f) and not '.' in f[len("nodal_loads_"):-len(".json")]):
#                     nodal_load_files.append(f)
        
#         # Load all required load cases
#         load_case_data = {}
#         for filename in nodal_load_files:
#             filepath = os.path.join(load_data_dir, filename)
#             with open(filepath, 'r') as f:
#                 data = json.load(f)
                
#                 if re.match(r"nodal_loads_(.+)\.(\d+)\.json", filename):
#                     load_case = re.match(r"nodal_loads_(.+)\.(\d+)\.json", filename).group(1)
#                     if load_case not in load_case_data:
#                         load_case_data[load_case] = {}
#                     load_case_data[load_case].update(data)
#                 else:
#                     load_case = filename[len("nodal_loads_"):-len(".json")]
#                     if load_case not in load_case_data:
#                         load_case_data[load_case] = {}
#                     load_case_data[load_case].update(data)
        
#         # Process each combination
#         for combo_name, combo_load_cases in load_combinations.items():
#             combined_loads = {}
#             all_nodes_in_combo = set()
            
#             # Initialize all nodes with zero loads
#             for load_case, _ in combo_load_cases:
#                 for node_key in load_case_data[load_case].keys():
#                     node_name = load_case_data[load_case][node_key][0]
#                     all_nodes_in_combo.add(node_name)
            
#             for node_name in all_nodes_in_combo:
#                 combined_loads[node_name] = {
#                     "Fx": 0.0, "Fy": 0.0, "Fz": 0.0,
#                     "Mx": 0.0, "My": 0.0, "Mz": 0.0
#                 }
            
#             # Apply load factors
#             for load_case, factor in combo_load_cases:
#                 for node_key, load_values in load_case_data[load_case].items():
#                     node_name = load_values[0]
#                     if node_name not in combined_loads: continue
                    
#                     combined_loads[node_name]["Fx"] += load_values[2] * factor
#                     combined_loads[node_name]["Fy"] += load_values[3] * factor
#                     combined_loads[node_name]["Fz"] += load_values[4] * factor
#                     combined_loads[node_name]["Mx"] += load_values[5] * factor
#                     combined_loads[node_name]["My"] += load_values[6] * factor
#                     combined_loads[node_name]["Mz"] += load_values[7] * factor
            
#             # Format output
#             output_data = {}
#             for node_name, loads in combined_loads.items():
#                 if combo_name.lower() == "mass":
#                     output_data[node_name] = [
#                         node_name_to_id.get(node_name),
#                         combo_name,
#                         loads["Fz"]/32.17,  # Convert to mass
#                         loads["Fz"]/32.17,
#                         0.0, 0.0, 0.0, 0.0
#                     ]
#                 else:
#                     output_data[node_name] = [
#                         node_name_to_id.get(node_name),
#                         combo_name,
#                         loads["Fx"], loads["Fy"], loads["Fz"],
#                         loads["Mx"], loads["My"], loads["Mz"]
#                     ]
            
#             # Save combination
#             output_file = f"nodal_loads_{combo_name}.json"
#             output_path = os.path.join(load_data_dir, output_file)
#             with open(output_path, 'w') as f:
#                 json.dump(output_data, f, indent=4)
    
#     # 5. Create combined loads from element loads
#     if load_combinations and all_element_loads:
#         combined_loads = {}
#         all_loads = {}
        
#         # Combine all element loads
#         for load_dict in all_element_loads:
#             for elem_ids, load_list in load_dict.items():
#                 if elem_ids not in all_loads:
#                     all_loads[elem_ids] = []
#                 all_loads[elem_ids].extend(load_list)
        
#         # Process each combination
#         for comb_name, case_factors in load_combinations.items():
#             combined_loads[comb_name] = {}
            
#             # Initialize structure
#             for elem_ids in all_loads.keys():
#                 combined_loads[comb_name][elem_ids] = [{
#                     "LoadCase": [comb_name],
#                     "uniform": {"x": 0.0, "y": 0.0, "z": 0.0},
#                     "point": {"x": 0.0, "y": 0.0, "z": 0.0, "location": 0.0},
#                     "temperature_points": []
#                 }]
            
#             # Apply factors
#             for case_name, factor in case_factors:
#                 for elem_ids, load_list in all_loads.items():
#                     for load_config in load_list:
#                         if case_name in load_config["LoadCase"]:
#                             target = combined_loads[comb_name][elem_ids][0]
                            
#                             # Uniform loads
#                             if "uniform" in load_config:
#                                 for dim in ["x", "y", "z"]:
#                                     target["uniform"][dim] += load_config["uniform"].get(dim, 0.0) * factor
                            
#                             # Point loads
#                             if "point" in load_config:
#                                 for dim in ["x", "y", "z"]:
#                                     target["point"][dim] += load_config["point"].get(dim, 0.0) * factor
#                                 if target["point"]["location"] == 0.0 and load_config["point"].get("location", 0.0) != 0.0:
#                                     target["point"]["location"] = load_config["point"]["location"]
                            
#                             # Temperature loads
#                             if "temperature_points" in load_config:
#                                 if not target["temperature_points"]:
#                                     target["temperature_points"] = [
#                                         {"temp": 0.0, "y": tp["y"]} 
#                                         for tp in load_config["temperature_points"]
#                                     ]
#                                 for i, tp in enumerate(load_config["temperature_points"]):
#                                     if i < len(target["temperature_points"]):
#                                         target["temperature_points"][i]["temp"] += tp["temp"] * factor
            
#             # Save each combination
#             output_data = {"load_case": comb_name, "load_data": combined_loads[comb_name]}
#             output_file = f"member_loads_{comb_name}.json"
#             output_path = os.path.join(JSON_FOLDER, "load_data", output_file)
#             with open(output_path, 'w') as f:
#                 json.dump(output_data, f, indent=4)
    
#     # 6. Convert member loads to nodal loads
#     if load_combinations:
#         # Load structure data
#         nodes_dict, members_dict, structure = load_structure(JSON_FOLDER)
#         combined_data = create_combined_structure_json(JSON_FOLDER)
#         combined_data = merge_structures(structure, combined_data)
        
#         # Create node and member mappings
#         node_coords = {node["id"]: {"x": node.get("x", 0), "y": node.get("y", 0), "z": node.get("z", 0)}
#                       for node in combined_data["nodes"]}
#         member_nodes = {member["id"]: {"start_node": member["start_node_id"], "end_node": member["end_node_id"]}
#                        for member in combined_data["members"]}
#         node_members = defaultdict(list)
#         for member_id, nodes in member_nodes.items():
#             node_members[nodes["start_node"]].append(member_id)
#             node_members[nodes["end_node"]].append(member_id)
        
#         # Process each combination
#         for combo_name in load_combinations.keys():
#             # Load member loads
#             member_loads_file = os.path.join(JSON_FOLDER, "load_data", f"member_loads_{combo_name}.json")
#             if not os.path.exists(member_loads_file): continue
            
#             with open(member_loads_file, 'r') as f:
#                 load_data = json.load(f)
            
#             nodal_loads = defaultdict(lambda: {"x": 0.0, "y": 0.0, "z": 0.0})
            
#             # Process each member group
#             for member_ids, load_configs in load_data.get("load_data", {}).items():
#                 try:
#                     member_id_list = eval(member_ids) if isinstance(member_ids, str) else member_ids
#                 except:
#                     member_id_list = [member_ids]
                
#                 for member_id in member_id_list:
#                     if member_id not in member_nodes: continue
                    
#                     start_node = member_nodes[member_id]["start_node"]
#                     end_node = member_nodes[member_id]["end_node"]
                    
#                     # Get coordinates and length
#                     start_coords = node_coords.get(start_node, {"x": 0, "y": 0, "z": 0})
#                     end_coords = node_coords.get(end_node, {"x": 0, "y": 0, "z": 0})
#                     dx = end_coords["x"] - start_coords["x"]
#                     dy = end_coords["y"] - start_coords["y"]
#                     dz = end_coords["z"] - start_coords["z"]
#                     length = math.sqrt(dx**2 + dy**2 + dz**2)
#                     if length == 0: continue
                    
#                     # Process each load config
#                     for load_config in load_configs:
#                         # Uniform loads
#                         if "uniform" in load_config:
#                             ux = load_config["uniform"].get("x", 0.0)
#                             uy = load_config["uniform"].get("y", 0.0)
#                             uz = load_config["uniform"].get("z", 0.0)
                            
#                             total_x = ux * length / 2
#                             total_y = uy * length / 2
#                             total_z = uz * length / 2
                            
#                             nodal_loads[start_node]["x"] += total_x
#                             nodal_loads[start_node]["y"] += total_y
#                             nodal_loads[start_node]["z"] += total_z
#                             nodal_loads[end_node]["x"] += total_x
#                             nodal_loads[end_node]["y"] += total_y
#                             nodal_loads[end_node]["z"] += total_z
                        
#                         # Point loads
#                         if "point" in load_config:
#                             px = load_config["point"].get("x", 0.0)
#                             py = load_config["point"].get("y", 0.0)
#                             pz = load_config["point"].get("z", 0.0)
#                             location = load_config["point"].get("location", 0.0)
                            
#                             if location < 0 or location > length: continue
                            
#                             a = location
#                             b = length - location
                            
#                             nodal_loads[start_node]["x"] += px * b / length
#                             nodal_loads[end_node]["x"] += px * a / length
#                             nodal_loads[start_node]["y"] += py * b / length
#                             nodal_loads[end_node]["y"] += py * a / length
#                             nodal_loads[start_node]["z"] += pz * b / length
#                             nodal_loads[end_node]["z"] += pz * a / length
            
#             # Save converted loads
#             node_id_to_name = {node["id"]: node["name"] for node in combined_data["nodes"]}
#             output_data = {}
#             for node_id, loads in nodal_loads.items():
#                 node_name = node_id_to_name.get(node_id, f"N{node_id}")
#                 output_data[node_name] = [
#                     node_id, combo_name,
#                     loads["x"], loads["y"], loads["z"],
#                     0.0, 0.0, 0.0
#                 ]
            
#             output_file = f"member_loads_converted_into_nodal_loads_{combo_name}.json"
#             output_path = os.path.join(JSON_FOLDER, "load_data", output_file)
#             with open(output_path, 'w') as f:
#                 json.dump(output_data, f, indent=4)
    
#     # 7. Create nodal mass combinations
#     if load_combinations:
#         for combo_name in load_combinations.keys():
#             # Initialize combined loads
#             combined_nodal_loads = {}
            
#             # Load direct nodal loads
#             member_loads_path = os.path.join(JSON_FOLDER, "load_data", f"nodal_loads_{combo_name}.json")
#             if os.path.exists(member_loads_path):
#                 with open(member_loads_path, 'r') as f:
#                     member_data = json.load(f)
#                     for node_name, loads in member_data.items():
#                         combined_nodal_loads[node_name] = loads.copy()
            
#             # Load converted nodal loads
#             converted_loads_path = os.path.join(JSON_FOLDER, "load_data", 
#                                               f"member_loads_converted_into_nodal_loads_{combo_name}.json")
#             if os.path.exists(converted_loads_path):
#                 with open(converted_loads_path, 'r') as f:
#                     converted_data = json.load(f)
#                     for node_name, loads in converted_data.items():
#                         if node_name not in combined_nodal_loads:
#                             if isinstance(loads, list):
#                                 combined_nodal_loads[node_name] = loads.copy()
#                             elif isinstance(loads, dict):
#                                 node_id = int(node_name.replace('N', '').replace('n', ''))
#                                 combined_nodal_loads[node_name] = [
#                                     node_id, combo_name,
#                                     loads.get("x", 0.0), loads.get("y", 0.0), loads.get("z", 0.0),
#                                     0.0, 0.0, 0.0
#                                 ]
#                         else:
#                             if isinstance(loads, list):
#                                 for i in range(2, 8):
#                                     combined_nodal_loads[node_name][i] += loads[i]
#                             elif isinstance(loads, dict):
#                                 combined_nodal_loads[node_name][2] += loads.get("x", 0.0)
#                                 combined_nodal_loads[node_name][3] += loads.get("y", 0.0)
#                                 combined_nodal_loads[node_name][4] += loads.get("z", 0.0)
            
#             # Save combined loads
#             if combined_nodal_loads:
#                 output_file = f"nodal_mass_loads_{combo_name}.json"
#                 output_path = os.path.join(JSON_FOLDER, "load_data", output_file)
#                 with open(output_path, 'w') as f:
#                     json.dump(combined_nodal_loads, f, indent=4)
    
#     print("\nOperation completed successfully!")
#     print(f"Ran the following parts: {run_parts}")


















def process_loads(JSON_FOLDER, surface_configurations=None, nodal_load_entries=None, 
                 load_combinations=None, all_element_loads=None, run_parts=None):
    """
    Consolidated function to handle all load processing with optimized repeated operations.
    """
    
    # ===== 1. Define repeated paths and operations once =====
    
    # Define load data directory path once
    load_data_dir = os.path.join(JSON_FOLDER, "load_data")
    os.makedirs(load_data_dir, exist_ok=True)  # Create directory once
    
    # Load structure data once (used in multiple places)
    nodes_dict, members_dict, structure = load_structure(JSON_FOLDER)
    combined_data = create_combined_structure_json(JSON_FOLDER)
    combined_data = merge_structures(structure, combined_data)
    # print(f'combined_data={combined_data}')
    # Create node mappings once
    node_id_to_name = {node["id"]: node["name"] for node in combined_data["nodes"]}
    node_name_to_id = {node["name"]: node["id"] for node in combined_data["nodes"]}
    
    # Create node coordinate mapping once
    node_coords = {node["id"]: {"x": node.get("x", 0), "y": node.get("y", 0), "z": node.get("z", 0)}
                  for node in combined_data["nodes"]}
    # print(f'node_coords={node_coords}')
    formatted_node_coords = {
        node["id"]: {
            "x": node["x"],
            "y": node["y"],
            "z": node["z"]
        }
        for node in combined_data["nodes"]
    }
    # Create member-node mapping once
    member_nodes = {member["id"]: {"start_node": member["start_node_id"], "end_node": member["end_node_id"]}
                   for member in combined_data["members"]}
    all_load_case_names = set()
    # Add "self_weight" as a default load case
    all_load_case_names.add("self_weight")  # <-- This ensures it's always included
    if surface_configurations is not None and surface_configurations:
        for config_data in surface_configurations.values():
            if "load_case_names" in config_data:
                all_load_case_names.update(config_data["load_case_names"])
    print(f'all_load_case_names={all_load_case_names}')

    # ===== 2. Convert shell element load to nodal load =====
    # ===== from f'mesh_data_with_predefined_points{numbering}.json' to f"nodal_loads_{load_case_name} =====
    if surface_configurations is not None and surface_configurations:
        for config_name, config_data in surface_configurations.items():           
            numbering = list(surface_configurations.keys()).index(config_name) + 1
            
            mesh_file = os.path.join(JSON_FOLDER, 'wall', f'mesh_data_with_predefined_points{numbering}.json')
            if not os.path.exists(mesh_file):
                raise FileNotFoundError(f"Mesh file not found at: {mesh_file}")
            
            with open(mesh_file, 'r') as f:
                mesh_data = json.load(f)
            
            # Use pre-loaded structure data
            all_nodes = {node["name"]: {"x": node["x"], "y": node["y"], "z": node["z"]} 
                        for node in combined_data["nodes"]}
            
            output_data_all = {}
            element_areas = {}
            
            # First pass to calculate element areas
            for elem_name, elem_data in mesh_data["elements"].items():
                node_names = elem_data["nodes"]
                if len(node_names) not in [3, 4]: continue
                
                node_coords = [np.array([all_nodes[name]['x'], all_nodes[name]['y'], all_nodes[name]['z']])
                             for name in node_names]
                result = surface_element_calculator(node_coords, pressure=1.0)
                element_areas[elem_name] = {"area": result['area'], "load_cases": {}}
            
            # Process each load case
            for load_case_name, pressure in zip(config_data["load_case_names"], config_data["pressures"]):
                nodal_loads = {}
                total_force = [0.0, 0.0, 0.0]
                
                for elem_name, elem_data in mesh_data["elements"].items():
                    node_names = elem_data["nodes"]
                    if len(node_names) not in [3, 4]: continue
                    
                    node_coords = [np.array([all_nodes[name]['x'], all_nodes[name]['y'], all_nodes[name]['z']])
                                 for name in node_names]
                    result = surface_element_calculator(node_coords, pressure=pressure)
                    
                    element_areas[elem_name]["load_cases"][load_case_name] = {
                        "pressure": pressure,
                        "total_force": -pressure * element_areas[elem_name]["area"]
                    }
                    
                    for i, node_name in enumerate(node_names):
                        if node_name not in nodal_loads:
                            nodal_loads[node_name] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        
                        nodal_loads[node_name][0] += result['forces'][i*3]
                        nodal_loads[node_name][1] += result['forces'][i*3 + 1]
                        nodal_loads[node_name][2] += result['forces'][i*3 + 2]
                        total_force[0] += result['forces'][i*3]
                        total_force[1] += result['forces'][i*3 + 1]
                        total_force[2] += result['forces'][i*3 + 2]
                
                # Save load case
                output_data = {}
                for node_name, forces in nodal_loads.items():
                    output_data[node_name] = [
                        node_name, load_case_name, 
                        forces[0], forces[1], forces[2], 
                        forces[3], forces[4], forces[5]
                    ]
                
                output_path = os.path.join(load_data_dir, f"nodal_loads_{load_case_name}.{numbering}.json")
                with open(output_path, 'w') as f:
                    json.dump(output_data, f, indent=4)
                
                output_data_all.update(output_data)
            
            # Save element areas
            areas_path = os.path.join(load_data_dir, f"element_areas_{numbering}.json")
            with open(areas_path, 'w') as f:
                json.dump(element_areas, f, indent=4)
    
    # ===== 3. Process self-weight =====
    # ===== from "element_data.json" to f"nodal_loads_{'_'.join(load_case_names)}.json"  =====
    load_case_names = "self_weight"
    if isinstance(load_case_names, str):
        load_case_names = [load_case_names]
    
    element_data_path = os.path.join(JSON_FOLDER, "element_data.json")
    if not os.path.exists(element_data_path):
        raise FileNotFoundError(f"Element data file not found at: {element_data_path}")
    
    with open(element_data_path, 'r') as f:
        element_data = json.load(f)
    
    # Calculate self-weight using pre-loaded node mappings
    nodal_loads = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    total_weight = 0.0
    
    for element in element_data["elements"]:
        if "unit_weight" not in element or element["unit_weight"] is None:
            continue
        
        unit_weight = element["unit_weight"]
        area = element["area"]
        length = element["length"]
        weight = unit_weight * area * length

        print(f'{unit_weight} * {area} * {length} = {weight}')
        total_weight += weight
        half_weight = weight / 2
        
        nodal_loads[element["node_i_id"]][2] -= half_weight
        nodal_loads[element["node_j_id"]][2] -= half_weight
    
    # Save self-weight loads
    output_data = {}
    for load_case in load_case_names:
        for node_id, forces in nodal_loads.items():
            node_name = node_id_to_name.get(node_id, str(node_id))
            output_data[f"{node_name}"] = [
                node_name, load_case, 
                forces[0], forces[1], forces[2], 
                forces[3], forces[4], forces[5]
            ]
    
    output_path = os.path.join(load_data_dir, f"nodal_loads_{'_'.join(load_case_names)}.json")
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=4)
    
    # ===== 4. nodal_load_entries nodal load  =====
    # ===== from nodal_load_entries tof"nodal_loads_{load_case}.json" =====
    if nodal_load_entries:
        load_cases = defaultdict(list)
        for entry in nodal_load_entries:
            load_cases[entry[1]].append(entry)
        
        for load_case, entries in load_cases.items():
            output_data = {entry[0]: entry for entry in entries}
            output_path = os.path.join(load_data_dir, f"nodal_loads_{load_case}.json")
            with open(output_path, 'w') as f:
                json.dump(output_data, f, indent=4)

    # ===== 5. Convert nodal load to load combinations =====
    # =====from f.startswith(f"nodal_loads_{load_case}") and f.endswith(".json") to f"nodal_loads_{combo_name}.json" =====

    if load_combinations:
        # Find all nodal load files
        nodal_load_files = []
        for f in os.listdir(load_data_dir):
            for load_case in all_load_case_names:
                if f.startswith(f"nodal_loads_{load_case}") and f.endswith(".json"):
                    nodal_load_files.append(f)

        # print(f'nodal_load_file31={nodal_load_files}')
        # Load all required load cases
        load_case_data = {}
        for filename in nodal_load_files:
            filepath = os.path.join(load_data_dir, filename)
            with open(filepath, 'r') as f:
                data = json.load(f)


            for load_case in all_load_case_names:
                if filename.startswith(f"nodal_loads_{load_case}"):
                    if load_case not in load_case_data:
                        load_case_data[load_case] = {}
                    
                    for key, new_values in data.items():
                        if key not in load_case_data[load_case]:
                            # If key doesn't exist, just store the new values
                            load_case_data[load_case][key] = new_values
                        else:
                            # If key exists, sum the numerical values
                            current_values = load_case_data[load_case][key]
                            combined_values = []
                            
                            # Assuming the structure is [str, str, float, float, float, float, float, float]
                            for i in range(len(current_values)):
                                if isinstance(current_values[i], (int, float)) and isinstance(new_values[i], (int, float)):
                                    # Sum numerical values
                                    combined_values.append(current_values[i] + new_values[i])
                                else:
                                    # Keep the first occurrence of non-numerical values
                                    combined_values.append(current_values[i])
                            
                            load_case_data[load_case][key] = combined_values
                    break
        print(f'load_case_data={load_case_data}')
        # Process each combination
        for combo_name, combo_load_cases in load_combinations.items():
            combined_loads = {}
            all_nodes_in_combo = set()
            
            # Initialize all nodes with zero loads
            for load_case, _ in combo_load_cases:
                for node_key in load_case_data[load_case].keys():
                    node_name = load_case_data[load_case][node_key][0]
                    all_nodes_in_combo.add(node_name)
            
            for node_name in all_nodes_in_combo:
                combined_loads[node_name] = {
                    "Fx": 0.0, "Fy": 0.0, "Fz": 0.0,
                    "Mx": 0.0, "My": 0.0, "Mz": 0.0
                }

            if combo_name=="Comb2":
                print(f'combo_name, combo_load_cases = {combo_name}, {combo_load_cases}')

                # Print like 1.2*DL + 1.6*LL
                expr = ' + '.join([f'{factor}*{load}' for load, factor in combo_load_cases])
                print(f'{expr} = ?')

                print(f'combo_load_cases = {combo_load_cases}')
            # print(f'load_case_data={load_case_data}')
            # Apply load factors
            for load_case, factor in combo_load_cases:
                for node_key, load_values in load_case_data[load_case].items():
                    node_name = load_values[0]
                    if combo_name=="Comb2":
                        if node_name == "N10006":
                            print(f'node_key : load_values : factor = {node_key} : {load_values} : {factor}')
                    if node_name not in combined_loads:
                        continue

                    combined_loads[node_name]["Fx"] += load_values[2] * factor
                    combined_loads[node_name]["Fy"] += load_values[3] * factor
                    combined_loads[node_name]["Fz"] += load_values[4] * factor
                    combined_loads[node_name]["Mx"] += load_values[5] * factor
                    combined_loads[node_name]["My"] += load_values[6] * factor
                    combined_loads[node_name]["Mz"] += load_values[7] * factor

                    # if node_name == "N10006":
                    #     print(f'combined_loads = {combined_loads[node_name]}')

            # Format output

            output_data = {}
            for node_name, loads in combined_loads.items():
                output_data[node_name] = [
                        node_name_to_id.get(node_name),
                        combo_name,
                        loads["Fx"], loads["Fy"], loads["Fz"],
                        loads["Mx"], loads["My"], loads["Mz"]
                    ]

            
            # Save combination
            output_path = os.path.join(load_data_dir, f"nodal_loads_{combo_name}.json")
            with open(output_path, 'w') as f:
                json.dump(output_data, f, indent=4)
    
    
    # ===== 6. Convert member load to member load combination =====
    # =====from all_element_loads to f"member_loads_{comb_name}.json" =====
    if load_combinations and all_element_loads:
        combined_loads = {}
        all_loads = {}
        
        # Combine all element loads
        for load_dict in all_element_loads:
            for elem_ids, load_list in load_dict.items():
                if elem_ids not in all_loads:
                    all_loads[elem_ids] = []
                all_loads[elem_ids].extend(load_list)
        
        # Process each combination
        for comb_name, case_factors in load_combinations.items():
            combined_loads[comb_name] = {}
            
            # Initialize structure
            for elem_ids in all_loads.keys():
                combined_loads[comb_name][elem_ids] = [{
                    "LoadCase": [comb_name],
                    "uniform": {"x": 0.0, "y": 0.0, "z": 0.0},
                    "point": {"x": 0.0, "y": 0.0, "z": 0.0, "location": 0.0},
                    "temperature_points": []
                }]
            
            # Apply factors
            for case_name, factor in case_factors:
                for elem_ids, load_list in all_loads.items():
                    for load_config in load_list:
                        if case_name in load_config["LoadCase"]:
                            target = combined_loads[comb_name][elem_ids][0]
                            
                            # Uniform loads
                            if "uniform" in load_config:
                                for dim in ["x", "y", "z"]:
                                    target["uniform"][dim] += load_config["uniform"].get(dim, 0.0) * factor
                            
                            # Point loads
                            if "point" in load_config:
                                for dim in ["x", "y", "z"]:
                                    target["point"][dim] += load_config["point"].get(dim, 0.0) * factor
                                if target["point"]["location"] == 0.0 and load_config["point"].get("location", 0.0) != 0.0:
                                    target["point"]["location"] = load_config["point"]["location"]
                            
                            # Temperature loads
                            if "temperature_points" in load_config:
                                if not target["temperature_points"]:
                                    target["temperature_points"] = [
                                        {"temp": 0.0, "y": tp["y"]} 
                                        for tp in load_config["temperature_points"]
                                    ]
                                for i, tp in enumerate(load_config["temperature_points"]):
                                    if i < len(target["temperature_points"]):
                                        target["temperature_points"][i]["temp"] += tp["temp"] * factor
            
            # Save each combination
            output_data = {"load_case": comb_name, "load_data": combined_loads[comb_name]}
            output_path = os.path.join(load_data_dir, f"member_loads_{comb_name}.json")
            with open(output_path, 'w') as f:
                json.dump(output_data, f, indent=4)
    
    # ===== 7. Convert member loads to nodal loads =====
    # ===== f"member_loads_converted_into_nodal_loads_{combo_name}.json" =====
    if load_combinations:
        # Create node to members mapping
        node_members = defaultdict(list)
        for member_id, nodes in member_nodes.items():
            node_members[nodes["start_node"]].append(member_id)
            node_members[nodes["end_node"]].append(member_id)
        
        # Process each combination
        for combo_name in load_combinations.keys():
            # Load member loads
            member_loads_path = os.path.join(load_data_dir, f"member_loads_{combo_name}.json")
            if not os.path.exists(member_loads_path): continue
            
            with open(member_loads_path, 'r') as f:
                load_data = json.load(f)
            
            nodal_loads = defaultdict(lambda: {"x": 0.0, "y": 0.0, "z": 0.0})
            
            # Process each member group
            for member_ids, load_configs in load_data.get("load_data", {}).items():
                try:
                    member_id_list = eval(member_ids) if isinstance(member_ids, str) else member_ids
                except:
                    member_id_list = [member_ids]
                
                for member_id in member_id_list:
                    if member_id not in member_nodes:
                        continue

                    start_node = member_nodes[member_id]["start_node"]
                    end_node = member_nodes[member_id]["end_node"]
                    start_coords = formatted_node_coords.get(start_node)
                    end_coords = formatted_node_coords.get(end_node)

                    dx = end_coords["x"] - start_coords["x"]
                    dy = end_coords["y"] - start_coords["y"]
                    dz = end_coords["z"] - start_coords["z"]

                    length = math.sqrt(dx**2 + dy**2 + dz**2)

                    if length == 0: continue
                    
                    # Process each load config
                    for load_config in load_configs:
                        # Uniform loads
                        if "uniform" in load_config:
                            ux = load_config["uniform"].get("x", 0.0)
                            uy = load_config["uniform"].get("y", 0.0)
                            uz = load_config["uniform"].get("z", 0.0)
                            
                            total_x = ux * length / 2
                            total_y = uy * length / 2
                            total_z = uz * length / 2
                            
                            nodal_loads[start_node]["x"] += total_x
                            nodal_loads[start_node]["y"] += total_y
                            nodal_loads[start_node]["z"] += total_z
                            nodal_loads[end_node]["x"] += total_x
                            nodal_loads[end_node]["y"] += total_y
                            nodal_loads[end_node]["z"] += total_z
                        
                        # Point loads
                        if "point" in load_config:
                            px = load_config["point"].get("x", 0.0)
                            py = load_config["point"].get("y", 0.0)
                            pz = load_config["point"].get("z", 0.0)
                            location = load_config["point"].get("location", 0.0)
                            
                            if location < 0 or location > length: continue
                            
                            a = location
                            b = length - location
                            
                            nodal_loads[start_node]["x"] += px * b / length
                            nodal_loads[end_node]["x"] += px * a / length
                            nodal_loads[start_node]["y"] += py * b / length
                            nodal_loads[end_node]["y"] += py * a / length
                            nodal_loads[start_node]["z"] += pz * b / length
                            nodal_loads[end_node]["z"] += pz * a / length
            
            # Save converted loads
            output_data = {}
            for node_id, loads in nodal_loads.items():
                node_name = node_id_to_name.get(node_id, f"N{node_id}")
                output_data[node_name] = [
                    node_id, combo_name,
                    loads["x"], loads["y"], loads["z"],
                    0.0, 0.0, 0.0
                ]
            
            output_path = os.path.join(load_data_dir, f"member_loads_converted_into_nodal_loads_{combo_name}.json")
            with open(output_path, 'w') as f:
                json.dump(output_data, f, indent=4)
    
    # ===== 8. Create nodal mass combinations =====
    # if load_combinations:
    #     for combo_name in load_combinations.keys():
    #         combined_nodal_loads = {}
            
    #         # Load direct nodal loads
    #         member_loads_path = os.path.join(load_data_dir, f"nodal_loads_{combo_name}.json")
    #         if os.path.exists(member_loads_path):
    #             with open(member_loads_path, 'r') as f:
    #                 member_data = json.load(f)
    #                 for node_name, loads in member_data.items():
    #                     combined_nodal_loads[node_name] = loads.copy()
            
    #         # Load converted nodal loads
    #         converted_loads_path = os.path.join(load_data_dir, 
    #                                           f"member_loads_converted_into_nodal_loads_{combo_name}.json")
    #         if os.path.exists(converted_loads_path):
    #             with open(converted_loads_path, 'r') as f:
    #                 converted_data = json.load(f)
    #                 for node_name, loads in converted_data.items():
    #                     if node_name not in combined_nodal_loads:
    #                         if isinstance(loads, list):
    #                             combined_nodal_loads[node_name] = loads.copy()
    #                         elif isinstance(loads, dict):
    #                             node_id = int(node_name.replace('N', '').replace('n', ''))
    #                             combined_nodal_loads[node_name] = [
    #                                 node_id, combo_name,
    #                                 loads.get("x", 0.0), loads.get("y", 0.0), loads.get("z", 0.0),
    #                                 0.0, 0.0, 0.0
    #                             ]
    #                     else:
    #                         if isinstance(loads, list):
    #                             for i in range(2, 8):
    #                                 combined_nodal_loads[node_name][i] += loads[i]
    #                         elif isinstance(loads, dict):
    #                             combined_nodal_loads[node_name][2] += loads.get("x", 0.0)
    #                             combined_nodal_loads[node_name][3] += loads.get("y", 0.0)
    #                             combined_nodal_loads[node_name][4] += loads.get("z", 0.0)
            
    #         # Save combined loads
    #         if combined_nodal_loads:
    #             output_path = os.path.join(load_data_dir, f"center_mass_loads_{combo_name}.json")
    #             with open(output_path, 'w') as f:
    #                 json.dump(combined_nodal_loads, f, indent=4)
    
    if load_combinations:
        for combo_name in load_combinations.keys():
            combined_nodal_loads = {}
            
            # Load direct nodal loads
            member_loads_path = os.path.join(load_data_dir, f"nodal_loads_{combo_name}.json")
            if os.path.exists(member_loads_path):
                with open(member_loads_path, 'r') as f:
                    member_data = json.load(f)
                    for node_name, loads in member_data.items():
                        combined_nodal_loads[node_name] = loads.copy()
            
            # Load converted nodal loads
            converted_loads_path = os.path.join(load_data_dir, 
                                            f"member_loads_converted_into_nodal_loads_{combo_name}.json")
            if os.path.exists(converted_loads_path):
                with open(converted_loads_path, 'r') as f:
                    converted_data = json.load(f)
                    for node_name, loads in converted_data.items():
                        if node_name not in combined_nodal_loads:
                            if isinstance(loads, list):
                                combined_nodal_loads[node_name] = loads.copy()
                        else:
                            if isinstance(loads, list):
                                for i in range(2, 8):
                                    combined_nodal_loads[node_name][i] += loads[i]
            
            # Convert to the new format with separate coordinates
            formatted_loads = {}
            for node_name, loads in combined_nodal_loads.items():
                # Find the node coordinates
                coords = {}
                for node in combined_data['nodes']:
                    if node['name'] == node_name:
                        coords = {
                            "x": node['x'],
                            "y": node['y'],
                            "z": node['z']
                        }
                        break
                
                # Create the new format
                formatted_loads[node_name] = {
                    "load_data": loads[:8],  # Keep only the first 8 elements (ID, combo, and 6 load components)
                    "coordinates": coords
                }
            
            # Save formatted loads
            if formatted_loads:
                output_path = os.path.join(load_data_dir, f"center_mass_loads_{combo_name}.json")
                with open(output_path, 'w') as f:
                    json.dump(formatted_loads, f, indent=4)
    
    print("\nOperation completed successfully!")
    print(f"Ran the following parts: {run_parts}")



# # Example usage focusing on the main operations
# if __name__ == "__main__":
#     JSON_FOLDER = "model_data"  # Path to JSON files
#     load_combinations = {
#         "Comb2": [("DL", 1.2), ("LL", 1.6)],
#         "Comb1": [("DL", 1.4)]
#     }
    
#     # Example 1: Apply nodal loads
#     applied_nodal_loads = process_structure_loads(
#         JSON_FOLDER, 
#         operation_type="nodal_loads",
#         load_combinations=load_combinations,
#         numbering=1,
#         load_combination="Comb2"
#     )
#     print("Applied nodal loads:", applied_nodal_loads)
    
#     # Example 2: Apply member loads (uniform and point loads)
#     member_loads = [
#         # Uniform load
#         {
#             'element_names': ["main_beam1", "main_beam2"],
#             'load_type': 'beamUniform',
#             'load_case': 'DL',
#             'wy': -10.0,  # Vertical load (negative = downward)
#             'wz': -5.0   # Lateral load
#         },
#         # Point load
#         {
#             'element_names': ["cantilever_beam"],
#             'load_type': 'beamPoint',
#             'Py': -15.0,  # Vertical point load
#             'x': 0.7       # 70% along the beam length
#         }
#     ]
#     process_structure_loads(
#         JSON_FOLDER,
#         operation_type="member_loads",
#         load_pattern_tag=1,
#         loads=member_loads,
#         load_combinations=load_combinations,
#         load_combination="Comb2"
#     )
    
#     # Example 3: Apply masses using node names
#     process_structure_loads(
#         JSON_FOLDER,
#         operation_type="masses",
#         loaded_nodes=["base_node1", "base_node2", "roof_node1", "roof_node2"],
#         m_1=1500.0,  # Total mass (kg)
#         mass_distribution=[0.4, 0.4, 0.1, 0.1]  # More mass at base nodes
#     )
    
#     # Example 4: Thermal and mechanical combined loads
#     combined_loads = [
#         # Thermal gradient load
#         {
#             'element_names': ["exposed_beam1", "exposed_beam2"],
#             'load_type': 'beamThermal',
#             'temp_points': [20.0, -0.1, 25.0, 0.0, 30.0, 0.1]  # Temp (C) and y-coord pairs
#         },
#         # Simultaneous uniform load
#         {
#             'element_names': ["exposed_beam1", "exposed_beam2"],
#             'load_type': 'beamUniform',
#             'wy': -3.0  # Snow load
#         }
#     ]
#     process_structure_loads(
#         JSON_FOLDER,
#         operation_type="member_loads",
#         load_pattern_tag=2,  # Different pattern for thermal effects
#         loads=combined_loads
#     )



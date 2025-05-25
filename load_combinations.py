# from build_model import build_model
from import_ import *
from Grid_and_structure_creation import load_structure
# from test import create_combined_structure_json
from wall_meshing import create_combined_structure_json
from units import *



    
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

# =============================================
# 12. apply_shell_load
# =============================================


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

def apply_shell_load(JSON_FOLDER, load_case_names, pressures, numbering=1):
    """
    Apply shell load on each element of mesh_data_with_predefined_points{numbering}.json for multiple load cases,
    convert the surface_load to nodal load, and save them all in one file.
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
    # combined_file = os.path.join(JSON_FOLDER, "combined_structure_data.json")
    
    # Create a dictionary of all nodes with their coordinates
    all_nodes = {}
    # print(f'combined_data["nodes"]={combined_data["nodes"]}')
    for node in combined_data["nodes"]:
        all_nodes[node["name"]] = {
            "x": node["x"],
            "y": node["y"],
            "z": node["z"]
        }
    
    # Initialize output data structure that will contain all load cases
    output_data = {}
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
        
        # Store the load case data in the output structure
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
        
        print(f"Processed load case: {load_case_name}")
    
    # Save all nodal loads to a single JSON file
    output_file = f"nodal_loads.{numbering}.json"
    output_path = os.path.join(JSON_FOLDER, "load_data", output_file)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=4)

    # Save element areas with pressure information
    areas_file = f"element_areas.{numbering}.json"
    areas_path = os.path.join(JSON_FOLDER, "load_data", areas_file)
    with open(areas_path, 'w') as f:
        json.dump(element_areas, f, indent=4)
    
    print(f"\nSuccessfully created combined nodal loads file at: {output_path}")
    print(f"Element areas with pressure information saved to: {areas_path}")
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
    
    # 6. Save the output to a JSON file
    output_file = "nodal_loads.self_weight.json"
    output_path = os.path.join(JSON_FOLDER, "load_data", output_file)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=4)
    
    print(f"Successfully saved self-weight nodal loads to: {output_path}")
    return True

def process_structure_loads(JSON_FOLDER, operation_type, load_cases=None, load_combinations=None, numbering=1, 
                          load_combination="Comb2", load_pattern_tag=None, time_series_tag=None, 
                          loads=None, loaded_nodes=None, m_1=None, mass_distribution=None):
    """
    Unified function to handle all structure loading operations including:
    - Node name to ID mapping
    - Element name to ID mapping
    - Nodal loads application
    - Member loads application
    - Mass application
    
    Parameters:
    -----------
    JSON_FOLDER : str
        Path to the folder containing JSON files
    operation_type : str
        Type of operation to perform ('nodal_loads', 'member_loads', 'masses')
    load_combinations : dict, optional
        Dictionary of load combinations and their factors
    numbering : int, optional
        The numbering identifier for the JSON file (default: 1)
    load_combination : str, optional
        The load combination to use (default: "Comb2")
    load_pattern_tag : int, optional
        Tag for the load pattern (required for member loads)
    time_series_tag : int, optional
        Tag for the time series (None for static loads)
    loads : list of dict, optional
        List of load specifications (required for member loads)
    loaded_nodes : list, optional
        List of node tags to apply mass to (required for masses)
    m_1 : float, optional
        Total mass to distribute (required for masses)
    mass_distribution : list, optional
        List of factors for each node (default: equal distribution)
        
    Returns:
    --------
    Varies by operation type:
    - For nodal loads: list of tuples (node_id, load_values)
    - For member loads: None
    - For masses: None
    """
    print(f'operation_type={operation_type}')
    def get_node_name_to_id_mapping():
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
    
    def get_element_name_to_id_mapping():
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
    
    def apply_member_load(elements, load_type, **load_params):
        """Helper function to apply member load with existence check"""
        # Create parameter list
        params = []
        for k, v in load_params.items():
            params += [f'-{k}', v]
        
        if isinstance(elements, (list, tuple)):
            for el in elements:
                try:
                    # Check if element exists by trying to get its nodes
                    ops.eleNodes(el)
                    ops.eleLoad('-ele', el, '-type', load_type, *params)
                except:
                    print(f"Warning: Element {el} does not exist - skipping")
                    continue
        else:
            try:
                ops.eleNodes(elements)
                ops.eleLoad('-ele', elements, '-type', load_type, *params)
            except:
                print(f"Warning: Element {elements} does not exist - skipping")
    
   
    # if operation_type == "nodal_loads":
    #     nodal_loads_file = os.path.join(JSON_FOLDER, "load_data", f"nodal_loads.{numbering}.json")
    #     print(f"Debug: Looking for nodal loads file at {nodal_loads_file}")
        
    #     if load_combination not in load_combinations:
    #         raise ValueError(f"Load combination '{load_combination}' not found")
        
    #     selected_combo = load_combinations[load_combination]
    #     print(f"Debug: Selected load combination - {selected_combo}")
        
    #     # Create load factors mapping including the special "load_cases" container
    #     load_factors = {}
    #     for load_case, factor in selected_combo:
    #         print("load_case100:", load_case)
    #         load_factors[load_case] = factor
        
    #     print(f"Debug: Effective load factors - {load_factors}")
        
    #     node_name_to_id = get_node_name_to_id_mapping()
    #     if not node_name_to_id:
    #         print("Warning: Could not obtain node name to ID mapping. Using node names directly.")
    #     else:
    #         print(f"Debug: Node name to ID mapping contains {len(node_name_to_id)} entries")
        
    #     try:
    #         with open(nodal_loads_file, 'r') as f:
    #             nodal_loads_data = json.load(f)
    #         print(f"Debug: Successfully loaded nodal loads data with {len(nodal_loads_data)} entries")
    #     except FileNotFoundError:
    #         print(f"Nodal loads file not found: {nodal_loads_file}")
    #         return []
    #     except json.JSONDecodeError:
    #         print(f"Error parsing JSON from file: {nodal_loads_file}")
    #         return []
        
    #     combined_loads = {}
    #     print("Debug: Starting to process nodal loads")
        
    #     for key, load_info in nodal_loads_data.items():
    #         print(f"\nDebug: Processing key: {key}")
    #         print(f"Debug: Full load_info: {load_info}")  # Print the entire load_info tuple
            
            

    #         node_name = load_info[0]
    #         load_case_identifier = load_info[1]
    #         loads = load_info[2:]
            
    #         print(f"Debug: Node: {node_name}, Case: {load_case_identifier}, Loads: {loads}")
    #         # print(f"Debug: Current node_name_to_id mapping: {node_name_to_id}")  # Show current mapping
            
    #         if node_name in node_name_to_id:
    #             node_tag = node_name_to_id[node_name]
    #             print(f"Debug: Found node {node_name} with tag {node_tag}")
    #         else:
    #             print(f"Warning: Node {node_name} not found in mapping (available nodes: {list(node_name_to_id.keys())}), skipping")
    #             continue
            

    #         if load_case_identifier in load_factors:
    #             factor = load_factors[load_case_identifier]
    #             print(f"Debug: Found factor {factor} for case {load_case_identifier}")
                
    #             factored_loads = [load * factor for load in loads]
    #             print(f"Debug: Factored loads: {factored_loads}")
                
    #             if node_tag in combined_loads:
    #                 print(f"Debug: Existing loads for node {node_tag}: {combined_loads[node_tag]}")
    #                 for i in range(len(combined_loads[node_tag])):
    #                     combined_loads[node_tag][i] += factored_loads[i]
    #                 print(f"Debug: Updated loads for node {node_tag}: {combined_loads[node_tag]}")
    #             else:
    #                 combined_loads[node_tag] = factored_loads
    #                 print(f"Debug: New loads for node {node_tag}: {combined_loads[node_tag]}")
    #         else:
    #             print(f"Warning: Load case identifier {load_case_identifier} not found in load_factors (available cases: {list(load_factors.keys())})")
        
    #     print(f"Debug: Combined loads prepared for {len(combined_loads)} nodes")

    #     for node_tag, load_values in combined_loads.items():
    #         try:
    #             print(f"Debug: Attempting to apply load to node {node_tag}: {load_values}")
    #             print(node_tag, *load_values)
    #             ops.load(node_tag, *load_values)

    #             print(f"Debug: Successfully applied load to node {node_tag}")
    #         except Exception as e:
    #             print(f"Error applying load to node {node_tag}: {e}")
        

    #     return True
    
    if operation_type == "nodal_loads":
        # Load primary nodal loads
        nodal_loads_file = os.path.join(JSON_FOLDER, "load_data", f"nodal_loads.{numbering}.json")
        print(f"Debug: Looking for nodal loads file at {nodal_loads_file}")
        
        if not os.path.exists(nodal_loads_file):
            print(f"Nodal loads file not found: {nodal_loads_file}")
            return False
        
        try:
            with open(nodal_loads_file, 'r') as f:
                nodal_loads_data = json.load(f)
            print(f"Debug: Loaded {len(nodal_loads_data)} nodal load entries")
        except Exception as e:
            print(f"Error loading nodal loads: {str(e)}")
            return False

        # Load self-weight if available
        self_weight_file = os.path.join(JSON_FOLDER, "load_data", "nodal_loads.self_weight.json")
        if os.path.exists(self_weight_file):
            try:
                with open(self_weight_file, 'r') as f:
                    self_weight_data = json.load(f)
                print(f"Debug: Loaded {len(self_weight_data)} self-weight entries")
            except Exception as e:
                print(f"Error loading self-weight: {str(e)}")
                self_weight_data = {}
        else:
            print("Debug: No self-weight file found")
            self_weight_data = {}

        # Process load combinations
        if load_combinations:
            if load_combination not in load_combinations:
                raise ValueError(f"Load combination '{load_combination}' not found")
            load_factors = dict(load_combinations[load_combination])
            print(f"Debug: Using load factors: {load_factors}")
        else:
            load_factors = {}
            print("Debug: No load combinations specified")

        # Get node mapping
        node_name_to_id = get_node_name_to_id_mapping()
        if not node_name_to_id:
            print("Warning: Empty node name to ID mapping")

        # Combine all loads
        combined_loads = {}
        
        # Process primary nodal loads
        for load_info in nodal_loads_data.values():
            node_name = load_info[0]
            load_case = load_info[1]
            loads = load_info[2:]
            
            if node_name not in node_name_to_id:
                print(f"Warning: Node {node_name} not found in mapping")
                continue
                
            node_tag = node_name_to_id[node_name]
            factor = load_factors.get(load_case, 1.0)
            factored_loads = [load * factor for load in loads]
            
            if node_tag in combined_loads:
                combined_loads[node_tag] = [a+b for a,b in zip(combined_loads[node_tag], factored_loads)]
            else:
                combined_loads[node_tag] = factored_loads

        # Process self-weight loads
        for load_info in self_weight_data.values():
            node_name = load_info[0]
            load_case = load_info[1]  # "self_weight"
            loads = load_info[2:]
            
            if node_name not in node_name_to_id:
                print(f"Warning: Node {node_name} not found in mapping (self-weight)")
                continue
                
            node_tag = node_name_to_id[node_name]
            factor = load_factors.get(load_case, 1.0)
            factored_loads = [load * factor for load in loads]
            
            if node_tag in combined_loads:
                combined_loads[node_tag] = [a+b for a,b in zip(combined_loads[node_tag], factored_loads)]
            else:
                combined_loads[node_tag] = factored_loads

        # Apply all combined loads
        for node_tag, load_values in combined_loads.items():
            try:
                ops.load(node_tag, *load_values)
                print(f"Debug: Applied load to node {node_tag}: {load_values}")
            except Exception as e:
                print(f"Error applying load to node {node_tag}: {str(e)}")

 

    elif operation_type == "member_loads":
        if not load_pattern_tag:
            raise ValueError("load_pattern_tag parameter required for member loads")
        if not loads:
            raise ValueError("loads parameter required for member loads")
            
        pattern_type = 'Plain'
        if time_series_tag is not None:
            ops.pattern(pattern_type, load_pattern_tag, time_series_tag)
        else:
            ops.pattern(pattern_type, load_pattern_tag, 'Linear')
            
        element_name_to_id = get_element_name_to_id_mapping()
        if not element_name_to_id:
            print("Warning: Could not obtain element name to ID mapping. Using names directly may cause errors.")
            
        if load_combinations:
            if load_combination not in load_combinations:
                raise ValueError(f"Load combination '{load_combination}' not found")
            
            selected_combo = load_combinations[load_combination]
            load_factors = {load_case: factor for load_case, factor in selected_combo}
            
            combined_loads = {}
            
            for load in loads:
                load_copy = load.copy()
                load_case = load_copy.get('load_case')
                if load_case not in load_factors:
                    continue
                
                factor = load_factors[load_case]
                element_names = load_copy['element_names']
                if isinstance(element_names, str):
                    element_names = [element_names]
                
                element_ids = []
                for name in element_names:
                    if name in element_name_to_id:
                        element_ids.append(element_name_to_id[name])
                    else:
                        print(f"Warning: Element name '{name}' not found in mapping, skipping")
                        continue
                
                if not element_ids:
                    continue
                    
                load_type = load_copy['load_type']
                
                if len(element_ids) > 1:
                    key = (tuple(sorted(element_ids)), load_type)
                else:
                    key = (element_ids[0], load_type)
                
                load_copy.pop('element_names')
                load_copy.pop('load_type')
                load_copy.pop('load_case', None)
                
                for k, v in load_copy.items():
                    if isinstance(v, (int, float)):
                        load_copy[k] = v * factor
                
                if key in combined_loads:
                    existing_load = combined_loads[key]
                    for k, v in load_copy.items():
                        if k in existing_load:
                            existing_load[k] += v
                        else:
                            existing_load[k] = v
                else:
                    combined_loads[key] = load_copy
            
            for (elements, load_type), load_params in combined_loads.items():
                try:
                    apply_member_load(elements, load_type, **load_params)
                except Exception as e:
                    print(f"Error applying load to elements {elements}: {e}")
        
        else:
            for load in loads:
                load_copy = load.copy()
                element_names = load_copy.pop('element_names')
                if isinstance(element_names, str):
                    element_names = [element_names]
                
                element_ids = []
                for name in element_names:
                    if name in element_name_to_id:
                        element_ids.append(element_name_to_id[name])
                    else:
                        print(f"Warning: Element name '{name}' not found in mapping, skipping")
                        continue
                
                if not element_ids:
                    continue
                    
                load_type = load_copy.pop('load_type')
                try:
                    apply_member_load(element_ids, load_type, **load_copy)
                except Exception as e:
                    print(f"Error applying load to elements {element_ids}: {e}")
    
    elif operation_type == "masses":
        # Validate inputs
        if loaded_nodes is None:
            raise ValueError("loaded_nodes parameter must be provided (list of node names)")
        if m_1 is None:
            raise ValueError("m_1 parameter (total mass) must be provided")
        if not isinstance(loaded_nodes, list):
            raise ValueError("loaded_nodes must be a list of node names")
        
        # Get node name to ID mapping
        node_name_to_id = get_node_name_to_id_mapping()
        if not node_name_to_id:
            raise ValueError("Could not load node name to ID mapping")
        
        # Convert node names to IDs
        node_ids = []
        invalid_nodes = []
        for node_name in loaded_nodes:
            if node_name in node_name_to_id:
                node_ids.append(node_name_to_id[node_name])
            else:
                invalid_nodes.append(node_name)
        
        if invalid_nodes:
            print(f"Warning: These nodes not found in mapping and will be skipped: {invalid_nodes}")
        
        if not node_ids:
            raise ValueError("No valid nodes found after name mapping")
        
        # Set default equal distribution if not provided
        if mass_distribution is None:
            mass_distribution = [1.0/len(node_ids)] * len(node_ids)
        
        # Validate mass distribution
        if len(mass_distribution) != len(node_ids):
            raise ValueError("mass_distribution length must match number of valid nodes")
        if abs(sum(mass_distribution) - 1.0) > 1e-6:
            raise ValueError("mass_distribution factors must sum to 1.0")
        
        # Apply masses
        for node_id, factor in zip(node_ids, mass_distribution):
            try:
                ops.mass(node_id, m_1*factor, m_1*factor, 0, 0, 0, 0)
                print(f"Applied mass {m_1*factor:.2f} to node ID {node_id}")
            except Exception as e:
                print(f"Error applying mass to node ID {node_id}: {str(e)}")


    else:
        raise ValueError(f"Unknown operation type: {operation_type}")

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



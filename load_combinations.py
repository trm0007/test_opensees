# from build_model import build_model
from import_ import *
from Grid_and_structure_creation import load_structure
# from test import create_combined_structure_json
# from wall_meshing import create_combined_structure_json
from units import *




def create_combined_structure_json(JSON_FOLDER):
    wall_dir = os.path.join(JSON_FOLDER, 'wall')

    # Load nodes data
    with open(os.path.join(JSON_FOLDER, 'node_input_data.json'), 'r') as f:
        nodes_data = json.load(f)

    # Load members data
    with open(os.path.join(JSON_FOLDER, 'members_input_data.json'), 'r') as f:
        members_data = json.load(f)

    # Load and combine all mesh data files
    mesh_data = {'elements': {}}
    mesh_files = []

    if os.path.exists(wall_dir):
        for filename in os.listdir(wall_dir):
            if filename.endswith('.json'):
                mesh_files.append(filename)

        for filename in mesh_files:
            mesh_data_path = os.path.join(wall_dir, filename)
            with open(mesh_data_path, 'r') as f:
                single_mesh = json.load(f)
                mesh_data['elements'].update(single_mesh.get('elements', {}))

    print("Number of nodes:", len(nodes_data))
    print("First 3 nodes:", nodes_data[:3])

    print("Number of members:", len(members_data))
    print("First 3 members:", members_data[:3])

    print("Number of mesh elements:", len(mesh_data['elements']))
    first_3_keys = list(mesh_data['elements'].keys())[:3]
    print("First 3 mesh elements:", {k: mesh_data['elements'][k] for k in first_3_keys})

    print("Mesh JSON files:", mesh_files)

    # Print data for R10001, R10002, R10003
    for rid in ['R10001', 'R10002', 'R10003']:
        if rid in mesh_data['elements']:
            print(f"Data for {rid}:", mesh_data['elements'][rid])

    return nodes_data, members_data, mesh_data

def get_node_name_to_id_mapping(JSON_FOLDER):
    """Helper function to get node name to ID mapping (structural nodes only)"""
    nodes_data, _, _ = create_combined_structure_json(JSON_FOLDER)  # Ignore members and mesh data
    
    node_name_to_id = {}
    for node in nodes_data:
        node_name_to_id[node["name"]] = node["id"]
    
    print(f"Loaded {len(node_name_to_id)} structural nodes from input data")
    return node_name_to_id

def get_element_name_to_id_mapping(JSON_FOLDER):
    """Helper function to get element name to ID mapping (structural members only)"""
    _, members_data, _ = create_combined_structure_json(JSON_FOLDER)  # Ignore nodes and mesh data
    
    element_name_to_id = {}
    for member in members_data:
        element_name_to_id[member["name"]] = member["id"]
    
    print(f"Loaded {len(element_name_to_id)} structural members from input data")
    return element_name_to_id



def parse_member_ids(member_ids):
    """Parse member IDs from various input formats into a list of integers."""
    if isinstance(member_ids, str):
        if member_ids.startswith('[') and member_ids.endswith(']'):
            try:
                return [int(id_str.strip()) for id_str in member_ids[1:-1].split(',') if id_str.strip()]
            except ValueError as e:
                raise ValueError(f"Invalid member ID format in string list: {member_ids} - {e}")
        try:
            return [int(member_ids)]
        except ValueError as e:
            raise ValueError(f"Invalid member ID string format: {member_ids} - {e}")
    elif isinstance(member_ids, (list, tuple)):
        try:
            return [int(mid) for mid in member_ids]
        except ValueError as e:
            raise ValueError(f"Invalid member ID in list: {member_ids} - {e}")
    elif isinstance(member_ids, int):
        return [member_ids]
    else:
        raise ValueError(f"Invalid member ID format. Expected str, list, or int, got {type(member_ids)}: {member_ids}")

def save_element_loads(all_element_loads, json_folder):
    """Save element loads to a JSON file in the specified folder."""
    if not isinstance(all_element_loads, (list, dict)):
        raise ValueError("all_element_loads must be a list or dictionary")
    
    # Create output directory for nodal loads
    nodal_loads_dir = os.path.join(json_folder, "load_data")
    os.makedirs(nodal_loads_dir, exist_ok=True)
    
    member_loads_file = os.path.join(nodal_loads_dir, "member_loads.json")
    try:
        with open(member_loads_file, 'w') as f:
            json.dump(all_element_loads, f, indent=4)
        print(f"Saved element loads to {member_loads_file}")
    except (IOError, TypeError) as e:
        raise IOError(f"Failed to save element loads: {e}")

def convert_member_loads_to_nodal(all_element_loads, json_folder):
    """
    Converts member loads (uniform and point) into equivalent nodal loads.
    Uses "LoadCase" from input as the key in output JSON.
    Output format: { "LoadCase": { "node_name": [Fx, Fy, Fz, Mx, My, Mz], ... } }
    """
    # Validate inputs
    if not os.path.exists(json_folder):
        raise FileNotFoundError(f"JSON folder not found: {json_folder}")
    
    # Load structure data
    try:
        nodes_data, members_data, mesh_data = create_combined_structure_json(json_folder)
    except Exception as e:
        raise RuntimeError(f"Failed to load structure data: {e}")
    
    # Save element loads first
    try:
        save_element_loads(all_element_loads, json_folder)
    except Exception as e:
        print(f"Warning: Could not save element loads: {e}")

    # Create node ID to coordinates mapping
    node_coords = {}
    for node in nodes_data:
        node_id = node.get("id")
        if node_id is None:
            continue
        node_coords[node_id] = {
            "x": node.get("x"),
            "y": node.get("y"),
            "z": node.get("z")
        }
    
    # Create member ID to node IDs mapping
    member_nodes = {}
    for member in members_data:
        member_id = member.get("id")
        if member_id is None:
            continue
        member_nodes[member_id] = {
            "start_node": member.get("start_node_id"),
            "end_node": member.get("end_node_id")
        }
    
    # Create node to members mapping
    node_members = defaultdict(list)
    for member_id, nodes in member_nodes.items():
        if nodes["start_node"] is not None:
            node_members[nodes["start_node"]].append(member_id)
        if nodes["end_node"] is not None:
            node_members[nodes["end_node"]].append(member_id)
    
    # Initialize dictionary to store all nodal loads
    all_nodal_loads = {}
    
    # Process each load combination in the all_element_loads list
    for load_data in all_element_loads:
        if not isinstance(load_data, dict):
            print(f"Warning: Skipping invalid load data format: {type(load_data)}")
            continue
            
        # Initialize a dictionary to store loads by LoadCase
        load_case_nodal_loads = defaultdict(lambda: defaultdict(lambda: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
        
        # Process each member group in the load data
        for member_ids_str, load_configs in load_data.items():
            if not load_configs:
                continue
                
            # Parse member IDs
            try:
                member_id_list = parse_member_ids(member_ids_str)
            except ValueError as e:
                print(f"Warning: {e} - skipping")
                continue
            
            for member_id in member_id_list:
                if member_id not in member_nodes:
                    print(f"Warning: Member ID {member_id} not found in structure - skipping")
                    continue
                    
                start_node = member_nodes[member_id].get("start_node")
                end_node = member_nodes[member_id].get("end_node")
                
                if start_node is None or end_node is None:
                    print(f"Warning: Member {member_id} has missing node references - skipping")
                    continue
                
                # Get node coordinates
                start_coords = node_coords.get(start_node)
                end_coords = node_coords.get(end_node)
                
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
                    if not isinstance(load_config, dict):
                        print(f"Warning: Skipping invalid load config format: {type(load_config)}")
                        continue
                    
                    if "LoadCase" not in load_config:
                        print("Warning: Load configuration missing LoadCase - skipping")
                        continue
                    
                    # Get the load case name (handle both single value and list)
                    load_cases = load_config["LoadCase"]
                    if isinstance(load_cases, str):
                        load_cases = [load_cases]
                    elif not isinstance(load_cases, (list, tuple)):
                        print(f"Warning: Invalid LoadCase format: {type(load_cases)} - skipping")
                        continue
                    
                    for load_case in load_cases:
                        if not isinstance(load_case, str):
                            print(f"Warning: Invalid LoadCase name type: {type(load_case)} - skipping")
                            continue
                        
                        # Process uniform loads
                        if "uniform" in load_config and isinstance(load_config["uniform"], dict):
                            uniform_load = load_config["uniform"]
                            ux = uniform_load.get("x", 0.0)
                            uy = uniform_load.get("y", 0.0)
                            uz = uniform_load.get("z", 0.0)
                            
                            # For uniform load, each node gets half the total load
                            total_x = ux * length
                            total_y = uy * length
                            total_z = uz * length
                            
                            # Get node names
                            start_node_name = next((n["name"] for n in nodes_data if n.get("id") == start_node), f"n{start_node}")
                            end_node_name = next((n["name"] for n in nodes_data if n.get("id") == end_node), f"n{end_node}")
                            
                            # Update nodal loads for this LoadCase
                            load_case_nodal_loads[load_case][start_node_name][0] += total_x / 2  # Fx
                            load_case_nodal_loads[load_case][start_node_name][1] += total_y / 2  # Fy
                            load_case_nodal_loads[load_case][start_node_name][2] += total_z / 2  # Fz
                            
                            load_case_nodal_loads[load_case][end_node_name][0] += total_x / 2    # Fx
                            load_case_nodal_loads[load_case][end_node_name][1] += total_y / 2    # Fy
                            load_case_nodal_loads[load_case][end_node_name][2] += total_z / 2    # Fz
                        
                        # Process point loads
                        if "point" in load_config and isinstance(load_config["point"], dict):
                            point_load = load_config["point"]
                            px = point_load.get("x", 0.0)
                            py = point_load.get("y", 0.0)
                            pz = point_load.get("z", 0.0)
                            location = point_load.get("location", 0.0)
                            
                            if location < 0 or location > length:
                                print(f"Warning: Point load location {location} out of bounds for member {member_id} - skipping")
                                continue
                            
                            # Calculate distribution factors
                            a = location
                            b = length - location
                            
                            # Get node names
                            start_node_name = next((n["name"] for n in nodes_data if n.get("id") == start_node), f"n{start_node}")
                            end_node_name = next((n["name"] for n in nodes_data if n.get("id") == end_node), f"n{end_node}")
                            
                            # Distribute point load to nodes
                            if length > 0:
                                load_case_nodal_loads[load_case][start_node_name][0] += px * b / length  # Fx
                                load_case_nodal_loads[load_case][end_node_name][0] += px * a / length   # Fx
                                load_case_nodal_loads[load_case][start_node_name][1] += py * b / length  # Fy
                                load_case_nodal_loads[load_case][end_node_name][1] += py * a / length   # Fy
                                load_case_nodal_loads[load_case][start_node_name][2] += pz * b / length  # Fz
                                load_case_nodal_loads[load_case][end_node_name][2] += pz * a / length   # Fz
        
        # Convert defaultdict to regular dict for JSON serialization
        for load_case in load_case_nodal_loads:
            all_nodal_loads[load_case] = dict(load_case_nodal_loads[load_case])
    
    # Create output directory for nodal loads
    nodal_loads_dir = os.path.join(json_folder, "load_data")
    os.makedirs(nodal_loads_dir, exist_ok=True)
    
    # Save ALL load cases to a SINGLE file
    output_file = "member_loads_converted_into_nodal_loads.json"
    output_path = os.path.join(nodal_loads_dir, output_file)
    
    try:
        with open(output_path, 'w') as f:
            json.dump(all_nodal_loads, f, indent=4)
        print(f"Saved ALL nodal loads to {output_path}")
    except (IOError, TypeError) as e:
        print(f"Error saving nodal loads: {e}")
    
    return all_nodal_loads


def surface_element_calculator(coords):
    """
    Placeholder function - you need to implement this based on your requirements
    Should return a dictionary with at least 'area' key
    """
    if len(coords) == 3:  # Triangle
        v1 = coords[1] - coords[0]
        v2 = coords[2] - coords[0]
        area = 0.5 * np.linalg.norm(np.cross(v1, v2))
    elif len(coords) == 4:  # Quadrilateral (approximate)
        v1 = coords[1] - coords[0]
        v2 = coords[2] - coords[0]
        v3 = coords[2] - coords[0]
        v4 = coords[3] - coords[0]
        area1 = 0.5 * np.linalg.norm(np.cross(v1, v2))
        area2 = 0.5 * np.linalg.norm(np.cross(v3, v4))
        area = area1 + area2
    else:
        raise ValueError(f"Unsupported number of nodes: {len(coords)}")
    
    return {'area': area}



def calculate_nodal_loads_from_mesh(nodes_data, mesh_data, nodal_load_entries, JSON_FOLDER):
    tolerance = 1e-6

    # Validate inputs
    if not nodes_data or not mesh_data:
        raise ValueError("nodes_data and mesh_data cannot be empty")
    
    if not os.path.exists(JSON_FOLDER):
        raise ValueError(f"JSON_FOLDER does not exist: {JSON_FOLDER}")

    # Additional nodal load entries to be applied

    # Create node name to ID mapping
    node_name_to_id = {node['name']: node['id'] for node in nodes_data}
    # Create node coordinates map
    node_coords_map = {node['name']: np.array([node['x'], node['y'], node['z']]) for node in nodes_data}
    
    # Initialize nodal loads data structure
    nodal_loads = defaultdict(lambda: defaultdict(lambda: np.zeros(6)))
    total_applied_force_per_load_case = defaultdict(float)
    total_member_self_weight = 0.0  # Total self-weight from all member elements
    total_mesh_self_weight = 0.0    # Total self-weight from all mesh elements

    # Validate mesh_data structure
    if 'elements' not in mesh_data:
        raise KeyError("mesh_data must contain 'elements' key")

    # Process mesh elements
    for elem_id, elem in mesh_data['elements'].items():
        required_keys = ['nodes', 'load_case_names', 'pressures']
        for key in required_keys:
            if key not in elem:
                print(f"Warning: Element {elem_id} missing key '{key}', skipping")
                continue
        
        element_nodes = elem['nodes']
        load_cases = elem['load_case_names']
        pressures = elem['pressures']
        
        
        if len(load_cases) != len(pressures):
            print(f"Warning: Element {elem_id} has mismatched load_cases and pressures, skipping")
            continue
        
        try:
            coords = []
            for node_name in element_nodes:
                # print(f'element_nodes={node_name}')
                if node_name not in node_coords_map:
                    raise KeyError(f"Node {node_name} not found in coordinates")
                coords.append(node_coords_map[node_name])
        except KeyError as e:
            print(f"Warning: {e}, skipping element {elem_id}")
            continue
            
        try:
            res = surface_element_calculator(coords)
            area = res['area']
        except Exception as e:
            print(f"Warning: Failed to calculate area for element {elem_id}: {e}")
            continue

        num_nodes = len(coords)
        if num_nodes == 0:
            continue

        for i, load_case in enumerate(load_cases):
            pressure = pressures[i]
            total_force = pressure * area
            total_applied_force_per_load_case[load_case] += total_force
            
            # Track mesh self-weight separately
            if load_case == "self_weight":
                total_mesh_self_weight += abs(total_force)

            force_per_node = total_force / num_nodes
            for node in element_nodes:
                if node not in node_name_to_id:
                    print(f"Warning: Node {node} not found in node_name_to_id mapping")
                    continue
                    
                node_id = node_name_to_id[node]
                nodal_loads[load_case][node_id] += np.array([0.0, 0.0, force_per_node, 0.0, 0.0, 0.0])



    # Process member self-weights from element_data.json
    element_data_path = os.path.join(JSON_FOLDER, "element_data.json")
    if os.path.exists(element_data_path):
        try:
            with open(element_data_path, 'r') as f:
                element_data = json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            print(f"Warning: Failed to load element_data.json: {e}")
            element_data = None
            
        if element_data and "elements" in element_data:
            # Track processed elements to avoid double counting
            processed_elements = set()
            
            for i, element in enumerate(element_data["elements"]):
                # Create unique element identifier
                element_id = element.get('id', i)
                
                # Skip if already processed (duplicate detection)
                if element_id in processed_elements:
                    print(f"Warning: Duplicate element detected (ID: {element_id}), skipping")
                    continue
                
                processed_elements.add(element_id)
                
                # Check required fields
                required_fields = ["unit_weight", "area", "length", "node_i_id", "node_j_id"]
                missing_fields = [field for field in required_fields if field not in element]
                
                if missing_fields:
                    print(f"Warning: Element {element_id} missing fields {missing_fields}, skipping")
                    continue
                    
                if element["unit_weight"] is None or element["area"] is None or element["length"] is None:
                    print(f"Warning: Element {element_id} has None values, skipping")
                    continue
                
                # Calculate total weight of the element
                unit_weight = float(element["unit_weight"])
                area = float(element["area"])
                length = float(element["length"])
                
                weight = unit_weight * area * length
                
                # Skip if weight is essentially zero
                if abs(weight) < tolerance:
                    continue
                
                total_member_self_weight += abs(weight)  # Track absolute value for verification
                half_weight = weight / 2
                
                # Get node IDs
                node_i_id = element["node_i_id"]
                node_j_id = element["node_j_id"]
                
                # Validate that nodes exist
                if node_i_id not in [node['id'] for node in nodes_data]:
                    print(f"Warning: Node ID {node_i_id} not found in nodes_data")
                    continue
                if node_j_id not in [node['id'] for node in nodes_data]:
                    print(f"Warning: Node ID {node_j_id} not found in nodes_data")
                    continue
                
                # Apply self-weight (downward, so negative Z)
                # Use abs(half_weight) to ensure consistent downward direction
                nodal_loads["self_weight"][node_i_id][2] -= abs(half_weight)
                nodal_loads["self_weight"][node_j_id][2] -= abs(half_weight)

    # Verify mesh element force distribution (excluding self_weight for now)
    for load_case, total_applied_force in total_applied_force_per_load_case.items():
        if load_case != "self_weight":
            total_nodal_force = sum(force[2] for force in nodal_loads[load_case].values())
            print(f"Mesh Load case: {load_case}, Applied force: {total_applied_force:.6f}, Nodal force sum: {total_nodal_force:.6f}")
            if abs(total_applied_force - total_nodal_force) >= tolerance:
                raise AssertionError(
                    f"Verification failed for mesh load case {load_case}: "
                    f"applied force {total_applied_force:.6f}, nodal force sum {total_nodal_force:.6f}, "
                    f"difference: {abs(total_applied_force - total_nodal_force):.6f}"
                )
    # Process additional nodal load entries
    print("\n" + "="*50)
    print("PROCESSING ADDITIONAL NODAL LOADS")
    print("="*50)
    
    processed_point_loads = defaultdict(lambda: defaultdict(float))
    
    for entry in nodal_load_entries:
        node_name, load_case, fx, fy, fz, mx, my, mz = entry
        print(f'node_name={node_name}, load_case={load_case}, fx={fx}, fy={fy}, fz={fz}, mx={mx}, my={my}, mz={mz}')
        # Check if node exists in the model
        if node_name not in node_name_to_id:
            print(f"Warning: Node {node_name} not found in model, skipping load entry")
            continue
            
        node_id = node_name_to_id[node_name]
        
        # Add the loads to the nodal_loads structure
        load_vector = np.array([fx, fy, fz, mx, my, mz])
        nodal_loads[load_case][node_id] += load_vector
        


    
    print("="*60)

    # Prepare output data
    output = {}
    for load_case, node_loads in nodal_loads.items():
        output[load_case] = {}
        for node_id, loads in node_loads.items():
            node_name = next((node['name'] for node in nodes_data if node['id'] == node_id), None)
            if node_name:
                output[load_case][node_name] = loads.tolist()
            else:
                print(f"Warning: Could not find node name for ID {node_id}")

    print("\nTotal load vector (6 DOF) per load case:")
    for load_case, node_loads in nodal_loads.items():
        total_load_vector = np.zeros(6)
        for load in node_loads.values():
            total_load_vector += load
        print(f"Total load for load case '{load_case}' = {total_load_vector}")


    
    # Save the combined loads
    save_path = os.path.join(JSON_FOLDER, "load_data", "nodal_loads_from_slab.json")
    try:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        with open(save_path, 'w') as f:
            json.dump(output, f, indent=2)
        print(f"Successfully saved nodal loads to: {save_path}")
    except (IOError, OSError) as e:
        print(f"Warning: Failed to save nodal loads: {e}")

    return output



def create_nodal_load_combinations(nodal_loads, load_combinations, JSON_FOLDER):
    """
    Creates load combinations from individual load cases and saves them to JSON files.
    
    Args:
        nodal_loads (dict): Dictionary of nodal loads from calculate_nodal_loads_from_mesh
        load_combinations (dict): Dictionary of load combinations with factors
        JSON_FOLDER (str): Path to the JSON output folder
    """
    tolerance = 1e-6
    
    # Validate inputs
    if not nodal_loads:
        raise ValueError("nodal_loads cannot be empty")
    
    if not os.path.exists(JSON_FOLDER):
        raise ValueError(f"JSON_FOLDER does not exist: {JSON_FOLDER}")
    
    # Process each combination
    for combo_name, components in load_combinations.items():
        combined_loads = defaultdict(lambda: np.zeros(6))
        
        # Combine load cases with factors
        for load_case, factor in components:
            if load_case not in nodal_loads:
                print(f"Warning: Load case '{load_case}' not found in nodal_loads, skipping for combination '{combo_name}'")
                continue
                
            if abs(factor) < tolerance:
                continue  # Skip zero factors
                
            for node_id, loads in nodal_loads[load_case].items():
                # Convert loads to numpy array if it's a list
                loads_array = np.array(loads) if isinstance(loads, (list, tuple)) else loads
                combined_loads[node_id] += loads_array * factor
        
        # Prepare output data
        output = {}
        for node_id, loads in combined_loads.items():
            # Convert numpy array back to list for JSON serialization
            output[node_id] = loads.tolist()
        
        # Save the combination
        save_dir = os.path.join(JSON_FOLDER, "load_data")
        os.makedirs(save_dir, exist_ok=True)
        
        # Create filename based on combination name
        filename = f"nodal_loads_{combo_name.lower()}.json"
        save_path = os.path.join(save_dir, filename)
        
        try:
            with open(save_path, 'w') as f:
                json.dump(output, f, indent=2)
            print(f"Successfully saved combination '{combo_name}' to: {save_path}")
        except (IOError, OSError) as e:
            print(f"Warning: Failed to save combination '{combo_name}': {e}")
    
    return True





def process_element_loads_with_node_coordinates(all_element_loads, nodes_data, members_data, JSON_FOLDER):
    """
    Combines member-node coordinate mapping and element load processing into a single function.
    Now includes calculation of point load coordinates based on their location along the member.
    
    Args:
        all_element_loads: List of element load dictionaries
        nodes_data: List of node data dictionaries
        members_data: List of member data dictionaries
        
    Returns:
        Dictionary containing combined element loads with node coordinates and point load coordinates
    """
    wall_dir = os.path.join(JSON_FOLDER, 'load_data')
    os.makedirs(wall_dir, exist_ok=True)

    # 1. First create the node coordinate mapping
    node_coords = {node['id']: [node['x'], node['y'], node['z']] for node in nodes_data}
    
    # 2. Create member name to node coordinates mapping and name to ID mapping
    member_nodes = {}
    name_to_id = {}
    for member in members_data:
        member_nodes[member['name']] = {
            'id': member['id'],
            'start_node': node_coords.get(member['start_node_id']),
            'end_node': node_coords.get(member['end_node_id'])
        }
        name_to_id[member['name']] = member['id']
    
    # 3. Process element loads with the coordinate data
    final_json = {"element_loads": []}
    existing_entries = {}

    for element_load in all_element_loads:
        for key, loads in element_load.items():
            member_names = eval(key)  # Convert string representation of list to actual list
            for load in loads:
                loadcase = tuple(load["LoadCase"])
                uniform = load["uniform"]
                point = load["point"]
                temp_points = [[tp["temp"], tp["y"]] for tp in load["temperature_points"]]

                for member_name in member_names:
                    # Get coordinates and ID for this member
                    member_data = member_nodes.get(member_name, {
                        'id': -1,
                        'start_node': [0, 0, 0],
                        'end_node': [0, 0, 0]
                    })
                    member_id = member_data['id']
                    coords = {
                        'start_node': member_data['start_node'],
                        'end_node': member_data['end_node']
                    }
                    
                    # Calculate point load coordinates
                    location = point["location"]  # Fractional location along member (0-1)
                    start_coords = coords['start_node']
                    end_coords = coords['end_node']
                    
                    # Linear interpolation between start and end coordinates
                    point_coords = [
                        start_coords[0] + location * (end_coords[0] - start_coords[0]),
                        start_coords[1] + location * (end_coords[1] - start_coords[1]),
                        start_coords[2] + location * (end_coords[2] - start_coords[2])
                    ]

                    # Calculate the center point
                    center = [(start_coords[0] + end_coords[0]) / 2,
                            (start_coords[1] + end_coords[1]) / 2,
                            (start_coords[2] + end_coords[2]) / 2]
                    
                    # Calculate the length of the element
                    member_length = math.sqrt(
                        (end_coords[0] - start_coords[0])**2 +
                        (end_coords[1] - start_coords[1])**2 +
                        (end_coords[2] - start_coords[2])**2
                    )

                    entry_key = (member_id, loadcase)
                    
                    if entry_key in existing_entries:
                        entry = existing_entries[entry_key]
                        entry["uniform"][0] += uniform["x"]
                        entry["uniform"][1] += uniform["y"]
                        entry["uniform"][2] += uniform["z"]
                        entry["point"][0] += point["x"]
                        entry["point"][1] += point["y"]
                        entry["point"][2] += point["z"]
                        # Update point coordinates by weighted average
                        if "point_coords" in entry:
                            entry["point_coords"] = point_coords

                        if "center" in entry:
                            entry["center"] = center

                        if "member_length" in entry:
                            entry["member_length"] = member_length
                    else:
                        new_entry = {
                            "member": member_id,
                            "member_name": member_name,  # Add member name to output
                            "load_case": list(loadcase),
                            "uniform": [uniform["x"], uniform["y"], uniform["z"]],
                            "point": [point["x"], point["y"], point["z"], point["location"]],
                            "temperature_points": temp_points,
                            "start_node_coords": coords['start_node'],
                            "end_node_coords": coords['end_node'],
                            "point_coords": point_coords,  # Add calculated point coordinates
                            "center": center,  # Add calculated center
                            "member_length": member_length,
                        }
                        existing_entries[entry_key] = new_entry
                        final_json["element_loads"].append(new_entry)

    # Create filename and save
    members_load_combination_filename = f"member_loads.json"
    output_path = os.path.join(wall_dir, members_load_combination_filename)
    
    with open(output_path, 'w') as f:
        json.dump(final_json, f, indent=2)
    
    print(f"Saved load combination to {output_path}")
        # Sum total loads for each loadcase
    loadcase_totals = {}

    for entry in final_json["element_loads"]:
        lc = tuple(entry["load_case"])
        uniform = entry["uniform"]
        point = entry["point"]

        if lc not in loadcase_totals:
            loadcase_totals[lc] = {
                "uniform_total": [0.0, 0.0, 0.0],
                "point_total": [0.0, 0.0, 0.0]
            }

        loadcase_totals[lc]["uniform_total"][0] += uniform[0]
        loadcase_totals[lc]["uniform_total"][1] += uniform[1]
        loadcase_totals[lc]["uniform_total"][2] += uniform[2]

        loadcase_totals[lc]["point_total"][0] += point[0]
        loadcase_totals[lc]["point_total"][1] += point[1]
        loadcase_totals[lc]["point_total"][2] += point[2]

    for lc, totals in loadcase_totals.items():
        print(f"sum total load for {lc} for uniform load = {totals['uniform_total']}")
        print(f"sum total load for {lc} for point load = {totals['point_total']}")

    return final_json



def apply_member_load_combinations(final_json_with_coords, load_combinations, JSON_FOLDER):
    """Apply load combinations to member loads and save combined results to files"""
    wall_dir = os.path.join(JSON_FOLDER, 'load_data')
    os.makedirs(wall_dir, exist_ok=True)

    for combo_name, factors in load_combinations.items():
        combined_results = {}
        combo_loads = defaultdict(list)

        total_uniform = [0.0, 0.0, 0.0]
        total_point = [0.0, 0.0, 0.0]

        for member_load in final_json_with_coords["element_loads"]:
            member_id = member_load["member"]
            combined_load = {
                "member": member_id,
                "load_case": combo_name,
                "uniform": [0.0, 0.0, 0.0],
                "point": [0.0, 0.0, 0.0, 0.0],  # x,y,z,location
                "temperature_points": [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
                "start_node_coords": member_load["start_node_coords"],
                "end_node_coords": member_load["end_node_coords"],
                "point_coords": member_load.get("point_coords", [0, 0, 0]),
                "center": member_load.get("center", [0, 0, 0]),
                "member_length": member_load.get("member_length", [0, 0, 0]),

            }

            for load_case, factor in factors:
                if load_case in member_load["load_case"]:
                    combined_load["uniform"][0] += member_load["uniform"][0] * factor
                    combined_load["uniform"][1] += member_load["uniform"][1] * factor
                    combined_load["uniform"][2] += member_load["uniform"][2] * factor

                    combined_load["point"][0] += member_load["point"][0] * factor
                    combined_load["point"][1] += member_load["point"][1] * factor
                    combined_load["point"][2] += member_load["point"][2] * factor
                    combined_load["point"][3] = member_load["point"][3]

                    for i, (temp, y) in enumerate(member_load["temperature_points"]):
                        combined_load["temperature_points"][i][0] += temp * factor
                        combined_load["temperature_points"][i][1] = y

            total_uniform[0] += combined_load["uniform"][0]
            total_uniform[1] += combined_load["uniform"][1]
            total_uniform[2] += combined_load["uniform"][2]

            total_point[0] += combined_load["point"][0]
            total_point[1] += combined_load["point"][1]
            total_point[2] += combined_load["point"][2]

            combo_loads[member_id].append(combined_load)

        output_data = {"element_loads": []}
        for member_id, loads in combo_loads.items():
            output_data["element_loads"].extend(loads)

        members_load_combination_filename = f"member_load_{combo_name}.json"
        output_path = os.path.join(wall_dir, members_load_combination_filename)

        with open(output_path, 'w') as f:
            json.dump(output_data, f, indent=2)

        print(f"Saved load combination '{combo_name}' to {output_path}")
        print(f"sum total load for {combo_name} for uniform load = {total_uniform}")
        print(f"sum total load for {combo_name} for point load = {total_point}")

    return True



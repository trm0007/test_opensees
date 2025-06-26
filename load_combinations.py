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

    # print("Number of mesh elements:", len(mesh_data['elements']))
    # first_3_keys = list(mesh_data['elements'].keys())[:3]
    # print("First 3 mesh elements:", {k: mesh_data['elements'][k] for k in first_3_keys})

    # print("Mesh JSON files:", mesh_files)

    # # Print data for R10001, R10002, R10003
    # for rid in ['R10001', 'R10002', 'R10003']:
    #     if rid in mesh_data['elements']:
    #         print(f"Data for {rid}:", mesh_data['elements'][rid])

    return nodes_data, members_data, mesh_data



def parse_member_ids(member_names, JSON_FOLDER):
    """Parse member IDs from various input formats into a list of integers."""
    nodes_data, members_data, mesh_data = create_combined_structure_json(JSON_FOLDER)
    print("Member Names Input:", member_names)
    print("All Members in Data:")
    
    matched_ids = []
    for member in members_data:
        # print(member)
        if member['name'] in member_names:
            # print("Matched Member:", member)
            matched_ids.append(member['id'])

    if not matched_ids:
        print("No matching members found.")
        return None
    return matched_ids
  

def surface_element_calculator(coords, pressure):
    """
    Convert a uniformly distributed load (UDL) to equivalent nodal loads
    for a surface element with 6 DOF per node (3 translations + 3 rotations).
    
    Args:
        coords: List of nodal coordinates (3D points) that define the surface element
        pressure: Pressure magnitude (positive is outward normal direction)
    
    Returns:
        Dictionary containing:
        - 'nodal_loads': List of load vectors (6 DOF) for each node [Fx, Fy, Fz, Mx, My, Mz]
        - 'area': Total area of the surface element
        - 'centroid': [x, y, z] coordinates of the element centroid
    """
    coords = np.array(coords)
    num_nodes = len(coords)
    nodal_loads = [np.zeros(6) for _ in range(num_nodes)]  # 6 DOF per node
    
    # Calculate element centroid
    centroid = np.mean(coords, axis=0)
    area = 0.0
    
    # Gaussian integration setup
    if num_nodes == 3:  # Triangular element
        # Use single point integration for triangular elements
        gauss_points = [np.array([1/3, 1/3])]  # Natural coordinates (xi, eta)
        weights = [0.5]  # Weight for triangular element
    elif num_nodes == 4:  # Quadrilateral element
        g = 1/np.sqrt(3)
        gauss_points = [np.array([-g, -g]), np.array([g, -g]), 
                        np.array([g, g]), np.array([-g, g])]
        weights = [1.0, 1.0, 1.0, 1.0]
    else:
        raise ValueError(f"Unsupported number of nodes: {num_nodes}")
    
    for gp, w in zip(gauss_points, weights):
        xi, eta = gp[0], gp[1]
        
        # Calculate shape functions and derivatives
        if num_nodes == 3:  # Triangular element
            N = np.array([1 - xi - eta, xi, eta])
            # Derivatives with respect to natural coordinates
            dN_dxi = np.array([-1, 1, 0])
            dN_deta = np.array([-1, 0, 1])
        else:  # Quadrilateral element
            N = 0.25 * np.array([(1-xi)*(1-eta), (1+xi)*(1-eta), 
                                (1+xi)*(1+eta), (1-xi)*(1+eta)])
            # Derivatives with respect to natural coordinates
            dN_dxi = 0.25 * np.array([-(1-eta), (1-eta), (1+eta), -(1+eta)])
            dN_deta = 0.25 * np.array([-(1-xi), -(1+xi), (1+xi), (1-xi)])
        
        # Calculate Jacobian matrix
        dx_dxi = np.dot(dN_dxi, coords[:, 0])
        dy_dxi = np.dot(dN_dxi, coords[:, 1])
        dz_dxi = np.dot(dN_dxi, coords[:, 2])
        
        dx_deta = np.dot(dN_deta, coords[:, 0])
        dy_deta = np.dot(dN_deta, coords[:, 1])
        dz_deta = np.dot(dN_deta, coords[:, 2])
        
        # Tangent vectors
        t1 = np.array([dx_dxi, dy_dxi, dz_dxi])
        t2 = np.array([dx_deta, dy_deta, dz_deta])
        
        # Normal vector (cross product of tangent vectors)
        normal = np.cross(t1, t2)
        detJ = np.linalg.norm(normal)
        
        if detJ < 1e-10:
            print(f"Warning: Very small Jacobian determinant: {detJ}")
            continue
            
        # Normalize normal vector
        normal = normal / detJ
        
        # Accumulate area
        area += w * detJ
        
        # Pressure contribution - preserve sign convention
        p_vec = pressure * normal * w * detJ
        
        # Position of current integration point
        point_pos = np.dot(N, coords)
        
        # Add to nodal forces and moments
        for i in range(num_nodes):
            # Translational forces (first 3 DOFs)
            nodal_loads[i][:3] += N[i] * p_vec  # Fx, Fy, Fz
            moment_arm = [0.0,0.0,0.0]
            nodal_loads[i][3:] += N[i] * np.cross(moment_arm, p_vec)  # Mx, My, Mz
    
            nodal_loads[i] = np.where(np.abs(nodal_loads[i]) < 1e-6, 0.0, nodal_loads[i])
    
    return {
        'nodal_loads': nodal_loads,
        'area': area,
        'centroid': centroid.tolist()
    }


def calculate_nodal_loads_from_mesh(nodes_data, mesh_data, nodal_load_entries, JSON_FOLDER):
    """
    Calculate nodal loads with proper sign handling:
    1. Self-weight calculated as negative Z (downward)
    2. Pressure loads from calculator used as-is
    3. Apply final sign correction (1,1,-1,1,1,1) to all automatic loads 
    4. Add manual nodal loads with their original signs
    """
    tolerance = 1e-6
    
    # Final sign correction matrix (fx,fy,fz,mx,my,mz)
    SIGN_CORRECTION = np.array([1, 1, 1, 1, 1, 1])
    
    # Validate inputs
    if not nodes_data or not mesh_data:
        raise ValueError("nodes_data and mesh_data cannot be empty")
    if not os.path.exists(JSON_FOLDER):
        raise ValueError(f"JSON_FOLDER does not exist: {JSON_FOLDER}")

    # Create node mappings
    node_name_to_id = {node['name']: node['id'] for node in nodes_data}
    node_coords_map = {node['name']: np.array([node['x'], node['y'], node['z']]) for node in nodes_data}
    
    # Initialize nodal loads
    nodal_loads = defaultdict(lambda: defaultdict(lambda: np.zeros(6)))

    # Process mesh elements (pressure loads)
    if 'elements' not in mesh_data:
        raise KeyError("mesh_data must contain 'elements' key")

    for elem_id, elem in mesh_data['elements'].items():
        if not all(k in elem for k in ['nodes', 'load_case_names', 'pressures']):
            continue
        if len(elem['load_case_names']) != len(elem['pressures']):
            continue
            
        coords = [node_coords_map[n] for n in elem['nodes'] if n in node_coords_map]
        if len(coords) != len(elem['nodes']):
            continue
            
        for lc, pressure in zip(elem['load_case_names'], elem['pressures']):
            result = surface_element_calculator(coords, abs(pressure))
            for node_name, load in zip(elem['nodes'], result['nodal_loads']):
                if node_name in node_name_to_id:
                    nodal_loads[lc][node_name_to_id[node_name]] += np.array(load)

    # Process member self-weights (always negative Z)
    element_data_path = os.path.join(JSON_FOLDER, "element_data.json")
    if os.path.exists(element_data_path):
        with open(element_data_path) as f:
            element_data = json.load(f)
            
        for elem in element_data.get("elements", []):
            # if elem.get("area") is None:
            #     print("Element with None area:", elem)
            # if elem.get("unit_weight") is None:
            #     print("Element with None unit_weight:", elem)
            # if elem.get("length") is None:
            #     print("Element with None length:", elem)

            
            if not all(f in elem for f in ["unit_weight", "area", "length", "node_i_id", "node_j_id"]):
                continue
                
            weight = abs(float(elem["unit_weight"]) * 
                     float(elem["area"]) * 
                     float(elem["length"]))
            
            if weight < tolerance:
                continue
                
            half_weight = weight / 2
            nodal_loads["self_weight"][elem["node_i_id"]][2] -= half_weight
            nodal_loads["self_weight"][elem["node_j_id"]][2] -= half_weight

    # Apply final sign correction to automatic loads
    for lc in nodal_loads:
        for node_id in nodal_loads[lc]:
            nodal_loads[lc][node_id] *= SIGN_CORRECTION

    # Add manual nodal loads (no sign correction)
    for node_name, lc, fx, fy, fz, mx, my, mz in nodal_load_entries:
        if node_name in node_name_to_id:
            nodal_loads[lc][node_name_to_id[node_name]] += np.array([fx, fy, fz, mx, my, mz])

    # Prepare and save output
    output = {
        lc: {
            next((n['name'] for n in nodes_data if n['id'] == id), str(id)): load.tolist()
            for id, load in node_loads.items()
        }
        for lc, node_loads in nodal_loads.items()
    }
    
    save_path = os.path.join(JSON_FOLDER, "load_data", "nodal_loads_from_slab.json")
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    with open(save_path, 'w') as f:
        json.dump(output, f, indent=2)

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
                
            if (factor) < tolerance:
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


# def process_element_loads_with_node_coordinates(all_element_loads, nodes_data, members_data, JSON_FOLDER):
#     """
#     Combines member-node coordinate mapping and element load processing into a single function.
#     Now includes calculation of point load coordinates based on their location along the member.
    
#     Args:
#         all_element_loads: List of element load dictionaries
#         nodes_data: List of node data dictionaries
#         members_data: List of member data dictionaries
        
#     Returns:
#         Dictionary containing combined element loads with node coordinates and point load coordinates
#     """
#     wall_dir = os.path.join(JSON_FOLDER, 'load_data')
#     os.makedirs(wall_dir, exist_ok=True)

#     # 1. First create the node coordinate mapping
#     node_coords = {node['id']: [node['x'], node['y'], node['z']] for node in nodes_data}
    
#     # 2. Create member name to node coordinates mapping and name to ID mapping
#     member_nodes = {}
#     name_to_id = {}
#     for member in members_data:
#         member_nodes[member['name']] = {
#             'id': member['id'],
#             'start_node': node_coords.get(member['start_node_id']),
#             'end_node': node_coords.get(member['end_node_id'])
#         }
#         name_to_id[member['name']] = member['id']
    
#     # 3. Process element loads with the coordinate data
#     final_json = {"element_loads": []}
#     existing_entries = {}

#     for element_load in all_element_loads:
#         for key, loads in element_load.items():
#             member_names = eval(key)  # Convert string representation of list to actual list
#             for load in loads:
#                 loadcase = tuple(load["LoadCase"])
#                 uniform = load["uniform"]
#                 point = load["point"]
#                 temp_points = [[tp["temp"], tp["y"]] for tp in load["temperature_points"]]

#                 for member_name in member_names:
#                     # Get coordinates and ID for this member
#                     member_data = member_nodes.get(member_name, {
#                         'id': -1,
#                         'start_node': [0, 0, 0],
#                         'end_node': [0, 0, 0]
#                     })
#                     member_id = member_data['id']
#                     coords = {
#                         'start_node': member_data['start_node'],
#                         'end_node': member_data['end_node']
#                     }
                    
#                     # Calculate point load coordinates
#                     location = point["location"]  # Fractional location along member (0-1)
#                     start_coords = coords['start_node']
#                     end_coords = coords['end_node']
                    
#                     # Linear interpolation between start and end coordinates
#                     point_coords = [
#                         start_coords[0] + location * (end_coords[0] - start_coords[0]),
#                         start_coords[1] + location * (end_coords[1] - start_coords[1]),
#                         start_coords[2] + location * (end_coords[2] - start_coords[2])
#                     ]

#                     # Calculate the center point
#                     center = [(start_coords[0] + end_coords[0]) / 2,
#                             (start_coords[1] + end_coords[1]) / 2,
#                             (start_coords[2] + end_coords[2]) / 2]
                    
#                     # Calculate the length of the element
#                     member_length = math.sqrt(
#                         (end_coords[0] - start_coords[0])**2 +
#                         (end_coords[1] - start_coords[1])**2 +
#                         (end_coords[2] - start_coords[2])**2
#                     )

#                     entry_key = (member_id, loadcase)
                    
#                     if entry_key in existing_entries:
#                         entry = existing_entries[entry_key]
#                         entry["uniform"][0] += uniform["x"]
#                         entry["uniform"][1] += uniform["y"]
#                         entry["uniform"][2] += uniform["z"]
#                         entry["point"][0] += point["x"]
#                         entry["point"][1] += point["y"]
#                         entry["point"][2] += point["z"]
#                         # Update point coordinates by weighted average
#                         if "point_coords" in entry:
#                             entry["point_coords"] = point_coords

#                         if "center" in entry:
#                             entry["center"] = center

#                         if "member_length" in entry:
#                             entry["member_length"] = member_length
#                     else:
#                         new_entry = {
#                             "member": member_id,
#                             "member_name": member_name,  # Add member name to output
#                             "load_case": list(loadcase),
#                             "uniform": [uniform["x"], uniform["y"], uniform["z"]],
#                             "point": [point["x"], point["y"], point["z"], point["location"]],
#                             "temperature_points": temp_points,
#                             "start_node_coords": coords['start_node'],
#                             "end_node_coords": coords['end_node'],
#                             "point_coords": point_coords,  # Add calculated point coordinates
#                             "center": center,  # Add calculated center
#                             "member_length": member_length,
#                         }
#                         existing_entries[entry_key] = new_entry
#                         final_json["element_loads"].append(new_entry)

#     # Create filename and save
#     members_load_combination_filename = f"member_loads.json"
#     output_path = os.path.join(wall_dir, members_load_combination_filename)
    
#     with open(output_path, 'w') as f:
#         json.dump(final_json, f, indent=2)
    
#     print(f"Saved load combination to {output_path}")


#     return final_json


def process_element_loads_with_node_coordinates(loading_mapping, nodes_data, members_data, JSON_FOLDER):
    print("\n=== DEBUG START ===")
    print("1. Checking loading_mapping structure:")
    for i, (members, loads) in enumerate(loading_mapping):
        print(f"  Group {i}:")
        print(f"  Members type: {type(members)}, first member: {members[0] if members else 'empty'}")
        print(f"  Loads count: {len(loads)}")

    print("\n2. Inspecting members_data structure:")
    if members_data:
        first_member = members_data[0]
        print(f"  First member keys: {list(first_member.keys())}")
        print(f"  First member content: {first_member}")
    
    print("\n3. Building member_nodes dictionary:")
    member_nodes = {}
    for member in members_data:
        # print(f"  Processing member: {member.get('name', 'NO_NAME')} (type: {type(member.get('name', 'NO_NAME'))})")
        member_name = str(member.get('name', '')) if member.get('name') is not None else ''
        
        # Try to determine the correct coordinate keys
        coord_keys = []
        for key in member.keys():
            if any(coord in key.lower() for coord in ['x', 'y', 'z', 'coord']):
                coord_keys.append(key)
        
        # print(f"    Available coordinate-related keys: {coord_keys}")
        
        # Attempt to extract coordinates based on common patterns
        start_coords = [0, 0, 0]  # Default values
        end_coords = [0, 0, 0]    # Default values
        
        # Try different possible key patterns
        coord_patterns = [
            # Pattern 1: x1, y1, z1, x2, y2, z2
            (['x1', 'y1', 'z1'], ['x2', 'y2', 'z2']),
            # Pattern 2: start_x, start_y, start_z, end_x, end_y, end_z
            (['start_x', 'start_y', 'start_z'], ['end_x', 'end_y', 'end_z']),
            # Pattern 3: xi, yi, zi, xj, yj, zj
            (['xi', 'yi', 'zi'], ['xj', 'yj', 'zj']),
            # Pattern 4: x_start, y_start, z_start, x_end, y_end, z_end
            (['x_start', 'y_start', 'z_start'], ['x_end', 'y_end', 'z_end']),
            # Pattern 5: node_i_x, node_i_y, node_i_z, node_j_x, node_j_y, node_j_z
            (['node_i_x', 'node_i_y', 'node_i_z'], ['node_j_x', 'node_j_y', 'node_j_z']),
        ]
        
        found_pattern = False
        for start_keys, end_keys in coord_patterns:
            if all(key in member for key in start_keys + end_keys):
                start_coords = [member[key] for key in start_keys]
                end_coords = [member[key] for key in end_keys]
                print(f"    Using pattern: {start_keys} -> {end_keys}")
                found_pattern = True
                break
        
        if not found_pattern:
            # print(f"    WARNING: No coordinate pattern found for member {member_name}")
            # Try to get coordinates from node references if available
            if 'start_node_id' in member and 'end_node_id' in member:
                start_node_id = member['start_node_id']
                end_node_id = member['end_node_id']
                
                # Look up coordinates in nodes_data
                for node in nodes_data:
                    if node.get('id') == start_node_id or node.get('name') == start_node_id:
                        start_coords = [node.get('x', 0), node.get('y', 0), node.get('z', 0)]
                    if node.get('id') == end_node_id or node.get('name') == end_node_id:
                        end_coords = [node.get('x', 0), node.get('y', 0), node.get('z', 0)]
                # print(f"    Using node references: start_node={start_node_id}, end_node={end_node_id}")
        
        member_nodes[member_name] = {
            'id': member.get('id', -1),
            'start_node': start_coords,
            'end_node': end_coords
        }
        # print(f"    Stored: {member_nodes[member_name]}")

    print("\n4. Processing element loads:")
    final_json = {"element_loads": []}
    
    for member_list, loads in loading_mapping:
        # print(f"\nProcessing member list: {type(member_list)} with {len(member_list) if hasattr(member_list, '__len__') else 'unknown'} members")
        
        # Handle different types of member_list
        if isinstance(member_list, (list, tuple)):
            members_to_process = member_list
        else:
            # If it's a single member, wrap it in a list
            members_to_process = [member_list]
        
        for member_name in members_to_process:
            # print(f"  Current member: {member_name} (type: {type(member_name)})")
            
            # Convert to string if needed
            lookup_key = str(member_name) if not isinstance(member_name, str) else member_name
            # print(f"  Lookup key: {lookup_key}")
            
            member_data = member_nodes.get(lookup_key, {
                'id': -1,
                'start_node': [0, 0, 0],
                'end_node': [0, 0, 0]
            })
            # print(f"  Found member data: ID={member_data['id']}, Start={member_data['start_node']}, End={member_data['end_node']}")

            for load in loads:
                # print(f"    Processing load case: {load.get('LoadCase', 'UNKNOWN')}")
                
                # Create the load entry
                load_entry = {
                    "member_name": lookup_key,
                    "member_id": member_data['id'],
                    "start_node_coords": member_data['start_node'],
                    "end_node_coords": member_data['end_node'],
                    "load_case": load.get('LoadCase', 'UNKNOWN'),
                    "load_data": load
                }
                
                final_json["element_loads"].append(load_entry)

    print(f"\n5. Final result: {len(final_json['element_loads'])} load entries created")
    print("=== DEBUG END ===\n")
    
    # Optionally save to file
    try:
        import json
        import os
        output_file = os.path.join(JSON_FOLDER, "element_loads_with_coordinates.json")
        with open(output_file, 'w') as f:
            json.dump(final_json, f, indent=2)
        print(f"Results saved to: {output_file}")
    except Exception as e:
        print(f"Could not save to file: {e}")
    
    return final_json

# def apply_member_load_combinations(final_json_with_coords, load_combinations, JSON_FOLDER):
#     """Apply load combinations to member loads and save combined results to files"""
#     wall_dir = os.path.join(JSON_FOLDER, 'load_data')
#     os.makedirs(wall_dir, exist_ok=True)

#     for combo_name, factors in load_combinations.items():
#         combined_results = {}
#         combo_loads = defaultdict(list)

#         total_uniform = [0.0, 0.0, 0.0]
#         total_point = [0.0, 0.0, 0.0]

#         for member_load in final_json_with_coords["element_loads"]:
#             member_id = member_load["member"]
#             member_name = member_load["member_name"]
#             combined_load = {
#                 "member_name": member_name,
#                 "member": member_id,
#                 "load_case": combo_name,
#                 "uniform": [0.0, 0.0, 0.0],
#                 "point": [0.0, 0.0, 0.0, 0.0],  # x,y,z,location
#                 "temperature_points": [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
#                 "start_node_coords": member_load["start_node_coords"],
#                 "end_node_coords": member_load["end_node_coords"],
#                 "point_coords": member_load.get("point_coords", [0, 0, 0]),
#                 "center": member_load.get("center", [0, 0, 0]),
#                 "member_length": member_load.get("member_length", [0, 0, 0]),

#             }

#             for load_case, factor in factors:
#                 if load_case in member_load["load_case"]:
#                     combined_load["uniform"][0] += member_load["uniform"][0] * factor
#                     combined_load["uniform"][1] += member_load["uniform"][1] * factor
#                     combined_load["uniform"][2] += member_load["uniform"][2] * factor

#                     combined_load["point"][0] += member_load["point"][0] * factor
#                     combined_load["point"][1] += member_load["point"][1] * factor
#                     combined_load["point"][2] += member_load["point"][2] * factor
#                     combined_load["point"][3] = member_load["point"][3]

#                     for i, (temp, y) in enumerate(member_load["temperature_points"]):
#                         combined_load["temperature_points"][i][0] += temp * factor
#                         combined_load["temperature_points"][i][1] = y

#             total_uniform[0] += combined_load["uniform"][0]
#             total_uniform[1] += combined_load["uniform"][1]
#             total_uniform[2] += combined_load["uniform"][2]

#             total_point[0] += combined_load["point"][0]
#             total_point[1] += combined_load["point"][1]
#             total_point[2] += combined_load["point"][2]

#             combo_loads[member_id].append(combined_load)

#         output_data = {"element_loads": []}
#         for member_id, loads in combo_loads.items():
#             output_data["element_loads"].extend(loads)

#         members_load_combination_filename = f"member_load_{combo_name}.json"
#         output_path = os.path.join(wall_dir, members_load_combination_filename)

#         with open(output_path, 'w') as f:
#             json.dump(output_data, f, indent=2)

#         print(f"Saved load combination '{combo_name}' to {output_path}")
#         print(f"sum total load for {combo_name} for uniform load = {total_uniform}")
#         print(f"sum total load for {combo_name} for point load = {total_point}")

#     return True


def apply_member_load_combinations(final_json_with_coords, load_combinations, JSON_FOLDER):
    """Apply load combinations to member loads and save combined results to files"""
    
    # First, let's debug the structure of final_json_with_coords
    print("\n=== DEBUGGING apply_member_load_combinations ===")
    print(f"final_json_with_coords keys: {list(final_json_with_coords.keys())}")
    
    if "element_loads" in final_json_with_coords and final_json_with_coords["element_loads"]:
        first_load = final_json_with_coords["element_loads"][0]
        print(f"First element_load keys: {list(first_load.keys())}")
        print(f"First element_load sample: {first_load}")
    
    wall_dir = os.path.join(JSON_FOLDER, 'load_data')
    os.makedirs(wall_dir, exist_ok=True)
    
    for combo_name, factors in load_combinations.items():
        print(f"\nProcessing load combination: {combo_name}")
        combined_results = {}
        combo_loads = defaultdict(list)
        
        total_uniform = [0.0, 0.0, 0.0]
        total_point = [0.0, 0.0, 0.0]
        
        for member_load in final_json_with_coords["element_loads"]:
            # Fix: Use the correct key names based on your data structure
            member_id = member_load.get("member_id", member_load.get("member", -1))
            member_name = member_load.get("member_name", f"member_{member_id}")
            
            # Extract load data from the nested structure
            load_data = member_load.get("load_data", {})
            load_case = member_load.get("load_case", "")
            
            # Initialize combined load structure
            combined_load = {
                "member_name": member_name,
                "member": member_id,
                "load_case": combo_name,
                "uniform": [0.0, 0.0, 0.0],
                "point": [0.0, 0.0, 0.0, 0.0],  # x,y,z,location
                "temperature_points": [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
                "start_node_coords": member_load.get("start_node_coords", [0, 0, 0]),
                "end_node_coords": member_load.get("end_node_coords", [0, 0, 0]),
                "point_coords": member_load.get("point_coords", [0, 0, 0]),
                "center": member_load.get("center", [0, 0, 0]),
                "member_length": member_load.get("member_length", [0, 0, 0]),
            }
            
            # Apply load combination factors
            for load_case_name, factor in factors:
                if load_case_name == load_case or load_case_name in str(load_case):
                    # Extract load values from load_data structure
                    # Adjust these based on your actual load_data structure
                    uniform_loads = [0.0, 0.0, 0.0]
                    point_loads = [0.0, 0.0, 0.0, 0.0]
                    temp_points = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
                    
                    # Try to extract loads from different possible structures
                    if isinstance(load_data, dict):
                        # Case 1: Direct load values
                        if "uniform" in load_data:
                            if isinstance(load_data["uniform"], (list, tuple)):
                                uniform_loads = list(load_data["uniform"])[:3] + [0.0] * (3 - len(load_data["uniform"]))
                        
                        if "point" in load_data:
                            if isinstance(load_data["point"], (list, tuple)):
                                point_loads = list(load_data["point"])[:4] + [0.0] * (4 - len(load_data["point"]))
                        
                        if "temperature_points" in load_data:
                            temp_points = load_data["temperature_points"]
                        
                        # Case 2: Load values in nested structure
                        for key, value in load_data.items():
                            try:
                                if key.lower() in ['fx', 'fy', 'fz', 'force_x', 'force_y', 'force_z']:
                                    if key.lower() in ['fx', 'force_x']:
                                        uniform_loads[0] = float(value)
                                    elif key.lower() in ['fy', 'force_y']:
                                        uniform_loads[1] = float(value)
                                    elif key.lower() in ['fz', 'force_z']:
                                        uniform_loads[2] = float(value)
                                
                                # Handle point loads
                                if key.lower() in ['px', 'py', 'pz', 'point_x', 'point_y', 'point_z']:
                                    if key.lower() in ['px', 'point_x']:
                                        point_loads[0] = float(value)
                                    elif key.lower() in ['py', 'point_y']:
                                        point_loads[1] = float(value)
                                    elif key.lower() in ['pz', 'point_z']:
                                        point_loads[2] = float(value)
                            except (ValueError, TypeError):
                                # Skip non-numeric values
                                print(f"Warning: Could not convert {key}={value} to float")
                                continue
                    
                    # Apply factors to loads with type checking
                    try:
                        combined_load["uniform"][0] += float(uniform_loads[0]) * factor
                        combined_load["uniform"][1] += float(uniform_loads[1]) * factor
                        combined_load["uniform"][2] += float(uniform_loads[2]) * factor
                        
                        combined_load["point"][0] += float(point_loads[0]) * factor
                        combined_load["point"][1] += float(point_loads[1]) * factor
                        combined_load["point"][2] += float(point_loads[2]) * factor
                        if len(point_loads) > 3:
                            combined_load["point"][3] = float(point_loads[3])
                    except (ValueError, TypeError) as e:
                        print(f"Warning: Error applying load factors: {e}")
                        print(f"  uniform_loads: {uniform_loads}")
                        print(f"  point_loads: {point_loads}")
                        print(f"  factor: {factor}")
                        continue
                    
                    # Apply factors to temperature points
                    for i, temp_point in enumerate(temp_points):
                        if i < len(combined_load["temperature_points"]):
                            # Handle different temperature point formats
                            if isinstance(temp_point, (list, tuple)) and len(temp_point) >= 2:
                                temp, y = temp_point[0], temp_point[1]
                            elif isinstance(temp_point, (int, float)):
                                temp, y = temp_point, 0.0
                            else:
                                temp, y = 0.0, 0.0
                            
                            # Ensure temp is numeric
                            try:
                                temp_val = float(temp) if temp is not None else 0.0
                                y_val = float(y) if y is not None else 0.0
                                combined_load["temperature_points"][i][0] += temp_val * factor
                                combined_load["temperature_points"][i][1] = y_val
                            except (ValueError, TypeError):
                                # If conversion fails, skip this temperature point
                                print(f"Warning: Could not convert temperature point {temp_point} to numeric values")
                                continue
            
            # Add to totals
            total_uniform[0] += combined_load["uniform"][0]
            total_uniform[1] += combined_load["uniform"][1]
            total_uniform[2] += combined_load["uniform"][2]
            
            total_point[0] += combined_load["point"][0]
            total_point[1] += combined_load["point"][1]
            total_point[2] += combined_load["point"][2]
            
            combo_loads[member_id].append(combined_load)
        
        # Create output data structure
        output_data = {"element_loads": []}
        for member_id, loads in combo_loads.items():
            output_data["element_loads"].extend(loads)
        
        # Save to file
        members_load_combination_filename = f"member_load_{combo_name}.json"
        output_path = os.path.join(wall_dir, members_load_combination_filename)
        
        with open(output_path, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        print(f"Saved load combination '{combo_name}' to {output_path}")
        print(f"Sum total load for {combo_name} for uniform load = {total_uniform}")
        print(f"Sum total load for {combo_name} for point load = {total_point}")
    
    print("=== END DEBUGGING ===")
    return True



import os
import json
import numpy as np
from typing import Dict, List, Union

def filter_nodes_at_z_levels(
    nodes_data: List[Dict],
    members_data: List[Dict],
    z_spacing: List[float],
    
) -> Dict:
    all_nodes = {}

    # Map full_nodes_data by id for lookup
    node_lookup = {node['id']: node for node in nodes_data}

    for member in members_data:
        for key in ['start_node_id', 'end_node_id']:
            node_id = member.get(key)
            node_name_key = 'start_node_name' if 'start' in key else 'end_node_name'
            if node_id not in all_nodes and node_id in node_lookup:
                node = node_lookup[node_id]
                all_nodes[node_id] = {
                    'id': node_id,
                    'name': member.get(node_name_key, ''),
                    'x': node.get('x'),
                    'y': node.get('y'),
                    'z': node.get('z')
                }

    nodes_data = list(all_nodes.values())
    # print("nodes_data")
    # print(nodes_data)

    cumulative_z_levels = np.cumsum([0] + z_spacing).tolist()[1:]

    results = {}
    for z_level in cumulative_z_levels:
        level_nodes = [node for node in nodes_data if node['z'] == z_level]
        results[z_level] = {
            'z_level': z_level,
            'nodes': level_nodes
        }

    return results



def calculate_center_of_mass(
    JSON_FOLDER: str,
    z_spacing: List[float],
    mass_from_positive_loads: bool = False
) -> Dict:
    """
    Enhanced center of mass calculation with additional features.
    Assigns unique names and IDs to center of mass points, checking against existing node data.
    
    Args:
        JSON_FOLDER: Path to folder containing structural data
        z_spacing: List of z-spacing values
        mass_from_positive_loads: If True, treats positive Z loads as mass
        
    Returns:
        Dictionary with center of mass results for each level, including COM node info
    """
    # Load structural data
    nodes_data, members_data, mesh_data = create_combined_structure_json(JSON_FOLDER)
    
    # Get filtered level data
    level_data = filter_nodes_at_z_levels(nodes_data, members_data, z_spacing)
    print(f'level_data={level_data}')
    # Save filtered data
    filtered_path = os.path.join(JSON_FOLDER, "filtered_nodes_members_at_z_levels.json")
    with open(filtered_path, 'w') as f:
        json.dump(level_data, f, indent=2)
    
    # Load load data with error handling
    try:
        with open(os.path.join(JSON_FOLDER, "load_data", "nodal_loads_from_slab.json")) as f:
            nodal_loads = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        nodal_loads = {}
    
    try:
        with open(os.path.join(JSON_FOLDER, "load_data", "member_loads.json")) as f:
            member_loads = json.load(f)
            if isinstance(member_loads, list):
                member_loads = {"element_loads": member_loads}
    except (FileNotFoundError, json.JSONDecodeError):
        member_loads = {"element_loads": []}
    
    # Create optimized node mappings
    node_id_to_coords = {n['id']: np.array([n['x'], n['y'], n['z']]) for n in nodes_data}
    node_name_to_id = {n['name']: n['id'] for n in nodes_data}
    
    # Find the next available ID for COM nodes (starting from 100000 to avoid conflicts)
    existing_ids = {n['id'] for n in nodes_data}
    com_id_start = 100000
    while com_id_start in existing_ids:
        com_id_start += 1
    
    results = {}
    tol = 0.01
    
    for z_level in level_data.keys():
        print(f"\nProcessing z_level: {z_level:.2f}")
        
        total_mass = 0.0
        weighted_sum = np.zeros(3)  # For x, y, z coordinates
        
        # Process nodal loads
        for load_case, nodes in nodal_loads.items():
            for node_name, load in nodes.items():
                if node_name not in node_name_to_id:
                    continue
                    
                node_coords = node_id_to_coords.get(node_name_to_id[node_name])
                if node_coords is None or abs(node_coords[2] - z_level) >= tol:
                    continue
                
                # Handle mass calculation direction
                mass = load[2] if mass_from_positive_loads else abs(load[2])
                if mass <= 0:
                    continue
                    
                total_mass += mass
                weighted_sum += mass * node_coords
        
        # Process member loads
        for member_load in member_loads.get("element_loads", []):
            start = np.array(member_load.get("start_node_coords", [0, 0, 0]))
            end = np.array(member_load.get("end_node_coords", [0, 0, 0]))
            
            # Skip members not at this level
            if all(abs(coord[2] - z_level) >= tol for coord in [start, end]):
                continue
                
            # Calculate member properties
            length = np.linalg.norm(end - start)
            if length < 1e-6:
                continue
                
            # Get uniform load mass
            uniform = member_load.get("uniform", [0, 0, 0])
            uniform_z = uniform.get("z", 0) if isinstance(uniform, dict) else uniform[2]
            mass_per_length = uniform_z if mass_from_positive_loads else abs(uniform_z)
            
            if mass_per_length > 0:
                total_mass += mass_per_length * length
                center = np.array(member_load.get("center", (start + end) / 2))
                weighted_sum += mass_per_length * length * center
            
            # Process point loads
            point_load = member_load.get("point", [0, 0, 0, 0])
            if isinstance(point_load, dict):
                point_mass = point_load.get("z", 0) if mass_from_positive_loads else abs(point_load.get("z", 0))
                point_coords = np.array(member_load.get("point_coords", [0, 0, 0]))
            else:
                point_mass = point_load[2] if mass_from_positive_loads else abs(point_load[2])
                point_coords = np.array(member_load.get("point_coords", [0, 0, 0]))
            
            if point_mass > 0 and abs(point_coords[2] - z_level) < tol:
                total_mass += point_mass
                weighted_sum += point_mass * point_coords
        
        # Calculate center of mass
        com = weighted_sum / total_mass if total_mass > 0 else np.array([0, 0, z_level])
        
        # Generate unique name and ID for COM point
        com_name = f"COM_{z_level:.2f}"
        com_id = com_id_start
        com_id_start += 1
        
        # Ensure name is unique
        while com_name in node_name_to_id:
            com_name = f"COM_{z_level:.2f}_{com_id}"
        
        print(f"Total mass: {total_mass:.2f} kg")
        print(f"Center of mass: ({com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f})")
        print(f"Assigned COM node: {com_name} (ID: {com_id})")
        
        # Create COM node data
        com_node = {
            'id': com_id,
            'name': com_name,
            'x': float(com[0]),
            'y': float(com[1]),
            'z': float(com[2]),
            'is_com_node': True,
            'mass': total_mass
        }
        
        results[z_level] = {
            'total_mass': total_mass,
            'center_of_mass': com.tolist(),
            'z_level': z_level,
            'com_node': com_node,
            'nodes_at_level': level_data[z_level]['nodes'],

        }
    
    # Save results
    output_path = os.path.join(JSON_FOLDER, "center_of_mass_results.json")
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    
    return results


def create_rigid_diaphragms(JSON_FOLDER: str, ndm: int = 3, ndf: int = 6) -> bool:
    """
    Creates rigid diaphragms at each level using center of mass nodes as master nodes.
    Properly initializes OpenSees model dimensions before creating nodes and constraints.

    Args:
        JSON_FOLDER: Path to folder containing structural data
        ndm: Number of model dimensions (default=3 for 3D models)
        ndf: Number of degrees of freedom (default=6 for 3D models)

    Returns:
        bool: True if all rigid diaphragms were successfully created, False otherwise
    """
    # Load center of mass results
    output_path = os.path.join(JSON_FOLDER, "center_of_mass_results.json")
    try:
        with open(output_path) as f:
            com_results = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error loading center of mass results: {e}")
        return False

    # Initialize model dimensions
    ops.model('BasicBuilder', '-ndm', ndm, '-ndf', ndf)
    print(f"Initialized model with ndm={ndm}, ndf={ndf}")

    success = True

    for z_level, level_data in com_results.items():
        try:
            # Extract COM node and level nodes
            com_node = level_data['com_node']
            nodes_at_level = level_data['nodes_at_level']

            # Create master node in OpenSees
            rNodeTag = int(com_node['id'])
            x, y, z = float(com_node['x']), float(com_node['y']), float(com_node['z'])
            ops.node(rNodeTag, x, y, z)
            print(f"Created master node {rNodeTag} at ({x:.2f}, {y:.2f}, {z:.2f})")

            # Prepare constrained nodes (must exist in model already)
            cNodeTags = [int(node['id']) for node in nodes_at_level]

            # Create rigid diaphragm (perpDirn=3 for XY plane constraint)
            ops.rigidDiaphragm(3, rNodeTag, *cNodeTags)

            print(f"Created rigid diaphragm at z-level {z_level}:")
            print(f"  Master node: {rNodeTag}")
            print(f"  Constrained nodes: {len(cNodeTags)} nodes\n")

        except Exception as e:
            print(f"\nERROR creating diaphragm at z-level {z_level}: {str(e)}")
            print("Possible issues:")
            print("- Constrained nodes don't exist in model")
            print("- Duplicate node definitions")
            print("- Invalid node coordinates")
            success = False

    return success


def create_zero_length_elements(nodes, K):
    import openseespy.opensees as ops

    # Define soil materials for 3 directions
    ops.uniaxialMaterial("ENT", 1, K)         # vertical (Z)
    ops.uniaxialMaterial("ENT", 2, K * 0.5)   # horizontal X
    ops.uniaxialMaterial("ENT", 3, K * 0.3)   # horizontal Y

    spring_tags = []
    spring_node_start = 10001
    element_tag_start = 20001

    for i, node in enumerate(nodes):
        node_id = node["id"]
        x, y, z = node["coord"]

        # Foundation node
        ops.node(node_id, x, y, z)

        # Spring node
        spring_node_id = node_id + spring_node_start
        ops.node(spring_node_id, x, y, z)
        ops.fix(spring_node_id, 1, 1, 1, 1, 1, 1)

        # zeroLength element with 3 directions: X (1), Y (2), Z (3)
        element_id = element_tag_start + i
        ops.element("zeroLength", element_id, node_id, spring_node_id,
                    "-mat", 2, 3, 1,
                    "-dir", 1, 2, 3)

        print(f"element zeroLength {element_id} {node_id} {spring_node_id} -mat 2 3 1 -dir 1 2 3")

        spring_tags.append({
            "spring_node": spring_node_id,
            "element_id": element_id
        })

    return spring_tags

K = 1e8

nodes = [
    {"id": 1, "name": "n1", "coord": [0.0, 0.0, 0.0]},
    {"id": 2, "name": "n2", "coord": [1.0, 0.0, 0.0]},
    {"id": 3, "name": "n3", "coord": [1.0, 1.0, 0.0]},
    {"id": 4, "name": "n4", "coord": [0.0, 1.0, 0.0]},
]

# spring_tags = create_zero_length_elements(nodes, K)


# def calculate_center_of_mass(JSON_FOLDER, z_spacing):
#     # Load all required data
#     nodes_data, members_data, mesh_data = create_combined_structure_json(JSON_FOLDER)
    
#     # Load nodal loads
#     nodal_loads_path = os.path.join(JSON_FOLDER, "load_data", "nodal_loads_from_slab.json")
#     with open(nodal_loads_path, 'r') as f:
#         nodal_loads = json.load(f)
    
#     # Load member loads - handle both list and dict formats
#     member_loads_path = os.path.join(JSON_FOLDER, "load_data", "member_loads.json")
#     with open(member_loads_path, 'r') as f:
#         member_loads = json.load(f)
    
#     # Convert member_loads to consistent format if it's a list
#     if isinstance(member_loads, list):
#         member_loads = {"element_loads": member_loads}
    
#     # Create node mappings
#     node_id_to_coords = {node['id']: np.array([node['x'], node['y'], node['z']]) for node in nodes_data}
#     node_name_to_id = {node['name']: node['id'] for node in nodes_data}
    
#     # Initialize results dictionary
#     com_results = {}
    
#     # Calculate cumulative z-levels
#     cumulative_z_levels = []
#     current_z = 0.0
#     for spacing in z_spacing:
#         current_z += spacing
#         cumulative_z_levels.append(current_z)
    
#     for z_level in cumulative_z_levels:
#         print(f"\nCalculating center of mass for cumulative z_level = {z_level}")
        
#         # Filter nodes at this z_level (with tolerance)
#         tol = 0.01  # Small tolerance for floating point comparison
#         level_nodes = [node for node in nodes_data if abs(node['z'] - z_level) < tol]
        
#         if not level_nodes:
#             print(f"No nodes found at z_level = {z_level}")
#             continue
        
#         # Initialize total mass and weighted coordinates
#         total_mass = 0.0
#         weighted_x = 0.0
#         weighted_y = 0.0
#         weighted_z = 0.0
        
#         # Process nodal loads (point masses)
#         for load_case, nodes in nodal_loads.items():
#             for node_name, load in nodes.items():
#                 node_id = node_name_to_id.get(node_name)
#                 if node_id is None:
#                     continue
                
#                 node_coords = node_id_to_coords.get(node_id)
#                 if node_coords is None:
#                     continue
                
#                 # Only consider nodes at this z_level
#                 if abs(node_coords[2] - z_level) >= tol:
#                     continue
                
#                 # Fz represents the mass (assuming gravity acts in -z direction)
#                 mass = abs(load[2])  # Take absolute value of z-force as mass
#                 total_mass += mass
#                 weighted_x += mass * node_coords[0]
#                 weighted_y += mass * node_coords[1]
#                 weighted_z += mass * node_coords[2]
        
#         # Process member loads (distributed masses)
#         if "element_loads" in member_loads:
#             for member_load in member_loads["element_loads"]:
#                 # Get member coordinates
#                 start_coords = np.array(member_load.get("start_node_coords", [0, 0, 0]))
#                 end_coords = np.array(member_load.get("end_node_coords", [0, 0, 0]))
                
#                 # Check if member is at this z_level (either start or end is at level)
#                 if (abs(start_coords[2] - z_level) >= tol and 
#                     abs(end_coords[2] - z_level) >= tol):
#                     continue
                
#                 # Calculate member length
#                 member_length = np.linalg.norm(end_coords - start_coords)
#                 if member_length < 1e-6:
#                     continue
                
#                 # Get uniform load (assume z-component represents mass per unit length)
#                 uniform_load = member_load.get("uniform", [0, 0, 0])
#                 if isinstance(uniform_load, dict):  # Handle different formats
#                     uniform_z = abs(uniform_load.get("z", 0))
#                 else:
#                     uniform_z = abs(uniform_load[2])
                
#                 total_mass_member = uniform_z * member_length
                
#                 # Use center point for member mass location
#                 center = np.array(member_load.get("center", [0, 0, 0]))
                
#                 total_mass += total_mass_member
#                 weighted_x += total_mass_member * center[0]
#                 weighted_y += total_mass_member * center[1]
#                 weighted_z += total_mass_member * center[2]
                
#                 # Process point loads on members
#                 point_load = member_load.get("point", [0, 0, 0, 0])
#                 if isinstance(point_load, dict):  # Handle different formats
#                     point_mass = abs(point_load.get("z", 0))
#                     point_coords = np.array(member_load.get("point_coords", [0, 0, 0]))
#                 else:
#                     point_mass = abs(point_load[2])
#                     point_coords = np.array(member_load.get("point_coords", [0, 0, 0]))
                
#                 if point_mass > 0 and abs(point_coords[2] - z_level) < tol:
#                     total_mass += point_mass
#                     weighted_x += point_mass * point_coords[0]
#                     weighted_y += point_mass * point_coords[1]
#                     weighted_z += point_mass * point_coords[2]
        
#         # Calculate center of mass
#         if total_mass > 0:
#             com_x = weighted_x / total_mass
#             com_y = weighted_y / total_mass
#             com_z = weighted_z / total_mass
#         else:
#             com_x, com_y, com_z = 0, 0, z_level
        
#         print(f"Total mass at cumulative z_level {z_level}: {total_mass:.2f} units")
#         print(f"Center of mass: ({com_x:.2f}, {com_y:.2f}, {com_z:.2f})")
        
#         com_results[z_level] = {
#             'total_mass': total_mass,
#             'center_of_mass': [com_x, com_y, com_z],
#             'z_level': z_level
#         }
    
#     # Save results to JSON file
#     output_path = os.path.join(JSON_FOLDER, "center_of_mass_results.json")
#     with open(output_path, 'w') as f:
#         json.dump(com_results, f, indent=2)

#     return com_results












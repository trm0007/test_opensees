# from build_model import build_model
from import_ import *
from Grid_and_structure_creation import load_structure
# from test import create_combined_structure_json
# from wall_meshing import create_combined_structure_json
from units import *
import os
import json
import numpy as np
from typing import Dict, List, Union


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

    # Create output directory if it doesn't exist
    save_dir = os.path.join(JSON_FOLDER, "load_data")
    os.makedirs(save_dir, exist_ok=True)
    
    all_success = True
    
    # Process each combination
    for combo_name, components in load_combinations.items():
        combined_loads = defaultdict(lambda: np.zeros(6))
        
        try:
            # Combine load cases with factors
            for load_case, factor in components:
                if load_case not in nodal_loads:
                    raise KeyError(f"Load case '{load_case}' not found in nodal_loads")
                
                for node_id, loads in nodal_loads[load_case].items():
                    loads_array = np.array(loads) if isinstance(loads, (list, tuple)) else loads
                    combined_loads[node_id] += loads_array * factor
            
            # Prepare output data
            output = {
                node_id: loads.tolist()
                for node_id, loads in combined_loads.items()
            }
            
            # Create filename based on combination name
            filename = f"nodal_loads_{combo_name.lower()}.json"
            save_path = os.path.join(save_dir, filename)
            
            with open(save_path, 'w') as f:
                json.dump(output, f, indent=2)
            print(f"Successfully saved combination '{combo_name}' to: {save_path}")
            
        except (IOError, OSError, KeyError) as e:
            all_success = False
            print(f"Warning: Failed to process combination '{combo_name}': {e}")
    
    return all_success


def process_element_loads_with_node_coordinates(loading_mapping, nodes_data, members_data, JSON_FOLDER):

    # Create node lookup dictionaries
    node_by_id = {node['id']: node for node in nodes_data}
    node_by_name = {node['name']: node for node in nodes_data}

    # Build member_nodes dictionary
    member_nodes = {}
    for member in members_data:
        member_name = member['name']
        member_id = member['id']
        
        # Get start node
        start_node = (node_by_id.get(member['start_node_id']) or 
                     node_by_name.get(member['start_node_name']))
        
        # Get end node
        end_node = (node_by_id.get(member['end_node_id']) or 
                   node_by_name.get(member['end_node_name']))

        member_nodes[member_name] = {
            'id': member_id,
            'start_node': [start_node['x'], start_node['y'], start_node['z']],
            'end_node': [end_node['x'], end_node['y'], end_node['z']]
        }

    # Process loads
    final_json = {"element_loads": []}
    
    for member_name_list, loads in loading_mapping:
        print("member_name_list, loads")
        print(member_name_list, loads)

        # Process each individual member
        for member_name_str in member_name_list:
            print("member_name_str")
            print(member_name_str)
            member_name_str = str(member_name_str).strip()
            print("member_name_str")
            print(member_name_str)
            member_data = member_nodes[member_name_str]
            
            for load in loads:
                final_json["element_loads"].append({
                    "member_name": member_name_str,
                    "member_id": member_data['id'],
                    "start_node_coords": member_data['start_node'],
                    "end_node_coords": member_data['end_node'],  
                    "load_case": load['LoadCase'],
                    "load_data": load
                })
            print("final_json")
            print(final_json)

    # Save to file
    output_file = os.path.join(JSON_FOLDER, "element_loads_with_coordinates.json")
    with open(output_file, 'w') as f:
        json.dump(final_json, f, indent=2)
    
    return final_json


def apply_member_load_combinations(load_combinations, JSON_FOLDER):
    """Apply load combinations to member loads and save combined results to files"""
    
    # Load the element loads with coordinates file
    input_file = os.path.join(JSON_FOLDER, "element_loads_with_coordinates.json")
    try:
        with open(input_file, 'r') as f:
            final_json_with_coords = json.load(f)
    except Exception as e:
        print(f"Error loading element loads file: {e}")
        return False

    if "element_loads" in final_json_with_coords and final_json_with_coords["element_loads"]:
        first_load = final_json_with_coords["element_loads"][0]
        print(f"First element_load keys: {list(first_load.keys())}")
        print(f"First element_load sample: {first_load}")
    
    wall_dir = os.path.join(JSON_FOLDER, 'load_data')
    os.makedirs(wall_dir, exist_ok=True)
    
    for combo_name, factors in load_combinations.items():
        print(f"\nProcessing load combination: {combo_name}")
        print(f"Load combination factors: {factors}")
        
        combo_loads = defaultdict(list)
        total_uniform = [0.0, 0.0, 0.0]
        total_point = [0.0, 0.0, 0.0]
        
        for member_load in final_json_with_coords["element_loads"]:
            member_id = member_load.get("member_id", -1)
            member_name = member_load.get("member_name", f"member_{member_id}")
            
            # Handle load case (could be string or list)
            load_case = member_load.get("load_case", "")
            if isinstance(load_case, list) and len(load_case) > 0:
                load_case = load_case[0]
            
            load_data = member_load.get("load_data", {})
            
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
            }

            # Extract uniform loads - handle both dict and list formats
            uniform_loads = [0.0, 0.0, 0.0]
            if 'uniform' in load_data:
                if isinstance(load_data['uniform'], dict):
                    uniform_loads = [
                        load_data['uniform'].get('x', 0.0),
                        load_data['uniform'].get('y', 0.0),
                        load_data['uniform'].get('z', 0.0)
                    ]
                elif isinstance(load_data['uniform'], list):
                    uniform_loads = list(load_data['uniform'])[:3]
            
            # Extract point loads - handle both dict and list formats
            point_loads = [0.0, 0.0, 0.0, 0.0]
            if 'point' in load_data:
                if isinstance(load_data['point'], dict):
                    point_loads = [
                        load_data['point'].get('x', 0.0),
                        load_data['point'].get('y', 0.0),
                        load_data['point'].get('z', 0.0),
                        load_data['point'].get('location', 0.0)
                    ]
                elif isinstance(load_data['point'], list):
                    point_loads = list(load_data['point'])[:4]
            
            # Extract temperature points
            temp_points = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
            if 'temperature_points' in load_data:
                if isinstance(load_data['temperature_points'], list):
                    for i, tp in enumerate(load_data['temperature_points'][:3]):
                        if isinstance(tp, dict):
                            temp_points[i] = [tp.get('temp', 0.0), tp.get('y', 0.0)]
                        elif isinstance(tp, (list, tuple)):
                            temp_points[i] = list(tp)[:2]
            
            # Apply load combination factors
            for load_case_name, factor in factors:
                if str(load_case_name).lower() == str(load_case).lower():
                    print(f"Applying factor {factor} to load case '{load_case}' for member {member_name}")
                    
                    # Apply factors to uniform loads
                    combined_load["uniform"][0] += float(uniform_loads[0]) * factor
                    combined_load["uniform"][1] += float(uniform_loads[1]) * factor
                    combined_load["uniform"][2] += float(uniform_loads[2]) * factor
                    
                    # Apply factors to point loads (except location)
                    combined_load["point"][0] += float(point_loads[0]) * factor
                    combined_load["point"][1] += float(point_loads[1]) * factor
                    combined_load["point"][2] += float(point_loads[2]) * factor
                    combined_load["point"][3] = float(point_loads[3])  # Location doesn't get factored
                    
                    # Apply factors to temperature points
                    for i in range(3):
                        combined_load["temperature_points"][i][0] += float(temp_points[i][0]) * factor
                        combined_load["temperature_points"][i][1] = float(temp_points[i][1])
            
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
    Only uses x and y coordinates from nodes for COM calculation, z comes from level data.
    
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
    # print(f'level_data={level_data}')
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
    com_id_start = 10000000
    while com_id_start in existing_ids:
        com_id_start += 1
    
    results = {}
    tol = 0.01
    
    for z_level in level_data.keys():
        print(f"\nProcessing z_level: {z_level:.2f}")
        
        total_mass = 0.0
        weighted_sum_x = 0.0  # For x coordinate only
        weighted_sum_y = 0.0  # For y coordinate only
        
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
                weighted_sum_x += mass * node_coords[0]  # Only x coordinate
                weighted_sum_y += mass * node_coords[1]  # Only y coordinate
        
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
                weighted_sum_x += mass_per_length * length * center[0]  # Only x
                weighted_sum_y += mass_per_length * length * center[1]  # Only y
            
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
                weighted_sum_x += point_mass * point_coords[0]  # Only x
                weighted_sum_y += point_mass * point_coords[1]  # Only y
        
        # Calculate center of mass (x and y only)
        if total_mass > 0:
            com_x = weighted_sum_x / total_mass
            com_y = weighted_sum_y / total_mass
        else:
            com_x, com_y = 0.0, 0.0
        
        # Use the z_level for the z-coordinate
        com_z = z_level
        
        # Generate unique name and ID for COM point
        com_name = f"COM_{z_level:.2f}"
        com_id = com_id_start
        com_id_start += 1
        
        # Ensure name is unique
        while com_name in node_name_to_id:
            com_name = f"COM_{z_level:.2f}_{com_id}"
        
        print(f"Total mass: {total_mass:.2f} kg")
        print(f"Center of mass: ({com_x:.2f}, {com_y:.2f}, {com_z:.2f})")
        print(f"Assigned COM node: {com_name} (ID: {com_id})")
        
        # Create COM node data
        com_node = {
            'id': com_id,
            'name': com_name,
            'x': float(com_x),
            'y': float(com_y),
            'z': float(com_z),
            'is_com_node': True,
            'mass': total_mass
        }
        
        results[z_level] = {
            'total_mass': total_mass,
            'center_of_mass': [com_x, com_y, com_z],
            'z_level': z_level,
            'com_node': com_node,
            'nodes_at_level': level_data[z_level]['nodes'],
        }
    
    # Save results
    output_path = os.path.join(JSON_FOLDER, "center_of_mass_results.json")
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    return results

        
def create_zero_length_elements(JSON_FOLDER, zero_length_nodes):

    nodes_data, _, _ = create_combined_structure_json(JSON_FOLDER)

    name_to_node = {n["name"]: n for n in nodes_data}

    spring_tags = []
    spring_node_start = 10001
    element_tag_start = 20001
    mat_tag_start = 1

    for i, zl_node in enumerate(zero_length_nodes):
        name = zl_node["name"]
        node = name_to_node[name]
        node_id = node["id"]
        x, y, z = node["x"], node["y"], node["z"]
        Kx = float(zl_node["Kx"])
        Ky = float(zl_node["Ky"])
        Kz = float(zl_node["Kz"])

        ops.node(node_id, x, y, z)

        spring_node_id = node_id + spring_node_start
        ops.node(spring_node_id, x, y, z)
        ops.fix(spring_node_id, 1, 1, 1, 1, 1, 1)

        mat_x = mat_tag_start + i * 3
        mat_y = mat_tag_start + i * 3 + 1
        mat_z = mat_tag_start + i * 3 + 2
        ops.uniaxialMaterial("ENT", mat_x, Kx)
        ops.uniaxialMaterial("ENT", mat_y, Ky)
        ops.uniaxialMaterial("ENT", mat_z, Kz)

        element_id = element_tag_start + i
        ops.element("zeroLength", element_id, node_id, spring_node_id,
                    "-mat", mat_x, mat_y, mat_z,
                    "-dir", 1, 2, 3)

        print(f"element zeroLength {element_id} {node_id} {spring_node_id} -mat {mat_x} {mat_y} {mat_z} -dir 1 2 3")

        spring_tags.append({
            "spring_node": spring_node_id,
            "element_id": element_id
        })

    return spring_tags


load_combinationsiiii = [

    # Modified Load Combinations (Seismic Category B)
    "1.4*DL",
    "1.2*DL + 1.6*LL + 0.5*Lr",
    "1.2*DL + 1.6*Lr + LL",
    "1.2*DL + 1.6*Lr + 0.8*Wx",
    "1.2*DL + 1.6*Lr - 0.8*Wx",
    "1.2*DL + 1.6*Lr + 0.8*Wy",
    "1.2*DL + 1.6*Lr - 0.8*Wy",
    "1.2*DL + 1.6*W + LL + 0.5*Lr",
    "1.2*DL - 1.6*Wx + LL + 0.5*Lr",
    "1.2*DL + 1.6*Wy + LL + 0.5*Lr",
    "1.2*DL - 1.6*Wy + LL + 0.5*Lr",
    "1.2*DL + EQx + LL",
    "1.2*DL - EQx + LL",
    "1.2*DL + EQy + LL",
    "1.2*DL - EQy + LL",
    "0.9*DL + 1.6*Wx",
    "0.9*DL - 1.6*Wx",
    "0.9*DL + 1.6*Wy",
    "0.9*DL - 1.6*Wy",
    "0.9*DL + EQx",
    "0.9*DL - EQx",
    "0.9*DL + EQy",
    "0.9*DL - EQy",

    # Seismic Load Combinations for Category C(v) and D (include vertical EQ Ev = 0.11*DL)
    "1.2*DL + EQx + 0.3*EQy + LL + 0.11*DL",
    "1.2*DL - EQx + 0.3*EQy + LL + 0.11*DL",
    "1.2*DL + EQx - 0.3*EQy + LL + 0.11*DL",
    "1.2*DL - EQx - 0.3*EQy + LL + 0.11*DL",
    "1.2*DL + 0.3*EQx + EQy + LL + 0.11*DL",
    "1.2*DL - 0.3*EQx + EQy + LL + 0.11*DL",
    "1.2*DL + 0.3*EQx - EQy + LL + 0.11*DL",
    "1.2*DL - 0.3*EQx - EQy + LL + 0.11*DL",
    "0.9*DL + EQx + 0.3*EQy + LL + 0.11*DL",
    "0.9*DL - EQx + 0.3*EQy + LL + 0.11*DL",
    "0.9*DL + EQx - 0.3*EQy + LL + 0.11*DL",
    "0.9*DL - EQx - 0.3*EQy + LL + 0.11*DL",
    "0.9*DL + 0.3*EQx + EQy + LL + 0.11*DL",
    "0.9*DL - 0.3*EQx + EQy + LL + 0.11*DL",
    "0.9*DL + 0.3*EQx - EQy + LL + 0.11*DL",
    "0.9*DL - 0.3*EQx - EQy + LL + 0.11*DL"
]


load_combinations = {
    "mass": [("DL", 1.0), ("LL", 0.25), ("self_weight", 1.0)],
    "Comb2": [("DL", 1.2), ("LL", 1.6), ("self_weight", 1.2)],
    "Comb1": [("DL", 1.4), ("self_weight", 1.4)],
    "unfactored_load": [("DL", 1.0), ("LL", 1.0),  ("self_weight", 1.0)],
}














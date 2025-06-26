from import_ import *
from units import *


# =============================================
# Grid and Structure Creation
# =============================================
def create_grid(x_spacing, y_spacing, z_spacing):
    """Ensure spacing values are clean floats"""
    def validate_spacing(spacing_list):
        return [float(round(s, 6)) for s in spacing_list]  # Round to avoid tiny offsets
    
    x_spacing = validate_spacing(x_spacing)
    y_spacing = validate_spacing(y_spacing)
    z_spacing = validate_spacing(z_spacing)
    
    def generate_points(spacing_list):
        points = [0.0]
        for spacing in spacing_list:
            points.append(points[-1] + spacing)
        return points
    
    return generate_points(x_spacing), generate_points(y_spacing), generate_points(z_spacing)

def create_structure(x_spacing, y_spacing, z_spacing):
    """Create nodes and members with clean numeric representation"""
    x_points, y_points, z_points = create_grid(x_spacing, y_spacing, z_spacing)
    
    # Node creation with consistent coordinates
    nodes = {}
    node_counter = 1
    for z in z_points:
        for y in y_points:
            for x in x_points:
                nodes[f"n{node_counter}"] = (float(x), float(y), float(z))
                node_counter += 1
    
    # Members creation with robust indexing
    members = {}
    member_id = 1
    num_x = len(x_points)
    num_y = len(y_points)
    num_z = len(z_points)
    
    # Column members along z
    for yi in range(num_y):
        for xi in range(num_x):
            for zi in range(num_z - 1):
                n1 = f"n{zi*num_y*num_x + yi*num_x + xi + 1}"
                n2 = f"n{(zi+1)*num_y*num_x + yi*num_x + xi + 1}"
                members[f"cz{member_id}"] = (n1, n2, member_id)
                member_id += 1
    
    # Beam members along x
    for zi in range(num_z):
        for yi in range(num_y):
            for xi in range(num_x - 1):
                n1 = f"n{zi*num_y*num_x + yi*num_x + xi + 1}"
                n2 = f"n{zi*num_y*num_x + yi*num_x + xi + 2}"
                members[f"bx{member_id}"] = (n1, n2, member_id)
                member_id += 1
    
    # Beam members along y
    for zi in range(num_z):
        for yi in range(num_y - 1):
            for xi in range(num_x):
                n1 = f"n{zi*num_y*num_x + yi*num_x + xi + 1}"
                n2 = f"n{zi*num_y*num_x + (yi+1)*num_x + xi + 1}"
                members[f"by{member_id}"] = (n1, n2, member_id)
                member_id += 1
    
    
    return nodes, members, x_points, y_points, z_points

# =============================================
# JSON Handling Functions (unchanged)
# =============================================


def save_structure_to_json(nodes, members):
    """Save structure data to JSON file with consistent formatting"""
    
    
    # Convert nodes to JSON format
    nodes_json = []
    node_id_map = {}
    for idx, (name, coords) in enumerate(nodes.items(), start=1):
        nodes_json.append({
            "id": idx,
            "name": name,
            "x": coords[0],
            "y": coords[1],
            "z": coords[2]
        })
        node_id_map[name] = idx
    
    # Convert members to JSON format
    members_json = []
    for name, (n1, n2, mid) in members.items():
        members_json.append({
            "id": mid,
            "name": name,
            "start_node_id": node_id_map[n1],
            "end_node_id": node_id_map[n2],
            "start_node_name": n1,
            "end_node_name": n2
        })
    structure_data = {"nodes": nodes_json, "members": members_json}

    return structure_data

def update_structure(OUTPUT_FOLDER, x_spacing, y_spacing, z_spacing, create_new_nodes={}, create_new_members={}, delete_nodes=[], delete_members=[], tolerance=1e-6):
    """
    Update a structural model by adding new nodes, adding new members, and deleting nodes and members.
    Prevents duplication based on both names and geometry.
    
    Parameters:
    - json_filepath: Path to the JSON file containing the structure data
    - create_new_nodes: Dictionary of new nodes {node_name: (x, y, z)}
    - create_new_members: Dictionary of new members {member_name: (start_node_name, end_node_name)}
    - delete_nodes: List of node names to delete
    - delete_members: List of member names to delete
    - tolerance: Floating point tolerance for coordinate comparison
    
    Returns:
    - Updated structure data
    """

    
    # full_path = os.path.join(OUTPUT_FOLDER, json_filepath)
    
    
    # if os.path.exists(full_path):
    #     with open(full_path, "r") as f:
    #         structure_data = json.load(f)
    # else:
    #     structure_data = {"nodes": [], "members": []}

    nodes, members, x_points, y_points, z_points = create_structure(x_spacing, y_spacing, z_spacing)

    structure_data = save_structure_to_json(nodes, members)

    # Create mapping of existing node names to IDs
    node_name_to_id = {}
    max_node_id = 0
    for node in structure_data["nodes"]:
        node_name_to_id[node["name"]] = node["id"]
        max_node_id = max(max_node_id, node["id"])
    
    # Create coordinate mapping to detect geometric duplicates
    coordinate_to_node = {}
    for node in structure_data["nodes"]:
        coord_key = (round(node["x"]), 
                     round(node["y"]), 
                     round(node["z"]))
        coordinate_to_node[coord_key] = node["name"]
    
    # Delete specified nodes
    if delete_nodes:
        # First find all members connected to nodes being deleted
        connected_members = []
        for member in structure_data["members"]:
            if member["start_node_name"] in delete_nodes or member["end_node_name"] in delete_nodes:
                connected_members.append(member["name"])
        
        # Add these members to the delete_members list if not already there
        delete_members = list(set(delete_members + connected_members))
        
        # Remove nodes
        structure_data["nodes"] = [node for node in structure_data["nodes"] 
                                 if node["name"] not in delete_nodes]
        # Update node_name_to_id mapping
        for node_name in delete_nodes:
            node_name_to_id.pop(node_name, None)
            
        # Update coordinate mapping
        coordinate_to_node = {}
        for node in structure_data["nodes"]:
            coord_key = (round(node["x"]/tolerance)*tolerance, 
                         round(node["y"]/tolerance)*tolerance, 
                         round(node["z"]/tolerance)*tolerance)
            coordinate_to_node[coord_key] = node["name"]
    
    # Delete specified members
    if delete_members:
        structure_data["members"] = [member for member in structure_data["members"] 
                                   if member["name"] not in delete_members]
    
    # Create a set to track node pairs for member duplication checks
    member_node_pairs = set()
    for member in structure_data["members"]:
        # Store both directions to handle undirected members
        member_node_pairs.add((member["start_node_name"], member["end_node_name"]))
        member_node_pairs.add((member["end_node_name"], member["start_node_name"]))
    
    # Add new nodes (with geometric duplicate check)
    nodes_added = {}
    for node_name, coordinates in create_new_nodes.items():
        if node_name in node_name_to_id:
            print(f"Node {node_name} already exists. Skipping.")
            continue
            
        # Check for geometric duplicates
        x, y, z = coordinates
        coord_key = (round(x/tolerance)*tolerance, 
                     round(y/tolerance)*tolerance, 
                     round(z/tolerance)*tolerance)
        
        if coord_key in coordinate_to_node:
            existing_node = coordinate_to_node[coord_key]
            print(f"Node {node_name} at ({x}, {y}, {z}) is a geometric duplicate of existing node {existing_node}. Using existing node.")
            nodes_added[node_name] = existing_node
            continue
            
        max_node_id += 1
        new_node = {
            "id": max_node_id,
            "name": node_name,
            "x": x,
            "y": y,
            "z": z
        }
        structure_data["nodes"].append(new_node)
        node_name_to_id[node_name] = max_node_id  # Update the mapping
        coordinate_to_node[coord_key] = node_name  # Update coordinate mapping
    
    # Find max member ID
    max_member_id = 0
    existing_member_names = set()
    for member in structure_data["members"]:
        max_member_id = max(max_member_id, member["id"])
        existing_member_names.add(member["name"])
    
    # Add new members (with duplicate node pair check)
    for member_name, (start_node_name, end_node_name) in create_new_members.items():
        if member_name in existing_member_names:
            print(f"Member {member_name} already exists. Skipping.")
            continue
        
        # Replace with existing geometric duplicate nodes if applicable
        if start_node_name in nodes_added:
            start_node_name = nodes_added[start_node_name]
        if end_node_name in nodes_added:
            end_node_name = nodes_added[end_node_name]
            
        if start_node_name not in node_name_to_id:
            raise ValueError(f"Start node '{start_node_name}' for member '{member_name}' does not exist.")
        if end_node_name not in node_name_to_id:
            raise ValueError(f"End node '{end_node_name}' for member '{member_name}' does not exist.")
        
        # Check for duplicate node pairs (members connecting same nodes)
        if (start_node_name, end_node_name) in member_node_pairs or (end_node_name, start_node_name) in member_node_pairs:
            print(f"Member {member_name} connecting {start_node_name} and {end_node_name} is a duplicate of an existing member. Skipping.")
            continue
        
        max_member_id += 1
        new_member = {
            "id": max_member_id,
            "name": member_name,
            "start_node_id": node_name_to_id[start_node_name],
            "end_node_id": node_name_to_id[end_node_name],
            "start_node_name": start_node_name,
            "end_node_name": end_node_name
        }
        structure_data["members"].append(new_member)
        # Add node pair to set to prevent future duplicates
        member_node_pairs.add((start_node_name, end_node_name))
        member_node_pairs.add((end_node_name, start_node_name))
    
    nodes_full_path = os.path.join(OUTPUT_FOLDER, "node_input_data.json")
    members_full_path = os.path.join(OUTPUT_FOLDER, "members_input_data.json")
    
    # Save nodes data
    with open(nodes_full_path, "w") as f:
        json.dump(structure_data["nodes"], f, indent=4)

    # Save members data
    with open(members_full_path, "w") as f:
        json.dump(structure_data["members"], f, indent=4)
    
        
    return structure_data




# =============================================
# Filtering Functions
# =============================================

def load_structure(OUTPUT_FOLDER):
    """
    Load structure data from separate node and member files and return in specific format.
    
    Parameters:
    - OUTPUT_FOLDER: Path to the folder containing the node and member JSON files
    
    Returns:
    - nodes: Dictionary {node_name: (x, y, z)}
    - members: Dictionary {member_name: (start_node_name, end_node_name, member_id)}
    - structure: Complete structure dictionary with original format
    """
    
    nodes_path = os.path.join(OUTPUT_FOLDER, "node_input_data.json")
    members_path = os.path.join(OUTPUT_FOLDER, "members_input_data.json")
    
    # Initialize return values
    nodes_dict = {}
    members_dict = {}
    structure = {"nodes": [], "members": []}
    
    # Load nodes
    try:
        with open(nodes_path, "r") as f:
            nodes_data = json.load(f)
            structure["nodes"] = nodes_data
            nodes_dict = {node['name']: (node['x'], node['y'], node['z']) for node in nodes_data}
    except FileNotFoundError:
        print(f"Warning: Node file not found at {nodes_path}")
    except json.JSONDecodeError:
        print(f"Error: Could not parse node file at {nodes_path}")
    
    # Load members
    try:
        with open(members_path, "r") as f:
            members_data = json.load(f)
            structure["members"] = members_data
            
            # Create mapping from node_id to node_name
            node_id_to_name = {node['id']: node['name'] for node in structure["nodes"]}
            
            for member in members_data:
                # Get node names (using ID mapping for robustness)
                start_name = node_id_to_name.get(member['start_node_id'], member['start_node_name'])
                end_name = node_id_to_name.get(member['end_node_id'], member['end_node_name'])
                
                members_dict[member['name']] = (
                    start_name,
                    end_name,
                    member['id']
                )
    except FileNotFoundError:
        print(f"Warning: Member file not found at {members_path}")
    except json.JSONDecodeError:
        print(f"Error: Could not parse member file at {members_path}")
    
    return nodes_dict, members_dict, structure



# def load_structure_from_json( OUTPUT_FOLDER,json_file):
#     """Load structure data from JSON file"""
#     full_path = os.path.join( OUTPUT_FOLDER, json_file)
#     with open(full_path, 'r') as f:
#         data = json.load(f)
    
#     nodes = {}
#     for node in data['nodes']:
#         nodes[node['name']] = (node['x'], node['y'], node['z'])
    
#     members = {}
#     node_id_to_name = {node['id']: node['name'] for node in data['nodes']}
#     for member in data['members']:
#         members[member['name']] = (
#             member['start_node_name'],
#             member['end_node_name'],
#             member['id']
#         )
    
#     return nodes, members

def get_z_levels(nodes):
    """Extract unique Z levels from nodes"""
    z_levels = sorted({z for _, (_, _, z) in nodes.items()})
    return z_levels

def plot_structure(nodes, members, OUTPUT_FOLDER, title="Full Structure", save_path=None):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    for name, (x, y, z) in nodes.items():
        ax.scatter(x, y, z, color='blue', s=50)
        ax.text(x, y, z, name, fontsize=10)
    
    for name, (n1, n2, mid) in members.items():
        x_vals = [nodes[n1][0], nodes[n2][0]]
        y_vals = [nodes[n1][1], nodes[n2][1]]
        z_vals = [nodes[n1][2], nodes[n2][2]]
        ax.plot(x_vals, y_vals, z_vals, color='black', linewidth=2)
        
        x_mid = (x_vals[0] + x_vals[1]) / 2
        y_mid = (y_vals[0] + y_vals[1]) / 2
        z_mid = (z_vals[0] + z_vals[1]) / 2
        ax.text(x_mid, y_mid, z_mid, name, fontsize=10, color='red')
    
    ax.set_xlabel('X (m)', fontsize=12)
    ax.set_ylabel('Y (m)', fontsize=12)
    ax.set_zlabel('Z (m)', fontsize=12)
    ax.set_title(title, fontsize=14)
    plt.tight_layout()
    
    if save_path:
        full_path = os.path.join( OUTPUT_FOLDER, save_path)
        os.makedirs(os.path.dirname(full_path), exist_ok=True)
        plt.savefig(full_path)
        print(f"Saved structure plot to {full_path}")
    plt.close()


def plot_filtered_elements(filtered_items, nodes, output_folder, title_prefix="Filtered", save_path=None):
    """
    Plot filtered structural elements (nodes or members) grouped by Z-level.
    Only plots the filtered items without showing all nodes in the background.
    
    Args:
        filtered_items (dict): Dictionary of filtered items grouped by Z-level
        nodes (dict): Complete dictionary of all nodes in the structure (used for coordinates)
        output_folder (str): Root output directory for saving images
        title_prefix (str): Prefix for the plot title
        save_path (str): Relative path to save the image (without extension)
    """
    for z_level, items in filtered_items.items():
        # Create figure and 3D axis
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Check if we have items to plot
        if isinstance(items, dict) and items:
            first_item = next(iter(items.values()))
            
            # Plot members if items are member data
            if isinstance(first_item, tuple) and len(first_item) == 3 and isinstance(first_item[0], str):
                for member_name, (n1, n2, mid) in items.items():
                    # Get node coordinates
                    x1, y1, z1 = nodes[n1]
                    x2, y2, z2 = nodes[n2]
                    
                    # Plot member line
                    ax.plot([x1, x2], [y1, y2], [z1, z2], 
                           color='red', linewidth=2)
                    
                    # Calculate midpoint for label
                    x_mid = (x1 + x2) / 2
                    y_mid = (y1 + y2) / 2
                    z_mid = (z1 + z2) / 2
                    
                    # Add member label
                    ax.text(x_mid, y_mid, z_mid, member_name, 
                           fontsize=16, color='black')
                    
                    # Highlight connected nodes
                    ax.scatter(x1, y1, z1, color='blue', s=50)
                    ax.scatter(x2, y2, z2, color='blue', s=50)
                    ax.text(x1, y1, z1, n1, fontsize=16)
                    ax.text(x2, y2, z2, n2, fontsize=16)
            
            # Plot nodes if items are node data
            else:
                for node_name, (x, y, z) in items.items():
                    ax.scatter(x, y, z, color='blue', s=50)
                    ax.text(x, y, z, node_name, fontsize=16)
        
        # Set axis labels and title
        ax.set_xlabel('X (m)', fontsize=20)
        ax.set_ylabel('Y (m)', fontsize=20)
        ax.set_zlabel('Z (m)', fontsize=20)
        ax.set_title(f"{title_prefix} at Z={z_level}m", fontsize=14)
        plt.tight_layout()
        
        # Save plot if requested
        if save_path:
            # Create full path by joining output_folder and save_path
            full_save_path = os.path.join(output_folder, save_path)
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(full_save_path), exist_ok=True)
            # Save the figure
            plt.savefig(full_save_path, dpi=300, bbox_inches='tight')
            print(f"Saved filtered plot to {full_save_path}")
        
        # Close the figure to free memory
        plt.close()

def filter_nodes_by_z(nodes, z_points):
    """Filter and group nodes by Z level"""
    grouped_nodes = defaultdict(dict)
    
    for node_name, (x, y, z) in nodes.items():
        for z_level in z_points:
            if abs(z - z_level) < 1e-6:
                grouped_nodes[z_level][node_name] = (x, y, z)
                break
    
    return dict(grouped_nodes)

def filter_column_members_by_z(members, nodes, z_points):
    """Filter and group column members by starting Z level"""
    grouped_members = defaultdict(dict)
    
    column_members = {name: data for name, data in members.items() if name.startswith('cz')}
    
    for member_name, (n1, n2, mid) in column_members.items():
        n1_z = nodes[n1][2]
        for z_level in z_points[:-1]:  # No columns start at the top level
            if abs(n1_z - z_level) < 1e-6:
                grouped_members[z_level][member_name] = (n1, n2, mid)
                break
    
    return dict(grouped_members)

def filter_beam_x_by_z(members, nodes, z_points):
    """Filter and group beam members along X by Z level"""
    grouped_members = defaultdict(dict)
    
    beam_x_members = {name: data for name, data in members.items() if name.startswith('bx')}
    
    for member_name, (n1, n2, mid) in beam_x_members.items():
        n1_z = nodes[n1][2]
        for z_level in z_points:
            if abs(n1_z - z_level) < 1e-6:
                grouped_members[z_level][member_name] = (n1, n2, mid)
                break
    
    return dict(grouped_members)

def filter_beam_y_by_z(members, nodes, z_points):
    """Filter and group beam members along Y by Z level"""
    grouped_members = defaultdict(dict)
    
    beam_y_members = {name: data for name, data in members.items() if name.startswith('by')}
    
    for member_name, (n1, n2, mid) in beam_y_members.items():
        n1_z = nodes[n1][2]
        for z_level in z_points:
            if abs(n1_z - z_level) < 1e-6:
                grouped_members[z_level][member_name] = (n1, n2, mid)
                break
    
    return dict(grouped_members)


from collections import defaultdict

def create_member_name_lists(filtered_nodes, filtered_columns, filtered_beams_x, filtered_beams_y):
    grouped_members = defaultdict(list)
    grouped_nodes = defaultdict(list)
    
    grouped_members['columns'] = [name for group in filtered_columns.values() for name in group]
    grouped_members['beamx'] = [name for group in filtered_beams_x.values() for name in group]
    grouped_members['beamy'] = [name for group in filtered_beams_y.values() for name in group]
    grouped_nodes['nodes'] = [name for group in filtered_nodes.values() for name in group]
    
    # Option 1: Return a merged dictionary
    return {**grouped_members, **grouped_nodes}



def create_member_name_lists_with_rotation(filtered_nodes, filtered_columns, filtered_beams_x, filtered_beams_y):
    grouped_members = defaultdict(list)
    grouped_nodes = defaultdict(list)

    grouped_members['columns'] = [name for group in filtered_columns.values() for name in group]
    grouped_members['beamx'] = [name for group in filtered_beams_x.values() for name in group]
    grouped_members['beamy'] = [name for group in filtered_beams_y.values() for name in group]
    grouped_nodes['nodes'] = [name for group in filtered_nodes.values() for name in group]

    member_dict = {}
    for key, names in {**grouped_members, **grouped_nodes}.items():
        member_dict[key] = [{"name": name, "rotation": 0} for name in names]

    return member_dict

# def create_section_mapping_from_filtered(filtered_columns_path, filtered_beams_x_path, filtered_beams_y_path, output_file="section_mapping.json"):
#     """
#     Create a section mapping organized by Z-levels using pre-filtered JSON files.
    
#     Args:
#         filtered_columns_path: Path to filtered_columns.json
#         filtered_beams_x_path: Path to beamX.json
#         filtered_beams_y_path: Path to beamY.json
#         output_file: Path to save the JSON output
#     """
#     # Load filtered data
#     with open(filtered_columns_path, 'r') as f:
#         filtered_columns = json.load(f)
    
#     with open(filtered_beams_x_path, 'r') as f:
#         filtered_beams_x = json.load(f)
    
#     with open(filtered_beams_y_path, 'r') as f:
#         filtered_beams_y = json.load(f)
    
#     # Initialize the mapping structure
#     section_mapping = {
#         "columns": defaultdict(dict),
#         "beams_x": defaultdict(dict),
#         "beams_y": defaultdict(dict)
#     }
    
#     # Process columns (already grouped by starting Z-level)
#     for z_level, columns in filtered_columns.items():
#         for member_name, member_data in columns.items():
#             # The key format should be "at_z_X.X" where X.X is the z_level
#             z_key = f"at_z_{float(z_level):.1f}"
#             section_mapping["columns"][z_key][member_name] = "section1"
    
#     # Process x-beams (already grouped by Z-level)
#     for z_level, beams in filtered_beams_x.items():
#         for member_name, member_data in beams.items():
#             z_key = f"at_z_{float(z_level):.1f}"
#             section_mapping["beams_x"][z_key][member_name] = "section2"
    
#     # Process y-beams (already grouped by Z-level)
#     for z_level, beams in filtered_beams_y.items():
#         for member_name, member_data in beams.items():
#             z_key = f"at_z_{float(z_level):.1f}"
#             section_mapping["beams_y"][z_key][member_name] = "section3"
    
#     # Convert defaultdict to regular dict for JSON serialization
#     section_mapping = {
#         "columns": dict(section_mapping["columns"]),
#         "beams_x": dict(section_mapping["beams_x"]),
#         "beams_y": dict(section_mapping["beams_y"])
#     }
    
#     # Save to JSON file
#     with open(output_file, 'w') as f:
#         json.dump(section_mapping, f, indent=4)
    
#     print(f"Section mapping saved to {output_file}")
#     return section_mapping

# Run the function

def save_filtered_data(filtered_data, output_folder, filename):
    """Save filtered data to JSON file in the specified output folder"""
    full_path = os.path.join(output_folder, filename)
    os.makedirs(os.path.dirname(full_path), exist_ok=True)
    with open(full_path, 'w') as f:
        json.dump(filtered_data, f, indent=4)
    # print(f"Saved filtered data to {full_path}")

# def create_member_section_mapping(JSON_FOLDER, output_file, member_section_mapping):
#     """
#     Convert a section->members dict to member->section mapping and save to a JSON file.
    
#     Args:
#         JSON_FOLDER (str): folder to save JSON file
#         output_file (str): JSON filename
#         member_section_mapping (dict): dict with sections as keys and member lists as values
    
#     Returns:
#         dict: member->section mapping
#     """
#     output_path = os.path.join(JSON_FOLDER, output_file)

#     # Convert section->members to member->section mapping
#     member_to_section = {}
#     for section, members in member_section_mapping.items():
#         for member in members:
#             member_to_section[member] = section

#     # Save the member->section mapping to JSON
#     os.makedirs(os.path.dirname(output_path), exist_ok=True)
#     with open(output_path, 'w') as file:
#         json.dump(member_to_section, file, indent=4)

#     print(f"Member section mapping manually created and saved to {output_path}")
#     print(f"Mapped {len(member_to_section)} members to sections")

#     return member_to_section

def create_member_section_mapping(JSON_FOLDER, output_file, member_section_mapping):
    """
    Convert a section->members dict to member->section mapping and save to a JSON file.
    
    Args:
        JSON_FOLDER (str): folder to save JSON file
        output_file (str): JSON filename
        member_section_mapping (dict): dict with sections as keys and member dicts as values
    
    Returns:
        dict: member->section mapping (with rotation info preserved)
    """
    output_path = os.path.join(JSON_FOLDER, output_file)

    # Convert section->members to member->section mapping while preserving rotation info
    member_to_section = {}
    for section, members in member_section_mapping.items():
        for member_dict in members:
            member_name = member_dict["name"]
            # Store both section and rotation information
            member_to_section[member_name] = {
                "section": section,
                "rotation": member_dict["rotation"]
            }

    # Save the member->section mapping to JSON
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as file:
        json.dump(member_to_section, file, indent=4)

    print(f"Member section mapping manually created and saved to {output_path}")
    print(f"Mapped {len(member_to_section)} members to sections")

    return member_to_section


    
    

# =============================================
# Main Execution
# =============================================

# if __name__ == "__main__":

#     # Create the structure
#     x_spacing = [2*m, 3*m]
#     y_spacing = [4*m, 3*m]
#     z_spacing = [5*m, 3*m]

#     nodes, members, x_points, y_points, z_points = create_structure(x_spacing, y_spacing, z_spacing)

     
#     # # Save initial structure to JSON
#     save_structure_to_json(nodes, members, "lamb/structure_data_json.json")
    
#     # Example of updating the structure
#     json_filepath = "lamb/structure_data_json.json"
    
#     create_new_nodes = {
#         # "n28": (1, 1, 0),
#         # "n10": (2, 0, 10),
#         # "n11": (0, 4, 10),
#         # "n12": (2, 4, 10)
#     }
    
#     create_new_members = {

#         # "bx38": ("n28", "n10"),
#         # "bx39": ("n2", "n15"),
#         # "bx40": ("n2", "n18"),
#         # "by41": ("n2", "n6"),
#     }
#     # Nodes and members to delete
#     # Using the same variable names as in your function parameters
#     # delete_nodes = ["n5", "n6"]  # Nodes to delete
#     # delete_members = ["bx38", "bx39"]  # Members to delete
#     delete_nodes = []
#     delete_members = []
    
#     updated_structure = update_structure(
#         json_filepath,
#         create_new_nodes,
#         create_new_members,
#         delete_nodes,
#         delete_members
#     )
    
#     print(f"Structure updated successfully with {len(create_new_nodes)} new nodes and {len(create_new_members)} new members.")
    
#     # Load structure from JSON file
#     json_file = "lamb/structure_data_json.json"
#     nodes, members = load_structure_from_json(json_file)
    
#     # Get unique Z levels
#     z_points = get_z_levels(nodes)
#     print(f"Found Z levels: {z_points}")
    
#     # Plot full structure
#     plot_structure(nodes, members, save_path="lamb/images/structure_plot.png")
    
#     # Filter elements by Z level
#     filtered_nodes = filter_nodes_by_z(nodes, z_points)
#     filtered_columns = filter_column_members_by_z(members, nodes, z_points)
#     filtered_beams_x = filter_beam_x_by_z(members, nodes, z_points)
#     filtered_beams_y = filter_beam_y_by_z(members, nodes, z_points)
    
#     # Plot filtered elements with specific filenames
#     plot_filtered_elements(filtered_nodes, nodes, "Nodes at Z Level", save_path="lamb/images/nodes.png")
#     plot_filtered_elements(filtered_columns, nodes, "Columns starting at Z Level", save_path="lamb/images/columns.png")
#     plot_filtered_elements(filtered_beams_x, nodes, "X Beams at Z Level", save_path="lamb/images/beams_x.png")
#     plot_filtered_elements(filtered_beams_y, nodes, "Y Beams at Z Level", save_path="lamb/images/beams_y.png")
    
#     # Save filtered data to JSON files
#     save_filtered_data(filtered_nodes, "lamb/filtered_nodes.json")
#     save_filtered_data(filtered_columns, "lamb/filtered_columns.json")
#     save_filtered_data(filtered_beams_x, "lamb/beamX.json")
#     save_filtered_data(filtered_beams_y, "lamb/beamY.json")

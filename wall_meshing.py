from import_ import *
from load_combinations import create_combined_structure_json
from units import *

def add_new_shells(mesh_elements, node_names, add_shell):
    """
    Add new shells to the mesh based on predefined points
    
    Parameters:
    -----------
    mesh_elements: dict
        Dictionary containing the current mesh elements
    node_names: dict
        Dictionary mapping node IDs to their coordinates
    add_shell: dict
        Dictionary with shell names as keys and lists of point names as values
        
    Returns:
    --------
    dict: Updated mesh elements with new shells added
    """
    # Create a copy of the input mesh_elements to avoid modifying the original
    updated_mesh = mesh_elements.copy()
    
    # Determine starting mesh ID
    mesh_counter = 1 if not mesh_elements else max(elem_data['id'] for elem_data in mesh_elements.values()) + 1
    
    # Process each shell to be added
    for shell_name, points in add_shell.items():
        # Determine element type based on number of points
        if len(points) == 3:  # Triangle
            element_id = mesh_counter
            mesh_name = f"T{element_id}"
            shell_type = 'triangle'
        elif len(points) == 4:  # Quadrilateral
            element_id = mesh_counter
            mesh_name = f"R{element_id}"
            shell_type = 'rectangle'
        else:
            print(f"Warning: Shell {shell_name} has {len(points)} points, only triangles (3) and quads (4) are supported")
            continue
        
        # Extract node IDs and coordinates
        nodes = []
        coords = []
        
        # Find node IDs corresponding to predefined points
        for point_name in points:
            # Direct dictionary lookup instead of looping
            if point_name in node_names:
                nodes.append(point_name)
                coords.append(node_names[point_name])
            else:
                print(f"Warning: Point {point_name} not found in mesh nodes")
        
        # Only create shell if we have all required points
        if len(nodes) == len(points):
            updated_mesh[shell_name] = {
                'type': shell_type,
                'nodes': nodes,
                'coordinates': coords,
                'id': element_id
            }
            mesh_counter += 1


    
    return updated_mesh


def remove_shells(mesh_elements, remove_shell):
    """
    Remove specified shells from the mesh
    
    Parameters:
    -----------
    mesh_elements: dict
        Dictionary containing the current mesh elements
    remove_shell: list
        List of shell names to be removed
        
    Returns:
    --------
    dict: Updated mesh elements with specified shells removed
    """
    # Create a copy of the input mesh_elements to avoid modifying the original
    final_mesh = mesh_elements.copy()
    
    for shell_name in remove_shell:
        if shell_name in final_mesh:
            del final_mesh[shell_name]
            print(f"Removed shell {shell_name}")
        else:
            print(f"Warning: Shell {shell_name} not found in mesh")
    
    return final_mesh


def sort_nodes_anticlockwise(nodes, coords):
    """
    Sort nodes in anti-clockwise order around their centroid.
    
    Args:
        nodes: List of node IDs
        coords: List of corresponding 3D coordinates
        
    Returns:
        Tuple of (sorted_nodes, sorted_coords)
    """
    if len(nodes) <= 3:  # Triangles are already planar
        return nodes, coords
        
    # Calculate centroid
    centroid = np.mean(coords, axis=0)
    
    # Center the points
    centered = np.array(coords) - centroid
    
    # Compute normal vector using first 3 points
    normal = np.cross(centered[1] - centered[0], centered[2] - centered[0])
    normal = normal / np.linalg.norm(normal)
    
    # Create basis vectors for projection
    if abs(normal[0]) > 0.1 or abs(normal[1]) > 0.1:
        u = np.array([normal[1], -normal[0], 0])  # Orthogonal to normal in XY plane
    else:
        u = np.array([0, normal[2], -normal[1]])  # Orthogonal to normal in YZ plane
    u = u / np.linalg.norm(u)
    v = np.cross(normal, u)
    
    # Project points to 2D plane
    projected = []
    for point in centered:
        x_proj = np.dot(point, u)
        y_proj = np.dot(point, v)
        projected.append((x_proj, y_proj))
    
    # Calculate angles from centroid and sort
    angles = np.arctan2([p[1] for p in projected], [p[0] for p in projected])
    positive_angles = np.where(angles < 0, angles + 2*np.pi, angles)
    sorted_indices = np.argsort(positive_angles)
    
    # Return sorted nodes and coordinates
    sorted_nodes = [nodes[i] for i in sorted_indices]
    sorted_coords = [coords[i] for i in sorted_indices]
    
    return sorted_nodes, sorted_coords



def extract_node_coordinates(JSON_FOLDER, point_names):
    json_file = "node_input_data.json"
    json_path = os.path.join(JSON_FOLDER, json_file)
    with open(json_path, 'r') as f:
        nodes = json.load(f)
    coord_map = {node['name']: [node['x'], node['y'], node['z']] for node in nodes}
    return [coord_map[name] for name in point_names if name in coord_map]

def convert_surface_configurations(JSON_FOLDER, surface_configurations):
    new_configurations = {}
    for key, config in surface_configurations.items():
        point_names = config.get("points", [])
        if isinstance(point_names[0], list):
            continue  # skip if already converted to coordinates

        points = extract_node_coordinates(JSON_FOLDER, point_names)
        updated_points = [[x, y, z] for x, y, z in points]
        
        predefined_names = list(config.get("predefined_points", []))
        predefined_coords = extract_node_coordinates(JSON_FOLDER, predefined_names)
        updated_predefined = {
            name: np.array([x, y, z])
            for name, (x, y, z) in zip(predefined_names, predefined_coords)
        }

        new_configurations[key] = {
            "points": updated_points,
            "section_name": config["section_name"],
            "add_shell": config["add_shell"],
            "remove_shell": config["remove_shell"],
            "predefined_points": updated_predefined,
            "num_x_div": config["num_x_div"],
            "num_y_div": config["num_y_div"],
            "load_case_names": config["load_case_names"],
            "pressures": config["pressures"],
            "thickness": config["thickness"],
        }
    return new_configurations

def find_existing_node_id(new_coord, node_names, tolerance=1e-6):
    """Check if a node with the same coordinates already exists."""
    for node_id, coord in node_names.items():
        if np.linalg.norm(np.array(new_coord) - np.array(coord)) < tolerance:
            return node_id
    return None

def plot_mesh_data_separately(json_folder, location_name, IMAGE_FOLDER):
    # Paths to the files
    wall_dir = os.path.join(json_folder, 'wall')
    mesh_data_path = os.path.join(wall_dir, f'{location_name}.json')
    nodes_data_path = os.path.join(json_folder, 'node_input_data.json')

    # Load the data
    with open(mesh_data_path, 'r') as f:
        mesh_data = json.load(f)

    with open(nodes_data_path, 'r') as f:
        nodes_data = json.load(f)

    elements = mesh_data["elements"]  # ✅ correctly access elements
    nodes = nodes_data               # assuming it's a list of dicts

    node_dict = {n["name"]: (n["x"], n["y"]) for n in nodes}

    plt.figure(figsize=(10, 10))

    for element_name, element_data in elements.items():
        node_names = element_data["nodes"]
        try:
            coords = [node_dict[n] for n in node_names]
        except KeyError:
            continue
        coords.append(coords[0])
        xs, ys = zip(*coords)
        plt.plot(xs, ys, 'b-')
        for n in node_names:
            if n in node_dict:
                x, y = node_dict[n]
                plt.text(x, y, n, fontsize=14, color='green')
        centroid_x = sum(x for x, y in coords[:-1]) / len(coords[:-1])
        centroid_y = sum(y for x, y in coords[:-1]) / len(coords[:-1])
        plt.text(centroid_x, centroid_y, element_name, fontsize=10, color='red', ha='center')

    plt.axis("equal")
    plt.grid(True)
    plt.title("Mesh Plot with Node and Element Names")

    os.makedirs(IMAGE_FOLDER, exist_ok=True)
    image_filename = f"{location_name}.png"
    image_path = os.path.join(IMAGE_FOLDER, image_filename)
    plt.savefig(image_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved plot to: {image_path}")

def plot_mesh_data_z_axis(json_folder, IMAGE_FOLDER):
    """
    Plot all meshes from all location files with the same z-coordinate in a single image.
    Creates one image per unique z-value across all locations.
    """
    import os
    import json
    import matplotlib.pyplot as plt
    from collections import defaultdict
    import glob
    
    # Paths
    wall_dir = os.path.join(json_folder, 'wall')
    nodes_data_path = os.path.join(json_folder, 'node_input_data.json')
    
    # Load nodes data once
    with open(nodes_data_path, 'r') as f:
        nodes_data = json.load(f)
    
    # Create dictionary with x, y, z coordinates
    node_dict = {n["name"]: (n["x"], n["y"], n.get("z", 0)) for n in nodes_data}
    
    # Find all JSON files in wall directory
    json_files = glob.glob(os.path.join(wall_dir, '*.json'))

    
    # Group ALL elements by z-coordinate across all files
    z_groups = defaultdict(list)
    total_elements = 0
    skipped_elements = []
    
    for json_file in json_files:
        location_name = os.path.splitext(os.path.basename(json_file))[0]
        
        with open(json_file, 'r') as f:
            mesh_data = json.load(f)

        elements = mesh_data.get("elements", {})
        
        for element_name, element_data in elements.items():
            node_names = element_data.get("nodes", [])
            if not node_names:
                skipped_elements.append(f"{location_name}.{element_name} (no nodes)")
                continue
            
            # Get z-coordinate from first node
            first_node = node_names[0]
            if first_node not in node_dict:
                skipped_elements.append(f"{location_name}.{element_name} (missing node: {first_node})")
                continue
            
            z = node_dict[first_node][2]
            
            # Verify all nodes exist and have same z-coordinate
            missing_nodes = [n for n in node_names if n not in node_dict]
            if missing_nodes:
                skipped_elements.append(f"{location_name}.{element_name} (missing nodes: {missing_nodes})")
                continue
                
            z_coords = [node_dict[n][2] for n in node_names]
            if not all(z_coord == z for z_coord in z_coords):
                skipped_elements.append(f"{location_name}.{element_name} (mixed z-coordinates)")
                continue
            
            # Add to z-group with location info
            z_groups[z].append((f"{location_name}.{element_name}", element_data, location_name))
            total_elements += 1
    
    # Report skipped elements
    if skipped_elements:
        print(f"⚠️  Skipped {len(skipped_elements)} elements:")
        for skip_msg in skipped_elements[:5]:  # Show first 5
            print(f"   - {skip_msg}")
        if len(skipped_elements) > 5:
            print(f"   ... and {len(skipped_elements) - 5} more")
    
    if not z_groups:
        print("❌ No valid elements found to plot")
        return
    
    print(f"📊 Total valid elements: {total_elements}")
    print(f"🎯 Found {len(z_groups)} unique z-levels: {sorted(z_groups.keys())}")
    
    # Create output directory
    os.makedirs(IMAGE_FOLDER, exist_ok=True)
    
    # Create one combined image per unique z-value
    for z in sorted(z_groups.keys()):
        elements_group = z_groups[z]
        
        plt.figure(figsize=(16, 12))  # Larger for combining multiple locations
        
        print(f"🎨 Creating combined image for z = {z} with {len(elements_group)} elements")
        
        # Use single color for all elements
        plotted_nodes = set()  # Track plotted nodes to avoid duplicates
        
        plotted_count = 0
        for element_name, element_data, location_name in elements_group:
            node_names = element_data.get("nodes", [])
            
            try:
                coords = [node_dict[n][:2] for n in node_names]
            except KeyError:
                continue
            
            if len(coords) < 2:
                continue
            
            # Close the polygon if needed
            if coords[0] != coords[-1]:
                coords.append(coords[0])
            
            xs, ys = zip(*coords)
            plt.plot(xs, ys, 'b-', linewidth=1.5, alpha=0.8)
            
            # Label nodes (only once per node to avoid duplicates)
            for n in node_names:
                if n in node_dict and n not in plotted_nodes:
                    x, y = node_dict[n][:2]
                    plt.text(x, y, n, fontsize=6, color='green', fontweight='bold')
                    plotted_nodes.add(n)
            
            # Label element at centroid (with larger font size)
            if len(coords) > 1:
                centroid_x = sum(x for x, _ in coords[:-1]) / len(coords[:-1])
                centroid_y = sum(y for _, y in coords[:-1]) / len(coords[:-1])
                plt.text(centroid_x, centroid_y, element_name.split('.')[-1], 
                        fontsize=6, color='red', ha='center', va='center', fontweight='bold',
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))
            
            plotted_count += 1
        
        plt.axis("equal")
        plt.grid(True, alpha=0.3)
        plt.title(f"All Meshes at z = {z} ({plotted_count} elements)", 
                 fontsize=16, fontweight='bold')
        plt.xlabel("X Coordinate", fontsize=14)
        plt.ylabel("Y Coordinate", fontsize=14)
        
        # Save single combined image for this z-value
        z_str = f"{z:.3f}".rstrip('0').rstrip('.')
        image_path = os.path.join(IMAGE_FOLDER, f"combined_all_meshes_z_{z_str}.png")
        
        plt.savefig(image_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"✅ Saved combined z={z} image with {plotted_count} elements: {image_path}")
    
    print(f"🎉 Generated {len(z_groups)} combined images for z-levels: {sorted(z_groups.keys())}")


def plot_mesh_data_y_axis(json_folder, IMAGE_FOLDER):
    """
    Plot all meshes from all location files with the same y-coordinate in a single image.
    Creates one image per unique y-value across all locations (for vertical walls).
    """
    import os
    import json
    import matplotlib.pyplot as plt
    from collections import defaultdict
    import glob
    
    # Paths
    wall_dir = os.path.join(json_folder, 'wall')
    nodes_data_path = os.path.join(json_folder, 'node_input_data.json')
    
    # Load nodes data once
    with open(nodes_data_path, 'r') as f:
        nodes_data = json.load(f)
    
    # Create dictionary with x, y, z coordinates
    node_dict = {n["name"]: (n["x"], n["y"], n.get("z", 0)) for n in nodes_data}
    
    # Find all JSON files in wall directory
    json_files = glob.glob(os.path.join(wall_dir, '*.json'))

    
    # Group ALL elements by y-coordinate across all files
    y_groups = defaultdict(list)
    total_elements = 0
    skipped_elements = []
    
    for json_file in json_files:
        location_name = os.path.splitext(os.path.basename(json_file))[0]
        
        with open(json_file, 'r') as f:
            mesh_data = json.load(f)

        elements = mesh_data.get("elements", {})
        
        for element_name, element_data in elements.items():
            node_names = element_data.get("nodes", [])
            if not node_names:
                skipped_elements.append(f"{location_name}.{element_name} (no nodes)")
                continue
            
            # Get y-coordinate from first node
            first_node = node_names[0]
            if first_node not in node_dict:
                skipped_elements.append(f"{location_name}.{element_name} (missing node: {first_node})")
                continue
            
            y = node_dict[first_node][1]
            
            # Verify all nodes exist and have same y-coordinate
            missing_nodes = [n for n in node_names if n not in node_dict]
            if missing_nodes:
                skipped_elements.append(f"{location_name}.{element_name} (missing nodes: {missing_nodes})")
                continue
                
            y_coords = [node_dict[n][1] for n in node_names]
            if not all(y_coord == y for y_coord in y_coords):
                skipped_elements.append(f"{location_name}.{element_name} (mixed y-coordinates)")
                continue
            
            # Add to y-group with location info
            y_groups[y].append((f"{location_name}.{element_name}", element_data, location_name))
            total_elements += 1
    
    # Report skipped elements
    if skipped_elements:
        print(f"⚠️  Skipped {len(skipped_elements)} elements:")
        for skip_msg in skipped_elements[:5]:  # Show first 5
            print(f"   - {skip_msg}")
        if len(skipped_elements) > 5:
            print(f"   ... and {len(skipped_elements) - 5} more")
    
    if not y_groups:
        print("❌ No valid elements found to plot")
        return
    
    print(f"📊 Total valid elements: {total_elements}")
    print(f"🎯 Found {len(y_groups)} unique y-levels: {sorted(y_groups.keys())}")
    
    # Create output directory
    os.makedirs(IMAGE_FOLDER, exist_ok=True)
    
    # Create one combined image per unique y-value
    for y in sorted(y_groups.keys()):
        elements_group = y_groups[y]
        
        plt.figure(figsize=(16, 12))  # Larger for combining multiple locations
        
        print(f"🎨 Creating combined image for y = {y} with {len(elements_group)} elements")
        
        # Use single color for all elements
        plotted_nodes = set()  # Track plotted nodes to avoid duplicates
        
        plotted_count = 0
        for element_name, element_data, location_name in elements_group:
            node_names = element_data.get("nodes", [])
            
            try:
                # For vertical walls, plot x-z coordinates (x as horizontal, z as vertical)
                coords = [(node_dict[n][0], node_dict[n][2]) for n in node_names]
            except KeyError:
                continue
            
            if len(coords) < 2:
                continue
            
            # Close the polygon if needed
            if coords[0] != coords[-1]:
                coords.append(coords[0])
            
            xs, zs = zip(*coords)
            plt.plot(xs, zs, 'b-', linewidth=1.5, alpha=0.8)
            
            # Label nodes (only once per node to avoid duplicates)
            for n in node_names:
                if n in node_dict and n not in plotted_nodes:
                    x, z = node_dict[n][0], node_dict[n][2]
                    plt.text(x, z, n, fontsize=10, color='green', fontweight='bold')
                    plotted_nodes.add(n)
            
            # Label element at centroid (with larger font size)
            if len(coords) > 1:
                centroid_x = sum(x for x, _ in coords[:-1]) / len(coords[:-1])
                centroid_z = sum(z for _, z in coords[:-1]) / len(coords[:-1])
                plt.text(centroid_x, centroid_z, element_name.split('.')[-1], 
                        fontsize=10, color='red', ha='center', va='center', fontweight='bold',
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))
            
            plotted_count += 1
        
        plt.axis("equal")
        plt.grid(True, alpha=0.3)
        plt.title(f"All Vertical Wall Meshes at y = {y} ({plotted_count} elements)", 
                 fontsize=16, fontweight='bold')
        plt.xlabel("X Coordinate", fontsize=14)
        plt.ylabel("Z Coordinate", fontsize=14)
        
        # Save single combined image for this y-value
        y_str = f"{y:.3f}".rstrip('0').rstrip('.')
        image_path = os.path.join(IMAGE_FOLDER, f"combined_vertical_wall_meshes_y_{y_str}.png")
        
        plt.savefig(image_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"✅ Saved combined y={y} image with {plotted_count} elements: {image_path}")
    
    print(f"🎉 Generated {len(y_groups)} combined images for y-levels: {sorted(y_groups.keys())}")

def plot_mesh_data_x_axis(json_folder, IMAGE_FOLDER):
    """
    Plot all meshes from all location files with the same x-coordinate in a single image.
    Creates one image per unique x-value across all locations (for vertical walls).
    """
    import os
    import json
    import matplotlib.pyplot as plt
    from collections import defaultdict
    import glob
    
    # Paths
    wall_dir = os.path.join(json_folder, 'wall')
    nodes_data_path = os.path.join(json_folder, 'node_input_data.json')
    
    # Load nodes data once
    with open(nodes_data_path, 'r') as f:
        nodes_data = json.load(f)
    
    # Create dictionary with x, y, z coordinates
    node_dict = {n["name"]: (n["x"], n["y"], n.get("z", 0)) for n in nodes_data}
    
    # Find all JSON files in wall directory
    json_files = glob.glob(os.path.join(wall_dir, '*.json'))

    
    # Group ALL elements by x-coordinate across all files
    x_groups = defaultdict(list)
    total_elements = 0
    skipped_elements = []
    
    for json_file in json_files:
        location_name = os.path.splitext(os.path.basename(json_file))[0]
        
        with open(json_file, 'r') as f:
            mesh_data = json.load(f)

        elements = mesh_data.get("elements", {})
        
        for element_name, element_data in elements.items():
            node_names = element_data.get("nodes", [])
            if not node_names:
                skipped_elements.append(f"{location_name}.{element_name} (no nodes)")
                continue
            
            # Get x-coordinate from first node
            first_node = node_names[0]
            if first_node not in node_dict:
                skipped_elements.append(f"{location_name}.{element_name} (missing node: {first_node})")
                continue
            
            x = node_dict[first_node][0]
            
            # Verify all nodes exist and have same x-coordinate
            missing_nodes = [n for n in node_names if n not in node_dict]
            if missing_nodes:
                skipped_elements.append(f"{location_name}.{element_name} (missing nodes: {missing_nodes})")
                continue
                
            x_coords = [node_dict[n][0] for n in node_names]
            if not all(x_coord == x for x_coord in x_coords):
                skipped_elements.append(f"{location_name}.{element_name} (mixed x-coordinates)")
                continue
            
            # Add to x-group with location info
            x_groups[x].append((f"{location_name}.{element_name}", element_data, location_name))
            total_elements += 1
    
    # Report skipped elements
    if skipped_elements:
        print(f"⚠️  Skipped {len(skipped_elements)} elements:")
        for skip_msg in skipped_elements[:5]:  # Show first 5
            print(f"   - {skip_msg}")
        if len(skipped_elements) > 5:
            print(f"   ... and {len(skipped_elements) - 5} more")
    
    if not x_groups:
        print("❌ No valid elements found to plot")
        return
    
    print(f"📊 Total valid elements: {total_elements}")
    print(f"🎯 Found {len(x_groups)} unique x-levels: {sorted(x_groups.keys())}")
    
    # Create output directory
    os.makedirs(IMAGE_FOLDER, exist_ok=True)
    
    # Create one combined image per unique x-value
    for x in sorted(x_groups.keys()):
        elements_group = x_groups[x]
        
        plt.figure(figsize=(16, 12))  # Larger for combining multiple locations
        
        print(f"🎨 Creating combined image for x = {x} with {len(elements_group)} elements")
        
        # Use single color for all elements
        plotted_nodes = set()  # Track plotted nodes to avoid duplicates
        
        plotted_count = 0
        for element_name, element_data, location_name in elements_group:
            node_names = element_data.get("nodes", [])
            
            try:
                # For vertical walls, plot y-z coordinates (y as horizontal, z as vertical)
                coords = [(node_dict[n][1], node_dict[n][2]) for n in node_names]
            except KeyError:
                continue
            
            if len(coords) < 2:
                continue
            
            # Close the polygon if needed
            if coords[0] != coords[-1]:
                coords.append(coords[0])
            
            ys, zs = zip(*coords)
            plt.plot(ys, zs, 'b-', linewidth=1.5, alpha=0.8)
            
            # Label nodes (only once per node to avoid duplicates)
            for n in node_names:
                if n in node_dict and n not in plotted_nodes:
                    y, z = node_dict[n][1], node_dict[n][2]
                    plt.text(y, z, n, fontsize=10, color='green', fontweight='bold')
                    plotted_nodes.add(n)
            
            # Label element at centroid (with larger font size)
            if len(coords) > 1:
                centroid_y = sum(y for y, _ in coords[:-1]) / len(coords[:-1])
                centroid_z = sum(z for _, z in coords[:-1]) / len(coords[:-1])
                plt.text(centroid_y, centroid_z, element_name.split('.')[-1], 
                        fontsize=10, color='red', ha='center', va='center', fontweight='bold',
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))
            
            plotted_count += 1
        
        plt.axis("equal")
        plt.grid(True, alpha=0.3)
        plt.title(f"All Vertical Wall Meshes at x = {x} ({plotted_count} elements)", 
                 fontsize=16, fontweight='bold')
        plt.xlabel("Y Coordinate", fontsize=14)
        plt.ylabel("Z Coordinate", fontsize=14)
        
        # Save single combined image for this x-value
        x_str = f"{x:.3f}".rstrip('0').rstrip('.')
        image_path = os.path.join(IMAGE_FOLDER, f"combined_vertical_wall_meshes_x_{x_str}.png")
        
        plt.savefig(image_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"✅ Saved combined x={x} image with {plotted_count} elements: {image_path}")
    
    print(f"🎉 Generated {len(x_groups)} combined images for x-levels: {sorted(x_groups.keys())}")

def plot_mesh_data(json_folder, IMAGE_FOLDER, z_spacing):
    # plot_mesh_data_separately(json_folder, location_name, IMAGE_FOLDER)
    plot_mesh_data_z_axis(json_folder, IMAGE_FOLDER)
    # plot_mesh_data_y_axis(json_folder, IMAGE_FOLDER)
    # plot_mesh_data_x_axis(json_folder, IMAGE_FOLDER)

def create_proper_mesh_for_closed_area_3d1(points, predefined_points, shell_section_name, JSON_FOLDER, IMAGE_FOLDER, thickness, load_case_names, pressures, num_x_div, num_y_div, numbering, add_shell, remove_shell, location_name):

    # Calculate ID offsets based on numbering parameter
    node_id_offset = 10000 + (numbering - 1) * 1000
    element_id_offset = 10000 + (numbering - 1) * 1000
    
    # Calculate the plane equation ax + by + cz + d = 0
    p0, p1, p2 = np.array(points[0]), np.array(points[1]), np.array(points[2])
    v1 = p1 - p0
    v2 = p2 - p0
    normal = np.cross(v1, v2)
    a, b, c = normal
    d = -np.dot(normal, p0)
    
    # Find two orthogonal vectors in the plane (basis vectors)
    if abs(a) > 0.1 or abs(b) > 0.1:
        u = np.array([b, -a, 0])  # Orthogonal to normal in XY plane
    else:
        u = np.array([0, c, -b])  # Orthogonal to normal in YZ plane
    u = u / np.linalg.norm(u)
    v = np.cross(normal, u)
    v = v / np.linalg.norm(v)
    
    # Function to project 3D points to 2D plane coordinates
    def project_to_plane(points_3d):
        return [(np.dot(p - p0, u), np.dot(p - p0, v)) for p in points_3d]
    
    # Project original points to 2D plane coordinates
    points_2d = project_to_plane(points)
    main_poly = ShapelyPolygon(points_2d)
    
    # Get bounding box of the polygon in plane coordinates
    min_x, min_y, max_x, max_y = main_poly.bounds
    
    # Calculate step sizes
    x_step = (max_x - min_x) / num_x_div
    y_step = (max_y - min_y) / num_y_div
    
    # Create dictionaries to store mesh and node information
    mesh_elements = {}
    node_positions = {}  # Stores {internal_id: (actual_node_id, coordinates)}
    node_counter = 1
    mesh_counter = 1
    
    # Function to add or find node
    def add_or_find_node(point_3d):
        nonlocal node_counter
        # Check if this 3D point already exists
        for internal_id, (existing_node_id, existing_point) in node_positions.items():
            if np.linalg.norm(point_3d - existing_point) < 1e-6:
                return existing_node_id, existing_point, False
        
        # Node doesn't exist, create new one
        node_id = node_counter + node_id_offset
        node_positions[node_counter] = (node_id, point_3d)
        node_counter += 1
        return node_id, point_3d, True
    
    # Function to process clipped polygon
    def process_polygon(poly, poly_type='clipped'):
        nonlocal mesh_counter
        
        if not isinstance(poly, ShapelyPolygon):
            return
            
        ext_coords = list(poly.exterior.coords)
        
        if len(ext_coords) < 3:  # Need at least 3 points
            return
            
        # Convert back to 3D coordinates
        node_indices = []
        coords_3d = []
        for coord in ext_coords[:-1]:  # Skip last point (same as first)
            x_proj, y_proj = coord
            point_3d = p0 + x_proj * u + y_proj * v
            node_id, coord_3d, _ = add_or_find_node(point_3d)
            node_indices.append(node_id)
            coords_3d.append(coord_3d)
        
        # Handle simple polygons (triangles/rectangles) vs complex ones
        if len(node_indices) <= 4:
            element_id = mesh_counter + element_id_offset
            mesh_name = f"R{element_id}" if len(node_indices) == 4 else f"T{element_id}"
            elem_type = 'rectangle' if len(node_indices) == 4 else 'triangle'
            mesh_elements[mesh_name] = {
                'type': elem_type,
                'nodes': node_indices,
                'coordinates': coords_3d,
                'id': element_id
            }
            mesh_counter += 1
        else:
            # Triangulate complex polygon
            triangles = triangulate(poly)
            for triangle in triangles:
                tri_coords = list(triangle.exterior.coords)
                tri_node_indices = []
                tri_coords_3d = []
                
                for coord in tri_coords[:-1]:
                    x_proj, y_proj = coord
                    point_3d = p0 + x_proj * u + y_proj * v
                    node_id, coord_3d, _ = add_or_find_node(point_3d)
                    tri_node_indices.append(node_id)
                    tri_coords_3d.append(coord_3d)
                
                element_id = mesh_counter + element_id_offset
                mesh_name = f"T{element_id}"
                mesh_elements[mesh_name] = {
                    'type': 'triangle',
                    'nodes': tri_node_indices,
                    'coordinates': tri_coords_3d,
                    'id': element_id
                }
                mesh_counter += 1
    
    # First pass: create mesh from grid
    for i in range(num_x_div):
        for j in range(num_y_div):
            x1 = min_x + i * x_step
            x2 = x1 + x_step
            y1 = min_y + j * y_step
            y2 = y1 + y_step
            
            # Create rectangle in plane coordinates and clip it
            rect = ShapelyPolygon([(x1, y1), (x2, y1), (x2, y2), (x1, y2)])
            clipped = rect.intersection(main_poly)
            
            if clipped.is_empty:
                continue
                
            if isinstance(clipped, MultiPolygon):
                for poly in clipped.geoms:
                    process_polygon(poly)
            else:
                process_polygon(clipped)
    
    # Second pass: triangulate remaining areas
    covered_area = ShapelyPolygon()
    for mesh in mesh_elements.values():
        projected = project_to_plane(mesh['coordinates'])
        covered_area = covered_area.union(ShapelyPolygon(projected))
    
    remaining_area = main_poly.difference(covered_area)
    
    if not remaining_area.is_empty:
        if isinstance(remaining_area, MultiPolygon):
            for poly in remaining_area.geoms:
                process_polygon(poly, 'remaining')
        else:
            process_polygon(remaining_area, 'remaining')
    
    # Create a mapping of all node IDs to their coordinates
    node_names = {node_id: coord for internal_id, (node_id, coord) in node_positions.items()}
    
    # Find closest mesh points to predefined points
    closest_mesh_points = {}
    
    for name, predefined_point in predefined_points.items():
        min_dist = float('inf')
        closest_node = None
        for node_id, mesh_point in node_names.items():
            dist = np.linalg.norm(np.array(predefined_point) - np.array(mesh_point))
            if dist < min_dist:
                min_dist = dist
                closest_node = node_id
        if closest_node is not None:
            closest_mesh_points[name] = closest_node
    
    # Replace closest mesh points with predefined points
    replaced_nodes = set()
    for p_name, node_id in closest_mesh_points.items():
        predefined_point = predefined_points[p_name]
        replaced_nodes.add(node_id)
        # Update node_names
        node_names[node_id] = predefined_point
        # Update coordinates in mesh elements
        for elem in mesh_elements.values():
            for i, n in enumerate(elem['nodes']):
                if n == node_id:
                    elem['coordinates'][i] = predefined_point
    
    # Add new shells if provided
    if add_shell:
        mesh_elements = add_new_shells(mesh_elements, node_names, add_shell)
    
    # Remove shells if provided
    if remove_shell:
        mesh_elements = remove_shells(mesh_elements, remove_shell)
        
    # Clean up unused nodes (except predefined points)
    used_nodes = set()
    for element in mesh_elements.values():
        used_nodes.update(element['nodes'])
    
    # Create a list of predefined point node IDs
    predefined_node_ids = set(closest_mesh_points.values())
    
    # Remove unused nodes (except those that are predefined points)
    unused_nodes = set(node_names.keys()) - used_nodes
    for node_id in unused_nodes:
        if node_id not in predefined_node_ids:
            del node_names[node_id]
            print(f"Removed unused node N{node_id}")
    


    

    
    # Sort nodes anti-clockwise for all elements before JSON creation
    for elem_name, elem_data in mesh_elements.items():
        nodes = elem_data['nodes']
        coords = elem_data['coordinates']
        
        if len(nodes) >= 3:
            sorted_nodes, sorted_coords = sort_nodes_anticlockwise(nodes, coords)
            elem_data['nodes'] = sorted_nodes
            elem_data['coordinates'] = sorted_coords

    # Create the output structure
    output_data = {
        "elements": {},
    }
    
    # Create a mapping of replaced node IDs to predefined point names
    replaced_nodes_mapping = {node_id: p_name for p_name, node_id in closest_mesh_points.items()}

    nodes_dir = os.path.join(JSON_FOLDER, 'node_input_data.json')

    with open(nodes_dir, 'r') as f:
        node_data = json.load(f)
        # print(f'node_data={node_data}')


    # Add all elements in the requested format with IDs
    for elem_name, elem_data in mesh_elements.items():
        node_names_in_elem = []
        for node_id in elem_data['nodes']:
            
            original_coord = node_names.get(node_id)

            matched = False
            if original_coord is not None:
                for nd in node_data:
                    nd_coords = [nd["x"], nd["y"], nd["z"]]
                    if all(abs(float(a) - float(b)) < 1e-6 for a, b in zip(nd_coords, original_coord)):
                        node_name = nd["name"]
                        node_id_new = nd["id"]
                        node_names_in_elem.append(node_name)
                        matched = True
                        break

            if not matched:
                node_names_in_elem.append(replaced_nodes_mapping.get(node_id, f"n{node_id}"))


        output_data["elements"][elem_name] = {
            "id": elem_data['id'],
            "type": elem_data['type'],
            "nodes": node_names_in_elem,
            "shell_section": shell_section_name,
            "thickness": thickness, "load_case_names": load_case_names, "pressures": pressures,
        }
        
    # Ensure directory exists and save JSON
    wall_dir = os.path.join(JSON_FOLDER, 'wall')
    os.makedirs(wall_dir, exist_ok=True)
    # full_path = os.path.join(wall_dir, f'mesh_data_with_predefined_points{numbering}.json')
    full_path = os.path.join(wall_dir, f'{location_name}.json')
    
    with open(full_path, 'w') as f:
        json.dump(output_data, f, indent=4)

    # Load existing node data
    if os.path.exists(nodes_dir):
        with open(nodes_dir, 'r') as f:
            existing_nodes = json.load(f)
    else:
        existing_nodes = []

    # Create a set of existing coordinate keys
    existing_keys = {
        tuple([round(float(nd["x"]), 6), round(float(nd["y"]), 6), round(float(nd["z"]), 6)])

        for nd in existing_nodes
    }

    # Prepare new nodes from element data
    for elem_data in mesh_elements.values():
        for node_id in elem_data['nodes']:
            coord = node_names.get(node_id)
            if coord is not None:
                key = tuple(round(float(c), 6) for c in coord)
                if key not in existing_keys:
                    existing_nodes.append({
                        "id": int(node_id),
                        "name": f"n{node_id}",
                        "x": float(coord[0]),
                        "y": float(coord[1]),
                        "z": float(coord[2])
                    })
                    existing_keys.add(key)

    # Write back to JSON
    with open(nodes_dir, 'w') as f:
        json.dump(existing_nodes, f, indent=4)
    
    

    return output_data



from collections import defaultdict
import json
import os

def filter_mesh_by_direction(mesh_data, nodes_data, direction='z', tolerance=1e-6):
    """
    Filter and group mesh elements by their position along a specified direction (X, Y, or Z).
    For Z direction: All nodes must have the same Z coordinate (horizontal elements like slabs)
    For Y direction: All nodes must have the same Y coordinate (vertical elements along X-Z plane)
    For X direction: All nodes must have the same X coordinate (vertical elements along Y-Z plane)
    
    Parameters:
    - mesh_data: Dictionary containing mesh elements data
    - nodes_data: List of node dictionaries with id, name, x, y, z coordinates
    - direction: 'x', 'y', or 'z' - the direction along which to filter
    - tolerance: Tolerance for comparing coordinate values
    
    Returns:
    - Dictionary grouped by coordinate levels: {level: {mesh_name: mesh_info}}
    """
    
    # Create node lookup dictionary for quick access
    node_lookup = {node['name']: (node['x'], node['y'], node['z']) for node in nodes_data}
    node_lookup.update({f"n{node['id']}": (node['x'], node['y'], node['z']) for node in nodes_data})
    
    # Direction mapping
    direction_index = {'x': 0, 'y': 1, 'z': 2}
    if direction.lower() not in direction_index:
        raise ValueError("Direction must be 'x', 'y', or 'z'")
    
    coord_index = direction_index[direction.lower()]
    
    # Group mesh elements by their coordinate level
    grouped_mesh = defaultdict(dict)
    
    # Get all mesh elements
    mesh_elements = mesh_data.get('elements', {})
    
    for mesh_name, mesh_info in mesh_elements.items():
        nodes_in_mesh = mesh_info.get('nodes', [])
        
        if not nodes_in_mesh:
            continue
            
        # Get coordinates for all nodes in this mesh element
        mesh_coordinates = []
        valid_mesh = True
        
        for node_name in nodes_in_mesh:
            if node_name in node_lookup:
                mesh_coordinates.append(node_lookup[node_name])
            else:
                print(f"Warning: Node {node_name} not found in node_lookup for mesh {mesh_name}")
                valid_mesh = False
                break
        
        if not valid_mesh or len(mesh_coordinates) == 0:
            continue
        
        # Check if ALL nodes have the same coordinate in the specified direction
        first_coord = mesh_coordinates[0][coord_index]
        all_same_coord = True
        
        for coord in mesh_coordinates[1:]:
            if abs(coord[coord_index] - first_coord) > tolerance:
                all_same_coord = False
                break
        
        # Only include mesh elements where all nodes are at the same level
        if all_same_coord:
            # Find which level this mesh belongs to by checking existing levels
            assigned_level = None
            for existing_level in grouped_mesh.keys():
                if abs(first_coord - existing_level) < tolerance:
                    assigned_level = existing_level
                    break
            
            # If no existing level found, create new level
            if assigned_level is None:
                assigned_level = round(first_coord, 6)  # Round to avoid floating point precision issues
            
            # Store mesh information with node coordinates
            grouped_mesh[assigned_level][mesh_name] = {
                'id': mesh_info.get('id'),
                'type': mesh_info.get('type'),
                'nodes': nodes_in_mesh,
                'node_coordinates': mesh_coordinates,
                'shell_section': mesh_info.get('shell_section'),
                'thickness': mesh_info.get('thickness'),
                'load_case_names': mesh_info.get('load_case_names'),
                'pressures': mesh_info.get('pressures'),
                f'{direction.lower()}_coord': first_coord
            }
    
    return dict(grouped_mesh)

def filter_mesh_by_z_levels(mesh_data, nodes_data, z_points, tolerance=1e-6):
    """
    Filter and group mesh elements by specific Z levels.
    Only includes mesh elements where ALL nodes have the same Z coordinate.
    
    Parameters:
    - mesh_data: Dictionary containing mesh elements data
    - nodes_data: List of node dictionaries
    - z_points: List of specific Z coordinates to filter by
    - tolerance: Tolerance for comparing Z values
    
    Returns:
    - Dictionary grouped by Z levels: {z_level: {mesh_name: mesh_info}}
    """
    
    # Create node lookup dictionary
    node_lookup = {node['name']: (node['x'], node['y'], node['z']) for node in nodes_data}
    node_lookup.update({f"n{node['id']}": (node['x'], node['y'], node['z']) for node in nodes_data})
    
    grouped_mesh = defaultdict(dict)
    
    # Get all mesh elements
    mesh_elements = mesh_data.get('elements', {})
    
    for mesh_name, mesh_info in mesh_elements.items():
        nodes_in_mesh = mesh_info.get('nodes', [])
        
        if not nodes_in_mesh:
            continue
            
        # Get coordinates for all nodes in this mesh element
        mesh_coordinates = []
        valid_mesh = True
        
        for node_name in nodes_in_mesh:
            if node_name in node_lookup:
                mesh_coordinates.append(node_lookup[node_name])
            else:
                print(f"Warning: Node {node_name} not found in node_lookup for mesh {mesh_name}")
                valid_mesh = False
                break
        
        if not valid_mesh or len(mesh_coordinates) == 0:
            continue
        
        # Check if ALL nodes have the same Z coordinate
        first_z = mesh_coordinates[0][2]
        all_same_z = True
        
        for coord in mesh_coordinates[1:]:
            if abs(coord[2] - first_z) > tolerance:
                all_same_z = False
                break
        
        # Only process if all nodes are at the same Z level
        if all_same_z:
            # Check which Z level this mesh belongs to
            for z_level in z_points:
                if abs(first_z - z_level) < tolerance:
                    grouped_mesh[z_level][mesh_name] = {
                        'id': mesh_info.get('id'),
                        'type': mesh_info.get('type'),
                        'nodes': nodes_in_mesh,
                        'node_coordinates': mesh_coordinates,
                        'shell_section': mesh_info.get('shell_section'),
                        'thickness': mesh_info.get('thickness'),
                        'load_case_names': mesh_info.get('load_case_names'),
                        'pressures': mesh_info.get('pressures'),
                        'z_coord': first_z
                    }
                    break
    
    return dict(grouped_mesh)

def print_mesh_summary(grouped_mesh, direction='z'):
    """
    Print a summary of the grouped mesh elements.
    
    Parameters:
    - grouped_mesh: Dictionary returned by filter_mesh_by_direction or filter_mesh_by_z_levels
    - direction: Direction used for filtering ('x', 'y', or 'z')
    """
    
    print(f"\n=== Mesh Elements Grouped by {direction.upper()} Direction ===")
    
    for level in sorted(grouped_mesh.keys()):
        mesh_dict = grouped_mesh[level]
        # print(f"\n{direction.upper()} Level: {level}")
        # print(f"Number of mesh elements: {len(mesh_dict)}")
        
        # Show first few mesh elements as examples
        mesh_names = list(mesh_dict.keys())[:3]
        for mesh_name in mesh_names:
            mesh_info = mesh_dict[mesh_name]
            # print(f"  {mesh_name}: {mesh_info['type']}, nodes: {mesh_info['nodes']}")
        
        if len(mesh_dict) > 3:
            print(f"  ... and {len(mesh_dict) - 3} more elements")

# Example usage function
def filter_mesh(JSON_FOLDER, z_spacing):
    """
    Example of how to use the mesh filtering functions and save results to JSON files.
    Based on spacing: x_spacing = [10.0, 12.0, 11.0] * foot
                      y_spacing = [10.0, 12.0, 11.0] * foot  
                      z_spacing = [10.0, 12.0, 11.0] * foot
    """
    
    # Load the combined structure data (assuming create_combined_structure_json exists)
    nodes_data, members_data, mesh_data = create_combined_structure_json(JSON_FOLDER)
    
    # Calculate expected Z levels (floors) dynamically
    z_levels = [0.0]
    for spacing in z_spacing:
        z_levels.append(z_levels[-1] + spacing)
    
    # Filter mesh by Z direction (horizontal elements like slabs at each floor)
    mesh_by_z = filter_mesh_by_direction(mesh_data, nodes_data, direction='z')
    print_mesh_summary(mesh_by_z, 'z')
    print(f"Expected Z levels: {z_levels}")
    print(f"Found Z levels: {sorted(mesh_by_z.keys())}")
    
    # Filter mesh by Y direction (vertical elements with same Y coordinate - walls along X-Z plane)
    mesh_by_y = filter_mesh_by_direction(mesh_data, nodes_data, direction='y')
    print_mesh_summary(mesh_by_y, 'y')
    
    # Filter mesh by X direction (vertical elements with same X coordinate - walls along Y-Z plane)
    mesh_by_x = filter_mesh_by_direction(mesh_data, nodes_data, direction='x')
    print_mesh_summary(mesh_by_x, 'x')
    

    
    # Create filtered_mesh directory inside wall directory
    wall_dir = os.path.join(JSON_FOLDER, 'wall')
    filtered_mesh_dir = os.path.join(wall_dir, 'filtered_mesh')
    os.makedirs(filtered_mesh_dir, exist_ok=True)
    
    # Save filtered mesh data to JSON files
    # Save mesh filtered by Z direction
    mesh_z_path = os.path.join(filtered_mesh_dir, 'mesh_filtered_by_z.json')
    with open(mesh_z_path, 'w') as f:
        json.dump(mesh_by_z, f, indent=4)
    print(f"Mesh filtered by Z direction saved to: {mesh_z_path}")
    
    # Save mesh filtered by Y direction
    mesh_y_path = os.path.join(filtered_mesh_dir, 'mesh_filtered_by_y.json')
    with open(mesh_y_path, 'w') as f:
        json.dump(mesh_by_y, f, indent=4)
    print(f"Mesh filtered by Y direction saved to: {mesh_y_path}")
    
    # Save mesh filtered by X direction
    mesh_x_path = os.path.join(filtered_mesh_dir, 'mesh_filtered_by_x.json')
    with open(mesh_x_path, 'w') as f:
        json.dump(mesh_by_x, f, indent=4)
    print(f"Mesh filtered by X direction saved to: {mesh_x_path}")
    
    
    # Save a summary file with all filtered data
    summary_data = {
        'summary': {
            'total_mesh_elements': len(mesh_data.get('elements', {})),
            'expected_z_levels': z_levels,
            'found_z_levels': sorted(mesh_by_z.keys()),
            'z_level_counts': {str(k): len(v) for k, v in mesh_by_z.items()},
            'y_level_counts': {str(k): len(v) for k, v in mesh_by_y.items()},
            'x_level_counts': {str(k): len(v) for k, v in mesh_by_x.items()},
            'filtering_criteria': {
                'z_direction': 'All nodes must have same Z coordinate (horizontal elements)',
                'y_direction': 'All nodes must have same Y coordinate (vertical walls along X-Z)',
                'x_direction': 'All nodes must have same X coordinate (vertical walls along Y-Z)'
            }
        },
    }
    
    summary_path = os.path.join(filtered_mesh_dir, 'mesh_filtering_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary_data, f, indent=4)
    print(f"Complete mesh filtering summary saved to: {summary_path}")
    
    return mesh_by_z, mesh_by_y, mesh_by_x



import math
import os
import json

def create_slabs(node_data, z_level):
    slabs = {}
    
    # Extract all nodes and sort them by y then x coordinates for proper ordering
    nodes = sorted(node_data.keys(), key=lambda n: (node_data[n][1], node_data[n][0]))
    
    # Group nodes by their y-coordinate to identify rows
    y_coords = sorted({node_data[n][1] for n in nodes})
    rows = {}
    for y in y_coords:
        rows[y] = sorted([n for n in nodes if node_data[n][1] == y], 
                         key=lambda n: node_data[n][0])
    
    # Create rectangular slabs (4-noded) where possible
    slab_count = 1
    
    for i in range(len(y_coords)-1):
        y1 = y_coords[i]
        y2 = y_coords[i+1]
        row1 = rows[y1]
        row2 = rows[y2]
        
        # Try to match nodes between rows to form rectangles
        j = 0
        while j < len(row1)-1 and j < len(row2)-1:
            n1 = row1[j]
            n2 = row1[j+1]
            n3 = row2[j+1]
            n4 = row2[j]
            
            # Check if these nodes form a proper quadrilateral
            if (node_data[n1][0] == node_data[n4][0] and  # x-coord matches vertically
                node_data[n2][0] == node_data[n3][0] and  # x-coord matches vertically
                node_data[n1][1] == node_data[n2][1] and  # y-coord matches horizontally
                node_data[n4][1] == node_data[n3][1]):   # y-coord matches horizontally
                
                # Create slab configuration
                slab_name = f"{z_level}_slab{slab_count}"
                slabs[slab_name] = {
                    "points": [n1, n2, n3, n4],  # Already ordered anticlockwise
                    "section_name": "floor_slab1",
                    "add_shell": {},
                    "remove_shell": [],
                    "predefined_points": {},
                    "num_x_div": 3,
                    "num_y_div": 3,
                    "load_case_names": ['DL', 'LL', 'self_weight'],
                    "pressures": [-0.02, -0.045, -0.075],  # in ksf
                    "thickness": 6  # in inches
                }
                slab_count += 1
                j += 1  # Move to next potential rectangle
            else:
                j += 1
    
    return slabs

def create_surface_configurations(JSON_FOLDER):
    """Main function to create surface configurations from node data"""
    try:
        # Load filtered nodes data directly
        filepath = os.path.join(JSON_FOLDER, "filtered_nodes.json")
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        
        with open(filepath, 'r') as f:
            filtered_nodes = json.load(f)
        
        surface_configurations = {}  # Dictionary to hold all surface configurations
        
        # Process each z-level
        for z_level, node_data in filtered_nodes.items():
            print(f"\nProcessing z-level: {z_level}")
            
            # Create slabs for this z-level
            slabs = create_slabs(node_data, z_level)
            
            # Add to surface configurations
            surface_configurations.update(slabs)
        
        # Save surface configurations to JSON file
        output_filename = "surface_configurations.json"
        output_filepath = os.path.join(JSON_FOLDER, output_filename)
        with open(output_filepath, 'w') as f:
            json.dump(surface_configurations, f, indent=4)
        print(f"\nSaved surface configurations to {output_filename}")
        
        return surface_configurations
    
    except Exception as e:
        print(f"Error in create_surface_configurations: {str(e)}")
        return None

def order_nodes_anticlockwise(nodes, node_data):
    """Order nodes in anticlockwise direction"""
    if len(nodes) != len(set(nodes)):  # Check for duplicate nodes
        return [n for n in nodes if nodes.count(n) == 1]  # Remove duplicates
    
    if len(nodes) == 3:
        # For triangles, calculate centroid
        x = sum(node_data[n][0] for n in nodes) / 3
        y = sum(node_data[n][1] for n in nodes) / 3
        
        # Calculate angles from centroid
        def angle_from_centroid(n):
            dx = node_data[n][0] - x
            dy = node_data[n][1] - y
            return math.atan2(dy, dx)
        
        return sorted(nodes, key=angle_from_centroid, reverse=True)
    elif len(nodes) == 4:
        # For quadrilaterals, order is already anticlockwise in our creation method
        return nodes
    return nodes

def create_slabs_nodes(JSON_FOLDER):
    """Main function to create slabs from node data"""
    try:
        # Load filtered nodes data directly
        filepath = os.path.join(JSON_FOLDER, "filtered_nodes.json")
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        
        with open(filepath, 'r') as f:
            filtered_nodes = json.load(f)
        
        all_slabs = {}  # Dictionary to hold slabs from all z-levels
        
        # Process each z-level
        for z_level, node_data in filtered_nodes.items():
            print(f"\nProcessing z-level: {z_level}")
            
            # Create slabs for this z-level
            slabs = create_slabs(node_data)
            
            # Store slabs with z-level prefix
            for slab_name, slab_data in slabs.items():
                # Add z-level information to each slab
                slab_data['z_level'] = z_level
                all_slabs[f"{z_level}_{slab_name}"] = slab_data
            
        # Save all slabs data to a single JSON file
        output_filename = "all_slabs.json"
        output_filepath = os.path.join(JSON_FOLDER, output_filename)
        with open(output_filepath, 'w') as f:
            json.dump(all_slabs, f, indent=4)
        print(f"\nSaved all slabs data to {output_filename}")
        
        return all_slabs
    
    except Exception as e:
        print(f"Error in create_slabs_nodes: {str(e)}")
        return None

# Example usage:
# create_slabs_nodes("path_to_your_json_folder")


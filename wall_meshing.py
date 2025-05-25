from import_ import *
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
            print(f"Added shell {shell_name} with ID {element_id}")
    
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


def create_proper_mesh_for_closed_area_3d1(points, predefined_points, shell_section_name, JSON_FOLDER, IMAGE_FOLDER, 
                                         num_x_div=4, num_y_div=4, numbering=1, add_shell=None, remove_shell=None):
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
    
    # 3D Plotting
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title('3D Mesh Elements with Nodes (excluding predefined points)')
    
    # Plot original shape
    original_verts = [points + [points[0]]]
    original_poly = Poly3DCollection(original_verts, alpha=0.3, 
                                   facecolors='red', linewidths=2, 
                                   edgecolors='red', linestyles='--')
    ax.add_collection3d(original_poly)
    
    # Plot mesh elements
    for name, data in mesh_elements.items():
        # Plot each element
        color = 'cyan' if data['type'] == 'rectangle' else 'lightgreen'
        alpha = 0.4 if data['type'] == 'rectangle' else 0.6
        
        verts = [data['coordinates'] + [data['coordinates'][0]]]
        poly = Poly3DCollection(verts, alpha=alpha, facecolors=color, 
                              edgecolors='blue', linewidths=1)
        ax.add_collection3d(poly)
        
        # Add element name at centroid
        centroid = np.mean(data['coordinates'], axis=0)
        ax.text(centroid[0], centroid[1], centroid[2], name, 
                ha='center', va='center', fontsize=8, weight='bold')
        
        # Plot nodes not replaced by predefined points
        for node_num in data['nodes']:
            if node_num not in replaced_nodes and node_num in node_names:
                coord = node_names[node_num]
                ax.scatter([coord[0]], [coord[1]], [coord[2]], c='red', s=50)
                ax.text(coord[0], coord[1], coord[2], f'N{node_num}', 
                        ha='right', va='bottom', fontsize=8, color='darkred')
    
    # Plot predefined points
    for p_name, p_coord in predefined_points.items():
        ax.scatter([p_coord[0]], [p_coord[1]], [p_coord[2]], c='blue', s=100, marker='*')
        ax.text(p_coord[0], p_coord[1], p_coord[2], p_name, 
                ha='left', va='top', fontsize=10, color='darkblue', weight='bold')
    
    # Set labels and legend
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    rect_patch = mpatches.Patch(color='cyan', label='Rectangular Elements')
    tri_patch = mpatches.Patch(color='lightgreen', label='Triangular Elements')
    node_patch = mpatches.Patch(color='red', label='Nodes (non-predefined)')
    predef_patch = mpatches.Patch(color='blue', label='Predefined Points')
    ax.legend(handles=[rect_patch, tri_patch, node_patch, predef_patch])
    
    plt.tight_layout()
    
    # Ensure the IMAGE_FOLDER exists and save plot
    os.makedirs(IMAGE_FOLDER, exist_ok=True)
    image_filename = f"mesh_plot_{numbering}.png"
    image_path = os.path.join(IMAGE_FOLDER, image_filename)
    plt.savefig(image_path, dpi=300, bbox_inches='tight')
    print(f"✅ Saved plot to: {image_path}")
    
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
        "nodes": {},
    }
    
    # Create a mapping of replaced node IDs to predefined point names
    replaced_nodes_mapping = {node_id: p_name for p_name, node_id in closest_mesh_points.items()}

    # Add nodes with their coordinates and IDs (only nodes that are either used or predefined)
    for node_id, coord in node_names.items():
        if node_id not in replaced_nodes_mapping:
            output_data["nodes"][f"N{node_id}"] = {
                "id": node_id,
                "coordinates": [float(coord[0]), float(coord[1]), float(coord[2])],
                "is_predefined": node_id in replaced_nodes
            }

    # Add all elements in the requested format with IDs
    for elem_name, elem_data in mesh_elements.items():
        # Replace node names with predefined point names where applicable
        node_names_in_elem = []
        for node_id in elem_data['nodes']:
            node_names_in_elem.append(replaced_nodes_mapping.get(node_id, f"N{node_id}"))
        
        output_data["elements"][elem_name] = {
            "id": elem_data['id'],
            "type": elem_data['type'],
            "nodes": node_names_in_elem,
            "shell_section": shell_section_name
        }
        
    # Ensure directory exists and save JSON
    wall_dir = os.path.join(JSON_FOLDER, 'wall')
    os.makedirs(wall_dir, exist_ok=True)
    full_path = os.path.join(wall_dir, f'mesh_data_with_predefined_points{numbering}.json')
    
    with open(full_path, 'w') as f:
        json.dump(output_data, f, indent=4)
        
    return output_data



def create_combined_structure_json(JSON_FOLDER):
    """
    Combines structure data and mesh data into a single JSON file.
    Also updates node_input_data.json with new shell nodes (if any).
    
    Returns: {'nodes': [], 'members': [], 'shell_elements': []}
    """
    # Initialize empty structure
    combined_data = {
        "nodes": [],
        "members": [],
        "shell_elements": []
    }

    # ==================================================================
    # 1. FIRST Load Mesh Data (if exists)
    # ==================================================================
    wall_dir = os.path.join(JSON_FOLDER, "wall")
    if os.path.exists(wall_dir):
        print(f"Found wall directory at: {wall_dir}")
        for filename in sorted(os.listdir(wall_dir)):
            if filename.startswith("mesh_data_with_predefined_points"):
                with open(os.path.join(wall_dir, filename), 'r') as f:
                    mesh_data = json.load(f)
                    print(f"Loaded mesh file: {filename} with {len(mesh_data.get('nodes', {}))} nodes")
                    
                    # Add mesh nodes to combined_data (in 'coordinates' format)
                    for node_name, node_info in mesh_data.get('nodes', {}).items():
                        combined_data["nodes"].append({
                            "id": node_info['id'],
                            "name": node_name,
                            "coordinates": node_info['coordinates']
                        })
                    
                    # Add shell elements to combined_data
                    for elem_name, elem_info in mesh_data.get('elements', {}).items():
                        combined_data["shell_elements"].append({
                            "id": elem_info['id'],
                            "name": elem_name,
                            "node_names": elem_info['nodes'],
                            "shell_section": elem_info.get('shell_section', 'default')
                        })
    else:
        print(f"No wall directory found at: {wall_dir}")

    # ==================================================================
    # 2. THEN Update node_input_data.json with New Shell Nodes (No Duplicates)
    # ==================================================================
    node_file = os.path.join(JSON_FOLDER, "node_input_data.json")
    
    # Load existing nodes (or initialize empty list)
    if os.path.exists(node_file):
        with open(node_file, 'r') as f:
            existing_nodes = json.load(f)
    else:
        existing_nodes = []
    
    existing_ids = {node["id"] for node in existing_nodes}
    updated_nodes = existing_nodes.copy()  # Start with existing nodes

    # Add new shell nodes (if they don't exist in node_input_data.json)
    for node in combined_data["nodes"]:
        if node["id"] not in existing_ids:
            updated_nodes.append({
                "id": node["id"],
                "name": node["name"],
                "x": node["coordinates"][0],
                "y": node["coordinates"][1],
                "z": node["coordinates"][2]
            })
            print(f"Added new node ID {node['id']} to node_input_data.json")

    # Save updated nodes back to node_input_data.json
    with open(node_file, 'w') as f:
        json.dump(updated_nodes, f, indent=4)
    print(f"Updated {node_file} with {len(updated_nodes)} nodes (original: {len(existing_nodes)})")


    return combined_data



# def add_new_shells(mesh_elements, node_names, add_shell):
#     """
#     Add new shells to the mesh based on predefined points
    
#     Parameters:
#     -----------
#     mesh_elements: dict
#         Dictionary containing the current mesh elements
#     node_names: dict
#         Dictionary mapping node IDs to their coordinates
#     add_shell: dict
#         Dictionary with shell names as keys and lists of point names as values
        
#     Returns:
#     --------
#     dict: Updated mesh elements with new shells added
#     """
#     # Create a copy of the input mesh_elements to avoid modifying the original
#     updated_mesh = mesh_elements.copy()
    
#     # Determine starting mesh ID
#     mesh_counter = 1 if not mesh_elements else max(elem_data['id'] for elem_data in mesh_elements.values()) + 1
    
#     # Process each shell to be added
#     for shell_name, points in add_shell.items():
#         # Determine element type based on number of points
#         if len(points) == 3:  # Triangle
#             element_id = mesh_counter
#             mesh_name = f"T{element_id}"
#             shell_type = 'triangle'
#         elif len(points) == 4:  # Quadrilateral
#             element_id = mesh_counter
#             mesh_name = f"R{element_id}"
#             shell_type = 'rectangle'
#         else:
#             print(f"Warning: Shell {shell_name} has {len(points)} points, only triangles (3) and quads (4) are supported")
#             continue
        
#         # Extract node IDs and coordinates
#         nodes = []
#         coords = []
        
#         # Find node IDs corresponding to predefined points
#         for point_name in points:
#             # Direct dictionary lookup instead of looping
#             if point_name in node_names:
#                 nodes.append(point_name)
#                 coords.append(node_names[point_name])
#             else:
#                 print(f"Warning: Point {point_name} not found in mesh nodes")
        
#         # Only create shell if we have all required points
#         if len(nodes) == len(points):
#             updated_mesh[shell_name] = {
#                 'type': shell_type,
#                 'nodes': nodes,
#                 'coordinates': coords,
#                 'id': element_id
#             }
#             mesh_counter += 1
#             print(f"Added shell {shell_name} with ID {element_id}")
    
#     return updated_mesh


# def remove_shells(mesh_elements, remove_shell):
#     """
#     Remove specified shells from the mesh
    
#     Parameters:
#     -----------
#     mesh_elements: dict
#         Dictionary containing the current mesh elements
#     remove_shell: list
#         List of shell names to be removed
        
#     Returns:
#     --------
#     dict: Updated mesh elements with specified shells removed
#     """
#     # Create a copy of the input mesh_elements to avoid modifying the original
#     final_mesh = mesh_elements.copy()
    
#     for shell_name in remove_shell:
#         if shell_name in final_mesh:
#             del final_mesh[shell_name]
#             print(f"Removed shell {shell_name}")
#         else:
#             print(f"Warning: Shell {shell_name} not found in mesh")
    
#     return final_mesh


# def sort_nodes_anticlockwise(nodes, coords):
#     """
#     Sort nodes in anti-clockwise order around their centroid.
    
#     Args:
#         nodes: List of node IDs
#         coords: List of corresponding 3D coordinates
        
#     Returns:
#         Tuple of (sorted_nodes, sorted_coords)
#     """
#     if len(nodes) <= 3:  # Triangles are already planar
#         return nodes, coords
        
#     # Calculate centroid
#     centroid = np.mean(coords, axis=0)
    
#     # Center the points
#     centered = np.array(coords) - centroid
    
#     # Compute normal vector using first 3 points
#     normal = np.cross(centered[1] - centered[0], centered[2] - centered[0])
#     normal = normal / np.linalg.norm(normal)
    
#     # Create basis vectors for projection
#     if abs(normal[0]) > 0.1 or abs(normal[1]) > 0.1:
#         u = np.array([normal[1], -normal[0], 0])  # Orthogonal to normal in XY plane
#     else:
#         u = np.array([0, normal[2], -normal[1]])  # Orthogonal to normal in YZ plane
#     u = u / np.linalg.norm(u)
#     v = np.cross(normal, u)
    
#     # Project points to 2D plane
#     projected = []
#     for point in centered:
#         x_proj = np.dot(point, u)
#         y_proj = np.dot(point, v)
#         projected.append((x_proj, y_proj))
    
#     # Calculate angles from centroid and sort
#     angles = np.arctan2([p[1] for p in projected], [p[0] for p in projected])
#     positive_angles = np.where(angles < 0, angles + 2*np.pi, angles)
#     sorted_indices = np.argsort(positive_angles)
    
#     # Return sorted nodes and coordinates
#     sorted_nodes = [nodes[i] for i in sorted_indices]
#     sorted_coords = [coords[i] for i in sorted_indices]
    
#     return sorted_nodes, sorted_coords


# def create_proper_mesh_for_closed_area_3d1(points, predefined_points, shell_section_name, JSON_FOLDER, IMAGE_FOLDER, 
#                                          num_x_div=4, num_y_div=4, numbering=1, add_shell=None, remove_shell=None):
#     # Calculate ID offsets based on numbering parameter
#     node_id_offset = 10000 + (numbering - 1) * 1000
#     element_id_offset = 10000 + (numbering - 1) * 1000
    
#     # Calculate the plane equation ax + by + cz + d = 0
#     p0, p1, p2 = np.array(points[0]), np.array(points[1]), np.array(points[2])
#     v1 = p1 - p0
#     v2 = p2 - p0
#     normal = np.cross(v1, v2)
#     a, b, c = normal
#     d = -np.dot(normal, p0)
    
#     # Find two orthogonal vectors in the plane (basis vectors)
#     if abs(a) > 0.1 or abs(b) > 0.1:
#         u = np.array([b, -a, 0])  # Orthogonal to normal in XY plane
#     else:
#         u = np.array([0, c, -b])  # Orthogonal to normal in YZ plane
#     u = u / np.linalg.norm(u)
#     v = np.cross(normal, u)
#     v = v / np.linalg.norm(v)
    
#     # Function to project 3D points to 2D plane coordinates
#     def project_to_plane(points_3d):
#         return [(np.dot(p - p0, u), np.dot(p - p0, v)) for p in points_3d]
    
#     # Project original points to 2D plane coordinates
#     points_2d = project_to_plane(points)
#     main_poly = ShapelyPolygon(points_2d)
    
#     # Get bounding box of the polygon in plane coordinates
#     min_x, min_y, max_x, max_y = main_poly.bounds
    
#     # Calculate step sizes
#     x_step = (max_x - min_x) / num_x_div
#     y_step = (max_y - min_y) / num_y_div
    
#     # Create dictionaries to store mesh and node information
#     mesh_elements = {}
#     node_positions = {}  # Stores {internal_id: (actual_node_id, coordinates)}
#     node_counter = 1
#     mesh_counter = 1
    
#     # Function to add or find node
#     def add_or_find_node(point_3d):
#         nonlocal node_counter
#         # Check if this 3D point already exists
#         for internal_id, (existing_node_id, existing_point) in node_positions.items():
#             if np.linalg.norm(point_3d - existing_point) < 1e-6:
#                 return existing_node_id, existing_point, False
        
#         # Node doesn't exist, create new one
#         node_id = node_counter + node_id_offset
#         node_positions[node_counter] = (node_id, point_3d)
#         node_counter += 1
#         return node_id, point_3d, True
    
#     # Function to process clipped polygon
#     def process_polygon(poly, poly_type='clipped'):
#         nonlocal mesh_counter
        
#         if not isinstance(poly, ShapelyPolygon):
#             return
            
#         ext_coords = list(poly.exterior.coords)
        
#         if len(ext_coords) < 3:  # Need at least 3 points
#             return
            
#         # Convert back to 3D coordinates
#         node_indices = []
#         coords_3d = []
#         for coord in ext_coords[:-1]:  # Skip last point (same as first)
#             x_proj, y_proj = coord
#             point_3d = p0 + x_proj * u + y_proj * v
#             node_id, coord_3d, _ = add_or_find_node(point_3d)
#             node_indices.append(node_id)
#             coords_3d.append(coord_3d)
        
#         # Handle simple polygons (triangles/rectangles) vs complex ones
#         if len(node_indices) <= 4:
#             element_id = mesh_counter + element_id_offset
#             mesh_name = f"R{element_id}" if len(node_indices) == 4 else f"T{element_id}"
#             elem_type = 'rectangle' if len(node_indices) == 4 else 'triangle'
#             mesh_elements[mesh_name] = {
#                 'type': elem_type,
#                 'nodes': node_indices,
#                 'coordinates': coords_3d,
#                 'id': element_id
#             }
#             mesh_counter += 1
#         else:
#             # Triangulate complex polygon
#             triangles = triangulate(poly)
#             for triangle in triangles:
#                 tri_coords = list(triangle.exterior.coords)
#                 tri_node_indices = []
#                 tri_coords_3d = []
                
#                 for coord in tri_coords[:-1]:
#                     x_proj, y_proj = coord
#                     point_3d = p0 + x_proj * u + y_proj * v
#                     node_id, coord_3d, _ = add_or_find_node(point_3d)
#                     tri_node_indices.append(node_id)
#                     tri_coords_3d.append(coord_3d)
                
#                 element_id = mesh_counter + element_id_offset
#                 mesh_name = f"T{element_id}"
#                 mesh_elements[mesh_name] = {
#                     'type': 'triangle',
#                     'nodes': tri_node_indices,
#                     'coordinates': tri_coords_3d,
#                     'id': element_id
#                 }
#                 mesh_counter += 1
    
#     # First pass: create mesh from grid
#     for i in range(num_x_div):
#         for j in range(num_y_div):
#             x1 = min_x + i * x_step
#             x2 = x1 + x_step
#             y1 = min_y + j * y_step
#             y2 = y1 + y_step
            
#             # Create rectangle in plane coordinates and clip it
#             rect = ShapelyPolygon([(x1, y1), (x2, y1), (x2, y2), (x1, y2)])
#             clipped = rect.intersection(main_poly)
            
#             if clipped.is_empty:
#                 continue
                
#             if isinstance(clipped, MultiPolygon):
#                 for poly in clipped.geoms:
#                     process_polygon(poly)
#             else:
#                 process_polygon(clipped)
    
#     # Second pass: triangulate remaining areas
#     covered_area = ShapelyPolygon()
#     for mesh in mesh_elements.values():
#         projected = project_to_plane(mesh['coordinates'])
#         covered_area = covered_area.union(ShapelyPolygon(projected))
    
#     remaining_area = main_poly.difference(covered_area)
    
#     if not remaining_area.is_empty:
#         if isinstance(remaining_area, MultiPolygon):
#             for poly in remaining_area.geoms:
#                 process_polygon(poly, 'remaining')
#         else:
#             process_polygon(remaining_area, 'remaining')
    
#     # Create a mapping of all node IDs to their coordinates
#     node_names = {node_id: coord for internal_id, (node_id, coord) in node_positions.items()}
    
#     # Find closest mesh points to predefined points
#     closest_mesh_points = {}
    
#     for name, predefined_point in predefined_points.items():
#         min_dist = float('inf')
#         closest_node = None
#         for node_id, mesh_point in node_names.items():
#             dist = np.linalg.norm(np.array(predefined_point) - np.array(mesh_point))
#             if dist < min_dist:
#                 min_dist = dist
#                 closest_node = node_id
#         if closest_node is not None:
#             closest_mesh_points[name] = closest_node
    
#     # Replace closest mesh points with predefined points
#     replaced_nodes = set()
#     for p_name, node_id in closest_mesh_points.items():
#         predefined_point = predefined_points[p_name]
#         replaced_nodes.add(node_id)
#         # Update node_names
#         node_names[node_id] = predefined_point
#         # Update coordinates in mesh elements
#         for elem in mesh_elements.values():
#             for i, n in enumerate(elem['nodes']):
#                 if n == node_id:
#                     elem['coordinates'][i] = predefined_point
    
#     # Add new shells if provided
#     if add_shell:
#         mesh_elements = add_new_shells(mesh_elements, node_names, add_shell)
    
#     # Remove shells if provided
#     if remove_shell:
#         mesh_elements = remove_shells(mesh_elements, remove_shell)
    
#     # 3D Plotting
#     fig = plt.figure(figsize=(12, 8))
#     ax = fig.add_subplot(111, projection='3d')
#     ax.set_title('3D Mesh Elements with Nodes (excluding predefined points)')
    
#     # Plot original shape
#     original_verts = [points + [points[0]]]
#     original_poly = Poly3DCollection(original_verts, alpha=0.3, 
#                                    facecolors='red', linewidths=2, 
#                                    edgecolors='red', linestyles='--')
#     ax.add_collection3d(original_poly)
    
#     # Plot mesh elements
#     for name, data in mesh_elements.items():
#         # Plot each element
#         color = 'cyan' if data['type'] == 'rectangle' else 'lightgreen'
#         alpha = 0.4 if data['type'] == 'rectangle' else 0.6
        
#         verts = [data['coordinates'] + [data['coordinates'][0]]]
#         poly = Poly3DCollection(verts, alpha=alpha, facecolors=color, 
#                               edgecolors='blue', linewidths=1)
#         ax.add_collection3d(poly)
        
#         # Add element name at centroid
#         centroid = np.mean(data['coordinates'], axis=0)
#         ax.text(centroid[0], centroid[1], centroid[2], name, 
#                 ha='center', va='center', fontsize=8, weight='bold')
        
#         # Plot nodes not replaced by predefined points
#         for node_num in data['nodes']:
#             if node_num not in replaced_nodes:
#                 coord = node_names[node_num]
#                 ax.scatter([coord[0]], [coord[1]], [coord[2]], c='red', s=50)
#                 ax.text(coord[0], coord[1], coord[2], f'N{node_num}', 
#                         ha='right', va='bottom', fontsize=8, color='darkred')
    
#     # Plot predefined points
#     for p_name, p_coord in predefined_points.items():
#         ax.scatter([p_coord[0]], [p_coord[1]], [p_coord[2]], c='blue', s=100, marker='*')
#         ax.text(p_coord[0], p_coord[1], p_coord[2], p_name, 
#                 ha='left', va='top', fontsize=10, color='darkblue', weight='bold')
    
#     # Set labels and legend
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
    
#     rect_patch = mpatches.Patch(color='cyan', label='Rectangular Elements')
#     tri_patch = mpatches.Patch(color='lightgreen', label='Triangular Elements')
#     node_patch = mpatches.Patch(color='red', label='Nodes (non-predefined)')
#     predef_patch = mpatches.Patch(color='blue', label='Predefined Points')
#     ax.legend(handles=[rect_patch, tri_patch, node_patch, predef_patch])
    
#     plt.tight_layout()
    
#     # Ensure the IMAGE_FOLDER exists and save plot
#     os.makedirs(IMAGE_FOLDER, exist_ok=True)
#     image_filename = f"mesh_plot_{numbering}.png"
#     image_path = os.path.join(IMAGE_FOLDER, image_filename)
#     plt.savefig(image_path, dpi=300, bbox_inches='tight')
#     print(f"✅ Saved plot to: {image_path}")
    
#     # Sort nodes anti-clockwise for all elements before JSON creation
#     for elem_name, elem_data in mesh_elements.items():
#         nodes = elem_data['nodes']
#         coords = elem_data['coordinates']
        
#         if len(nodes) >= 3:
#             sorted_nodes, sorted_coords = sort_nodes_anticlockwise(nodes, coords)
#             elem_data['nodes'] = sorted_nodes
#             elem_data['coordinates'] = sorted_coords

#     # Create the output structure
#     output_data = {
#         "elements": {},
#         "nodes": {},
#     }
    
#     # Create a mapping of replaced node IDs to predefined point names
#     replaced_nodes_mapping = {node_id: p_name for p_name, node_id in closest_mesh_points.items()}

#     # Add all nodes with their coordinates and IDs
#     for node_id, coord in node_names.items():
#         if node_id not in replaced_nodes_mapping:
#             output_data["nodes"][f"N{node_id}"] = {
#                 "id": node_id,
#                 "coordinates": [float(coord[0]), float(coord[1]), float(coord[2])],
#                 "is_predefined": node_id in replaced_nodes
#             }

#     # Add all elements in the requested format with IDs
#     for elem_name, elem_data in mesh_elements.items():
#         # Replace node names with predefined point names where applicable
#         node_names_in_elem = []
#         for node_id in elem_data['nodes']:
#             node_names_in_elem.append(replaced_nodes_mapping.get(node_id, f"N{node_id}"))
        
#         output_data["elements"][elem_name] = {
#             "id": elem_data['id'],
#             "type": elem_data['type'],
#             "nodes": node_names_in_elem,
#             "shell_section": shell_section_name
#         }
        
#     # Ensure directory exists and save JSON
#     wall_dir = os.path.join(JSON_FOLDER, 'wall')
#     os.makedirs(wall_dir, exist_ok=True)
#     full_path = os.path.join(wall_dir, f'mesh_data_with_predefined_points{numbering}.json')
    
#     with open(full_path, 'w') as f:
#         json.dump(output_data, f, indent=4)
        
#     return output_data

# def create_proper_mesh_for_closed_area_3d(points, predefined_points, num_x_div=4, num_y_div=4, numbering=1):
#     # Calculate ID offsets based on numbering parameter
#     node_id_offset = 10000 + (numbering - 1) * 1000
#     element_id_offset = 10000 + (numbering - 1) * 1000
    
#     # Calculate the plane equation ax + by + cz + d = 0
#     p0, p1, p2 = np.array(points[0]), np.array(points[1]), np.array(points[2])
#     v1 = p1 - p0
#     v2 = p2 - p0
#     normal = np.cross(v1, v2)
#     a, b, c = normal
#     d = -np.dot(normal, p0)
    
#     # Find two orthogonal vectors in the plane (basis vectors)
#     if abs(a) > 0.1 or abs(b) > 0.1:
#         u = np.array([b, -a, 0])  # Orthogonal to normal in XY plane
#     else:
#         u = np.array([0, c, -b])  # Orthogonal to normal in YZ plane
#     u = u / np.linalg.norm(u)
#     v = np.cross(normal, u)
#     v = v / np.linalg.norm(v)
    
#     # Function to project 3D points to 2D plane coordinates
#     def project_to_plane(points_3d):
#         projected = []
#         for p in points_3d:
#             vec = p - p0
#             x_proj = np.dot(vec, u)
#             y_proj = np.dot(vec, v)
#             projected.append((x_proj, y_proj))
#         return projected
    
#     # Function to ensure counter-clockwise ordering
#     def ensure_counter_clockwise(nodes, coords):
#         if len(nodes) < 3:
#             return nodes, coords
        
#         # Calculate normal vector for the polygon
#         if len(nodes) == 3:
#             # For triangles
#             v1 = np.array(coords[1]) - np.array(coords[0])
#             v2 = np.array(coords[2]) - np.array(coords[0])
#             cross = np.cross(v1, v2)
#             normal = cross  # For triangles, use the cross product as the normal
#         else:
#             # For polygons with more than 3 points
#             # Use Newell's method to compute normal
#             normal = np.zeros(3)
#             for i in range(len(coords)):
#                 current = np.array(coords[i])
#                 next_point = np.array(coords[(i+1)%len(coords)])
#                 normal[0] += (current[1] - next_point[1]) * (current[2] + next_point[2])
#                 normal[1] += (current[2] - next_point[2]) * (current[0] + next_point[0])
#                 normal[2] += (current[0] - next_point[0]) * (current[1] + next_point[1])
#             cross = normal
        
#         # Project onto plane normal to check winding
#         dot_product = np.dot(cross, normal)
        
#         # If winding is clockwise (dot product negative), reverse the order
#         if dot_product < 0:
#             nodes = nodes[::-1]
#             coords = coords[::-1]
        
#         return nodes, coords
    
#     # Project original points to 2D plane coordinates
#     points_2d = project_to_plane(points)
#     main_poly = ShapelyPolygon(points_2d)
    
#     # Get bounding box of the polygon in plane coordinates
#     min_x, min_y, max_x, max_y = main_poly.bounds
    
#     # Calculate step sizes
#     x_step = (max_x - min_x) / num_x_div
#     y_step = (max_y - min_y) / num_y_div
    
#     # Create dictionaries to store mesh and node information
#     mesh_elements = {}
#     node_positions = {}  # Stores {internal_id: (actual_node_id, coordinates)}
#     node_counter = 1
#     mesh_counter = 1
    
#     # First pass: create rectangular elements clipped to the polygon
#     for i in range(num_x_div):
#         for j in range(num_y_div):
#             x1 = min_x + i * x_step
#             x2 = x1 + x_step
#             y1 = min_y + j * y_step
#             y2 = y1 + y_step
            
#             # Create rectangle in plane coordinates and clip it
#             rect = ShapelyPolygon([(x1, y1), (x2, y1), (x2, y2), (x1, y2)])
#             clipped = rect.intersection(main_poly)
            
#             if clipped.is_empty or not isinstance(clipped, (ShapelyPolygon, MultiPolygon)):
#                 continue
                
#             if isinstance(clipped, MultiPolygon):
#                 polygons = list(clipped.geoms)
#             else:
#                 polygons = [clipped]
            
#             for poly in polygons:
#                 if not isinstance(poly, ShapelyPolygon):
#                     continue
                    
#                 ext_coords = list(poly.exterior.coords)
                
#                 if len(ext_coords) >= 3:  # At least 3 points needed for a polygon
#                     # Convert back to 3D coordinates
#                     node_indices = []
#                     coords_3d = []
#                     for coord in ext_coords[:-1]:
#                         x_proj, y_proj = coord
#                         point_3d = p0 + x_proj * u + y_proj * v
                        
#                         # Check if this 3D point already exists
#                         found = False
#                         for internal_id, (existing_node_id, existing_point) in node_positions.items():
#                             if np.linalg.norm(point_3d - existing_point) < 1e-6:
#                                 node_indices.append(existing_node_id)
#                                 coords_3d.append(existing_point)
#                                 found = True
#                                 break
                        
#                         if not found:
#                             node_id = node_counter + node_id_offset
#                             node_positions[node_counter] = (node_id, point_3d)
#                             node_indices.append(node_id)
#                             coords_3d.append(point_3d)
#                             node_counter += 1
                    
#                     # Ensure counter-clockwise ordering
#                     node_indices, coords_3d = ensure_counter_clockwise(node_indices, coords_3d)
                    
#                     # Handle polygons with more than 4 points by triangulating them
#                     if len(node_indices) > 4:
#                         # Convert to 2D coordinates for triangulation
#                         poly_2d = ShapelyPolygon(ext_coords)
#                         triangles = triangulate(poly_2d)
                        
#                         for triangle in triangles:
#                             tri_coords = list(triangle.exterior.coords)
#                             tri_node_indices = []
#                             tri_coords_3d = []
                            
#                             for coord in tri_coords[:-1]:
#                                 x_proj, y_proj = coord
#                                 point_3d = p0 + x_proj * u + y_proj * v
                                
#                                 # Find or create nodes for this triangle
#                                 found = False
#                                 for nid, coord_3d in zip(node_indices, coords_3d):
#                                     if np.linalg.norm(point_3d - coord_3d) < 1e-6:
#                                         tri_node_indices.append(nid)
#                                         tri_coords_3d.append(coord_3d)
#                                         found = True
#                                         break
                                
#                                 if not found:
#                                     node_id = node_counter + node_id_offset
#                                     node_positions[node_counter] = (node_id, point_3d)
#                                     tri_node_indices.append(node_id)
#                                     tri_coords_3d.append(point_3d)
#                                     node_counter += 1
                            
#                             # Ensure counter-clockwise ordering for triangles
#                             tri_node_indices, tri_coords_3d = ensure_counter_clockwise(tri_node_indices, tri_coords_3d)
                            
#                             element_id = mesh_counter + element_id_offset
#                             mesh_name = f"T{element_id}"
#                             mesh_elements[mesh_name] = {
#                                 'type': 'triangle',
#                                 'nodes': tri_node_indices,
#                                 'coordinates': tri_coords_3d,
#                                 'id': element_id
#                             }
#                             mesh_counter += 1
#                     else:
#                         element_id = mesh_counter + element_id_offset
#                         mesh_name = f"R{element_id}" if len(node_indices) == 4 else f"T{element_id}"
#                         elem_type = 'rectangle' if len(node_indices) == 4 else 'triangle'
#                         mesh_elements[mesh_name] = {
#                             'type': elem_type,
#                             'nodes': node_indices,
#                             'coordinates': coords_3d,
#                             'id': element_id
#                         }
#                         mesh_counter += 1
    
#     # Second pass: triangulate remaining areas
#     covered_area = ShapelyPolygon()
#     for mesh in mesh_elements.values():
#         projected = project_to_plane(mesh['coordinates'])
#         covered_area = covered_area.union(ShapelyPolygon(projected))
    
#     remaining_area = main_poly.difference(covered_area)
    
#     if not remaining_area.is_empty and isinstance(remaining_area, (ShapelyPolygon, MultiPolygon)):
#         if isinstance(remaining_area, MultiPolygon):
#             remaining_polys = list(remaining_area.geoms)
#         else:
#             remaining_polys = [remaining_area]
        
#         for poly in remaining_polys:
#             if not isinstance(poly, ShapelyPolygon):
#                 continue
                
#             ext_coords = list(poly.exterior.coords)
#             coords = ext_coords[:-1]
            
#             # Check if this is a simple polygon we can handle
#             if len(coords) <= 4:
#                 # Handle as either triangle or rectangle
#                 node_indices = []
#                 coords_3d = []
#                 for coord in coords:
#                     x_proj, y_proj = coord
#                     point_3d = p0 + x_proj * u + y_proj * v
                    
#                     found = False
#                     for internal_id, (existing_node_id, existing_point) in node_positions.items():
#                         if np.linalg.norm(point_3d - existing_point) < 1e-6:
#                             node_indices.append(existing_node_id)
#                             coords_3d.append(existing_point)
#                             found = True
#                             break
                    
#                     if not found:
#                         node_id = node_counter + node_id_offset
#                         node_positions[node_counter] = (node_id, point_3d)
#                         node_indices.append(node_id)
#                         coords_3d.append(point_3d)
#                         node_counter += 1
                
#                 # Ensure counter-clockwise ordering
#                 node_indices, coords_3d = ensure_counter_clockwise(node_indices, coords_3d)
                
#                 element_id = mesh_counter + element_id_offset
#                 mesh_name = f"R{element_id}" if len(node_indices) == 4 else f"T{element_id}"
#                 elem_type = 'rectangle' if len(node_indices) == 4 else 'triangle'
#                 mesh_elements[mesh_name] = {
#                     'type': elem_type,
#                     'nodes': node_indices,
#                     'coordinates': coords_3d,
#                     'id': element_id
#                 }
#                 mesh_counter += 1
#             else:
#                 # Complex polygon - triangulate it
#                 triangles = triangulate(poly)
#                 for triangle in triangles:
#                     tri_coords = list(triangle.exterior.coords)
#                     tri_node_indices = []
#                     tri_coords_3d = []
                    
#                     for coord in tri_coords[:-1]:
#                         x_proj, y_proj = coord
#                         point_3d = p0 + x_proj * u + y_proj * v
                        
#                         # Find or create nodes for this triangle
#                         found = False
#                         for internal_id, (existing_node_id, existing_point) in node_positions.items():
#                             if np.linalg.norm(point_3d - existing_point) < 1e-6:
#                                 tri_node_indices.append(existing_node_id)
#                                 tri_coords_3d.append(existing_point)
#                                 found = True
#                                 break
                        
#                         if not found:
#                             node_id = node_counter + node_id_offset
#                             node_positions[node_counter] = (node_id, point_3d)
#                             tri_node_indices.append(node_id)
#                             tri_coords_3d.append(point_3d)
#                             node_counter += 1
                    
#                     # Ensure counter-clockwise ordering for triangles
#                     tri_node_indices, tri_coords_3d = ensure_counter_clockwise(tri_node_indices, tri_coords_3d)
                    
#                     element_id = mesh_counter + element_id_offset
#                     mesh_name = f"T{element_id}"
#                     mesh_elements[mesh_name] = {
#                         'type': 'triangle',
#                         'nodes': tri_node_indices,
#                         'coordinates': tri_coords_3d,
#                         'id': element_id
#                     }
#                     mesh_counter += 1
    
#     # [Rest of the function remains the same...]
#     # (The rest of the function including plotting and JSON output is unchanged)

#     fig = plt.figure(figsize=(12, 8))
#     ax = fig.add_subplot(111, projection='3d')
#     ax.set_title('3D Mesh Elements with Nodes (excluding predefined points)')
    
#     # Plot original shape
#     original_verts = [points + [points[0]]]
#     original_poly = Poly3DCollection(original_verts, alpha=0.3, 
#                                    facecolors='red', linewidths=2, 
#                                    edgecolors='red', linestyles='--')
#     ax.add_collection3d(original_poly)
    
#     # Plot mesh elements
#     for name, data in mesh_elements.items():
#         if data['type'] == 'rectangle':
#             color = 'cyan'
#             alpha = 0.4
#         else:
#             color = 'lightgreen'
#             alpha = 0.6
        
#         verts = [data['coordinates'] + [data['coordinates'][0]]]
#         poly = Poly3DCollection(verts, alpha=alpha, facecolors=color, 
#                               edgecolors='blue', linewidths=1)
#         ax.add_collection3d(poly)
        
#         centroid = np.mean(data['coordinates'], axis=0)
#         ax.text(centroid[0], centroid[1], centroid[2], name, 
#                 ha='center', va='center', fontsize=8, weight='bold')
        
        
#         # Initialize replaced_nodes as an empty set at the beginning of the plotting section
#         # Create a mapping of all node IDs to their coordinates
#         node_names = {node_id: coord for internal_id, (node_id, coord) in node_positions.items()}

#         # Initialize replaced_nodes as an empty set
#         replaced_nodes = set()

#         # Find closest mesh points to predefined points
#         predefined_points_list = list(predefined_points.values())
#         closest_mesh_points = {}

#         for i, predefined_point in enumerate(predefined_points_list):
#             min_dist = float('inf')
#             closest_node = None
#             for node_id, mesh_point in node_names.items():
#                 dist = np.linalg.norm(predefined_point - mesh_point)
#                 if dist < min_dist:
#                     min_dist = dist
#                     closest_node = node_id
#             if closest_node is not None:
#                 closest_mesh_points[f'P{i+1}'] = closest_node

#         # Replace closest mesh points with predefined points
#         for p_name, node_id in closest_mesh_points.items():
#             predefined_point = predefined_points[p_name]
#             # Mark this node as replaced
#             replaced_nodes.add(node_id)
#             # Update node_names
#             node_names[node_id] = predefined_point
#             # Update coordinates in mesh elements
#             for elem in mesh_elements.values():
#                 for i, n in enumerate(elem['nodes']):
#                     if n == node_id:
#                         elem['coordinates'][i] = predefined_point





#         # Only plot nodes that weren't replaced by predefined points
#         for node_num in data['nodes']:
#             if node_num not in replaced_nodes:
#                 coord = node_names[node_num]
#                 ax.scatter([coord[0]], [coord[1]], [coord[2]], c='red', s=50)
#                 ax.text(coord[0], coord[1], coord[2], f'N{node_num}', 
#                         ha='right', va='bottom', fontsize=8, color='darkred')
    
#     # Plot predefined points
#     for p_name, p_coord in predefined_points.items():
#         ax.scatter([p_coord[0]], [p_coord[1]], [p_coord[2]], c='blue', s=100, marker='*')
#         ax.text(p_coord[0], p_coord[1], p_coord[2], p_name, 
#                 ha='left', va='top', fontsize=10, color='darkblue', weight='bold')
    
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
    
#     rect_patch = mpatches.Patch(color='cyan', label='Rectangular Elements')
#     tri_patch = mpatches.Patch(color='lightgreen', label='Triangular Elements')
#     node_patch = mpatches.Patch(color='red', label='Nodes (non-predefined)')
#     predef_patch = mpatches.Patch(color='blue', label='Predefined Points')
#     ax.legend(handles=[rect_patch, tri_patch, node_patch, predef_patch])
    
#     plt.tight_layout()
#     # plt.show()
    
#     # Create the output structure with elements and nodes including IDs
#     output_data = {
#         "elements": {},
#         "nodes": {},
#     }

#     # Add all nodes with their coordinates and IDs
#     for node_id, coord in node_names.items():
#         output_data["nodes"][f"N{node_id}"] = {
#             "id": node_id,
#             "coordinates": [float(coord[0]), float(coord[1]), float(coord[2])],
#             "is_predefined": node_id in replaced_nodes
#         }

#     # Create a mapping of replaced node IDs to predefined point names
#     replaced_nodes_mapping = {node_id: p_name for p_name, node_id in closest_mesh_points.items()}

#     # Add all elements in the requested format with IDs
#     for elem_name, elem_data in mesh_elements.items():
#         # Replace node names with predefined point names where applicable
#         node_names_in_elem = []
#         for node_id in elem_data['nodes']:
#             if node_id in replaced_nodes_mapping:
#                 node_names_in_elem.append(replaced_nodes_mapping[node_id])
#             else:
#                 node_names_in_elem.append(f"N{node_id}")
        
#         output_data["elements"][elem_name] = {
#             "id": elem_data['id'],
#             "type": elem_data['type'],
#             "nodes": node_names_in_elem,
#             # "nodes_coordinate": [[float(coord[0]), float(coord[1]), float(coord[2])] 
#             #                     for coord in elem_data['coordinates']]
#         }

#     # Save to JSON file
#     with open(f'lamb/mesh_data_with_predefined_points{numbering}.json', 'w') as f:
#         json.dump(output_data, f, indent=4)

#     return output_data



# def add_new_shells(mesh_elements, node_names, add_shell):
#     """
#     Add new shells to the mesh based on predefined points
    
#     Parameters:
#     -----------
#     mesh_elements: dict
#         Dictionary containing the current mesh elements
#     node_names: dict
#         Dictionary mapping node IDs to their coordinates
#     add_shell: dict
#         Dictionary with shell names as keys and lists of point names as values
        
#     Returns:
#     --------
#     dict: Updated mesh elements with new shells added
#     """
    
#     # Create a copy of the input mesh_elements to avoid modifying the original
#     updated_mesh = mesh_elements.copy()
#     # Handle empty case
#     if not mesh_elements:
#         mesh_counter = 1  # Start from 1 if no elements exist
#     else:
#         mesh_counter = max(elem_data['id'] for elem_data in mesh_elements.values()) + 1
    
#     # Process each shell to be added
#     for shell_name, points in add_shell.items():
#         if len(points) == 3:  # Triangle
#             element_id = mesh_counter
#             mesh_name = f"T{element_id}"
#             shell_type = 'triangle'
#         elif len(points) == 4:  # Quadrilateral
#             element_id = mesh_counter
#             mesh_name = f"R{element_id}"
#             shell_type = 'rectangle'
#         else:
#             print(f"Warning: Shell {shell_name} has {len(points)} points, only triangles (3) and quads (4) are supported")
#             continue
        
#         # Extract node IDs and coordinates
#         nodes = []
#         coords = []
        
#         # Find node IDs corresponding to predefined points
#         for point_name in points:
#             found = False
#             for node_id, node_info in node_names.items():
#                 if node_id == point_name:
#                     found = True
#                     nodes.append(node_id)
#                     coords.append(node_info)
#                     break
            
#             if not found:
#                 print(f"Warning: Point {point_name} not found in mesh nodes")
#                 # You could implement logic to create new nodes here if needed
        
#         # Only create shell if we have all required points
#         if len(nodes) == len(points):
#             updated_mesh[shell_name] = {
#                 'type': shell_type,
#                 'nodes': nodes,
#                 'coordinates': coords,
#                 'id': element_id
#             }
#             mesh_counter += 1
#             print(f"Added shell {shell_name} with ID {element_id}")
    
#     return updated_mesh


# def remove_shells(mesh_elements, remove_shell):
#     """
#     Remove specified shells from the mesh
    
#     Parameters:
#     -----------
#     mesh_elements: dict
#         Dictionary containing the current mesh elements
#     remove_shell: list
#         List of shell names to be removed
        
#     Returns:
#     --------
#     dict: Updated mesh elements with specified shells removed
#     """
#     # Create a copy of the input mesh_elements to avoid modifying the original
#     final_mesh = mesh_elements.copy()
    
#     for shell_name in remove_shell:
#         if shell_name in final_mesh:
#             del final_mesh[shell_name]
#             print(f"Removed shell {shell_name}")
#         else:
#             print(f"Warning: Shell {shell_name} not found in mesh")
    
#     return final_mesh



# def sort_nodes_anticlockwise(nodes, coords):
#     """
#     Sort nodes in anti-clockwise order around their centroid.
    
#     Args:
#         nodes: List of node IDs
#         coords: List of corresponding 3D coordinates
        
#     Returns:
#         Tuple of (sorted_nodes, sorted_coords)
#     """
#     if len(nodes) <= 3:  # Triangles are already planar
#         return nodes, coords
        
#     # Calculate centroid
#     centroid = np.mean(coords, axis=0)
    
#     # Center the points
#     centered = np.array(coords) - centroid
    
#     # Compute normal vector using first 3 points
#     normal = np.cross(centered[1] - centered[0], centered[2] - centered[0])
#     normal = normal / np.linalg.norm(normal)
    
#     # Create basis vectors for projection
#     if abs(normal[0]) > 0.1 or abs(normal[1]) > 0.1:
#         u = np.array([normal[1], -normal[0], 0])  # Orthogonal to normal in XY plane
#     else:
#         u = np.array([0, normal[2], -normal[1]])  # Orthogonal to normal in YZ plane
#     u = u / np.linalg.norm(u)
#     v = np.cross(normal, u)
    
#     # Project points to 2D plane
#     projected = []
#     for point in centered:
#         x_proj = np.dot(point, u)
#         y_proj = np.dot(point, v)
#         projected.append((x_proj, y_proj))
    
#     # Calculate angles from centroid
#     angles = np.arctan2([p[1] for p in projected], [p[0] for p in projected])
    
#     # Convert to positive angles [0, 2π] and sort
#     positive_angles = np.where(angles < 0, angles + 2*np.pi, angles)
#     sorted_indices = np.argsort(positive_angles)
    
#     # Return sorted nodes and coordinates
#     sorted_nodes = [nodes[i] for i in sorted_indices]
#     sorted_coords = [coords[i] for i in sorted_indices]
    
#     return sorted_nodes, sorted_coords

# def create_proper_mesh_for_closed_area_3d1(points, predefined_points , shell_section_name,  JSON_FOLDER, IMAGE_FOLDER, 
#                                          num_x_div=4, num_y_div=4, numbering=1,add_shell=None, remove_shell=None):

#     # Calculate ID offsets based on numbering parameter
#     node_id_offset = 10000 + (numbering - 1) * 1000
#     element_id_offset = 10000 + (numbering - 1) * 1000
    
#     # Calculate the plane equation ax + by + cz + d = 0
#     p0, p1, p2 = np.array(points[0]), np.array(points[1]), np.array(points[2])
#     v1 = p1 - p0
#     v2 = p2 - p0
#     normal = np.cross(v1, v2)
#     a, b, c = normal
#     d = -np.dot(normal, p0)
    
#     # Find two orthogonal vectors in the plane (basis vectors)
#     if abs(a) > 0.1 or abs(b) > 0.1:
#         u = np.array([b, -a, 0])  # Orthogonal to normal in XY plane
#     else:
#         u = np.array([0, c, -b])  # Orthogonal to normal in YZ plane
#     u = u / np.linalg.norm(u)
#     v = np.cross(normal, u)
#     v = v / np.linalg.norm(v)
    
#     # Function to project 3D points to 2D plane coordinates
#     def project_to_plane(points_3d):
#         projected = []
#         for p in points_3d:
#             vec = p - p0
#             x_proj = np.dot(vec, u)
#             y_proj = np.dot(vec, v)
#             projected.append((x_proj, y_proj))
#         return projected
    
#     # Project original points to 2D plane coordinates
#     points_2d = project_to_plane(points)
#     main_poly = ShapelyPolygon(points_2d)
    
#     # Get bounding box of the polygon in plane coordinates
#     min_x, min_y, max_x, max_y = main_poly.bounds
    
#     # Calculate step sizes
#     x_step = (max_x - min_x) / num_x_div
#     y_step = (max_y - min_y) / num_y_div
    
#     # Create dictionaries to store mesh and node information
#     mesh_elements = {}
#     node_positions = {}  # Stores {internal_id: (actual_node_id, coordinates)}
#     node_counter = 1
#     mesh_counter = 1
    
#     # First pass: create rectangular elements clipped to the polygon
#     for i in range(num_x_div):
#         for j in range(num_y_div):
#             x1 = min_x + i * x_step
#             x2 = x1 + x_step
#             y1 = min_y + j * y_step
#             y2 = y1 + y_step
            
#             # Create rectangle in plane coordinates and clip it
#             rect = ShapelyPolygon([(x1, y1), (x2, y1), (x2, y2), (x1, y2)])
#             clipped = rect.intersection(main_poly)
            
#             if clipped.is_empty or not isinstance(clipped, (ShapelyPolygon, MultiPolygon)):
#                 continue
                
#             if isinstance(clipped, MultiPolygon):
#                 polygons = list(clipped.geoms)
#             else:
#                 polygons = [clipped]
            
#             for poly in polygons:
#                 if not isinstance(poly, ShapelyPolygon):
#                     continue
                    
#                 ext_coords = list(poly.exterior.coords)
                
#                 if len(ext_coords) >= 3:  # At least 3 points needed for a polygon
#                     # Convert back to 3D coordinates
#                     node_indices = []
#                     coords_3d = []
#                     for coord in ext_coords[:-1]:
#                         x_proj, y_proj = coord
#                         point_3d = p0 + x_proj * u + y_proj * v
                        
#                         # Check if this 3D point already exists
#                         found = False
#                         for internal_id, (existing_node_id, existing_point) in node_positions.items():
#                             if np.linalg.norm(point_3d - existing_point) < 1e-6:
#                                 node_indices.append(existing_node_id)
#                                 coords_3d.append(existing_point)
#                                 found = True
#                                 break
                        
#                         if not found:
#                             node_id = node_counter + node_id_offset
#                             node_positions[node_counter] = (node_id, point_3d)
#                             node_indices.append(node_id)
#                             coords_3d.append(point_3d)
#                             node_counter += 1
                    
#                     # Handle polygons with more than 4 points by triangulating them
#                     if len(node_indices) > 4:
#                         # Convert to 2D coordinates for triangulation
#                         poly_2d = ShapelyPolygon(ext_coords)
#                         triangles = triangulate(poly_2d)
                        
#                         for triangle in triangles:
#                             tri_coords = list(triangle.exterior.coords)
#                             tri_node_indices = []
#                             tri_coords_3d = []
                            
#                             for coord in tri_coords[:-1]:
#                                 x_proj, y_proj = coord
#                                 point_3d = p0 + x_proj * u + y_proj * v
                                
#                                 # Find or create nodes for this triangle
#                                 found = False
#                                 for nid, coord_3d in zip(node_indices, coords_3d):
#                                     if np.linalg.norm(point_3d - coord_3d) < 1e-6:
#                                         tri_node_indices.append(nid)
#                                         tri_coords_3d.append(coord_3d)
#                                         found = True
#                                         break
                                
#                                 if not found:
#                                     node_id = node_counter + node_id_offset
#                                     node_positions[node_counter] = (node_id, point_3d)
#                                     tri_node_indices.append(node_id)
#                                     tri_coords_3d.append(point_3d)
#                                     node_counter += 1
                            
#                             element_id = mesh_counter + element_id_offset
#                             mesh_name = f"T{element_id}"
#                             mesh_elements[mesh_name] = {
#                                 'type': 'triangle',
#                                 'nodes': tri_node_indices,
#                                 'coordinates': tri_coords_3d,
#                                 'id': element_id
#                             }
#                             mesh_counter += 1
#                     else:
#                         element_id = mesh_counter + element_id_offset
#                         mesh_name = f"R{element_id}" if len(node_indices) == 4 else f"T{element_id}"
#                         elem_type = 'rectangle' if len(node_indices) == 4 else 'triangle'
#                         mesh_elements[mesh_name] = {
#                             'type': elem_type,
#                             'nodes': node_indices,
#                             'coordinates': coords_3d,
#                             'id': element_id
#                         }
#                         mesh_counter += 1
    
#     # Second pass: triangulate remaining areas
#     covered_area = ShapelyPolygon()
#     for mesh in mesh_elements.values():
#         projected = project_to_plane(mesh['coordinates'])
#         covered_area = covered_area.union(ShapelyPolygon(projected))
    
#     remaining_area = main_poly.difference(covered_area)
    
#     if not remaining_area.is_empty and isinstance(remaining_area, (ShapelyPolygon, MultiPolygon)):
#         if isinstance(remaining_area, MultiPolygon):
#             remaining_polys = list(remaining_area.geoms)
#         else:
#             remaining_polys = [remaining_area]
        
#         for poly in remaining_polys:
#             if not isinstance(poly, ShapelyPolygon):
#                 continue
                
#             ext_coords = list(poly.exterior.coords)
#             coords = ext_coords[:-1]
            
#             # Check if this is a simple polygon we can handle
#             if len(coords) <= 4:
#                 # Handle as either triangle or rectangle
#                 node_indices = []
#                 coords_3d = []
#                 for coord in coords:
#                     x_proj, y_proj = coord
#                     point_3d = p0 + x_proj * u + y_proj * v
                    
#                     found = False
#                     for internal_id, (existing_node_id, existing_point) in node_positions.items():
#                         if np.linalg.norm(point_3d - existing_point) < 1e-6:
#                             node_indices.append(existing_node_id)
#                             coords_3d.append(existing_point)
#                             found = True
#                             break
                    
#                     if not found:
#                         node_id = node_counter + node_id_offset
#                         node_positions[node_counter] = (node_id, point_3d)
#                         node_indices.append(node_id)
#                         coords_3d.append(point_3d)
#                         node_counter += 1
                
#                 element_id = mesh_counter + element_id_offset
#                 mesh_name = f"R{element_id}" if len(node_indices) == 4 else f"T{element_id}"
#                 elem_type = 'rectangle' if len(node_indices) == 4 else 'triangle'
#                 mesh_elements[mesh_name] = {
#                     'type': elem_type,
#                     'nodes': node_indices,
#                     'coordinates': coords_3d,
#                     'id': element_id
#                 }
#                 mesh_counter += 1
#             else:
#                 # Complex polygon - triangulate it
#                 triangles = triangulate(poly)
#                 for triangle in triangles:
#                     tri_coords = list(triangle.exterior.coords)
#                     tri_node_indices = []
#                     tri_coords_3d = []
                    
#                     for coord in tri_coords[:-1]:
#                         x_proj, y_proj = coord
#                         point_3d = p0 + x_proj * u + y_proj * v
                        
#                         # Find or create nodes for this triangle
#                         found = False
#                         for internal_id, (existing_node_id, existing_point) in node_positions.items():
#                             if np.linalg.norm(point_3d - existing_point) < 1e-6:
#                                 tri_node_indices.append(existing_node_id)
#                                 tri_coords_3d.append(existing_point)
#                                 found = True
#                                 break
                        
#                         if not found:
#                             node_id = node_counter + node_id_offset
#                             node_positions[node_counter] = (node_id, point_3d)
#                             tri_node_indices.append(node_id)
#                             tri_coords_3d.append(point_3d)
#                             node_counter += 1
                    
#                     element_id = mesh_counter + element_id_offset
#                     mesh_name = f"T{element_id}"
#                     mesh_elements[mesh_name] = {
#                         'type': 'triangle',
#                         'nodes': tri_node_indices,
#                         'coordinates': tri_coords_3d,
#                         'id': element_id
#                     }
#                     mesh_counter += 1
    
#     # Third pass: handle any remaining unconnected nodes
#     all_nodes = {node_id: coord for internal_id, (node_id, coord) in node_positions.items()}
#     used_nodes = set()
#     for elem in mesh_elements.values():
#         used_nodes.update(elem['nodes'])
    
#     unused_nodes = set(all_nodes.keys()) - used_nodes
#     if unused_nodes:
#         # Create a KDTree for spatial queries
#         from scipy.spatial import KDTree
#         coords = np.array([all_nodes[node_id] for node_id in all_nodes.keys()])
#         kdtree = KDTree(coords)
        
#         # Try to connect unused nodes to nearby elements
#         for node_id in unused_nodes:
#             point = all_nodes[node_id]
#             # Find the closest used node
#             dist, idx = kdtree.query(point, k=2)
#             closest_node_id = list(all_nodes.keys())[idx[1]]  # idx[0] is the point itself
            
#             # Find elements containing the closest node
#             connected_elements = []
#             for elem_name, elem_data in mesh_elements.items():
#                 if closest_node_id in elem_data['nodes']:
#                     connected_elements.append(elem_data)
            
#             # Try to add this node to one of the connected elements
#             for elem in connected_elements:
#                 if elem['type'] == 'triangle' and len(elem['nodes']) == 3:
#                     # Try to convert triangle to quadrilateral
#                     new_nodes = elem['nodes'] + [node_id]
#                     if len(new_nodes) == 4:
#                         # Create new quadrilateral
#                         element_id = mesh_counter + element_id_offset
#                         mesh_name = f"R{element_id}"
#                         new_coords = elem['coordinates'] + [all_nodes[node_id]]
#                         mesh_elements[mesh_name] = {
#                             'type': 'rectangle',
#                             'nodes': new_nodes,
#                             'coordinates': new_coords,
#                             'id': element_id
#                         }
#                         mesh_counter += 1
#                         used_nodes.add(node_id)
#                         break
#                 elif elem['type'] == 'rectangle' and len(elem['nodes']) == 4:
#                     # Split rectangle into two triangles
#                     element_id1 = mesh_counter + element_id_offset
#                     element_id2 = mesh_counter + 1 + element_id_offset
#                     mesh_name1 = f"T{element_id1}"
#                     mesh_name2 = f"T{element_id2}"
                    
#                     # Create two triangles (simple diagonal split)
#                     nodes = elem['nodes']
#                     coords = elem['coordinates']
                    
#                     # Triangle 1: nodes 0,1,2
#                     mesh_elements[mesh_name1] = {
#                         'type': 'triangle',
#                         'nodes': [nodes[0], nodes[1], nodes[2]],
#                         'coordinates': [coords[0], coords[1], coords[2]],
#                         'id': element_id1
#                     }
                    
#                     # Triangle 2: nodes 0,2,3
#                     mesh_elements[mesh_name2] = {
#                         'type': 'triangle',
#                         'nodes': [nodes[0], nodes[2], nodes[3]],
#                         'coordinates': [coords[0], coords[2], coords[3]],
#                         'id': element_id2
#                     }
                    
#                     mesh_counter += 2
#                     # Now try to add the unused node to one of these new triangles
#                     break
    
#     # Create a mapping of all node IDs to their coordinates
#     node_names = {node_id: coord for internal_id, (node_id, coord) in node_positions.items()}
    
#     # Find closest mesh points to predefined points
#     predefined_points_list = list(predefined_points.values())
#     predefined_points_names = list(predefined_points.keys())
#     closest_mesh_points = {}
    
#     for i, (name, predefined_point) in enumerate(predefined_points.items()):  # Use both name and point
#         min_dist = float('inf')
#         closest_node = None
#         for node_id, mesh_point in node_names.items():
#             dist = np.linalg.norm(predefined_point - mesh_point)
#             if dist < min_dist:
#                 min_dist = dist
#                 closest_node = node_id
#         if closest_node is not None:
#             closest_mesh_points[name] = closest_node  # Use the original name (n5, n6, etc.)

#     # Replace closest mesh points with predefined points
#     replaced_nodes = set()
#     for p_name, node_id in closest_mesh_points.items():
#         predefined_point = predefined_points[p_name]  # Now this will work as we're using correct keys
#         # Mark this node as replaced
#         replaced_nodes.add(node_id)
#         # Update node_names
#         node_names[node_id] = predefined_point
#         # Update coordinates in mesh elements
#         for elem in mesh_elements.values():
#             for i, n in enumerate(elem['nodes']):
#                 if n == node_id:
#                     elem['coordinates'][i] = predefined_point

#     # 2. Add new shells
#     updated_mesh = add_new_shells(mesh_elements, node_names, add_shell)
    
#     # 3. Remove shells (if they exist)
#     final_mesh = remove_shells(updated_mesh, remove_shell)
    
#     # Update mesh_elements to use the final mesh
#     mesh_elements = final_mesh
    
#     # [Rest of the function remains the same...]
#     # (Plotting and JSON output code from previous version)
#     # 3D Plotting - Skip plotting replaced nodes
#     fig = plt.figure(figsize=(12, 8))
#     ax = fig.add_subplot(111, projection='3d')
#     ax.set_title('3D Mesh Elements with Nodes (excluding predefined points)')
    
#     # Plot original shape
#     original_verts = [points + [points[0]]]
#     original_poly = Poly3DCollection(original_verts, alpha=0.3, 
#                                    facecolors='red', linewidths=2, 
#                                    edgecolors='red', linestyles='--')
#     ax.add_collection3d(original_poly)
    
#     # Plot mesh elements
#     for name, data in mesh_elements.items():
#         if data['type'] == 'rectangle':
#             color = 'cyan'
#             alpha = 0.4
#         else:
#             color = 'lightgreen'
#             alpha = 0.6
        
#         verts = [data['coordinates'] + [data['coordinates'][0]]]
#         poly = Poly3DCollection(verts, alpha=alpha, facecolors=color, 
#                               edgecolors='blue', linewidths=1)
#         ax.add_collection3d(poly)
        
#         centroid = np.mean(data['coordinates'], axis=0)
#         ax.text(centroid[0], centroid[1], centroid[2], name, 
#                 ha='center', va='center', fontsize=8, weight='bold')
        
#         # Only plot nodes that weren't replaced by predefined points
#         for node_num in data['nodes']:
#             if node_num not in replaced_nodes:
#                 coord = node_names[node_num]
#                 ax.scatter([coord[0]], [coord[1]], [coord[2]], c='red', s=50)
#                 ax.text(coord[0], coord[1], coord[2], f'N{node_num}', 
#                         ha='right', va='bottom', fontsize=8, color='darkred')
    
#     # Plot predefined points
#     for p_name, p_coord in predefined_points.items():
#         ax.scatter([p_coord[0]], [p_coord[1]], [p_coord[2]], c='blue', s=100, marker='*')
#         ax.text(p_coord[0], p_coord[1], p_coord[2], p_name, 
#                 ha='left', va='top', fontsize=10, color='darkblue', weight='bold')
    
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
    
#     rect_patch = mpatches.Patch(color='cyan', label='Rectangular Elements')
#     tri_patch = mpatches.Patch(color='lightgreen', label='Triangular Elements')
#     node_patch = mpatches.Patch(color='red', label='Nodes (non-predefined)')
#     predef_patch = mpatches.Patch(color='blue', label='Predefined Points')
#     ax.legend(handles=[rect_patch, tri_patch, node_patch, predef_patch])
    
#     plt.tight_layout()
#     # Ensure the IMAGE_FOLDER exists
#     os.makedirs(IMAGE_FOLDER, exist_ok=True)

#     # Define the image file path (e.g., "mesh_plot.png")
#     image_filename = f"mesh_plot_{numbering}.png"  # Include numbering if needed
#     image_path = os.path.join(IMAGE_FOLDER, image_filename)

#     # Save the figure
#     plt.savefig(image_path, dpi=300, bbox_inches='tight')
#     print(f"✅ Saved plot to: {image_path}")
#     # plt.show()
    
#     # Sort nodes anti-clockwise for ALL elements before JSON creation
#     for elem_name, elem_data in mesh_elements.items():
#         # Get current nodes and coordinates
#         nodes = elem_data['nodes']
#         coords = elem_data['coordinates']
        
#         # Sort them anti-clockwise (only if 3+ nodes)
#         if len(nodes) >= 3:
#             sorted_nodes, sorted_coords = sort_nodes_anticlockwise(nodes, coords)
#             elem_data['nodes'] = sorted_nodes
#             elem_data['coordinates'] = sorted_coords

#     # Create the output structure with elements and nodes including IDs
#     output_data = {
#         "elements": {},
#         "nodes": {},
#     }
#     # print(f'node_names.items()={node_names.items()}')
#     replaced_nodes_mapping = {node_id: p_name for p_name, node_id in closest_mesh_points.items()}

#     # Add all nodes with their coordinates and IDs
#     for node_id, coord in node_names.items():
#         if node_id in replaced_nodes_mapping:
#                 continue
#         else:
#             # node_names_in_elem.append(f"N{node_id}")
#             output_data["nodes"][f"N{node_id}"] = {
#                 "id": node_id,
#                 "coordinates": [float(coord[0]), float(coord[1]), float(coord[2])],
#                 "is_predefined": node_id in replaced_nodes
#             }

#     # Sort nodes anti-clockwise before creating element
    
#     # Create a mapping of replaced node IDs to predefined point names
    
#     # print(f'mesh_elements.items()={mesh_elements.items()}')
#     # Add all elements in the requested format with IDs
#     for elem_name, elem_data in mesh_elements.items():
#         # Replace node names with predefined point names where applicable
#         node_names_in_elem = []
#         for node_id in elem_data['nodes']:
#             if node_id in replaced_nodes_mapping:
#                 node_names_in_elem.append(replaced_nodes_mapping[node_id])
#             else:
#                 node_names_in_elem.append(f"N{node_id}")
        
#         output_data["elements"][elem_name] = {
#             "id": elem_data['id'],
#             "type": elem_data['type'],
#             "nodes": node_names_in_elem,
#             "shell_section": shell_section_name 
#             # "nodes_coordinate": [[float(coord[0]), float(coord[1]), float(coord[2])] 
#             #                     for coord in elem_data['coordinates']]
#         }
#     wall_dir = os.path.join(JSON_FOLDER, 'wall')
#     os.makedirs(wall_dir, exist_ok=True)

#     full_path = os.path.join(wall_dir, f'mesh_data_with_predefined_points{numbering}.json')

#     # Save to JSON file
#     with open(full_path, 'w') as f:
#         json.dump(output_data, f, indent=4)

#     return output_data




# def create_shell_mesh(JSON_FOLDER="output_folder", IMAGE_FOLDER=None):
#     """
#     Create shell elements from mesh data with consistent file paths
#     Aligned with create_proper_mesh_for_closed_area_3d1() file structure
    
#     Args:
#         JSON_FOLDER (str): Path to folder containing structure data and wall subfolder
#         IMAGE_FOLDER (str, optional): Path for saving images (not used in current implementation)
#     """
    
#     # =============================================
#     # 1. File Path Configuration and Validation
#     # =============================================
#     # Structure data path
#     structure_file = os.path.join(JSON_FOLDER, "structure_data_json.json")
    
#     # Mesh files path (in wall subdirectory)
#     wall_dir = os.path.join(JSON_FOLDER, "wall")

    
#     # Validate paths exist
#     if not os.path.exists(structure_file):
#         raise FileNotFoundError(f"Structure file not found at: {structure_file}")
#     if not os.path.exists(wall_dir):
#         raise FileNotFoundError(f"Wall directory not found at: {wall_dir}")

#     # =============================================
#     # 2. Material Definition
#     # =============================================
    
#     # Define material properties
#     E = 2.1e11       # Elastic modulus (Pa)
#     nu = 0.3         # Poisson's ratio
#     rho = 2500       # Density (kg/m^3)
#     t = 0.5          # Shell thickness (m)
    
#     # Define material and section
#     ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
#     section_tag = 100
#     ops.section("PlateFiber", section_tag, 1, t)

#     # =============================================
#     # 3. Load Structure Data
#     # =============================================
    
#     with open(structure_file, 'r') as f:
#         structure_data = json.load(f)
#         # Create mapping of structure nodes: {"nID": {node_data}}
#         structure_nodes = {f"n{node['id']}": node for node in structure_data['nodes']}

#     # =============================================
#     # 4. Find and Process Mesh Files
#     # =============================================
    
#     # Find all mesh files in wall directory
#     mesh_files = []
#     for filename in os.listdir(wall_dir):
#         match = re.match(r'mesh_data_with_predefined_points(\d+)\.json', filename)
#         if match:
#             numbering = int(match.group(1))
#             mesh_files.append((numbering, os.path.join(wall_dir, filename)))
    
#     if not mesh_files:
#         raise FileNotFoundError(f"No valid mesh files found in {wall_dir}")

#     # Process each mesh file
#     for numbering, mesh_file in sorted(mesh_files, key=lambda x: x[0]):
#         print(f"\nProcessing mesh {numbering} from {mesh_file}")
        
#         with open(mesh_file, 'r') as f:
#             mesh_data = json.load(f)
        
#         # =============================================
#         # 5. Create Nodes
#         # =============================================
#         print("\nCreating Nodes:")
#         mesh_node_ids = set()
#         for node_name, node_info in mesh_data['nodes'].items():
#             node_id = node_info['id']
#             x, y, z = node_info['coordinates']
#             ops.node(node_id, x, y, z)
#             mesh_node_ids.add(node_id)
#             print(f"  Created node {node_id} at ({x:.3f}, {y:.3f}, {z:.3f})")

#         # =============================================
#         # 6. Create Shell Elements
#         # =============================================
#         print("\nCreating Shell Elements:")
#         for elem_name, elem_info in mesh_data['elements'].items():
#             elem_id = elem_info['id']
#             node_tags = []
            
#             # Resolve node references
#             for node_ref in elem_info['nodes']:
#                 if node_ref.startswith('N'):  # Regular mesh node
#                     node_id = int(node_ref[1:])
#                     if node_id in mesh_node_ids:
#                         node_tags.append(node_id)
#                     else:
#                         print(f"  Warning: Mesh node {node_ref} not found")
#                 elif node_ref.startswith('P'):  # Predefined point
#                     # Find by coordinate matching
#                     target_coords = next(
#                         (coord for n, coord in zip(elem_info['nodes'], 
#                                                 elem_info.get('nodes_coordinate', [])) 
#                         if n == node_ref), None)
                    
#                     if target_coords:
#                         # Search in mesh nodes
#                         found = False
#                         for n_id, n_data in mesh_data['nodes'].items():
#                             if np.allclose(n_data['coordinates'], target_coords, atol=1e-6):
#                                 node_tags.append(n_data['id'])
#                                 found = True
#                                 break
                        
#                         # Search in structure nodes if not found
#                         if not found:
#                             for s_id, s_node in structure_nodes.items():
#                                 if np.allclose([s_node['x'], s_node['y'], s_node['z']], 
#                                               target_coords, atol=1e-6):
#                                     node_tags.append(int(s_id[1:]))
#                                     found = True
#                                     break
                        
#                         if not found:
#                             print(f"  Warning: Could not find node for {node_ref} at {target_coords}")
#                 elif node_ref.startswith('n'):  # Structure node
#                     node_id = int(node_ref[1:])
#                     if node_ref in structure_nodes:
#                         node_tags.append(node_id)
#                     else:
#                         print(f"  Warning: Structure node {node_ref} not found")

#             # Create element if all nodes were resolved
#             if len(node_tags) == len(elem_info['nodes']):
#                 if len(node_tags) == 4:
#                     ops.element("ShellMITC4", elem_id, *node_tags, section_tag)
#                     print(f"  Created ShellMITC4 element {elem_id} with nodes {node_tags}")
#                 elif len(node_tags) == 3:
#                     ops.element("ShellDKGT", elem_id, *node_tags, section_tag)
#                     print(f"  Created ShellDKGT element {elem_id} with nodes {node_tags}")
#                 else:
#                     print(f"  Warning: Element {elem_name} has unsupported number of nodes ({len(node_tags)})")
#             else:
#                 print(f"  Warning: Element {elem_name} is missing nodes (expected {len(elem_info['nodes'])}, found {len(node_tags)})")

#     print("\nShell mesh creation completed successfully")
#     return True




# Create dictionaries for adding and removing shells (with corrected spelling)
# add_shell = {
#     "Shell1": ["P1", "P2", "P3", "P4"],  # Quadrilateral shell
#     "Shell2": ["P5", "P6", "P7"],       # Triangular shell
#     "Shell3": ["P8", "P9", "P10", "P11"] # Another quadrilateral
# }

# remove_shell = ["Shell4", "Shell5"]  # Shells to be removed (also corrected spelling)

# predefined_points = {
#         # 'P1': np.array([0, 0, 0]),
#         # 'P2': np.array([2, 0, 1.5]),

#     }

# Example usage:
# horizontal_points = [
#     [0, 0, 0],
#     [2, 0, 0],
#     [2.5, 1.5, 0],
#     [1.5, 2.5, 0],
#     [0.5, 2, 0],
#     [-0.5, 1, 0]
# ]


# vertical_points = [
#     [0, 0, 0],    # Bottom-front
#     [3, 0, 0],    # Bottom-back
#     [3, 0, 2],    # Top-back
#     [0, 0, 2]     # Top-front
# ]



# inclined_points = [
#     [0, 0, 0],    # Base point 1
#     [3, 0, 0],    # Base point 2
#     [2, 2, 2],    # Top point 1
#     [0, 2, 2]     # Top point 2
# ]




# complex_inclined_points = [
#     [0, 0, 0],     # Base point
#     [4, 0, 1],     # Right point (slightly elevated)
#     [3, 3, 3],     # Top point
#     [1, 3, 2],     # Left point
#     [0, 2, 1.5]    # Front point
# ]




# # Horizontal element (works as before)
# mesh_elements = create_proper_mesh_for_closed_area_3d1(horizontal_points , predefined_points)

# # # Vertical element (now works correctly)
# # mesh_elements = create_proper_mesh_for_closed_area_3d(vertical_points)
# mesh_elements = create_proper_mesh_for_closed_area_3d1(vertical_points, predefined_points)

# # # Inclined element (now works correctly)
# mesh_elements = create_proper_mesh_for_closed_area_3d1(inclined_points, predefined_points)
# # print(mesh_elements)
# # # Complex inclined element
# mesh_elements = create_proper_mesh_for_closed_area_3d1(complex_inclined_points, predefined_points)



# create_shell_mesh()


# import os
# import numpy as np
# from shapely.geometry import Polygon as ShapelyPolygon
# from shapely.ops import triangulate
# from shapely.geometry import MultiPolygon
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
# import matplotlib.patches as mpatches
# import json
# import re
# import openseespy.opensees as ops

# # Set up directories
# OUTPUT_FOLDER = "output_folder"  # Main output directory
# os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist

# # Subdirectories
# JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
# IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
# os.makedirs(JSON_FOLDER, exist_ok=True)
# os.makedirs(IMAGE_FOLDER, exist_ok=True)

# # Define all input variables
# add_shell = {
#     # "Shell1": ["P1", "P2", "P3", "P4"],  # Quadrilateral shell
#     # "Shell2": ["P5", "P6", "P7"],       # Triangular shell
#     # "Shell3": ["P8", "P9", "P10", "P11"] # Another quadrilateral
# }

# # remove_shell = ["Shell4", "Shell5"]  # Shells to be removed
# remove_shell = []  # Shells to be removed

# predefined_points = {
#     # 'P1': np.array([0, 0, 0]),
#     # 'P2': np.array([2, 0, 1.5]),
#     # 'P3': np.array([2, 2, 1.5]),
#     # 'P4': np.array([0, 2, 0]),
#     # 'P5': np.array([1, 1, 1]),
#     # 'P6': np.array([1.5, 0.5, 1.2]),
#     # 'P7': np.array([0.5, 1.5, 0.8]),
#     # 'P8': np.array([0.5, 0, 0.5]),
#     # 'P9': np.array([1.5, 0, 1]),
#     # 'P10': np.array([1.5, 2, 1]),
#     # 'P11': np.array([0.5, 2, 0.5])
# }

# # Define all surface configurations
# surface_configurations = {
#     "horizontal": [
#         [0, 0, 0],
#         [2, 0, 0],
#         [2.5, 1.5, 0],
#         [1.5, 2.5, 0],
#         [0.5, 2, 0],
#         [-0.5, 1, 0]
#     ],
#     "vertical": [
#         [0, 0, 0],    # Bottom-front
#         [3, 0, 0],    # Bottom-back
#         [3, 0, 2],    # Top-back
#         [0, 0, 2]     # Top-front
#     ],
#     "inclined": [
#         [0, 0, 0],    # Base point 1
#         [3, 0, 0],    # Base point 2
#         [2, 2, 2],    # Top point 1
#         [0, 2, 2]     # Top point 2
#     ],
#     "complex_inclined": [
#         [0, 0, 0],     # Base point
#         [4, 0, 1],     # Right point (slightly elevated)
#         [3, 3, 3],     # Top point
#         [1, 3, 2],     # Left point
#         [0, 2, 1.5]    # Front point
#     ]
# }


# # Process each configuration
# for config_name, points in surface_configurations.items():
#     print(f"\nProcessing {config_name} configuration...")
    
#     # Call mesh creation function with shell parameters
#     mesh_data = create_proper_mesh_for_closed_area_3d1(
#         points=points,
#         predefined_points=predefined_points,
#         JSON_FOLDER=JSON_FOLDER,
#         IMAGE_FOLDER=IMAGE_FOLDER,
#         num_x_div=4,
#         num_y_div=4,
#         numbering=list(surface_configurations.keys()).index(config_name) + 1,
#         add_shell=add_shell,
#         remove_shell=remove_shell
#     )
    
#     # Call shell creation function
#     create_shell_mesh(
#         JSON_FOLDER=JSON_FOLDER,
#         IMAGE_FOLDER=IMAGE_FOLDER
#     )

# print("\nAll configurations processed successfully!")
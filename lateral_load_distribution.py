# import math
# import numpy as np
# from numpy.linalg import inv, solve

# def calculate_wall_rotation(vertices):
#     """Calculate wall rotation angle from vertices based on principal axis"""
#     # Extract x and y coordinates
#     x = vertices[:, 0]
#     y = vertices[:, 1]
    
#     # Calculate centroid
#     cx = np.mean(x)
#     cy = np.mean(y)
    
#     # Center coordinates
#     x_centered = x - cx
#     y_centered = y - cy
    
#     # Calculate covariance matrix
#     cov_matrix = np.cov(x_centered, y_centered)
    
#     # Calculate eigenvalues and eigenvectors
#     eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
    
#     # Get principal eigenvector (direction of maximum variance)
#     principal_eigenvector = eigenvectors[:, np.argmax(eigenvalues)]
    
#     # Calculate angle in radians (-π to π)
#     angle_rad = math.atan2(principal_eigenvector[1], principal_eigenvector[0])
    
#     # Convert to degrees (0-360)
#     angle_deg = math.degrees(angle_rad) % 360
    
#     # Adjust angle to be between -90 and 90 degrees (for walls)
#     if angle_deg > 90:
#         angle_deg -= 180
#     elif angle_deg < -90:
#         angle_deg += 180
    
#     return angle_deg

# def calculate_wall_properties(vertices):
#     """Calculate wall properties from vertices using shoelace formula"""
#     # Extract x and y coordinates (ignore z for 2D calculations)
#     x = vertices[:, 0]
#     y = vertices[:, 1]
    
#     # Close the polygon by adding the first point at the end
#     x_closed = np.append(x, x[0])
#     y_closed = np.append(y, y[0])
    
#     # Calculate area using shoelace formula
#     area = 0.5 * abs(sum(x_closed[i] * y_closed[i+1] - x_closed[i+1] * y_closed[i] 
#                          for i in range(len(x_closed)-1)))
    
#     # Calculate centroid
#     if area > 0:
#         cx = sum((x_closed[i] + x_closed[i+1]) * (x_closed[i] * y_closed[i+1] - x_closed[i+1] * y_closed[i]) 
#                  for i in range(len(x_closed)-1)) / (6 * area)
#         cy = sum((y_closed[i] + y_closed[i+1]) * (x_closed[i] * y_closed[i+1] - x_closed[i+1] * y_closed[i]) 
#                  for i in range(len(y_closed)-1)) / (6 * area)
#     else:
#         cx = cy = 0
    
#     # Calculate second moments of area (moments of inertia)
#     Ix = 0
#     Iy = 0
    
#     for i in range(len(x_closed)-1):
#         x1, y1 = x_closed[i], y_closed[i]
#         x2, y2 = x_closed[i+1], y_closed[i+1]
        
#         # Contribution from this edge to Ix and Iy
#         cross_product = x1 * y2 - x2 * y1
        
#         Ix += cross_product * (y1**2 + y1*y2 + y2**2)
#         Iy += cross_product * (x1**2 + x1*x2 + x2**2)
    
#     Ix = abs(Ix) / 12
#     Iy = abs(Iy) / 12
    
#     return area, Ix, Iy, cx, cy

# def transform_stiffness_matrix(Kx_local, Ky_local, rotation_deg):
#     """Transform local stiffness to global coordinates"""
#     theta = math.radians(rotation_deg)
#     cos_theta = math.cos(theta)
#     sin_theta = math.sin(theta)
    
#     # Full stiffness transformation
#     Kxx = Kx_local * cos_theta**2 + Ky_local * sin_theta**2
#     Kyy = Kx_local * sin_theta**2 + Ky_local * cos_theta**2
#     Kxy = (Kx_local - Ky_local) * sin_theta * cos_theta
    
#     return Kxx, Kyy, Kxy

# def calculate_mesh_properties(node_coordinates):
#     """
#     Calculate area and moments of inertia for a single mesh element
#     Supports both triangular (3 nodes) and quadrilateral (4 nodes) meshes
#     """
#     coords = np.array(node_coordinates)
#     num_nodes = len(coords)
    
#     if num_nodes == 3:
#         # Triangular mesh
#         x1, y1, z1 = coords[0]
#         x2, y2, z2 = coords[1]
#         x3, y3, z3 = coords[2]
        
#         # Area using cross product
#         area = 0.5 * abs((x2-x1)*(z3-z1) - (x3-x1)*(z2-z1))
        
#         # Centroid
#         cx = (x1 + x2 + x3) / 3
#         cz = (z1 + z2 + z3) / 3
        
#         # Moments of inertia about centroidal axes
#         # Moments of inertia about centroidal axes (corrected formula)
#         Ix = area * ((z1-cz)**2 + (z2-cz)**2 + (z3-cz)**2) / 6
#         Iy = area * ((x1-cx)**2 + (x2-cx)**2 + (x3-cx)**2) / 6
        
#     elif num_nodes == 4:
#         # Quadrilateral mesh - use shoelace formula
#         x = coords[:, 0]
#         z = coords[:, 2]  # Using z-coordinate as the "y" for 2D calculations
        
#         # Close the polygon
#         x_closed = np.append(x, x[0])
#         z_closed = np.append(z, z[0])
        
#         # Area using shoelace formula
#         area = 0.5 * abs(sum(x_closed[i] * z_closed[i+1] - x_closed[i+1] * z_closed[i] 
#                            for i in range(len(x_closed)-1)))
        
#         # Centroid
#         if area > 0:
#             cx = sum((x_closed[i] + x_closed[i+1]) * (x_closed[i] * z_closed[i+1] - x_closed[i+1] * z_closed[i]) 
#                      for i in range(len(x_closed)-1)) / (6 * area)
#             cz = sum((z_closed[i] + z_closed[i+1]) * (x_closed[i] * z_closed[i+1] - x_closed[i+1] * z_closed[i]) 
#                      for i in range(len(z_closed)-1)) / (6 * area)
#         else:
#             cx = np.mean(x)
#             cz = np.mean(z)
        
#         # Moments of inertia using shoelace-based formulas
#         Ix = 0
#         Iy = 0
        
#         for i in range(len(x_closed)-1):
#             x1, z1 = x_closed[i], z_closed[i]
#             x2, z2 = x_closed[i+1], z_closed[i+1]
            
#             cross_product = x1 * z2 - x2 * z1
            
#             Ix += cross_product * (z1**2 + z1*z2 + z2**2)
#             Iy += cross_product * (x1**2 + x1*x2 + x2**2)
        
#         Ix = abs(Ix) / 12
#         Iy = abs(Iy) / 12
        
#     else:
#         raise ValueError(f"Unsupported number of nodes: {num_nodes}. Only 3 or 4 nodes supported.")
    
#     return area, Ix, Iy, cx, cz

# def calculate_mesh_based_wall_stiffness(mesh_data, E, poisson_ratio, thickness, story_height=None):
#     """
#     Calculate shear wall stiffness using mesh-wise approach
    
#     Args:
#         mesh_data: Dictionary containing mesh information for a floor level
#         E: Elastic modulus (ksi)
#         poisson_ratio: Material Poisson's ratio
#         thickness: Wall thickness (inches)
#         story_height: Story height (inches) - optional, used only if mesh doesn't provide height info
    
#     Returns:
#         Dictionary with detailed stiffness calculations
#     """
    
#     total_area = 0
#     total_Ix = 0
#     total_Iy = 0
#     weighted_cx = 0
#     weighted_cz = 0
#     mesh_details = []
    
#     # Calculate wall dimensions dynamically from mesh coordinates
#     all_x_coords = []
#     all_z_coords = []
    
#     for mesh_name, mesh_info in mesh_data.items():
#         node_coords = mesh_info["node_coordinates"]
        
#         # Collect all coordinates to find wall bounds
#         for coord in node_coords:
#             all_x_coords.append(coord[0])
#             all_z_coords.append(coord[2])
        
#         # Calculate properties for this mesh
#         area, Ix, Iy, cx, cz = calculate_mesh_properties(node_coords)
        
#         # Store mesh details
#         mesh_details.append({
#             "name": mesh_name,
#             "area": area,
#             "Ix": Ix,
#             "Iy": Iy,
#             "centroid_x": cx,
#             "centroid_z": cz,
#             "node_coordinates": node_coords
#         })
        
#         # Accumulate total properties
#         total_area += area
#         total_Ix += Ix
#         total_Iy += Iy
#         weighted_cx += area * cx
#         weighted_cz += area * cz
    
#     # Calculate overall wall properties
#     if total_area > 0:
#         overall_cx = weighted_cx / total_area
#         overall_cz = weighted_cz / total_area
#     else:
#         overall_cx = np.mean(all_x_coords) if all_x_coords else 0
#         overall_cz = np.mean(all_z_coords) if all_z_coords else 0
    
#     # Calculate wall dimensions FROM MESH GEOMETRY
#     wall_width_x = max(all_x_coords) - min(all_x_coords) if all_x_coords else 0
#     wall_height_z = max(all_z_coords) - min(all_z_coords) if all_z_coords else 0
    
#     # CRITICAL FIX: Use actual wall dimensions, not story height
#     # For structural analysis, we need the actual wall dimensions
#     effective_height_x = wall_height_z  # Height for bending about X-axis
#     effective_height_z = wall_width_x   # Height for bending about Z-axis
    
#     # Handle edge cases
#     if effective_height_x <= 0 and story_height:
#         print(f"Warning: Wall height from mesh is {effective_height_x}, using story height {story_height}")
#         effective_height_x = story_height
    
#     if effective_height_z <= 0 and story_height:
#         print(f"Warning: Wall width from mesh is {effective_height_z}, using story height {story_height}")
#         effective_height_z = story_height
    
#     # Calculate effective dimensions for aspect ratio considerations
#     Dx = wall_width_x
#     Dz = wall_height_z
    
#     # Calculate stiffness with shear deformation effects
#     if effective_height_x > 0 and effective_height_z > 0:
#         # Aspect ratio factors (accounts for shear deformation)
#         # Use actual wall dimensions for aspect ratio calculations
#         denom_x = 1 + 0.6 * (1 + poisson_ratio) * (Dx**2 / effective_height_x**2)
#         denom_z = 1 + 0.6 * (1 + poisson_ratio) * (Dz**2 / effective_height_z**2)
        
#         # Flexural stiffness using ACTUAL wall dimensions
#         Kx_flexural = (3 * E * total_Ix) / (effective_height_x**3 * denom_x)
#         Kz_flexural = (3 * E * total_Iy) / (effective_height_z**3 * denom_z)
        
#     else:
#         Kx_flexural = Kz_flexural = 0
#         print("Warning: Invalid wall dimensions, setting stiffness to zero")
    
#     return {
#         "Kx": Kx_flexural,
#         "Kz": Kz_flexural,
#         "centroid_x": overall_cx,
#         "centroid_z": overall_cz,
#         "total_Ix": total_Ix,
#         "total_Iy": total_Iy,
#         "total_area": total_area,
#         "width_x": wall_width_x,
#         "height_z": wall_height_z,  # Renamed for clarity
#         "effective_height_x": effective_height_x,  # Added for transparency
#         "effective_height_z": effective_height_z,  # Added for transparency
#         "mesh_details": mesh_details
#     }

# def calculate_element_stiffness(element, E, poisson_ratio):
#     """Calculate element stiffness properties with corrected wall handling"""
#     x = element["x"]
#     y = element["y"]
#     story_height = element["height"]  # Renamed for clarity
#     element_type = element.get("type", "column")
#     rotation = element.get("rotation", 0.0)
    
#     # Calculate local stiffness based on element shape
#     if element_type == "wall":
#         if "meshes" in element:
#             # Calculate stiffness from meshes using ACTUAL mesh dimensions
#             mesh_data = element["meshes"]
#             thickness = element.get("thickness", 0.5)
            
#                         # Calculate rotation dynamically from first mesh if not specified
#             if "rotation" not in element:
#                 first_mesh = next(iter(mesh_data.values()))
#                 vertices = np.array(first_mesh["node_coordinates"])[:, [0, 2]]  # Use x-z coordinates
#                 rotation = calculate_wall_rotation(vertices)
#                 element["rotation"] = rotation  # Store calculated rotation

#             # Pass story_height as optional parameter, but prioritize mesh dimensions
#             mesh_results = calculate_mesh_based_wall_stiffness(
#                 mesh_data, E, poisson_ratio, thickness, story_height
#             )
            
#             Kx_local = mesh_results["Kx"]
#             Ky_local = mesh_results["Kz"]
#             cx = mesh_results["centroid_x"]
#             cz = mesh_results["centroid_z"]
            
#             # Update element position based on mesh centroid
#             element["x"] = x + cx
#             element["y"] = y + cz
            
#             # Store additional mesh information for debugging
#             element["mesh_info"] = {
#                 "actual_width": mesh_results["width_x"],
#                 "actual_height": mesh_results["height_z"],
#                 "effective_height_x": mesh_results["effective_height_x"],
#                 "effective_height_z": mesh_results["effective_height_z"],
#                 "story_height_used": story_height
#             }
            
#         else:
#             # Fall back to simple rectangular wall calculation
#             vertices = np.array(element["vertices"])
#             area, Ix, Iy, cx, cy = calculate_wall_properties(vertices)
            
#             # For vertex-based walls, use story height as intended
#             x_coords = vertices[:, 0]
#             y_coords = vertices[:, 1]
#             Dx = max(x_coords) - min(x_coords)
#             Dy = max(y_coords) - min(y_coords)
            
#             # Calculate local stiffness using story height
#             denom_x = 1 + 0.6 * (1 + poisson_ratio) * (Dx**2 / story_height**2)
#             denom_y = 1 + 0.6 * (1 + poisson_ratio) * (Dy**2 / story_height**2)
#             Kx_local = (3 * E * Ix) / (story_height**3 * denom_x)
#             Ky_local = (3 * E * Iy) / (story_height**3 * denom_y)
        
#     elif element_type == "circular_column":
#         diameter = element["diameter"]
#         I = math.pi * diameter**4 / 64
#         Kx_local = Ky_local = (12 * E * I) / (story_height**3)
        
#     elif element_type == "L_shaped_column":
#         a = element["a"]
#         b = element["b"]
#         t = element["t"]
        
#         # Calculate centroid and moments of inertia for L-shape
#         A1 = a * t
#         A2 = (b - t) * t
#         total_area = A1 + A2
        
#         # Centroid calculation
#         y1 = t / 2
#         x1 = a / 2
#         y2 = (b - t) / 2 + t
#         x2 = t / 2
        
#         xc = (A1 * x1 + A2 * x2) / total_area
#         yc = (A1 * y1 + A2 * y2) / total_area
        
#         # Moments of inertia about centroidal axes
#         Ix1 = (t * a**3)/12 + A1 * (yc - y1)**2
#         Ix2 = ((b - t) * t**3)/12 + A2 * (yc - y2)**2
#         Iy1 = (a * t**3)/12 + A1 * (xc - x1)**2
#         Iy2 = (t * (b - t)**3)/12 + A2 * (xc - x2)**2
        
#         Ix = Ix1 + Ix2
#         Iy = Iy1 + Iy2
        
#         Kx_local = (12 * E * Ix) / (story_height**3)
#         Ky_local = (12 * E * Iy) / (story_height**3)
        
#     else:  # Rectangular column
#         Dx = element["Dx"]
#         Dy = element["Dy"]
#         Ix = (Dy * Dx**3) / 12
#         Iy = (Dx * Dy**3) / 12
#         Kx_local = (12 * E * Ix) / (story_height**3)
#         Ky_local = (12 * E * Iy) / (story_height**3)
    
#     # Transform to global coordinates if rotated
#     if abs(rotation) > 1e-6:
#         Kxx, Kyy, Kxy = transform_stiffness_matrix(Kx_local, Ky_local, rotation)
#     else:
#         Kxx, Kyy, Kxy = Kx_local, Ky_local, 0.0
    
#     return Kxx, Kyy, Kxy

# def lateral_load_distribution_enhanced(Lx, Ly, E, poisson_ratio, fc, Vx, Vy, components, elements):
#     """
#     Enhanced lateral load distribution analysis accounting for:
#     - Rotated elements (walls/columns not aligned with X/Y axes)
#     - Full stiffness coupling (Kxy terms)
#     - Accurate center of rigidity calculation
#     - Displacement-based force distribution
    
#     Args:
#         Lx, Ly: Building dimensions (inches)
#         E: Elastic modulus (ksi)
#         poisson_ratio: Material Poisson's ratio
#         fc: Concrete compressive strength (ksi)
#         Vx, Vy: Applied lateral loads (kips)
#         components: List of mass components (dicts with 'wt', 'xi', 'yi')
#         elements: List of structural elements (walls/columns)
    
#     Returns:
#         Dictionary with full analysis results
#     """
    
#     # --- 1. Calculate Center of Mass (COM) ---
#     total_wt = sum(comp["wt"] for comp in components)
#     sum_wt_xi = sum(comp["wt"] * comp["xi"] for comp in components)
#     sum_wt_yi = sum(comp["wt"] * comp["yi"] for comp in components)
#     Xm = sum_wt_xi / total_wt if total_wt > 1e-10 else 0
#     Ym = sum_wt_yi / total_wt if total_wt > 1e-10 else 0
#     COM = (Xm, Ym)
    
#     # --- 2. Calculate Element Stiffness Properties ---
#     element_details = []
#     total_Kxx, total_Kyy, total_Kxy = 0, 0, 0
#     sum_Kxx_x, sum_Kyy_y, sum_Kxy_x, sum_Kxy_y = 0, 0, 0, 0
    
#     for element in elements:
#         # Get element properties
#         x = element["x"]
#         y = element["y"]
#         height = element["height"]
#         element_type = element.get("type", "column")
#         rotation = element.get("rotation", 0.0)
        
#         # Calculate stiffness
#         Kxx, Kyy, Kxy = calculate_element_stiffness(element, E, poisson_ratio)
        
#         # Store element data
#         element_details.append({
#             "name": element.get("name", "Unnamed"),
#             "Kxx": Kxx, "Kyy": Kyy, "Kxy": Kxy,
#             "x": element["x"], "y": element["y"],  # Updated position if using meshes
#             "rotation": rotation
#         })
        
#         # Accumulate global stiffness
#         total_Kxx += Kxx
#         total_Kyy += Kyy
#         total_Kxy += Kxy
#         sum_Kxx_x += Kxx * element["x"]
#         sum_Kyy_y += Kyy * element["y"]
#         sum_Kxy_x += Kxy * element["x"]
#         sum_Kxy_y += Kxy * element["y"]
    
#     # --- 3. Calculate Center of Rigidity (COR) ---
#     det = total_Kxx * total_Kyy - total_Kxy**2
#     # Standard method: weighted average of element positions by their stiffness
#     if total_Kxx > 1e-10:
#         Xr = sum_Kxx_x / total_Kxx
#     else:
#         Xr = 0
        
#     if total_Kyy > 1e-10:
#         Yr = sum_Kyy_y / total_Kyy
#     else:
#         Yr = 0
        
#     COR = (Xr, Yr)
    
#     # --- 4. Calculate Eccentricity ---
#     ex = Xm - Xr
#     ey = Ym - Yr
#     min_ex = 0.05 * Lx
#     min_ey = 0.05 * Ly
#     ex_total = ex + (-min_ex if ex < 0 else min_ex)
#     ey_total = ey + (-min_ey if ey < 0 else min_ey)
    
#     # --- 5. Torsional Moment ---
#     T = Vy * ex_total + Vx * ey_total
    
#     # --- 6. Torsional Rigidity ---
#     Jr_total = 0
#     Jr_terms = []
#     for elem in element_details:
#         dx, dy = elem["x"] - Xr, elem["y"] - Yr
#         Jr_xx = elem["Kxx"] * dy**2
#         Jr_yy = elem["Kyy"] * dx**2
#         Jr_xy = elem["Kxy"] * dx * dy
#         Jr_total += Jr_xx + Jr_yy + 2 * Jr_xy
#         Jr_terms.append((Jr_xx, Jr_yy, Jr_xy))
    
#     # --- 7. Displacement-Based Force Distribution ---
#     results = []
#     sum_Fx, sum_Fy = 0, 0
    
#     # Calculate global displacements
#     # Calculate global displacements (simplified approach)
#     delta_x = Vx / total_Kxx if total_Kxx > 1e-10 else 0
#     delta_y = Vy / total_Kyy if total_Kyy > 1e-10 else 0
    
#     theta_torsion = T / Jr_total if Jr_total > 1e-10 else 0
    
#     for i, elem in enumerate(element_details):
#         dx, dy = elem["x"] - Xr, elem["y"] - Yr
        
#         # Direct + coupling + torsion
#         Fx = (elem["Kxx"] * delta_x + elem["Kxy"] * delta_y - 
#               (elem["Kyy"] * dy + elem["Kxy"] * dx) * theta_torsion)
#         Fy = (elem["Kxy"] * delta_x + elem["Kyy"] * delta_y + 
#               (elem["Kxx"] * dx + elem["Kxy"] * dy) * theta_torsion)
        
#         # Breakdown components
#         Fx_direct = elem["Kxx"] * delta_x
#         Fy_direct = elem["Kyy"] * delta_y
#         Fx_coupling = elem["Kxy"] * delta_y
#         Fy_coupling = elem["Kxy"] * delta_x
#         Fx_torsion = -(elem["Kyy"] * dy + elem["Kxy"] * dx) * theta_torsion
#         Fy_torsion = (elem["Kxx"] * dx + elem["Kxy"] * dy) * theta_torsion
        
#         results.append({
#             "name": elem["name"],
#             "Fx": Fx, "Fy": Fy,
#             "Fx_direct": Fx_direct, "Fy_direct": Fy_direct,
#             "Fx_coupling": Fx_coupling, "Fy_coupling": Fy_coupling,
#             "Fx_torsion": Fx_torsion, "Fy_torsion": Fy_torsion,
#             "dx": dx, "dy": dy
#         })
#         sum_Fx += Fx
#         sum_Fy += Fy
    
#     # --- 8. Normalize Forces to Ensure Equilibrium ---
#     if abs(sum_Fy) > 1e-10 and abs(sum_Fy - Vy) > 1e-6:
#         scale = Vy / sum_Fy
#         for res in results:
#             res["Fy"] *= scale
#             res["Fy_direct"] *= scale
#             res["Fy_coupling"] *= scale
#             res["Fy_torsion"] *= scale
#         sum_Fy = Vy

#     if abs(sum_Fx) > 1e-10 and abs(sum_Fx - Vx) > 1e-6:
#         scale = Vx / sum_Fx
#         for res in results:
#             res["Fx"] *= scale
#             res["Fx_direct"] *= scale
#             res["Fx_coupling"] *= scale
#             res["Fx_torsion"] *= scale
#         sum_Fx = Vx
    
#     # --- 9. Return Results ---
#     return {
#         "center_of_mass": COM,
#         "center_of_rigidity": COR,
#         "total_stiffness": (total_Kxx, total_Kyy, total_Kxy),
#         "element_details": element_details,
#         "eccentricity": (ex, ey),
#         "eccentricity_total": (ex_total, ey_total),
#         "torsional_moment": T,
#         "torsional_rigidity": Jr_total,
#         "shear_distribution": results,
#         "total_shear": (sum_Fx, sum_Fy),
#         "displacements": (delta_x, delta_y, theta_torsion)
#     }

# # Example usage with mesh-based walls
# components = [
#     {"wt": 629.4652, "xi": 314.9608, "yi": 246.063},
#     {"wt": -18.88396, "xi": 551.1814, "yi": 462.598},
#     {"wt": 25.898, "xi": 5.905515, "yi": 236.2206},
#     {"wt": 25.898, "xi": 314.9608, "yi": 236.2206},
#     {"wt": 25.898, "xi": 623.9961, "yi": 236.2206},
#     {"wt": 25.898, "xi": 78.7402, "yi": 5.905515},
#     {"wt": 25.898, "xi": 551.1814, "yi": 5.905515},
#     {"wt": 25.898, "xi": 78.7402, "yi": 486.2206},
#     {"wt": 25.898, "xi": 551.1814, "yi": 427.1655},
# ]

# # Define mesh data for walls
# wall_mesh_data = {
#     "R34001": {
#         "node_coordinates": [
#             [10.0*12, 0.0*12, 6.666666666666666*12],
#             [6.666666666666666*12, 0.0*12, 6.666666666666666*12],
#             [6.666666666666666*12, 0.0*12, 10.0*12],
#             [10.0*12, 0.0*12, 10.0*12]
#         ]
#     },
#     "R34002": {
#         "node_coordinates": [
#             [10.0*12, 0.0, 3.3333333333333326*12],
#             [6.666666666666666*12, 0.0, 3.3333333333333326*12],
#             [6.666666666666666*12, 0.0, 6.666666666666666*12],
#             [10.0*12, 0.0, 6.666666666666666*12]
#         ]
#     },
#     "R34003": {
#         "node_coordinates": [
#             [10.0*12, 0.0, 0.0],
#             [6.666666666666666*12, 0.0, 0.0],
#             [6.666666666666666*12, 0.0, 3.3333333333333326*12],
#             [10.0*12, 0.0, 3.3333333333333326*12]
#         ]
#     }
# }
# Lx = 629.9216
# Ly = 492.126
# # Define elements with mesh-based walls
# elements = [
#     {
#         "name": "W2-1", 
#         "type": "wall", 
#         "meshes": wall_mesh_data,
#         "height": 0.0, 
#         "x": 78.7402,
#         "y": 5.905515,
#         "thickness": 0.5,
#         "rotation": 0.0
#     },
#     {
#         "name": "C1", 
#         "type": "circular_column", 
#         "diameter": 15.0, 
#         "height": 157.4804, 
#         "x": 200.0, 
#         "y": 200.0,
#         "rotation": 0.0
#     },
#     {
#         "name": "L1", 
#         "type": "L_shaped_column", 
#         "a": 12.0, 
#         "b": 12.0, 
#         "t": 4.0, 
#         "height": 157.4804, 
#         "x": 400.0, 
#         "y": 300.0,
#         "rotation": 45.0
#     },
#     {
#         "name": "R1", 
#         "type": "rectangular_column", 
#         "Dx": 18.0, 
#         "Dy": 12.0, 
#         "height": 157.4804, 
#         "x": 300.0, 
#         "y": 150.0,
#         "rotation": 22.5
#     },
#     # Additional rectangular columns (R2-R9)
#     {
#         "name": "R2",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": 50.0,
#         "y": 50.0,
#         "rotation": 0.0
#     },
#     {
#         "name": "R3",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": 50.0,
#         "y": Ly - 50.0,
#         "rotation": 0.0
#     },
#     {
#         "name": "R4",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": Lx - 50.0,
#         "y": 50.0,
#         "rotation": 0.0
#     },
#     {
#         "name": "R5",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": Lx - 50.0,
#         "y": Ly - 50.0,
#         "rotation": 0.0
#     },
#     {
#         "name": "R6",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": Lx/2,
#         "y": 50.0,
#         "rotation": 0.0
#     },
#     {
#         "name": "R7",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": Lx/2,
#         "y": Ly - 50.0,
#         "rotation": 0.0
#     },
#     {
#         "name": "R8",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": 50.0,
#         "y": Ly/2,
#         "rotation": 0.0
#     },
#     {
#         "name": "R9",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": Lx - 50.0,
#         "y": Ly/2,
#         "rotation": 0.0
#     }
# ]

# # Material and geometric properties

# fc = 4.35
# E = 57000 * math.sqrt(fc)
# poisson_ratio = 0.25
# Vx = -100.0
# Vy = -0.0

# # Run enhanced analysis with mesh-based walls
# results = lateral_load_distribution_enhanced(Lx, Ly, E, poisson_ratio, fc, Vx, Vy, components, elements)

# # Enhanced output with rotation effects
# print("\n" + "="*80)
# print("ENHANCED STRUCTURAL ANALYSIS WITH MESH-BASED WALLS".center(80))
# print("="*80)

# print("\n1. CENTER OF MASS & RIGIDITY")
# print("-"*80)
# print(f"{'Center of Mass (COM)':<40}: X = {results['center_of_mass'][0]:.3f} in, Y = {results['center_of_mass'][1]:.3f} in")
# print(f"{'Center of Rigidity (COR)':<40}: X = {results['center_of_rigidity'][0]:.3f} in, Y = {results['center_of_rigidity'][1]:.3f} in")

# print("\n2. ENHANCED STIFFNESS PROPERTIES")
# print("-"*80)
# print(f"{'Total Stiffness Kxx':<40}: {results['total_stiffness'][0]:.3f} kips/in")
# print(f"{'Total Stiffness Kyy':<40}: {results['total_stiffness'][1]:.3f} kips/in")
# print(f"{'Total Coupling Stiffness Kxy':<40}: {results['total_stiffness'][2]:.3f} kips/in")

# print("\n3. ELEMENT STIFFNESS DETAILS WITH ROTATION")
# print("-"*80)
# print(f"{'Element':<10}{'Kxx':<12}{'Kyy':<12}{'Kxy':<12}{'Rotation':<10}{'X (in)':<10}{'Y (in)':<10}")
# print("-"*80)
# for element in results['element_details']:
#     print(f"{element['name']:<10}{element['Kxx']:<20.3f}{element['Kyy']:<20.3f}"
#           f"{element['Kxy']:<20.3f}{element['rotation']:<10.1f}"
#           f"{element['x']:<10.3f}{element['y']:<10.3f}")

# print("\n4. ECCENTRICITY")
# print("-"*80)
# print(f"{'Natural Eccentricity (ex, ey)':<40}: ex = {results['eccentricity'][0]:.3f} in, ey = {results['eccentricity'][1]:.3f} in")
# print(f"{'Design Eccentricity (ex_total, ey_total)':<40}: ex = {results['eccentricity_total'][0]:.3f} in, ey = {results['eccentricity_total'][1]:.3f} in")

# print("\n5. ENHANCED TORSIONAL PROPERTIES")
# print("-"*80)
# print(f"{'Torsional Moment (T)':<40}: {results['torsional_moment']:.3f} kip·in")
# print(f"{'Total Torsional Rigidity Jr_total':<40}: {results['torsional_rigidity']:.3f} kip·in²")

# print("\n6. ENHANCED SHEAR DISTRIBUTION")
# print("-"*80)
# print(f"{'Element':<10}{'Total Fx':<12}{'Total Fy':<12}{'Direct Fx':<12}{'Direct Fy':<12}{'Torsion Fx':<12}{'Torsion Fy':<12}")
# print("-"*80)
# for result in results['shear_distribution']:
#     print(f"{result['name']:<10}{result['Fx']:<12.3f}{result['Fy']:<12.3f}"
#           f"{result['Fx_direct']:<12.3f}{result['Fy_direct']:<12.3f}"
#           f"{result['Fx_torsion']:<12.3f}{result['Fy_torsion']:<12.3f}")

# print("-"*80)
# print(f"{'TOTAL':<10}{results['total_shear'][0]:<12.3f}{results['total_shear'][1]:<12.3f}")

# print("\n7. COUPLING EFFECTS SUMMARY")
# print("-"*80)
# total_coupling = sum(abs(element['Kxy']) for element in results['element_details'])
# if total_coupling > 1e-6:
#     print(f"{'Total Coupling Stiffness':<40}: {total_coupling:.3f} kips/in")
#     print("WARNING: Significant coupling effects detected due to element rotations!")
#     print("The analysis includes full coupling between X and Y directions.")
# else:
#     print("No significant coupling effects - elements are primarily axis-aligned.")

# print("Applied Vx:", Vx, "kips")
# print("Applied Vy:", Vy, "kips")
# print("Total Distributed Fx:", results['total_shear'][0], "kips")
# print("Total Distributed Fy:", results['total_shear'][1], "kips")
# print("Force Balance Check - Fx:", abs(Vx - results['total_shear'][0]), "kips")
# print("Force Balance Check - Fy:", abs(Vy - results['total_shear'][1]), "kips")

# print("="*80)



from collections import defaultdict
import json
import math
import os
import numpy as np
from numpy.linalg import inv, solve

def calculate_wall_rotation(vertices):
    """Calculate wall rotation angle from vertices based on principal axis"""
    # Extract x and y coordinates
    x = vertices[:, 0]
    y = vertices[:, 1]

    # Calculate centroid
    cx = np.mean(x)
    cy = np.mean(y)

    # Center coordinates
    x_centered = x - cx
    y_centered = y - cy

    # Calculate covariance matrix
    cov_matrix = np.cov(x_centered, y_centered)

    # Get principal eigenvector (direction of maximum variance)
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
    principal_eigenvector = eigenvectors[:, np.argmax(eigenvalues)]

    # Calculate angle in radians (-π to π)
    angle_rad = math.atan2(principal_eigenvector[1], principal_eigenvector[0])

    # Convert to degrees (0-360)
    angle_deg = math.degrees(angle_rad) % 360

    # Adjust angle to be between -90 and 90 degrees (for walls)
    if angle_deg > 90:
        angle_deg -= 180
    elif angle_deg < -90:
        angle_deg += 180

    return angle_deg

def calculate_wall_properties(vertices):
    """
    Calculate wall properties (area, centroid, and centroidal moments of inertia)
    from vertices using shoelace formula and parallel axis theorem.
    Assumes vertices are in order (clockwise or counter-clockwise).
    """
    x = vertices[:, 0]
    y = vertices[:, 1]

    # Close the polygon by adding the first point at the end
    x_closed = np.append(x, x[0])
    y_closed = np.append(y, y[0])
    n = len(x_closed) - 1 # Number of segments/vertices

    # Calculate area using shoelace formula
    area_term = sum(x_closed[i] * y_closed[i+1] - x_closed[i+1] * y_closed[i] for i in range(n))
    area = 0.5 * area_term

    if abs(area) < 1e-9: # Handle degenerate cases (e.g., a line)
        cx = np.mean(x)
        cy = np.mean(y)
        return 0, 0, 0, cx, cy

    # Calculate centroid (Cx, Cy)
    cx = sum((x_closed[i] + x_closed[i+1]) * (x_closed[i] * y_closed[i+1] - x_closed[i+1] * y_closed[i])
             for i in range(n)) / (6 * area_term)
    cy = sum((y_closed[i] + y_closed[i+1]) * (x_closed[i] * y_closed[i+1] - x_closed[i+1] * y_closed[i])
             for i in range(n)) / (6 * area_term)

    # Calculate moments of inertia about the origin (Ixo, Iyo)
    Ixo_sum = 0
    Iyo_sum = 0

    for i in range(n):
        x1, y1 = x_closed[i], y_closed[i]
        x2, y2 = x_closed[i+1], y_closed[i+1]

        cross_product = (x1 * y2 - x2 * y1)

        Ixo_sum += cross_product * (y1**2 + y1*y2 + y2**2)
        Iyo_sum += cross_product * (x1**2 + x1*x2 + x2**2)

    Ixo = (1/12) * Ixo_sum
    Iyo = (1/12) * Iyo_sum

    # Apply Parallel Axis Theorem to get centroidal moments of inertia
    effective_area = abs(area)
    Ix_centroidal = Ixo - effective_area * cy**2
    Iy_centroidal = Iyo - effective_area * cx**2

    return effective_area, Ix_centroidal, Iy_centroidal, cx, cy

def transform_stiffness_matrix(Kx_local, Ky_local, rotation_deg):
    """Transform local stiffness to global coordinates for planar element"""
    theta = math.radians(rotation_deg)
    cos_theta = math.cos(theta)
    sin_theta = math.sin(theta)

    Kxx = Kx_local * cos_theta**2 + Ky_local * sin_theta**2
    Kyy = Kx_local * sin_theta**2 + Ky_local * cos_theta**2
    Kxy = (Kx_local - Ky_local) * sin_theta * cos_theta

    return Kxx, Kyy, Kxy

def calculate_mesh_properties(node_coordinates):
    """
    Calculate area and moments of inertia for a single mesh element.
    Supports both triangular (3 nodes) and quadrilateral (4 nodes) meshes.
    Assumes coordinates are (x, y, z), and we are working in the X-Z plane for walls.
    So, x-coordinate is `x`, and z-coordinate is `y` in 2D context.
    """
    coords = np.array(node_coordinates)
    num_nodes = len(coords)

    x = coords[:, 0]
    z = coords[:, 2] # Using z-coordinate as the "y" for 2D calculations of the mesh cross-section

    if num_nodes == 3:
        x_closed = np.append(x, x[0])
        z_closed = np.append(z, z[0])
        n_tri = len(x_closed) - 1

        area_term = sum(x_closed[i] * z_closed[i+1] - x_closed[i+1] * z_closed[i] for i in range(n_tri))
        area = 0.5 * area_term

        if abs(area) < 1e-9:
            cx = np.mean(x)
            cz = np.mean(z)
            return 0, 0, 0, cx, cz

        cx = sum((x_closed[i] + x_closed[i+1]) * (x_closed[i] * z_closed[i+1] - x_closed[i+1] * z_closed[i])
                 for i in range(n_tri)) / (6 * area_term)
        cz = sum((z_closed[i] + z_closed[i+1]) * (x_closed[i] * z_closed[i+1] - x_closed[i+1] * z_closed[i])
                 for i in range(n_tri)) / (6 * area_term)

        Ixo_tri = 0
        Izo_tri = 0

        for i in range(n_tri):
            x_i, z_i = x_closed[i], z_closed[i]
            x_i1, z_i1 = x_closed[i+1], z_closed[i+1]

            cross_product = (x_i * z_i1 - x_i1 * z_i)

            Ixo_tri += cross_product * (z_i**2 + z_i*z_i1 + z_i1**2)
            Izo_tri += cross_product * (x_i**2 + x_i*x_i1 + x_i1**2)

        Ixo = (1/12) * Ixo_tri
        Izo = (1/12) * Izo_tri

        effective_area = abs(area)
        Ix_centroidal = Ixo - effective_area * cz**2
        Iy_centroidal = Izo - effective_area * cx**2

    elif num_nodes == 4:
        x_closed = np.append(x, x[0])
        z_closed = np.append(z, z[0])
        n_quad = len(x_closed) - 1

        area_term = sum(x_closed[i] * z_closed[i+1] - x_closed[i+1] * z_closed[i]
                        for i in range(n_quad))
        area = 0.5 * area_term

        if abs(area) < 1e-9:
            cx = np.mean(x)
            cz = np.mean(z)
            return 0, 0, 0, cx, cz

        cx = sum((x_closed[i] + x_closed[i+1]) * (x_closed[i] * z_closed[i+1] - x_closed[i+1] * z_closed[i])
                 for i in range(n_quad)) / (6 * area_term)
        cz = sum((z_closed[i] + z_closed[i+1]) * (x_closed[i] * z_closed[i+1] - x_closed[i+1] * z_closed[i])
                 for i in range(n_quad)) / (6 * area_term)

        Ixo_sum = 0
        Izo_sum = 0

        for i in range(n_quad):
            x1, z1 = x_closed[i], z_closed[i]
            x2, z2 = x_closed[i+1], z_closed[i+1]

            cross_product = x1 * z2 - x2 * z1

            Ixo_sum += cross_product * (z1**2 + z1*z2 + z2**2)
            Izo_sum += cross_product * (x1**2 + x1*x2 + x2**2)

        Ixo = abs(Ixo_sum) / 12
        Izo = abs(Izo_sum) / 12

        effective_area = abs(area)
        Ix_centroidal = Ixo - effective_area * cz**2
        Iy_centroidal = Izo - effective_area * cx**2

    else:
        raise ValueError(f"Unsupported number of nodes: {num_nodes}. Only 3 or 4 nodes supported.")

    Ix_centroidal = max(0, Ix_centroidal)
    Iy_centroidal = max(0, Iy_centroidal)

    return abs(area), Ix_centroidal, Iy_centroidal, cx, cz

def calculate_mesh_based_wall_stiffness(mesh_data, E, poisson_ratio, thickness, story_height):
    """
    Calculate shear wall stiffness using mesh-wise approach.
    """
    total_area = 0
    total_Ix = 0
    total_Iy = 0
    weighted_cx = 0
    weighted_cz = 0
    mesh_details = []

    all_x_coords = []
    all_z_coords = []

    for mesh_name, mesh_info in mesh_data.items():
        node_coords = mesh_info["node_coordinates"]

        for coord in node_coords:
            all_x_coords.append(coord[0])
            all_z_coords.append(coord[2])

        area, Ix_mesh, Iy_mesh, cx_mesh, cz_mesh = calculate_mesh_properties(node_coords)

        mesh_details.append({
            "name": mesh_name,
            "area": area,
            "Ix": Ix_mesh,
            "Iy": Iy_mesh,
            "centroid_x": cx_mesh,
            "centroid_z": cz_mesh,
            "node_coordinates": node_coords
        })

        total_area += area
        weighted_cx += area * cx_mesh
        weighted_cz += area * cz_mesh

    if total_area > 1e-9:
        overall_cx = weighted_cx / total_area
        overall_cz = weighted_cz / total_area
    else:
        overall_cx = np.mean(all_x_coords) if all_x_coords else 0
        overall_cz = np.mean(all_z_coords) if all_z_coords else 0
        print("Warning: Total area for mesh-based wall is zero. Centroid based on mean node coords.")

    for mesh_detail in mesh_details:
        dx = mesh_detail["centroid_x"] - overall_cx
        dz = mesh_detail["centroid_z"] - overall_cz
        total_Ix += mesh_detail["Ix"] + mesh_detail["area"] * dz**2
        total_Iy += mesh_detail["Iy"] + mesh_detail["area"] * dx**2

    wall_width_x = max(all_x_coords) - min(all_x_coords) if all_x_coords else 0
    wall_height_z = max(all_z_coords) - min(all_z_coords) if all_z_coords else 0

    H = story_height

    if H <= 1e-9:
        print(f"Warning: Story height for wall is zero or too small ({H}). Setting stiffness to zero.")
        Kx_flexural = 0
        Kz_flexural = 0
    else:
        denom_x = 1 + 0.6 * (1 + poisson_ratio) * (wall_width_x**2 / H**2)
        denom_z = 1 + 0.6 * (1 + poisson_ratio) * (wall_height_z**2 / H**2)

        Kx_flexural = (3 * E * total_Iy) / (H**3 * denom_x)
        Kz_flexural = (3 * E * total_Ix) / (H**3 * denom_z)

    return {
        "Kx": Kx_flexural,
        "Kz": Kz_flexural,
        "centroid_x": overall_cx,
        "centroid_z": overall_cz,
        "total_Ix": total_Ix,
        "total_Iy": total_Iy,
        "total_area": total_area,
        "width_x": wall_width_x,
        "height_z": wall_height_z,
        "effective_flexural_height": H,
        "mesh_details": mesh_details
    }


def create_components_and_elements(JSON_FOLDER, load_combinations):
    """
    Creates components and elements for each z_level using JSON files from the specified folder,
    including properly filtered load values from different load combinations.
    """
    # Define file paths
    filtered_nodes_path = os.path.join(JSON_FOLDER, "filtered_nodes.json")
    filtered_columns_path = os.path.join(JSON_FOLDER, "filtered_columns.json")
    element_data_path = os.path.join(JSON_FOLDER, "element_data.json")
    load_data_dir = os.path.join(JSON_FOLDER, "load_data")
    lateral_load_dir = os.path.join(JSON_FOLDER, "lateral_load_distribution")
    
    # Create output directories
    os.makedirs(lateral_load_dir, exist_ok=True)

    # Load the JSON files
    try:
        with open(filtered_nodes_path, 'r') as f:
            filtered_nodes = json.load(f)
        with open(filtered_columns_path, 'r') as f:
            filtered_columns = json.load(f)
        with open(element_data_path, 'r') as f:
            element_data = json.load(f)
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Required JSON file not found: {e}")

    # Process each load combination
    all_results = {}
    for z_level, nodes in filtered_nodes.items():
        current_z = float(z_level)  # Convert z_level to float for comparison
        
        # Calculate Lx and Ly for this z_level
        x_coords = [coords[0] for coords in nodes.values()]
        y_coords = [coords[1] for coords in nodes.values()]
        Lx = max(x_coords) - min(x_coords)
        Ly = max(y_coords) - min(y_coords)
        
        for combo_name, combo_factors in load_combinations.items():
            components = []
            elements = []
            
            # Process nodal loads (existing code remains same)
            nodal_load_file = os.path.join(load_data_dir, f"nodal_loads_{combo_name.lower()}.json")
            if os.path.exists(nodal_load_file):
                with open(nodal_load_file, 'r') as f:
                    nodal_loads = json.load(f)
                    for node_id, loads in nodal_loads.items():
                        if node_id in nodes:  # Only process nodes at this z_level
                            components.append({
                                "wt": abs(loads[2]),
                                "xi": nodes[node_id][0],
                                "yi": nodes[node_id][1]
                            })
            
            # Process member loads with z-level filtering
            member_load_file = os.path.join(load_data_dir, f"member_loads_{combo_name.lower()}.json")
            if os.path.exists(member_load_file):
                with open(member_load_file, 'r') as f:
                    member_loads = json.load(f).get("element_loads", [])
                    for member in member_loads:
                        # Check if member exists at current z_level
                        start_z = member.get("start_node_coords", [0,0,0])[2]
                        end_z = member.get("end_node_coords", [0,0,0])[2]
                        
                        # Skip if member doesn't match current z_level
                        if not (abs(start_z - current_z) < 0.001 and abs(end_z - current_z) < 0.001):

                            continue
                            
                        # Process uniform loads
                        if member.get("uniform"):
                            uniform_z_load = member["uniform"][2]
                            member_length = member["member_length"]
                            equivalent_load = abs(uniform_z_load * member_length)
                            center = member.get("center", [0, 0, 0])
                            
                            # Verify center is at current z_level
                            if abs(center[2] - current_z) < 0.001:
                                components.append({
                                    "wt": equivalent_load,
                                    "xi": center[0],
                                    "yi": center[1]
                                })
                        
                        # Process point loads
                        if member.get("point") and member.get("point")[2] != 0:
                            point_z_load = member["point"][2]
                            point_coords = member.get("point_coords", [0, 0, 0])
                            
                            # Verify point is at current z_level
                            if abs(point_coords[2] - current_z) < 0.001:
                                components.append({
                                    "wt": abs(point_z_load),
                                    "xi": point_coords[0],
                                    "yi": point_coords[1]
                                })
            
            # Create elements (columns)
            for column_id, column_info in filtered_columns.get(z_level, {}).items():
                ele_tag = column_info[2]
                element = next(
                    (ele for ele in element_data["elements"] if ele["eleTag"] == ele_tag),
                    None
                )
                
                if element:
                    node_i_id = element["node_i_id"]
                    node_key = f"n{node_i_id}"
                    x, y, _ = nodes.get(node_key, [0, 0, 0])
                    
                    elem_type = element["type"].lower()
                    
                    if elem_type == "circular":
                        diameter = (element["area"] * 4 / 3.1416) ** 0.5
                        elements.append({
                            "name": column_id,
                            "type": "circular_column",
                            "diameter": round(diameter, 1),
                            "height": element["length"],
                            "x": x,
                            "y": y,
                            "rotation": element.get("rotation", 0.0)
                        })
                    elif elem_type == "rectangular":
                        elements.append({
                            "name": column_id,
                            "type": "rectangular_column",
                            "Dx": element.get("Dx", 12.0),
                            "Dy": element.get("Dy", 18.0),
                            "height": element["length"],
                            "x": x,
                            "y": y,
                            "rotation": element.get("rotation", 0.0)
                        })
                    elif elem_type in ["l_shaped", "l-shaped"]:
                        elements.append({
                            "name": column_id,
                            "type": "L_shaped_column",
                            "a": element.get("a", 12.0),
                            "b": element.get("b", 12.0),
                            "t": element.get("t", 4.0),
                            "height": element["length"],
                            "x": x,
                            "y": y,
                            "rotation": element.get("rotation", 0.0)
                        })
            
            # Store results
            result_key = f"{z_level}_{combo_name}"
            all_results[result_key] = {
                "components": components,
                "elements": elements,
                "Lx": Lx,
                "Ly": Ly,
                "load_combination": combo_name,
                "load_factors": combo_factors
            }
            
            # Save individual file
            filename = f"lateral_load_distribution_z_level_{z_level.replace('.', '_')}_{combo_name}.json"
            save_path = os.path.join(lateral_load_dir, filename)
            with open(save_path, 'w') as f:
                json.dump(all_results[result_key], f, indent=2)
    
    print(f"Processing complete. Results saved in {lateral_load_dir}")
    return all_results


def calculate_element_stiffness(element, E, poisson_ratio):
    """Calculate element stiffness properties with corrected wall handling"""
    x = element["x"]
    y = element["y"]
    story_height = element["height"]
    element_type = element.get("type", "column")
    rotation = element.get("rotation", 0.0)

    if story_height <= 1e-9:
        print(f"Warning: Element '{element.get('name', 'Unnamed')}' has zero or negligible height. Stiffness set to 0.")
        return 0.0, 0.0, 0.0

    if element_type == "wall":
        if "meshes" in element:
            mesh_data = element["meshes"]
            thickness = element.get("thickness", 0.5)

            if "rotation" not in element or element["rotation"] is None:
                first_mesh = next(iter(mesh_data.values()))
                vertices_2d = np.array(first_mesh["node_coordinates"])[:, [0, 2]]
                rotation = calculate_wall_rotation(vertices_2d)
                element["rotation"] = rotation

            mesh_results = calculate_mesh_based_wall_stiffness(
                mesh_data, E, poisson_ratio, thickness, story_height
            )

            Kx_local = mesh_results["Kx"]
            Ky_local = mesh_results["Kz"]
            cx = mesh_results["centroid_x"]
            cz = mesh_results["centroid_z"]

            element["x"] = cx
            element["y"] = y
            element["z_coordinate"] = cz

            element["mesh_info"] = {
                "total_area": mesh_results["total_area"],
                "total_Ix": mesh_results["total_Ix"],
                "total_Iy": mesh_results["total_Iy"],
                "wall_width_x": mesh_results["width_x"],
                "wall_height_z": mesh_results["height_z"],
                "effective_flexural_height": mesh_results["effective_flexural_height"]
            }

        else:
            vertices = np.array(element["vertices"])
            area, Ix_wall, Iy_wall, cx_wall, cy_wall = calculate_wall_properties(vertices)

            x_coords = vertices[:, 0]
            y_coords = vertices[:, 1]
            Dx = max(x_coords) - min(x_coords)
            Dy = max(y_coords) - min(y_coords)

            denom_x = 1 + 0.6 * (1 + poisson_ratio) * (Dx**2 / story_height**2)
            denom_y = 1 + 0.6 * (1 + poisson_ratio) * (Dy**2 / story_height**2)

            Kx_local = (3 * E * Iy_wall) / (story_height**3 * denom_x)
            Ky_local = (3 * E * Ix_wall) / (story_height**3 * denom_y)
            
            element["x"] = x + cx_wall
            element["y"] = y + cy_wall

    elif element_type == "circular_column":
        diameter = element["diameter"]
        I = math.pi * diameter**4 / 64
        Kx_local = Ky_local = (12 * E * I) / (story_height**3)

    elif element_type == "L_shaped_column":
        a = element["a"]
        b = element["b"]
        t = element["t"]

        A1 = a * t
        A2 = (b - t) * t
        total_area = A1 + A2

        y1 = t / 2
        x1 = a / 2
        y2 = (b - t) / 2 + t
        x2 = t / 2

        xc = (A1 * x1 + A2 * x2) / total_area
        yc = (A1 * y1 + A2 * y2) / total_area

        Ix1_c = (a * t**3) / 12
        Iy1_c = (t * a**3) / 12

        Ix2_c = (t * (b - t)**3) / 12
        Iy2_c = ((b - t) * t**3) / 12

        Ix = Ix1_c + A1 * (yc - y1)**2 + Ix2_c + A2 * (yc - y2)**2
        Iy = Iy1_c + A1 * (xc - x1)**2 + Iy2_c + A2 * (xc - x2)**2

        Kx_local = (12 * E * Iy) / (story_height**3)
        Ky_local = (12 * E * Ix) / (story_height**3)

    else:  # Rectangular column
        Dx = element["Dx"]
        Dy = element["Dy"]
        Ix = (Dx * Dy**3) / 12
        Iy = (Dy * Dx**3) / 12

        Kx_local = (12 * E * Iy) / (story_height**3)
        Ky_local = (12 * E * Ix) / (story_height**3)

    if abs(rotation) > 1e-6:
        Kxx, Kyy, Kxy = transform_stiffness_matrix(Kx_local, Ky_local, rotation)
    else:
        Kxx, Kyy, Kxy = Kx_local, Ky_local, 0.0

    return Kxx, Kyy, Kxy

def lateral_load_distribution_enhanced(Lx, Ly, E, poisson_ratio, fc, Vx_total, Vy_total, components, elements):
    """
    Enhanced lateral load distribution analysis.
    """
    # --- 1. Calculate Center of Mass (COM) ---
    total_wt = sum(comp["wt"] for comp in components)
    sum_wt_xi = sum(comp["wt"] * comp["xi"] for comp in components)
    sum_wt_yi = sum(comp["wt"] * comp["yi"] for comp in components)
    Xm = sum_wt_xi / total_wt if total_wt > 1e-10 else 0
    Ym = sum_wt_yi / total_wt if total_wt > 1e-10 else 0
    COM = (Xm, Ym)

    # --- 2. Calculate Element Stiffness Properties ---
    element_details = []
    
    sum_Kxx = 0
    sum_Kyy = 0
    sum_Kxy = 0

    sum_Kxx_yi = 0
    sum_Kyy_xi = 0
    sum_Kxy_xi = 0
    sum_Kxy_yi = 0
    sum_Kxx_y_sq = 0
    sum_Kyy_x_sq = 0
    sum_Kxy_xy = 0

    for element in elements:
        Kxx_i, Kyy_i, Kxy_i = calculate_element_stiffness(element, E, poisson_ratio)
        
        x_i = element["x"]
        y_i = element["y"]

        element_details.append({
            "name": element.get("name", "Unnamed"),
            "Kxx": Kxx_i, "Kyy": Kyy_i, "Kxy": Kxy_i,
            "x": x_i, "y": y_i,
            "rotation": element.get("rotation", 0.0)
        })

        sum_Kxx += Kxx_i
        sum_Kyy += Kyy_i
        sum_Kxy += Kxy_i
        
        sum_Kxx_yi += Kxx_i * y_i
        sum_Kyy_xi += Kyy_i * x_i
        sum_Kxy_xi += Kxy_i * x_i
        sum_Kxy_yi += Kxy_i * y_i
        sum_Kxx_y_sq += Kxx_i * y_i**2
        sum_Kyy_x_sq += Kyy_i * x_i**2
        sum_Kxy_xy += Kxy_i * x_i * y_i
    
    # --- 3. Calculate Center of Rigidity (COR) ---
    A = np.array([
        [sum_Kyy, -sum_Kxy],
        [-sum_Kxy, sum_Kxx]
    ])
    B = np.array([
        sum_Kyy_xi - sum_Kxy_yi,
        sum_Kxx_yi - sum_Kxy_xi
    ])
    
    if abs(np.linalg.det(A)) > 1e-10:
        Xr, Yr = np.linalg.solve(A, B)
    else:
        print("Warning: Stiffness matrix for COR calculation is singular. COR may be inaccurate.")
        Xr = (sum(elem['Kxx'] * elem['x'] for elem in element_details) / sum_Kxx) if sum_Kxx > 1e-10 else 0
        Yr = (sum(elem['Kyy'] * elem['y'] for elem in element_details) / sum_Kyy) if sum_Kyy > 1e-10 else 0

    COR = (Xr, Yr)

    # --- 4. Calculate Eccentricity ---
    ex = Xm - Xr
    ey = Ym - Yr
    min_ex = 0.05 * Lx
    min_ey = 0.05 * Ly
    
    ex_total = ex + (min_ex if ex >= 0 else -min_ex)
    ey_total = ey + (min_ey if ey >= 0 else -min_ey)
    
    # --- 5. Torsional Moment applied at COR ---
    T_applied = Vx_total * ey_total - Vy_total * ex_total

    # --- 6. Global Stiffness Matrix [K] and Solution for Displacements ---
    K_global = np.zeros((3, 3))
    Jr_total = 0

    for elem in element_details:
        dx, dy = elem["x"] - Xr, elem["y"] - Yr

        K_global[0, 0] += elem["Kxx"]
        K_global[0, 1] += elem["Kxy"]
        K_global[1, 0] += elem["Kxy"]
        K_global[1, 1] += elem["Kyy"]

        K_global[0, 2] += (elem["Kxx"] * dy - elem["Kxy"] * dx)
        K_global[1, 2] += (elem["Kxy"] * dy - elem["Kyy"] * dx)
        
        Jr_elem = elem["Kxx"] * dy**2 + elem["Kyy"] * dx**2 - 2 * elem["Kxy"] * dx * dy
        Jr_total += Jr_elem

    K_global[2, 0] = K_global[0, 2]
    K_global[2, 1] = K_global[1, 2]
    K_global[2, 2] = Jr_total

    V_applied = np.array([Vx_total, Vy_total, T_applied])

    try:
        displacements = solve(K_global, V_applied)
        delta_x, delta_y, theta_torsion = displacements[0], displacements[1], displacements[2]
    except np.linalg.LinAlgError:
        print("Error: Global stiffness matrix is singular. Cannot solve for displacements.")
        delta_x, delta_y, theta_torsion = 0, 0, 0

    # --- 7. Displacement-Based Force Distribution ---
    results = []
    sum_Fx_check, sum_Fy_check = 0, 0

    for elem in element_details:
        dx, dy = elem["x"] - Xr, elem["y"] - Yr

        Fx_trans = elem["Kxx"] * delta_x + elem["Kxy"] * delta_y
        Fy_trans = elem["Kxy"] * delta_x + elem["Kyy"] * delta_y

        Fx_torsion = elem["Kxx"] * (-dy * theta_torsion) + elem["Kxy"] * (dx * theta_torsion)
        Fy_torsion = elem["Kxy"] * (-dy * theta_torsion) + elem["Kyy"] * (dx * theta_torsion)
        
        Fx = Fx_trans + Fx_torsion
        Fy = Fy_trans + Fy_torsion

        results.append({
            "name": elem["name"],
            "Fx": Fx, "Fy": Fy,
            "Fx_direct_trans": elem["Kxx"] * delta_x,
            "Fy_direct_trans": elem["Kyy"] * delta_y,
            "Fx_coupling_trans": elem["Kxy"] * delta_y,
            "Fy_coupling_trans": elem["Kxy"] * delta_x,
            "Fx_torsion": Fx_torsion,
            "Fy_torsion": Fy_torsion,
            "Fx_trans": Fx_trans, # Added Fx_trans
            "Fy_trans": Fy_trans, # Added Fy_trans
            "dx_from_cor": dx,
            "dy_from_cor": dy,
            "Kxx": elem["Kxx"], "Kyy": elem["Kyy"], "Kxy": elem["Kxy"],
            "x": elem["x"], "y": elem["y"],
            "rotation": elem["rotation"]
        })
        sum_Fx_check += Fx
        sum_Fy_check += Fy

    return {
        "center_of_mass": COM,
        "center_of_rigidity": COR,
        "total_stiffness_matrix": K_global,
        "element_details": element_details,
        "eccentricity": (ex, ey),
        "eccentricity_total": (ex_total, ey_total),
        "torsional_moment_applied": T_applied,
        "torsional_rigidity_Jr_total": Jr_total,
        "shear_distribution": results,
        "total_shear_distributed": (sum_Fx_check, sum_Fy_check),
        "displacements": (delta_x, delta_y, theta_torsion)
    }


# def save_lateral_load_results(JSON_FOLDER, load_combinations, E, poisson_ratio, fc, Vx_total, Vy_total):
#     """
#     Calculates lateral load distribution for each z-level and load combination,
#     then saves results to organized file structure.
    
#     Args:
#         JSON_FOLDER: Path to folder containing input JSON files
#         load_combinations: Dictionary of load combinations
#         E: Elastic modulus
#         poisson_ratio: Poisson's ratio
#         fc: Concrete strength
#         Vx_total: Total lateral load in x-direction
#         Vy_total: Total lateral load in y-direction
#     """
#     # Create output directory
#     results_dir = os.path.join(JSON_FOLDER, "lateral_load_results")
#     os.makedirs(results_dir, exist_ok=True)
    
#     # First create components and elements for each z-level and load combo
#     all_components = create_components_and_elements(JSON_FOLDER, load_combinations)
    
#     # Process each z-level and load combination
#     for result_key, data in all_components.items():
#         z_level, combo_name = result_key.split("_", 1)
        
#         # Get the components and elements for this z-level and load combo
#         components = data["components"]
#         elements = data["elements"]
#         Lx = data["Lx"]
#         Ly = data["Ly"]
        
#         # Calculate lateral load distribution
#         results = lateral_load_distribution_enhanced(
#             Lx, Ly, E, poisson_ratio, fc, Vx_total, Vy_total, 
#             components, elements
#         )
        
#         # Create subfolder for this z-level
#         z_level_dir = os.path.join(results_dir, f"z_level_{z_level.replace('.', '_')}")
#         os.makedirs(z_level_dir, exist_ok=True)
        
#         # Save results for this load combination
#         output_file = os.path.join(z_level_dir, f"results_{combo_name}.json")
        
#         # Convert numpy arrays to lists for JSON serialization
#         results_serializable = {
#             "center_of_mass": results["center_of_mass"],
#             "center_of_rigidity": results["center_of_rigidity"],
#             "total_stiffness_matrix": results["total_stiffness_matrix"].tolist(),
#             "element_details": results["element_details"],
#             "eccentricity": results["eccentricity"],
#             "eccentricity_total": results["eccentricity_total"],
#             "torsional_moment_applied": results["torsional_moment_applied"],
#             "torsional_rigidity_Jr_total": results["torsional_rigidity_Jr_total"],
#             "shear_distribution": results["shear_distribution"],
#             "total_shear_distributed": results["total_shear_distributed"],
#             "displacements": results["displacements"],
#             "load_combination": combo_name,
#             "z_level": z_level,
#             "Lx": Lx,
#             "Ly": Ly,
#             "Vx_total": Vx_total,
#             "Vy_total": Vy_total
#         }
        
#         with open(output_file, 'w') as f:
#             json.dump(results_serializable, f, indent=2)
    
#     print(f"All lateral load results saved in {results_dir}")
#     return results_dir


import os
import json

def save_lateral_load_results(JSON_FOLDER, load_combinations, E, poisson_ratio, fc, Vx_total, Vy_total):
    """
    Calculates lateral load distribution for each z-level and load combination,
    then saves results to organized file structure in TXT format.
    
    Args:
        JSON_FOLDER: Path to folder containing input JSON files
        load_combinations: Dictionary of load combinations
        E: Elastic modulus
        poisson_ratio: Poisson's ratio
        fc: Concrete strength
        Vx_total: Total lateral load in x-direction
        Vy_total: Total lateral load in y-direction
    """
    # Create output directory
    results_dir = os.path.join(JSON_FOLDER, "lateral_load_results")
    os.makedirs(results_dir, exist_ok=True)
    
    # First create components and elements for each z-level and load combo
    all_components = create_components_and_elements(JSON_FOLDER, load_combinations)
    
    # Process each z-level and load combination
    for result_key, data in all_components.items():
        z_level, combo_name = result_key.split("_", 1)
        
        # Get the components and elements for this z-level and load combo
        components = data["components"]
        elements = data["elements"]
        Lx = data["Lx"]
        Ly = data["Ly"]
        
        # Calculate lateral load distribution
        results = lateral_load_distribution_enhanced(
            Lx, Ly, E, poisson_ratio, fc, Vx_total, Vy_total, 
            components, elements
        )
        
        # Create subfolder for this z-level
        z_level_dir = os.path.join(results_dir, f"z_level_{z_level.replace('.', '_')}")
        os.makedirs(z_level_dir, exist_ok=True)
        
        # Save results for this load combination as TXT
        output_file = os.path.join(z_level_dir, f"results_{combo_name}.txt")
        
        # Format the results for TXT output
        with open(output_file, 'w') as f:
            # Write header information
            f.write(f"LATERAL LOAD ANALYSIS RESULTS\n")
            f.write(f"="*50 + "\n")
            f.write(f"Load Combination: {combo_name}\n")
            f.write(f"Z-Level: {z_level}\n")
            f.write(f"Building Dimensions: Lx = {Lx}, Ly = {Ly}\n")
            f.write(f"Total Loads: Vx = {Vx_total}, Vy = {Vy_total}\n\n")
            
            # Write centers information
            f.write("CENTERS\n")
            f.write("-"*50 + "\n")
            f.write(f"Center of Mass: X = {results['center_of_mass'][0]:.3f}, Y = {results['center_of_mass'][1]:.3f}\n")
            f.write(f"Center of Rigidity: X = {results['center_of_rigidity'][0]:.3f}, Y = {results['center_of_rigidity'][1]:.3f}\n\n")
            
            # Write eccentricities
            f.write("ECCENTRICITIES\n")
            f.write("-"*50 + "\n")
            f.write(f"Design Eccentricity: X = {results['eccentricity'][0]:.3f}, Y = {results['eccentricity'][1]:.3f}\n")
            f.write(f"Total Eccentricity: X = {results['eccentricity_total'][0]:.3f}, Y = {results['eccentricity_total'][1]:.3f}\n\n")
            
            # Write stiffness matrix
            f.write("TOTAL STIFFNESS MATRIX\n")
            f.write("-"*50 + "\n")
            stiffness_matrix = results["total_stiffness_matrix"]
            f.write(f"X-axis: {stiffness_matrix[0][0]:.2f}, {stiffness_matrix[0][1]:.2f}, {stiffness_matrix[0][2]:.2e}\n")
            f.write(f"Y-axis: {stiffness_matrix[1][0]:.2f}, {stiffness_matrix[1][1]:.2f}, {stiffness_matrix[1][2]:.2e}\n")
            f.write(f"Rotation: {stiffness_matrix[2][0]:.2e}, {stiffness_matrix[2][1]:.2e}, {stiffness_matrix[2][2]:.2f}\n\n")
            
            # Write torsional properties
            f.write("TORSIONAL PROPERTIES\n")
            f.write("-"*50 + "\n")
            f.write(f"Applied Torsional Moment: {results['torsional_moment_applied']:.2f}\n")
            f.write(f"Total Torsional Rigidity (Jr): {results['torsional_rigidity_Jr_total']:.2f}\n\n")
            
            # Write displacements
            f.write("DISPLACEMENTS\n")
            f.write("-"*50 + "\n")
            f.write(f"X translation: {results['displacements'][0]:.6f}\n")
            f.write(f"Y translation: {results['displacements'][1]:.6f}\n")
            f.write(f"Rotation: {results['displacements'][2]:.6e}\n\n")
            
            # Write total shear distribution
            f.write("TOTAL SHEAR DISTRIBUTION\n")
            f.write("-"*50 + "\n")
            f.write(f"X-direction: {results['total_shear_distributed'][0]:.2f}\n")
            f.write(f"Y-direction: {results['total_shear_distributed'][1]:.2f}\n\n")
            
            # Write element details
            f.write("ELEMENT DETAILS\n")
            f.write("-"*50 + "\n")
            for element in results["element_details"]:
                f.write(f"Element {element['name']}:\n")
                f.write(f"  Position: X = {element['x']}, Y = {element['y']}, Rotation = {element['rotation']}\n")
                f.write(f"  Stiffness: Kxx = {element['Kxx']:.2f}, Kyy = {element['Kyy']:.2f}, Kxy = {element['Kxy']:.2f}\n")
            f.write("\n")
            
            # Write shear distribution
            f.write("SHEAR DISTRIBUTION BY ELEMENT\n")
            f.write("-"*50 + "\n")
            for shear in results["shear_distribution"]:
                f.write(f"Element {shear['name']}:\n")
                f.write(f"  Total Forces: Fx = {shear['Fx']:.2f}, Fy = {shear['Fy']:.2f}\n")
                f.write(f"  Direct Translation: Fx = {shear['Fx_direct_trans']:.2f}, Fy = {shear['Fy_direct_trans']:.2f}\n")
                f.write(f"  Torsional Components: Fx = {shear['Fx_torsion']:.2f}, Fy = {shear['Fy_torsion']:.2f}\n")
                f.write(f"  Distance from COR: dx = {shear['dx_from_cor']:.2f}, dy = {shear['dy_from_cor']:.2f}\n\n")
    
    print(f"All lateral load results saved in TXT format in {results_dir}")
    return results_dir

# def create_components_and_elements(JSON_FOLDER, load_combinations):
#     """
#     Creates components and elements for each z_level using JSON files from the specified folder,
#     including load values from different load combinations, and saves separate files for each combination.
    
#     Args:
#         JSON_FOLDER (str): Path to the folder containing the required JSON files
#         load_combinations (dict): Dictionary of load combinations
#     Returns:
#         dict: A dictionary containing components and elements for each z_level and combination
#     """
#     # Define file paths based on the JSON_FOLDER
#     filtered_nodes_path = os.path.join(JSON_FOLDER, "filtered_nodes.json")
#     filtered_columns_path = os.path.join(JSON_FOLDER, "filtered_columns.json")
#     element_data_path = os.path.join(JSON_FOLDER, "element_data.json")
#     load_data_dir = os.path.join(JSON_FOLDER, "load_data")
#     lateral_load_dir = os.path.join(JSON_FOLDER, "lateral_load_distribution")
    
#     # Create output directories if they don't exist
#     os.makedirs(lateral_load_dir, exist_ok=True)
    
#     # Load the JSON files
#     try:
#         with open(filtered_nodes_path, 'r') as f:
#             filtered_nodes = json.load(f)
        
#         with open(filtered_columns_path, 'r') as f:
#             filtered_columns = json.load(f)
        
#         with open(element_data_path, 'r') as f:
#             element_data = json.load(f)
#     except FileNotFoundError as e:
#         raise FileNotFoundError(f"Required JSON file not found: {e}")

#     # Dictionary to store load values for each node and combination
#     node_loads = defaultdict(dict)
    
#     # Process each load combination to get nodal loads
#     for combo_name, _ in load_combinations.items():
#         load_file = os.path.join(load_data_dir, f"nodal_loads_{combo_name.lower()}.json")
#         try:
#             with open(load_file, 'r') as f:
#                 combo_loads = json.load(f)
#                 for node_id, loads in combo_loads.items():
#                     node_loads[node_id][combo_name] = abs(loads[2])  # Using absolute Z-load
#         except FileNotFoundError:
#             print(f"Warning: Load file not found for combination {combo_name}")
#             continue
#         member_load_file = os.path.join(load_data_dir, f"member_load_{combo_name.lower()}.json")
#         # try:
#         #     with open(member_load_file, 'r') as f:
#         #         member_load_file = json.load(f)
#         #         for node_id, loads in member_load_file.items():
#         #             node_loads[node_id][combo_name] = abs(loads[2])  # Using absolute Z-load
#         # except FileNotFoundError:
#         #     print(f"Warning: Load file not found for combination {combo_name}")
#         #     continue

#     # Process each z_level and combination
#     all_results = {}
#     for z_level, nodes in filtered_nodes.items():
#         # Calculate Lx and Ly for this z_level
#         x_coords = [coords[0] for coords in nodes.values()]
#         y_coords = [coords[1] for coords in nodes.values()]
#         Lx = max(x_coords) - min(x_coords)
#         Ly = max(y_coords) - min(y_coords)
        
#         # Process each load combination separately
#         for combo_name, combo_factors in load_combinations.items():
#             components = []
#             elements = []
            
#             # Create components with loads for this combination
#             for node_id, coords in nodes.items():
#                 wt = node_loads.get(node_id, {}).get(combo_name, 25.898)  # Default if not found
                
#                 components.append({
#                     "wt": wt,
#                     "xi": coords[0],
#                     "yi": coords[1]
#                 })
            
#             # Create elements (same for all combinations)
#             for column_id, column_info in filtered_columns.get(z_level, {}).items():
#                 ele_tag = column_info[2]
#                 element = next(
#                     (ele for ele in element_data["elements"] if ele["eleTag"] == ele_tag),
#                     None
#                 )
                
#                 if element:
#                     node_i_id = element["node_i_id"]
#                     node_key = f"n{node_i_id}"
#                     x, y, _ = nodes.get(node_key, [0, 0, 0])
                    
#                     elem_type = element["type"].lower()
                    
#                     if elem_type == "circular":
#                         diameter = (element["area"] * 4 / 3.1416) ** 0.5
#                         elements.append({
#                             "name": column_id,
#                             "type": "circular_column",
#                             "diameter": round(diameter, 1),
#                             "height": element["length"],
#                             "x": x,
#                             "y": y,
#                             "rotation": element.get("rotation", 0.0)
#                         })
#                     elif elem_type == "rectangular":
#                         elements.append({
#                             "name": column_id,
#                             "type": "rectangular_column",
#                             "Dx": element.get("Dx", 12.0),
#                             "Dy": element.get("Dy", 18.0),
#                             "height": element["length"],
#                             "x": x,
#                             "y": y,
#                             "rotation": element.get("rotation", 0.0)
#                         })
#                     elif elem_type in ["l_shaped", "l-shaped"]:
#                         elements.append({
#                             "name": column_id,
#                             "type": "L_shaped_column",
#                             "a": element.get("a", 12.0),
#                             "b": element.get("b", 12.0),
#                             "t": element.get("t", 4.0),
#                             "height": element["length"],
#                             "x": x,
#                             "y": y,
#                             "rotation": element.get("rotation", 0.0)
#                         })
            
#             # Store results for this z_level and combination
#             result_key = f"{z_level}_{combo_name}"
#             all_results[result_key] = {
#                 "components": components,
#                 "elements": elements,
#                 "Lx": Lx,
#                 "Ly": Ly,
#             }
            
#             # Save individual file for this z_level and combination
#             filename = f"lateral_load_distribution_z_level_{z_level.replace('.', '_')}_{combo_name}.json"
#             save_path = os.path.join(lateral_load_dir, filename)
            
#             with open(save_path, 'w') as f:
#                 json.dump(all_results[result_key], f, indent=2)
    
#     print(f"Processing complete. Results saved in {lateral_load_dir}")
#     return all_results


# # Example usage with mesh-based walls
# components = [
#     {"wt": 629.4652, "xi": 314.9608, "yi": 246.063},
#     {"wt": -18.88396, "xi": 551.1814, "yi": 462.598},
#     {"wt": 25.898, "xi": 5.905515, "yi": 236.2206},
#     {"wt": 25.898, "xi": 314.9608, "yi": 236.2206},
#     {"wt": 25.898, "xi": 623.9961, "yi": 236.2206},
#     {"wt": 25.898, "xi": 78.7402, "yi": 5.905515},
#     {"wt": 25.898, "xi": 551.1814, "yi": 5.905515},
#     {"wt": 25.898, "xi": 78.7402, "yi": 486.2206},
#     {"wt": 25.898, "xi": 551.1814, "yi": 427.1655},
# ]

# wall_mesh_data = {
#     "R34001": {
#         "node_coordinates": [
#             [10.0*12, 0.0*12, 6.666666666666666*12],
#             [6.666666666666666*12, 0.0*12, 6.666666666666666*12],
#             [6.666666666666666*12, 0.0*12, 10.0*12],
#             [10.0*12, 0.0*12, 10.0*12]
#         ]
#     },
#     "R34002": {
#         "node_coordinates": [
#             [10.0*12, 0.0, 3.3333333333333326*12],
#             [6.666666666666666*12, 0.0, 3.3333333333333326*12],
#             [6.666666666666666*12, 0.0, 6.666666666666666*12],
#             [10.0*12, 0.0, 6.666666666666666*12]
#         ]
#     },
#     "R34003": {
#         "node_coordinates": [
#             [10.0*12, 0.0, 0.0],
#             [6.666666666666666*12, 0.0, 0.0],
#             [6.666666666666666*12, 0.0, 3.3333333333333326*12],
#             [10.0*12, 0.0, 3.3333333333333326*12]
#         ]
#     }
# }
#     # {
#     #     "name": "W2-1",
#     #     "type": "wall",
#     #     "meshes": wall_mesh_data,
#     #     "height": 157.4804,
#     #     "x": 0.0,
#     #     "y": 0.0,
#     #     "thickness": 0.5,
#     #     "rotation": None
#     # },

# Lx = 629.9216
# Ly = 492.126


# elements = [

#     {
#         "name": "C1",
#         "type": "circular_column",
#         "diameter": 15.0,
#         "height": 157.4804,
#         "x": 200.0,
#         "y": 200.0,
#         "rotation": 0.0
#     },
#     {
#         "name": "L1",
#         "type": "L_shaped_column",
#         "a": 12.0,
#         "b": 12.0,
#         "t": 4.0,
#         "height": 157.4804,
#         "x": 400.0,
#         "y": 300.0,
#         "rotation": 45.0
#     },
#     {
#         "name": "R1",
#         "type": "rectangular_column",
#         "Dx": 18.0,
#         "Dy": 12.0,
#         "height": 157.4804,
#         "x": 300.0,
#         "y": 150.0,
#         "rotation": 22.5
#     },
#     {
#         "name": "R2",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": 50.0,
#         "y": 50.0,
#         "rotation": 0.0
#     },
#     {
#         "name": "R3",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": 50.0,
#         "y": Ly - 50.0,
#         "rotation": 0.0
#     },
#     {
#         "name": "R4",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": Lx - 50.0,
#         "y": 50.0,
#         "rotation": 0.0
#     },
#     {
#         "name": "R5",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": Lx - 50.0,
#         "y": Ly - 50.0,
#         "rotation": 0.0
#     },
#     {
#         "name": "R6",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": Lx/2,
#         "y": 50.0,
#         "rotation": 0.0
#     },
#     {
#         "name": "R7",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": Lx/2,
#         "y": Ly - 50.0,
#         "rotation": 0.0
#     },
#     {
#         "name": "R8",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": 50.0,
#         "y": Ly/2,
#         "rotation": 0.0
#     },
#     {
#         "name": "R9",
#         "type": "rectangular_column",
#         "Dx": 12.0,
#         "Dy": 18.0,
#         "height": 157.4804,
#         "x": Lx - 50.0,
#         "y": Ly/2,
#         "rotation": 0.0
#     }
# ]

# fc = 4.35
# E = 57000 * math.sqrt(fc)
# poisson_ratio = 0.25
# Vx = -100.0
# Vy = -0.0

# results = lateral_load_distribution_enhanced(Lx, Ly, E, poisson_ratio, fc, Vx, Vy, components, elements)

# print("Element   Kxx         Kyy         Kxy    Total Fx    Total Fy    Fx_Trans    Fy_Trans    Fx_Torsion  Fy_Torsion")
# print("------------------------------------------------------------------------------------------------------------------------------------------------------")
# for item in results['shear_distribution']:
#     print(f"{item['name']:<9} {item['Kxx']:<11.3f} {item['Kyy']:<11.3f} {item['Kxy']:<11.3f} {item['Fx']:<11.3f} {item['Fy']:<11.3f} {item['Fx_trans']:<11.3f} {item['Fy_trans']:<11.3f} {item['Fx_torsion']:<11.3f} {item['Fy_torsion']:<11.3f}")

# print("\n--- Summary ---")
# print(f"Center of Mass (COM): X = {results['center_of_mass'][0]:.3f} in, Y = {results['center_of_mass'][1]:.3f} in")
# print(f"Center of Rigidity (COR): X = {results['center_of_rigidity'][0]:.3f} in, Y = {results['center_of_rigidity'][1]:.3f} in")
# print(f"Applied Torsional Moment: {results['torsional_moment_applied']:.3f} kip·in")
# print(f"Global Displacements: Delta_x = {results['displacements'][0]:.6f} in, Delta_y = {results['displacements'][1]:.6f} in, Theta = {results['displacements'][2]:.8f} radians")
# print(f"Total Distributed Shear: Fx = {results['total_shear_distributed'][0]:.3f} kips, Fy = {results['total_shear_distributed'][1]:.3f} kips")
# print(f"Force Balance Check: Fx = {(Vx - results['total_shear_distributed'][0]):.6f}, Fy = {(Vy - results['total_shear_distributed'][1]):.6f}")

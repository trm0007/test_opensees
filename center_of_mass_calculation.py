
import math

def polygon_area(vertices):
    n = len(vertices)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += vertices[i][0] * vertices[j][1]
        area -= vertices[j][0] * vertices[i][1]
    return abs(area) / 2.0

def polygon_centroid(vertices):
    n = len(vertices)
    cx, cy = 0.0, 0.0
    area = polygon_area(vertices)
    if area == 0:
        return (0.0, 0.0)
    for i in range(n):
        j = (i + 1) % n
        factor = vertices[i][0] * vertices[j][1] - vertices[j][0] * vertices[i][1]
        cx += (vertices[i][0] + vertices[j][0]) * factor
        cy += (vertices[i][1] + vertices[j][1]) * factor
    cx /= (6 * area)
    cy /= (6 * area)
    return (cx, cy)

def calculate_center_of_mass(slab_vertices, slab_thickness, columns, walls, density, H):
    slab_area = polygon_area(slab_vertices)
    slab_volume = slab_area * slab_thickness
    slab_mass = slab_volume * density
    slab_cx, slab_cy = polygon_centroid(slab_vertices)

    column_cm_x, column_cm_y, total_column_mass = 0.0, 0.0, 0.0
    for col in columns:
        b, h = col['b'], col['h']
        volume = b * h * H
        mass = volume * density
        total_column_mass += mass
        column_cm_x += mass * col['x']
        column_cm_y += mass * col['y']

    wall_cm_x, wall_cm_y, total_wall_mass = 0.0, 0.0, 0.0
    for wall in walls:
        x1, y1, x2, y2 = wall['x1'], wall['y1'], wall['x2'], wall['y2']
        t = wall['t']
        length = math.hypot(x2 - x1, y2 - y1)
        volume = length * t * H
        mass = volume * density
        total_wall_mass += mass
        centroid_x = (x1 + x2) / 2
        centroid_y = (y1 + y2) / 2
        wall_cm_x += mass * centroid_x
        wall_cm_y += mass * centroid_y

    total_mass = slab_mass + total_column_mass + total_wall_mass
    if total_mass == 0:
        return (0.0, 0.0)
    com_x = (slab_mass * slab_cx + column_cm_x + wall_cm_x) / total_mass
    com_y = (slab_mass * slab_cy + column_cm_y + wall_cm_y) / total_mass

    
    return (com_x, com_y)

def calculate_center_of_rigidity(columns, walls, E, H):
    sum_ky, sum_ky_xi, sum_kx, sum_kx_yi = 0.0, 0.0, 0.0, 0.0

    for col in columns:
        x, y, b, h = col['x'], col['y'], col['b'], col['h']
        Ix = (b * h**3) / 12
        Iy = (h * b**3) / 12
        ky_i = (12 * E * Ix) / H**3
        kx_i = (12 * E * Iy) / H**3
        sum_ky += ky_i
        sum_ky_xi += ky_i * x
        sum_kx += kx_i
        sum_kx_yi += kx_i * y

    for wall in walls:
        x1, y1, x2, y2 = wall['x1'], wall['y1'], wall['x2'], wall['y2']
        t = wall['t']
        dx, dy = x2 - x1, y2 - y1
        length = math.hypot(dx, dy)
        if length == 0:
            continue
        theta = math.atan2(dy, dx)
        sin_t, cos_t = math.sin(theta), math.cos(theta)
        I_lx = (t * length**3) / 12
        I_ly = (length * t**3) / 12
        Ix_global = I_lx * sin_t**2 + I_ly * cos_t**2
        Iy_global = I_lx * cos_t**2 + I_ly * sin_t**2
        ky_i = (12 * E * Ix_global) / H**3
        kx_i = (12 * E * Iy_global) / H**3
        centroid_x = (x1 + x2) / 2
        centroid_y = (y1 + y2) / 2
        sum_ky += ky_i
        sum_ky_xi += ky_i * centroid_x
        sum_kx += kx_i
        sum_kx_yi += kx_i * centroid_y

    CRx = sum_ky_xi / sum_ky if sum_ky != 0 else 0.0
    CRy = sum_kx_yi / sum_kx if sum_kx != 0 else 0.0
    return (CRx, CRy)

# Example usage
if __name__ == "__main__":
    slab_vertices = [(0, 0), (10, 0), (10, 5), (5, 5), (0, 5)]
    slab_thickness = 0.2  # meters
    columns = [
        {'x': 2, 'y': 2, 'b': 0.4, 'h': 0.4},
        {'x': 8, 'y': 2, 'b': 0.4, 'h': 0.4},
        {'x': 5, 'y': 4, 'b': 0.4, 'h': 0.4},
    ]
    walls = [
        {'x1': 0, 'y1': 0, 'x2': 0, 'y2': 5, 't': 0.2},
        {'x1': 10, 'y1': 0, 'x2': 10, 'y2': 5, 't': 0.2},
    ]
    E = 25e9  # Pa (Modulus of Elasticity for concrete)
    density = 2400  # kg/m³
    H = 3  # meters

    com = calculate_center_of_mass(slab_vertices, slab_thickness, columns, walls, density, H)
    cor = calculate_center_of_rigidity(columns, walls, E, H)

    # print(f"Center of Mass (CoM): ({com[0]:.2f}, {com[1]:.2f}) meters")
    # print(f"Center of Rigidity (CoR): ({cor[0]:.2f}, {cor[1]:.2f}) meters")


import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def plot_structure(slab_vertices, columns, walls, com, cor):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_aspect('equal')
    
    # Plot slab
    slab_polygon = patches.Polygon(slab_vertices, closed=True, 
                                  edgecolor='black', facecolor='lightgray', 
                                  label='Slab')
    ax.add_patch(slab_polygon)
    
    # Plot columns
    for col in columns:
        x, y, b, h = col['x'], col['y'], col['b'], col['h']
        rect = patches.Rectangle((x - b/2, y - h/2), b, h,
                                edgecolor='blue', facecolor='skyblue', 
                                label='Column')
        ax.add_patch(rect)
    
    # Plot walls
    for wall in walls:
        x1, y1, x2, y2 = wall['x1'], wall['y1'], wall['x2'], wall['y2']
        t = wall['t']
        
        # Calculate wall polygon coordinates
        dx = x2 - x1
        dy = y2 - y1
        length = np.hypot(dx, dy)
        if length == 0:
            continue
            
        # Perpendicular direction components
        px = -dy * t / (2 * length)
        py = dx * t / (2 * length)
        
        # Create wall polygon
        wall_vertices = [
            (x1 + px, y1 + py),
            (x1 - px, y1 - py),
            (x2 - px, y2 - py),
            (x2 + px, y2 + py)
        ]
        
        wall_polygon = patches.Polygon(wall_vertices, closed=True,
                                      edgecolor='darkgreen', facecolor='limegreen',
                                      label='Wall')
        ax.add_patch(wall_polygon)
    
    # Plot centers
    ax.scatter(*com, c='red', s=100, marker='*', label='Center of Mass (CoM)')
    ax.scatter(*cor, c='orange', s=100, marker='X', label='Center of Rigidity (CoR)')
    
    # Annotate coordinates
    ax.annotate(f'CoM: ({com[0]:.2f}, {com[1]:.2f})', com, 
                textcoords="offset points", xytext=(15,5))
    ax.annotate(f'CoR: ({cor[0]:.2f}, {cor[1]:.2f})', cor, 
                textcoords="offset points", xytext=(15,-15))
    
    # Configure plot
    ax.set_xlabel('X Position (m)')
    ax.set_ylabel('Y Position (m)')
    ax.set_title('Structural Layout with CoM and CoR')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # Remove duplicate labels
    ax.legend(by_label.values(), by_label.keys(), loc='upper right', bbox_to_anchor=(1.15, 1.15))
    ax.grid(True)
    plt.autoscale()
    plt.show()

# Example usage with previous data
if __name__ == "__main__":
    slab_vertices = [(0, 0), (10, 0), (10, 5), (5, 5), (0, 5)]
    columns = [
        {'x': 2, 'y': 2, 'b': 0.4, 'h': 0.4},
        {'x': 8, 'y': 2, 'b': 0.4, 'h': 0.4},
        {'x': 5, 'y': 4, 'b': 0.4, 'h': 0.4},
    ]
    walls = [
        {'x1': 0, 'y1': 0, 'x2': 0, 'y2': 5, 't': 0.2},
        {'x1': 10, 'y1': 0, 'x2': 10, 'y2': 5, 't': 0.2},
    ]
    
    # Use previously calculated CoM and CoR
    com = (4.35, 2.81)  # Example values from previous calculation
    cor = (5.00, 2.50)  # Example values from previous calculation
    
    plot_structure(slab_vertices, columns, walls, com, cor)


    # ... [Keep all previous functions unchanged] ...

if __name__ == "__main__":
    # Common parameters for all examples
    E = 25e9  # Pa
    density = 2400  # kg/m³
    H = 3  # meters
    slab_thickness = 0.2  # meters

    # ========================================================================
    # Example 1: Asymmetric Hexagonal Structure
    # ========================================================================
    slab_vertices_1 = [
        (0, 0), (12, 0), (15, 5), 
        (12, 10), (5, 10), (0, 5)
    ]
    
    columns_1 = [
        {'x': 2, 'y': 2, 'b': 0.4, 'h': 0.4},
        {'x': 8, 'y': 2, 'b': 0.5, 'h': 0.6},  # Larger column
        {'x': 13, 'y': 5, 'b': 0.4, 'h': 0.4},
        {'x': 6, 'y': 9, 'b': 0.5, 'h': 0.5},
    ]
    
    walls_1 = [
        {'x1': 10, 'y1': 0, 'x2': 15, 'y2': 5, 't': 0.3},
        {'x1': 5, 'y1': 10, 'x2': 12, 'y2': 10, 't': 0.25},
        {'x1': 0, 'y1': 5, 'x2': 5, 'y2': 10, 't': 0.2},
    ]
    
    com1 = calculate_center_of_mass(slab_vertices_1, slab_thickness, 
                                   columns_1, walls_1, density, H)
    cor1 = calculate_center_of_rigidity(columns_1, walls_1, E, H)
    
    plot_structure(slab_vertices_1, columns_1, walls_1, com1, cor1)

    # ========================================================================
    # Example 2: Circular Structure with Radial Elements
    # ========================================================================
    # Generate circular slab (12-sided polygon)
    n_sides = 12
    radius = 8
    center = (10, 10)
    angles = [i*(2*np.pi/n_sides) for i in range(n_sides)]
    slab_vertices_2 = [(center[0] + radius*np.cos(a), 
                       center[1] + radius*np.sin(a)) 
                      for a in angles]
    
    columns_2 = [
        {'x': center[0] + 6*np.cos(angle), 
         'y': center[1] + 6*np.sin(angle), 
         'b': 0.4, 'h': 0.4}
        for angle in angles[::3]
    ] + [{'x': center[0], 'y': center[1], 'b': 1.0, 'h': 1.0}]
    
    walls_2 = [
        # Radial walls
        {'x1': center[0], 'y1': center[1], 
         'x2': center[0]+radius*np.cos(a), 
         'y2': center[1]+radius*np.sin(a), 't': 0.2}
        for a in angles[::4]
    ] + [
        # Chord walls
        {'x1': center[0]+7*np.cos(angles[i]), 
         'y1': center[1]+7*np.sin(angles[i]),
         'x2': center[0]+7*np.cos(angles[i+2]), 
         'y2': center[1]+7*np.sin(angles[i+2]), 't': 0.15}
        for i in range(0, n_sides, 3)
    ]
    
    com2 = calculate_center_of_mass(slab_vertices_2, slab_thickness,
                                   columns_2, walls_2, density, H)
    cor2 = calculate_center_of_rigidity(columns_2, walls_2, E, H)
    
    plot_structure(slab_vertices_2, columns_2, walls_2, com2, cor2)

    # ========================================================================
    # Example 3: L-Shaped Structure with Stiffness Concentration
    # ========================================================================
    slab_vertices_3 = [
        (0, 0), (10, 0), (10, 5), 
        (5, 5), (5, 15), (0, 15)
    ]
    
    columns_3 = [
        {'x': 2, 'y': 2, 'b': 0.4, 'h': 0.4},
        {'x': 2, 'y': 8, 'b': 0.6, 'h': 0.6},
        {'x': 2, 'y': 12, 'b': 0.5, 'h': 0.5},
        {'x': 8, 'y': 2, 'b': 0.4, 'h': 0.4},
        {'x': 8, 'y': 12, 'b': 0.4, 'h': 0.4},
    ]
    
    walls_3 = [
        {'x1': 10, 'y1': 0, 'x2': 10, 'y2': 5, 't': 0.3},
        {'x1': 5, 'y1': 5, 'x2': 5, 'y2': 15, 't': 0.25},
        {'x1': 0, 'y1': 15, 'x2': 5, 'y2': 15, 't': 0.2},
        {'x1': 7, 'y1': 10, 'x2': 3, 'y2': 14, 't': 0.2},
    ]
    
    com3 = calculate_center_of_mass(slab_vertices_3, slab_thickness,
                                   columns_3, walls_3, density, H)
    cor3 = calculate_center_of_rigidity(columns_3, walls_3, E, H)
    
    plot_structure(slab_vertices_3, columns_3, walls_3, com3, cor3)
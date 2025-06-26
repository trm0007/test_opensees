from import_ import *     
from input import *
from load_combinations import *
from Grid_and_structure_creation import *
from wall_meshing import *
from sections_function import *
from units import *
import matplotlib.pyplot as plt
import opsvis as opsv
from ipywidgets import widgets


def apply_fixity_conditions(fixity_data, JSON_FOLDER):
    """Applies boundary conditions, supporting both node IDs and names."""
    nodes_data, _, _ = create_combined_structure_json(JSON_FOLDER)
    
    # Create mappings for both IDs and names
    existing_ids = {node["id"] for node in nodes_data}  # Set of all valid IDs
    name_to_id = {node["name"]: node["id"] for node in nodes_data}  # Dict mapping names to IDs

    for group_name, group in fixity_data.items():
        for node_ref in group["nodes"]:
            # Case 1: Node reference is an ID (integer)
            if isinstance(node_ref, int):
                if node_ref in existing_ids:
                    ops.fix(node_ref, *group['vals'])
                    print(f"Fixed node ID {node_ref}")
                else:
                    print(f"Warning: Node ID {node_ref} not found (skipping)")
            
            # Case 2: Node reference is a name (string, e.g., "n1")
            elif isinstance(node_ref, str):
                if node_ref in name_to_id:
                    node_id = name_to_id[node_ref]
                    ops.fix(node_id, *group['vals'])
                    print(f"Fixed node '{node_ref}' (ID: {node_id})")
                else:
                    print(f"Warning: Node name '{node_ref}' not found (skipping)")
            
            # Invalid type
            else:
                print(f"Error: Invalid node reference {node_ref} (must be int or str)")
                
# Extruded shapes visualization - Fixed version

# Updated visualization function with better error handling

# =============================================
# Section Visualization Setup
# =============================================


def define_rc_section10_opensees(sec_tag, core_tag, cover_tag, steel_tag, H, B, cover_H, cover_B, offset,
                    n_bars_top, dia_top, n_bars_bot, dia_bot,
                    n_bars_secondary_top, dia_sec_top, n_bars_secondary_bot, dia_sec_bot,
                    n_bars_int, dia_int):
    """
    Draw a rectangular RC section with visualization.
    
    Parameters:
    -----------
    sec_tag : int
        Tag for the section
    core_tag : int
        Material tag for core concrete
    cover_tag : int
        Material tag for cover concrete
    steel_tag : int
        Material tag for reinforcement steel
    H : float
        Height of the section
    B : float
        Width of the section
    cover_H : float
        Cover thickness in H direction
    cover_B : float
        Cover thickness in B direction
    offset : float
        Offset for secondary reinforcement from primary reinforcement
    n_bars_top : int
        Number of primary reinforcement bars at top
    dia_top : float
        Diameter of primary reinforcement bars at top
    n_bars_bot : int
        Number of primary reinforcement bars at bottom
    dia_bot : float
        Diameter of primary reinforcement bars at bottom
    n_bars_secondary_top : int
        Number of secondary reinforcement bars at top
    dia_sec_top : float
        Diameter of secondary reinforcement bars at top
    n_bars_secondary_bot : int
        Number of secondary reinforcement bars at bottom
    dia_sec_bot : float
        Diameter of secondary reinforcement bars at bottom
    n_bars_int : int
        Number of intermediate reinforcement bars
    dia_int : float
        Diameter of intermediate reinforcement bars
    """
    # Calculate coordinates for core and cover
    core_y = H / 2 - cover_H
    core_z = B / 2 - cover_B
    cover_y = H / 2
    cover_z = B / 2

    # Create the OpenSees fiber section
    ops.section('Fiber', sec_tag, '-GJ', 1.0e9)
    
    # Generate a quadrilateral shaped patch (core patch)
    ops.patch('quad', core_tag,
          10, 10,  # nfIJ, nfJK
          *[-core_y, core_z],
          *[-core_y, -core_z],
          *[core_y, -core_z],
          *[core_y, core_z])
    
    # Define the four cover patches
    ops.patch('quad', cover_tag,
          2, 10,
          *[-cover_y, cover_z],
          *[-core_y, core_z],
          *[core_y, core_z],
          *[cover_y, cover_z])    
    ops.patch('quad', cover_tag,
          2, 10,
          *[-core_y, -core_z],
          *[-cover_y, -cover_z],
          *[cover_y, -cover_z],
          *[core_y, -core_z])    
    ops.patch('quad', cover_tag,
          10, 2,
          *[-cover_y, cover_z],
          *[-cover_y, -cover_z],
          *[-core_y, -core_z],
          *[-core_y, core_z])     
    ops.patch('quad', cover_tag,
          10, 2,
          *[core_y, core_z],
          *[core_y, -core_z],
          *[cover_y, -cover_z],
          *[cover_y, cover_z]) 

    # Construct straight lines of fibers for reinforcement
    # Top primary reinforcement
    ops.layer('straight', steel_tag, n_bars_top, 0.25*np.pi*dia_top**2,
          *[core_y, core_z], *[core_y, -core_z])
    
    # Bottom primary reinforcement
    ops.layer('straight', steel_tag, n_bars_bot, 0.25*np.pi*dia_bot**2,
      *[-core_y, core_z], *[-core_y, -core_z])


    # Secondary top reinforcement
    if n_bars_secondary_top > 0:
        ops.layer('straight', steel_tag, n_bars_secondary_top,
                0.25 * np.pi * dia_sec_top**2,
                *[core_y - offset, core_z],
                *[core_y - offset, -core_z])

    
        # Secondary bottom reinforcement
        ops.layer('straight', steel_tag, n_bars_secondary_bot, 0.25*np.pi*dia_sec_bot**2,
            *[-core_y + offset, core_z], *[-core_y + offset, -core_z])
    
    # Intermediate reinforcement (top and bottom)
    ops.layer('straight', steel_tag, n_bars_int // 2, 0.25*np.pi*dia_int**2,
          *[-core_y, core_z], *[core_y, core_z])
    ops.layer('straight', steel_tag, n_bars_int // 2, 0.25*np.pi*dia_int**2,
          *[-core_y, -core_z], *[core_y, -core_z])
    
def define_circular_rc_section10_opensees(sec_tag, core_tag, cover_tag, steel_tag,
                             D_Sec, cover_Sec, num_Bars_Sec, bar_dia_Sec,
                             ri=0.0, nf_Core_R=8, nf_Core_T=8, nf_Cover_R=4, nf_Cover_T=8):
    """
    Draw a circular RC section with visualization.
    
    Parameters:
    -----------
    sec_tag : int
        Tag for the section
    core_tag : int
        Material tag for core concrete
    cover_tag : int
        Material tag for cover concrete
    steel_tag : int
        Material tag for reinforcement steel
    D_Sec : float
        Outer diameter of the section
    cover_Sec : float
        Cover thickness to reinforcing steel
    num_Bars_Sec : int
        Number of longitudinal reinforcement bars
    bar_dia_Sec : float
        Diameter of longitudinal reinforcement bars
    ri : float, optional
        Inner radius (for hollow sections), default = 0.0
    nf_Core_R : int, optional
        Number of radial divisions in core, default = 8
    nf_Core_T : int, optional
        Number of theta divisions in core, default = 8
    nf_Cover_R : int, optional
        Number of radial divisions in cover, default = 4
    nf_Cover_T : int, optional
        Number of theta divisions in cover, default = 8
    """
    # Calculate key dimensions
    ro = D_Sec / 2  # outer radius
    rc = ro - cover_Sec  # core radius
    bar_area_Sec = 0.25 * np.pi * bar_dia_Sec**2  # bar area
    
    # constructs a FiberSection 
    ops.section('Fiber', sec_tag, '-GJ', 1.0e9)
    
    # Define the core patch (confined concrete)
    ops.patch('circ', core_tag, nf_Core_T, nf_Core_R, 0, 0, ri, rc, 0, 360)
    
    # Define the cover patch (unconfined concrete)
    ops.patch('circ', cover_tag, nf_Cover_T, nf_Cover_R, 0, 0, rc, ro, 0, 360)
    
    # Define the reinforcing layer
    theta = 360.0 / num_Bars_Sec  # angle increment between bars
    ops.layer('circ', steel_tag, num_Bars_Sec, bar_area_Sec, 0, 0, rc, 0, 360)


# def define_L_rc_section10_opensees(sec_tag, core_tag, cover_tag, steel_tag, H, B, t, cover_H, cover_B, offset,
#                                   n_bars_vertical, dia_vertical, n_bars_horizontal, dia_horizontal,
#                                   n_bars_corner, dia_corner):
#     """
#     Define an L-shaped RC section in OpenSees.
    
#     Parameters:
#     -----------
#     sec_tag : int
#         Tag for the section
#     core_tag : int
#         Material tag for core concrete
#     cover_tag : int
#         Material tag for cover concrete
#     steel_tag : int
#         Material tag for reinforcement steel
#     H : float
#         Total height of the L-section
#     B : float
#         Total width of the L-section
#     t : float
#         Thickness of the L-section legs
#     cover_H : float
#         Cover thickness in vertical direction
#     cover_B : float
#         Cover thickness in horizontal direction
#     offset : float
#         Offset for secondary reinforcement from primary reinforcement
#     n_bars_vertical : int
#         Number of vertical reinforcement bars
#     dia_vertical : float
#         Diameter of vertical reinforcement bars
#     n_bars_horizontal : int
#         Number of horizontal reinforcement bars
#     dia_horizontal : float
#         Diameter of horizontal reinforcement bars
#     n_bars_corner : int
#         Number of corner reinforcement bars
#     dia_corner : float
#         Diameter of corner reinforcement bars
#     """
#     # Calculate coordinates for core and cover
#     core_y1 = 0
#     core_y2 = H - cover_H
#     core_z1 = 0
#     core_z2 = t - cover_B  # Vertical leg
#     core_z3 = B - cover_B  # Horizontal leg
    
#     # Create the OpenSees fiber section
#     ops.section('Fiber', sec_tag, '-GJ', 1.0e9)
    
#     # Define core patches for vertical and horizontal legs
#     # Vertical leg core
#     ops.patch('quad', core_tag, 10, 10,
#               core_z1, core_y1,
#               core_z2, core_y1,
#               core_z2, core_y2,
#               core_z1, core_y2)
    
#     # Horizontal leg core
#     ops.patch('quad', core_tag, 10, 10,
#               core_z1, core_y1,
#               core_z3, core_y1,
#               core_z3, t - cover_H,
#               core_z1, t - cover_H)
    
#     # Define cover patches
#     # Vertical leg left cover
#     ops.patch('quad', cover_tag, 2, 10,
#               0, 0,
#               cover_B, 0,
#               cover_B, H,
#               0, H)
    
#     # Vertical leg right cover
#     ops.patch('quad', cover_tag, 2, 10,
#               t - cover_B, 0,
#               t, 0,
#               t, H,
#               t - cover_B, H)
    
#     # Horizontal leg bottom cover
#     ops.patch('quad', cover_tag, 10, 2,
#               0, 0,
#               B, 0,
#               B, cover_H,
#               0, cover_H)
    
#     # Horizontal leg top cover
#     ops.patch('quad', cover_tag, 10, 2,
#               0, t - cover_H,
#               B, t - cover_H,
#               B, t,
#               0, t)
    
#     # Define reinforcement layers
#     bar_area_vertical = 0.25 * np.pi * dia_vertical**2
#     bar_area_horizontal = 0.25 * np.pi * dia_horizontal**2
#     bar_area_corner = 0.25 * np.pi * dia_corner**2
    
#     # Vertical reinforcement
#     vert_bar_x = t - cover_B - offset
#     ops.layer('straight', steel_tag, n_bars_vertical, bar_area_vertical,
#               vert_bar_x, cover_H,
#               vert_bar_x, H - cover_H)
    
#     # Horizontal reinforcement
#     horiz_bar_y = t - cover_H - offset
#     ops.layer('straight', steel_tag, n_bars_horizontal, bar_area_horizontal,
#               cover_B, horiz_bar_y,
#               B - cover_B, horiz_bar_y)
    
#     # Corner reinforcement (placed at intersection of vertical and horizontal bars)
#     ops.layer('straight', steel_tag, n_bars_corner, bar_area_corner,
#               vert_bar_x, horiz_bar_y,
#               vert_bar_x, horiz_bar_y)
  
def define_L_rc_section10_opensees(sec_tag, core_tag, cover_tag, steel_tag, H, B, t, cover_H, cover_B,
                                  n_bars_vertical_outer, dia_vertical_outer, 
                                  n_bars_horizontal_outer, dia_horizontal_outer,
                                  n_bars_vertical_inner, dia_vertical_inner,
                                  n_bars_horizontal_inner, dia_horizontal_inner,
                                  corner_bar_dia):
    """
    Define an L-shaped RC section in OpenSees with detailed reinforcement layout.
    
    Parameters:
    -----------
    sec_tag : int
        Tag for the section
    core_tag : int
        Material tag for core concrete
    cover_tag : int
        Material tag for cover concrete
    steel_tag : int
        Material tag for reinforcement steel
    H : float
        Total height of the L-section
    B : float
        Total width of the L-section
    t : float
        Thickness of the L-section legs
    cover_H : float
        Cover thickness in vertical direction
    cover_B : float
        Cover thickness in horizontal direction
    n_bars_vertical_outer : int
        Number of outer vertical reinforcement bars
    dia_vertical_outer : float
        Diameter of outer vertical reinforcement bars
    n_bars_horizontal_outer : int
        Number of outer horizontal reinforcement bars
    dia_horizontal_outer : float
        Diameter of outer horizontal reinforcement bars
    n_bars_vertical_inner : int
        Number of inner vertical reinforcement bars
    dia_vertical_inner : float
        Diameter of inner vertical reinforcement bars
    n_bars_horizontal_inner : int
        Number of inner horizontal reinforcement bars
    dia_horizontal_inner : float
        Diameter of inner horizontal reinforcement bars
    corner_bar_dia : float
        Diameter of corner reinforcement bars
    """
    import numpy as np
    import openseespy.opensees as ops
    
    # Calculate coordinates for core and cover
    core_y1 = cover_H
    core_y2 = H - cover_H
    core_z1 = cover_B
    core_z2 = t - cover_B  # Vertical leg
    core_z3 = B - cover_B  # Horizontal leg
    
    # Create the OpenSees fiber section
    ops.section('Fiber', sec_tag, '-GJ', 1.0e9)
    
    # Define core patches for vertical and horizontal legs
    # Vertical leg core
    ops.patch('quad', core_tag, 10, 10,
              core_z1, core_y1,
              core_z2, core_y1,
              core_z2, core_y2,
              core_z1, core_y2)
    
    # Horizontal leg core (only the part that doesn't overlap with vertical leg)
    ops.patch('quad', core_tag, 10, 10,
              core_z2, core_y1,
              core_z3, core_y1,
              core_z3, t - cover_H,
              core_z2, t - cover_H)
    
    # Define cover patches
    # Vertical leg left cover
    ops.patch('quad', cover_tag, 2, 10,
              0, 0,
              cover_B, 0,
              cover_B, H,
              0, H)
    
    # Vertical leg right cover
    ops.patch('quad', cover_tag, 2, 10,
              t - cover_B, 0,
              t, 0,
              t, H,
              t - cover_B, H)
    
    # Vertical leg top cover
    ops.patch('quad', cover_tag, 10, 2,
              cover_B, H - cover_H,
              t - cover_B, H - cover_H,
              t - cover_B, H,
              cover_B, H)
    
    # Horizontal leg bottom cover
    ops.patch('quad', cover_tag, 10, 2,
              0, 0,
              B, 0,
              B, cover_H,
              0, cover_H)
    
    # Horizontal leg top cover
    ops.patch('quad', cover_tag, 10, 2,
              t, t - cover_H,
              B, t - cover_H,
              B, t,
              t, t)
    
    # Horizontal leg right cover
    ops.patch('quad', cover_tag, 2, 10,
              B - cover_B, cover_H,
              B, cover_H,
              B, t - cover_H,
              B - cover_B, t - cover_H)
    
    # Calculate bar areas
    bar_area_vertical_outer = 0.25 * np.pi * dia_vertical_outer**2
    bar_area_horizontal_outer = 0.25 * np.pi * dia_horizontal_outer**2
    bar_area_vertical_inner = 0.25 * np.pi * dia_vertical_inner**2
    bar_area_horizontal_inner = 0.25 * np.pi * dia_horizontal_inner**2
    bar_area_corner = 0.25 * np.pi * corner_bar_dia**2
    
    # Outer vertical reinforcement (along the outer edge of vertical leg)
    if n_bars_vertical_outer > 0:
        outer_vert_z = cover_B
        if n_bars_vertical_outer == 1:
            ops.fiber(outer_vert_z, (core_y1 + core_y2) / 2, bar_area_vertical_outer, steel_tag)
        else:
            for i in range(n_bars_vertical_outer):
                y_pos = core_y1 + i * (core_y2 - core_y1) / (n_bars_vertical_outer - 1)
                ops.fiber(outer_vert_z, y_pos, bar_area_vertical_outer, steel_tag)
    
    # Inner vertical reinforcement (along the inner edge of vertical leg)
    if n_bars_vertical_inner > 0:
        inner_vert_z = core_z2
        if n_bars_vertical_inner == 1:
            ops.fiber(inner_vert_z, (core_y1 + core_y2) / 2, bar_area_vertical_inner, steel_tag)
        else:
            for i in range(n_bars_vertical_inner):
                y_pos = core_y1 + i * (core_y2 - core_y1) / (n_bars_vertical_inner - 1)
                ops.fiber(inner_vert_z, y_pos, bar_area_vertical_inner, steel_tag)
    
    # Outer horizontal reinforcement (along the outer edge of horizontal leg)
    if n_bars_horizontal_outer > 0:
        outer_horiz_y = core_y1
        if n_bars_horizontal_outer == 1:
            ops.fiber((core_z1 + core_z3) / 2, outer_horiz_y, bar_area_horizontal_outer, steel_tag)
        else:
            for i in range(n_bars_horizontal_outer):
                z_pos = core_z1 + i * (core_z3 - core_z1) / (n_bars_horizontal_outer - 1)
                ops.fiber(z_pos, outer_horiz_y, bar_area_horizontal_outer, steel_tag)
    
    # Inner horizontal reinforcement (along the inner edge of horizontal leg)
    if n_bars_horizontal_inner > 0:
        inner_horiz_y = t - cover_H
        # Only place bars in the extended part of horizontal leg (beyond vertical leg)
        if n_bars_horizontal_inner == 1:
            z_pos = (core_z2 + core_z3) / 2
            ops.fiber(z_pos, inner_horiz_y, bar_area_horizontal_inner, steel_tag)
        else:
            for i in range(n_bars_horizontal_inner):
                z_pos = core_z2 + i * (core_z3 - core_z2) / (n_bars_horizontal_inner - 1)
                ops.fiber(z_pos, inner_horiz_y, bar_area_horizontal_inner, steel_tag)
    
    # Corner reinforcement (at the intersection of legs)
    ops.fiber(core_z2, core_y1, bar_area_corner, steel_tag)  # Bottom corner
    ops.fiber(core_z1, inner_horiz_y, bar_area_corner, steel_tag)  # Top corner

def assign_member_section_into_opensees(section_definitions, materials):
    # Extract material name → matTag mapping
    material_tag_map = {mat["name"]: mat["matTag"] for mat in materials}

    for category in section_definitions.values():
        for section_name, section_info in category.items():
            tag = section_info["section_tag"]

            core_mat_tag = material_tag_map[section_info["core_tag"]]
            cover_mat_tag = material_tag_map[section_info["cover_tag"]]
            steel_mat_tag = material_tag_map[section_info["steel_tag"]]

            if section_info["type"] == "rectangular":
                define_rc_section10_opensees(
                    sec_tag=tag,
                    core_tag=core_mat_tag,
                    cover_tag=cover_mat_tag,
                    steel_tag=steel_mat_tag,
                    H=section_info["H"],
                    B=section_info["B"],
                    cover_H=section_info["cover_H"],
                    cover_B=section_info["cover_B"],
                    offset=section_info["offset"],
                    n_bars_top=section_info["n_bars_top"],
                    dia_top=section_info["dia_top"],
                    n_bars_bot=section_info["n_bars_bot"],
                    dia_bot=section_info["dia_bot"],
                    n_bars_secondary_top=section_info["n_bars_secondary_top"],
                    dia_sec_top=section_info["dia_sec_top"],
                    n_bars_secondary_bot=section_info["n_bars_secondary_bot"],
                    dia_sec_bot=section_info["dia_sec_bot"],
                    n_bars_int=section_info["n_bars_int"],
                    dia_int=section_info["dia_int"]
                )

            elif section_info["type"] == "circular":
                define_circular_rc_section10_opensees(
                    sec_tag=tag,
                    core_tag=core_mat_tag,
                    cover_tag=cover_mat_tag,
                    steel_tag=steel_mat_tag,
                    D_Sec=section_info["D_Sec"],
                    cover_Sec=section_info["cover_Sec"],
                    num_Bars_Sec=section_info["num_Bars_Sec"],
                    bar_dia_Sec=section_info["bar_dia_Sec"],
                    ri=section_info.get("ri", 0.0),
                    nf_Core_R=section_info.get("nf_Core_R", 8),
                    nf_Core_T=section_info.get("nf_Core_T", 8),
                    nf_Cover_R=section_info.get("nf_Cover_R", 4),
                    nf_Cover_T=section_info.get("nf_Cover_T", 8)
                )

            elif section_info["type"] == "L":
                # define_L_rc_section10_opensees(
                #     sec_tag=tag,
                #     core_tag=core_mat_tag,
                #     cover_tag=cover_mat_tag,
                #     steel_tag=steel_mat_tag,
                #     H=section_info["H"],
                #     B=section_info["B"],
                #     t=section_info["t"],
                #     cover_H=section_info["cover_H"],
                #     cover_B=section_info["cover_B"],
                #     offset=section_info["offset"],
                #     n_bars_vertical=section_info["n_bars_vertical"],
                #     dia_vertical=section_info["dia_vertical"],
                #     n_bars_horizontal=section_info["n_bars_horizontal"],
                #     dia_horizontal=section_info["dia_horizontal"],
                #     n_bars_corner=section_info["n_bars_corner"],
                #     dia_corner=section_info["dia_corner"]
                # )
                # elif section_info["type"] == "L":
                    define_L_rc_section10_opensees(
                        sec_tag=tag,
                        core_tag=core_mat_tag,
                        cover_tag=cover_mat_tag,
                        steel_tag=steel_mat_tag,
                        H=section_info["H"],
                        B=section_info["B"],
                        t=section_info["t"],
                        cover_H=section_info["cover_H"],
                        cover_B=section_info["cover_B"],
                        n_bars_vertical_outer=section_info["n_bars_vertical_outer"],
                        dia_vertical_outer=section_info["dia_vertical_outer"],
                        n_bars_horizontal_outer=section_info["n_bars_horizontal_outer"],
                        dia_horizontal_outer=section_info["dia_horizontal_outer"],
                        n_bars_vertical_inner=section_info["n_bars_vertical_inner"],
                        dia_vertical_inner=section_info["dia_vertical_inner"],
                        n_bars_horizontal_inner=section_info["n_bars_horizontal_inner"],
                        dia_horizontal_inner=section_info["dia_horizontal_inner"],
                        corner_bar_dia=section_info["corner_bar_dia"]
                    )



def extract_materials(materials_config):
    all_materials = []
    for mat in materials_config["materials"]:
        mat_data = {
            "name": mat["name"],
            "type": mat["type"],
            "id": mat["id"],
            **mat["properties"]
        }
        all_materials.append(mat_data)
    return all_materials

def extract_sections(sections_config):
    all_sections = []
    for sec in sections_config["sections"]:
        sec_data = {
            "name": sec["name"],
            "type": sec["type"],
            "id": sec["id"],
            "material_id": sec["material_id"],
            **sec["properties"]
        }
        all_sections.append(sec_data)
    return all_sections




def assign_node_into_opensees(JSON_FOLDER):
    """
    Creates nodes in OpenSees directly from the combined JSON file.

    Args:
        JSON_FOLDER (str): Path to folder containing the combined JSON file
        json_filename (str): Name of the JSON file (default matches create_combined_structure_json output)

    Returns:
        dict: Mapping from node names to OpenSees node IDs
    """
    # Load the combined JSON file
    # Call the first function to load structure data

    # Call the second function to create combined structure JSON
    nodes_data, members_data, mesh_data = create_combined_structure_json(JSON_FOLDER)
    
    # 2. Load the combined structure data to get all node coordinates
    # Merge structural data


    # Process nodes and create OpenSees nodes
    nodes = nodes_data
    all_nodes = {}  # Dictionary to store node information

    for node in nodes:
        try:
            # Validate required fields exist
            node_id = node['id']
            node_name = node['name']
            
            # Require explicit coordinate values - no defaults
            if 'x' not in node or 'y' not in node or 'z' not in node:
                raise ValueError(f"Node {node_id} ({node_name}) is missing coordinate data")
                
            node_coordinates = {
                'x': node['x'],
                'y': node['y'], 
                'z': node['z']
            }
            
            # Store node information
            all_nodes[node_name] = {
                'id': node_id,
                'coordinates': node_coordinates
            }
            
            # Create node in OpenSees
            ops.node(node_id, node_coordinates['x'], node_coordinates['y'], node_coordinates['z'])
            
            # print(f"Created node - ID: {node_id}, Name: {node_name}, "
            #     f"Coordinates: ({node_coordinates['x']}, {node_coordinates['y']}, {node_coordinates['z']})")
        
        except KeyError as e:
            print(f"Error: Node {node.get('id', 'UNKNOWN')} is missing required field {str(e)}")
            continue
        except ValueError as e:
            print(f"Error: {str(e)}")
            continue


def define_materials(materials_list):
    for mat in materials_list:
        ops.uniaxialMaterial(*mat.values())





def assign_beam_column_into_opensees(opensees_element_json_file):
    with open(opensees_element_json_file, "r") as f:
        d = json.load(f)

    for ele in d["elements"]:
        ops.geomTransf(ele["transType"], ele["transTag"],
                       float(ele["vecxz"][0]), float(ele["vecxz"][1]), float(ele["vecxz"][2]))

        ops.element('nonlinearBeamColumn', ele["eleTag"], ele["node_i_id"], ele["node_j_id"],
                    ele["numIntgrPts"], ele["secTag"], ele["transTag"], '-integration', 'Lobatto')



def assign_shell_material_and_section_into_opensees(materials_config, sections_config):
    material_name_to_id = {}
    shell_section_name_to_id = {}

    # Create materials
    for mat in materials_config['materials']:
        mat_id = mat['id']
        mat_type = mat['type']
        props = mat['properties']
        mat_name = mat['name']

        if mat_type == "ENT":
            ops.uniaxialMaterial('ENT', mat_id, props['E'])
        elif mat_type == "Elastic":
            ops.uniaxialMaterial("Elastic", mat_id, props['E'])
        elif mat_type == "ElasticIsotropic":
            ops.nDMaterial('ElasticIsotropic', mat_id, props['E'], props['nu'], props['rho'])

        material_name_to_id[mat_name] = mat_id

    # Create sections
    for sec in sections_config['sections']:
        sec_id = sec['section_id']
        mat_name = sec['material_id']  # material_id is name in your dict
        mat_id = material_name_to_id.get(mat_name)
        props = sec['properties']
        sec_name = sec['name']

        if sec['type'] == "PlateFiber":
            # PlateFiber section command requires material tag(s), thickness, and possibly more parameters
            # Usually, PlateFiber section uses uniaxialMaterial tags of fiber materials, not nDMaterial directly.
            # If mat_id refers to nDMaterial, this could cause error.
            # You may need to create a fiber material or use appropriate material tag here.

            # For now, assume PlateFiber expects material tag and thickness
            # But if mat_id is nDMaterial, it will error. You must ensure PlateFiber material tags are uniaxialMaterial tags.

            # This is likely the cause of error "invalid tags"
            # So to fix: use correct uniaxial material tag for PlateFiber section, not nDMaterial id.

            # If you only have nDMaterial, consider replacing PlateFiber with shell or other section type compatible with nDMaterial.

            ops.section('PlateFiber', sec_id, mat_id, props['thickness'])
            # shell_section_name_to_id[sec_name] = sec_id

    return material_name_to_id, shell_section_name_to_id



def assign_shell_element_into_opensees(JSON_FOLDER, sections_config):
    nodes_data, members_data, mesh_data = create_combined_structure_json(JSON_FOLDER)
    # print("sections_config")
    # print(sections_config)
    # for sec in sections_config['sections']:
    #     print(f"Section Name: {sec['name']}, Section ID: {sec['section_id']}")

    # for node in nodes_data:
    #     print(f"Node Name: {node['name']}, Node ID: {node['id']}")

    section_name_to_id = {sec['name']: sec['section_id'] for sec in sections_config['sections']}

    node_name_to_id = {node['name']: node['id'] for node in nodes_data}

    for elem_key, shell_elem in mesh_data["elements"].items():
        # print(f"Processing element key: {elem_key}")

        elem_id = shell_elem.get("id")
        section_name = shell_elem.get("shell_section")
        node_names = shell_elem.get("nodes")
        section_id = section_name_to_id.get(section_name)

        # print(f"Element ID: {elem_id}")
        # print(f"Element nodes: {node_names}, section: {section_name}")
        # print(f"Resolved section ID: {section_id}")

        if elem_id is None or node_names is None or section_id is None:
            print(f"⚠️ Cannot process element {elem_key}. Skipping.")
            continue

        try:
            node_tags = [node_name_to_id[name] for name in node_names]
        except KeyError as e:
            print(f"❌ Node name {e.args[0]} not found. Skipping element {elem_key}.")
            continue

        if len(node_tags) == 4:
            ops.element("ShellMITC4", elem_id, *node_tags, section_id)
        elif len(node_tags) == 3:
            ops.element("ShellDKGT", elem_id, *node_tags, section_id)
        else:
            print(f"⚠️ Unsupported number of nodes ({len(node_tags)}) for element {elem_key}. Skipping.")




def apply_nodal_loads(JSON_FOLDER, load_combination="Comb2"):
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
    is_mass_load = load_combination == "mass"

    
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

def apply_member_loads(JSON_FOLDER, load_combination="Comb2"):
    load_file_path = os.path.join(JSON_FOLDER, "load_data", f"member_loads_{load_combination}.json")
    
    if not os.path.exists(load_file_path):
        print(f"Load file not found: {load_file_path}")
        return

    with open(load_file_path, 'r') as f:
        data = json.load(f)

    element_loads = data.get("element_loads", [])

    for load in element_loads:
        elem_id = load["member"]

        # Apply uniform load
        uniform = load.get("uniform", [0.0, 0.0, 0.0])
        if any(abs(v) > 1e-10 for v in uniform):
            print(f"Applying uniform load to element {elem_id}: {uniform}")
            ops.eleLoad('-ele', elem_id, '-type', '-beamUniform', uniform[1], uniform[2], uniform[0])

        # Apply point load
        point = load.get("point", [0.0, 0.0, 0.0, 0.0])
        if any(abs(v) > 1e-10 for i, v in enumerate(point) if i != 3):  # skip location in check
            print(f"Applying point load to element {elem_id}: {point}")
            ops.eleLoad('-ele', elem_id, '-type', '-beamPoint', point[1], point[2], point[3], point[0])

        # Apply temperature load
        temp_pts = load.get("temperature_points", [])
        flat_temp = [val for pair in temp_pts for val in pair]
        if any(abs(v) > 1e-10 for v in flat_temp):
            print(f"Applying thermal load to element {elem_id}: {flat_temp}")
            ops.eleLoad('-ele', elem_id, '-type', '-beamThermal', *flat_temp)


def assign_opensees_loads(JSON_FOLDER, combo_name="Comb2"):
    # 1. First get all structural data
    nodes_data, members_data, mesh_data = create_combined_structure_json(JSON_FOLDER)
    
    # 2. Create node name to ID mapping (n10002 -> actual ID)
    node_name_to_id = {node['name']: node['id'] for node in nodes_data}
    
    # 3. Assign nodal point loads
    load_data_dir = os.path.join(JSON_FOLDER, "load_data")
    nodal_file = os.path.join(load_data_dir, f"nodal_loads_{combo_name.lower()}.json")
    print(f'nodal_file={nodal_file}')
    
    with open(nodal_file) as f:
        for node_name, loads in json.load(f).items():
            node_id = node_name_to_id[node_name]
            ops.load(node_id, *[float(v) for v in loads])
            # print(f"Assigned nodal load to node {node_id}: {loads}")
    
    # 4. Assign member loads
    member_file = os.path.join(load_data_dir, f"member_load_{combo_name}.json")
    with open(member_file) as f:
        for elem in json.load(f)["element_loads"]:
            member_id = elem["member"]
            
            # Uniform loads - corrected syntax
            if any(elem["uniform"]):
                print(f"Assigned uniform load to member {member_id}: {elem['uniform']}")
                # Correct order: wx, wy, wz (no reordering needed)
                ops.eleLoad('-ele', member_id, '-type', '-BeamUniform', elem['uniform'])
                print("Successfully assigned member uniform load")
            
            # Point loads - corrected syntax for OpenSees beamPoint
            if any(elem["point"][:3]):  # Check if any of the force components are non-zero
                print(f"Assigned point load to member {member_id}: {elem['point']}")
                # OpenSees beamPoint syntax: Py, Pz, xDivL, Px
                # Your data format: [Px, Py, Pz, location]
                px, py, pz, location = elem["point"]
                ops.eleLoad('-ele', member_id, '-type', '-BeamPoint', py, pz, location, px)
                print("Successfully assigned member point load")

# Analysis Functions - Improved with better recorder handling
def ensure_output_dir():
    """Ensure output directory exists."""
    os.makedirs('FGU_RC3DF_files', exist_ok=True)
# =============================================
# 6. run_gravity
# =============================================
# def run_gravity(JSON_FOLDER, steps=10, comb_name=None):
#     """Run gravity analysis with dynamic recorder naming."""
#     ensure_output_dir()
#     reaction_file = os.path.join(JSON_FOLDER, "post_processing", f"Gravity_Reactions_{comb_name}.out" if comb_name else "Gravity_Reactions.out")
#     # reaction_file = f"Gravity_Reactions_{comb_name}.out" if comb_name else "Gravity_Reactions.out"
#     reaction_path = os.path.join('FGU_RC3DF_files', reaction_file)
    
#     # ops.recorder('Node', '-file', reaction_path,
#     #             '-time', '-node', '-dof', 'reaction')




#     ops.constraints('Transformation')
#     ops.numberer('RCM')
#     ops.system('BandGeneral')
#     # ops.test('NormDispIncr', 1.0e-6, 100, 0, 2)
#     ops.test('NormDispIncr', 1.0e-4, 200, 0, 2)
#     ops.algorithm('Newton')
#     # ops.algorithm('NewtonLineSearch') 
#     # ops.algorithm('KrylovNewton')
#     ops.integrator('LoadControl', 1/steps)
#     # ops.integrator('LoadControl', 1/(steps*2))
#     ops.analysis('Static')
#     # ops.analyze(steps)  # e.g., steps=10
    
#     ops.record()
    
#     ok = ops.analyze(steps)
#     ops.reactions()     # Must call to update reactions!
    
#     if ok == 0:
#         print(f"Gravity analysis {'for ' + comb_name if comb_name else ''} completed successfully")
#         return True
#     else:
#         print(f"Gravity analysis {'for ' + comb_name if comb_name else ''} failed")
#         return False

def run_gravity(JSON_FOLDER, steps=10, comb_name="Comb2"):
        
    """
    Runs gravity analysis.
    Note that the model should be built before
    calling this function.
    
    Keyword arguments:
    steps -- total number of analysis steps

    """
    
    ops.initialize()
    # Records the response of a number of nodes at every converged step
    ops.recorder('Node', '-file', 
                os.path.join('FGU_RC3DF_files','Gravity_Reactions.out'),
        '-time','-node', *list(range(1,5)), '-dof', *list(range(1,7)), 'reaction')

    # enforces the constraints using the transformation method
    ops.constraints('Transformation')

    # RCM numberer uses the reverse Cuthill-McKee scheme to order the matrix equations
    ops.numberer('RCM')

    # Construct a BandGeneralSOE linear system of equation
    ops.system('BandGeneral')

    # Uses the norm of the left hand side solution vector of the matrix equation to
    # determine if convergence has been reached
    ops.test('NormDispIncr', 1.0e-6, 100, 0, 2)

    # Uses the Newton-Raphson algorithm to solve the nonlinear residual equation
    ops.algorithm('Newton')

    # Uses LoadControl integrator object
    ops.integrator('LoadControl', 1/steps)

    # Constructs the Static Analysis object
    ops.analysis('Static')

    # Records the current state of the model
    ops.record()
    # Performs the analysis
    ops.analyze(steps)    
    
    print("Gravity analysis Done!")   


 
# =============================================
# 7. run_modal
# =============================================
def run_modal(n_evs=3, comb_name=None):
    """
    Runs Modal analysis with support for load combinations.
    
    Args:
        n_evs (int): Number of eigenvalues to compute
        comb_name (str): Name of the load combination (for output files)
        
    Returns:
        np.array: Array of eigenvalues
    """
    ensure_output_dir()
    
    # Create unique recorder names if combination is specified
    if comb_name:
        eigen_files = [
            os.path.join('FGU_RC3DF_files', f'ModalAnalysis_EigenVec1_{comb_name}.out'),
            os.path.join('FGU_RC3DF_files', f'ModalAnalysis_EigenVec2_{comb_name}.out'),
            os.path.join('FGU_RC3DF_files', f'ModalAnalysis_EigenVec3_{comb_name}.out')
        ]
        eigenval_file = os.path.join('FGU_RC3DF_files', f'ModalAnalysis_EigenVal_{comb_name}.out')
    else:
        eigen_files = [
            os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVec1.out'),
            os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVec2.out'),
            os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVec3.out')
        ]
        eigenval_file = os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVal.out')

    # Set up recorders for each mode
    for i, file in enumerate(eigen_files[:n_evs], 1):
        ops.recorder('Node', '-file', file,
                    '-node', *list(range(5,9)), '-dof', 1, 2, f'eigen {i}')

    # Analysis configuration
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('BandGen')
    ops.test('NormDispIncr', 1.0e-12, 25, 0, 2)
    ops.algorithm('Newton')
    ops.analysis('Transient')  # Needed for eigen commands
    
    # Compute eigenvalues
    lamda = np.array(ops.eigen(n_evs))
    
    # Write eigenvalues to file
    with open(eigenval_file, "w") as eig_file:
        eig_file.write("lambda omega period frequency\n")
        for l in lamda:
            omega = l**0.5
            period = 2*np.pi/omega
            freq = omega/(2*np.pi)
            eig_file.write(f"{l:2.6e} {omega:2.6e} {period:2.6e} {freq:2.6e}\n")

    ops.record()
    print(f"Modal analysis {'for ' + comb_name if comb_name else ''} completed")
    return lamda
# =============================================
# 8. run_pushover
# =============================================
def run_pushover(steps=10000, direction='X', comb_name=None):
    """Run pushover analysis with dynamic recorder naming."""
    ensure_output_dir()
    
    reaction_file = f"Pushover_Horizontal_Reactions{direction}_{comb_name}.out" if comb_name else f"Pushover_Horizontal_Reactions{direction}.out"
    disp_file = f"Pushover_Story_Displacement{direction}_{comb_name}.out" if comb_name else f"Pushover_Story_Displacement{direction}.out"
    
    reaction_path = os.path.join('FGU_RC3DF_files', reaction_file)
    disp_path = os.path.join('FGU_RC3DF_files', disp_file)

    d_o_f = 1 if direction == 'X' else 2
    phi = 1.0
    
    ops.recorder('Node', '-file', reaction_path,
                '-time', '-node', *list(range(1,5)), '-dof', d_o_f, 'reaction')
    ops.recorder('Node', '-file', disp_path,
                '-time', '-node', *list(range(5,9)), '-dof', d_o_f, 'disp')

    # Create lateral load pattern
    pattern_tag = 100 if comb_name else 2  # Use high tag for combination cases
    ops.pattern('Plain', pattern_tag, 1)
    step = 1.0e-05
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGen')
    ops.test('NormDispIncr', 0.000001, 100)
    ops.algorithm('NewtonLineSearch', True, 0.8, 1000, 0.1, 10.0)
    ops.integrator('DisplacementControl', 5, d_o_f, step)
    ops.analysis('Static')
    ops.record()

    ok = ops.analyze(steps)
    
    if ok == 0:
        print(f'Pushover Analysis in {direction} {"for " + comb_name if comb_name else ""} completed successfully')
        return True
    else:
        print(f'Pushover Analysis in {direction} {"for " + comb_name if comb_name else ""} failed')
        return False
# =============================================
# 9. run_time_history
# =============================================
def run_time_history(direction='X', g_motion_id=1, scaling_id=1,
                    lamda=1.0, acc_file='FGU_RC3DF_files/acc_1.txt',
                    comb_name=None, analysis_duration_ratio=0.29):
    """
    Runs Time history analysis with support for load combinations.
    
    Args:
        direction (str): Direction of excitation ('X' or 'Y')
        g_motion_id (int): Ground motion identifier
        scaling_id (int): Scaling factor identifier
        lamda (float): Scaling factor for ground motion
        acc_file (str): Path to acceleration file
        comb_name (str): Name of the load combination
        analysis_duration_ratio (float): Fraction of full duration to analyze (0-1)
        
    Returns:
        bool: True if analysis succeeded, False otherwise
    """
    ensure_output_dir()
    
    # Create output file names
    if comb_name:
        reaction_file = os.path.join('FGU_RC3DF_files', 
                                   f'TimeHistory_Reactions_{direction}_{g_motion_id}_{scaling_id}_{comb_name}.out')
        disp_file = os.path.join('FGU_RC3DF_files', 
                               f'TimeHistory_Displacement_{direction}_{g_motion_id}_{scaling_id}_{comb_name}.out')
        accel_file = os.path.join('FGU_RC3DF_files', 
                                 f'TimeHistory_Acceleration_{direction}_{g_motion_id}_{scaling_id}_{comb_name}.out')
    else:
        reaction_file = os.path.join('FGU_RC3DF_files', 
                                    f'TimeHistory_Reactions_{direction}_{g_motion_id}.out')
        disp_file = os.path.join('FGU_RC3DF_files', 
                               f'TimeHistory_Displacement_{direction}_{g_motion_id}.out')
        accel_file = os.path.join('FGU_RC3DF_files', 
                                f'TimeHistory_Acceleration_{direction}_{g_motion_id}.out')

    # Determine DOF and damping parameters
    dof = 1 if direction == 'X' else 2
    
    # Get modal properties (assuming first mode dominates)
    try:
        eigenvals = np.loadtxt(os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVal.out'), 
                             skiprows=1)
        omega = eigenvals[0,1] if direction == 'X' else eigenvals[2,1]
    except:
        print("Warning: Could not read eigenvalue file, using default omega=2π")
        omega = 2*np.pi  # Fallback value
    
    xi = 0.05  # Damping ratio (5%)
    alpha_M = 0.0       # Mass proportional damping
    beta_K = 2*xi/omega # Stiffness proportional damping
    
    # Set up recorders
    ops.recorder('Node', '-file', reaction_file,
                '-time', '-node', *list(range(1,5)), '-dof', dof, 'reaction')
    ops.recorder('Node', '-file', disp_file,
                '-time', '-node', *list(range(5,9)), '-dof', dof, 'disp')
    ops.recorder('Node', '-file', accel_file,
                '-time', '-node', *list(range(5,9)), '-dof', dof, 'accel')

    # Load acceleration time history
    accelerogram = np.loadtxt(acc_file)
    dt = 0.02  # Time step in acceleration file
    n_steps = len(accelerogram)
    
    # Analysis parameters
    tol = 1.0e-6
    max_iter = 500
    analysis_dt = 0.01  # Analysis time step (should be ≤ dt/2 for accuracy)

    # Define time series and pattern
    ops.timeSeries('Path', 2, '-dt', dt, '-values', *accelerogram, '-factor', lamda)
    ops.pattern('UniformExcitation', 3, dof, '-accel', 2)
    
    # Analysis configuration
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', tol, max_iter, 0, 2)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.rayleigh(alpha_M, beta_K, 0.0, 0.0)
    ops.analysis('Transient')

    # Run analysis
    print(f"Running Time-History analysis (λ={lamda}) {'for ' + comb_name if comb_name else ''}")
    start_time = time.time()
    
    ok = 0
    current_time = ops.getTime()
    final_time = n_steps * dt
    target_time = analysis_duration_ratio * final_time
    
    while ok == 0 and current_time < target_time:
        ok = ops.analyze(1, analysis_dt)
        current_time = ops.getTime()
        
        # Optional: Print progress
        if int(current_time/dt) % 100 == 0:
            print(f"Time: {current_time:.2f}s ({current_time/target_time:.1%})")
    
    elapsed_time = time.time() - start_time
    
    if ok == 0:
        print(f"Time-History completed in {elapsed_time:.2f}s")
        return True
    else:
        print(f"Time-History failed at {current_time:.2f}s")
        return False
# =============================================
# 10. reset_analysis
# =============================================
def reset_analysis():
    """Reset the analysis state."""
    ops.setTime(0.0)
    ops.loadConst()
    ops.remove('recorders')
    ops.wipeAnalysis()

     



# OUTPUT_FOLDER = "output_folder"  # Main output directory
# os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist

# Subdirectories
# JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
# IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
# os.makedirs(JSON_FOLDER, exist_ok=True)
# os.makedirs(IMAGE_FOLDER, exist_ok=True)
# opensees_element_json_file = os.path.join(JSON_FOLDER, "element_data.json")
# final_run(JSON_FOLDER, materials, section_definitions, materials_config, sections_config, opensees_element_json_file)






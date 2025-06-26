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
import json
import os
import numpy as np
import matplotlib.tri as mtri
from matplotlib import colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# =============================================
# Common Utility Functions
# =============================================
def create_output_directories():
    """Create necessary output directories"""
    OUTPUT_FOLDER = "postprocessing_folder"
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    
    # Subdirectories
    JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
    IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
    
    # Create subdirectories for different element types
    os.makedirs(os.path.join(JSON_FOLDER, "member"), exist_ok=True)
    os.makedirs(os.path.join(JSON_FOLDER, "shell_element"), exist_ok=True)
    os.makedirs(os.path.join(JSON_FOLDER, "slab_reinforcement"), exist_ok=True)
    
    os.makedirs(os.path.join(IMAGE_FOLDER, "member"), exist_ok=True)
    os.makedirs(os.path.join(IMAGE_FOLDER, "shell_element"), exist_ok=True)
    os.makedirs(os.path.join(IMAGE_FOLDER, "slab_reinforcement"), exist_ok=True)
    os.makedirs(os.path.join(IMAGE_FOLDER, "structural_model"), exist_ok=True)
    
    return JSON_FOLDER, IMAGE_FOLDER

def save_to_json(data, filepath):
    """Save data to JSON file"""
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=4)
    print(f"Saved results to: {filepath}")

# =============================================
# Beam Element Functions
# =============================================
def create_beam_response_json(load_combination="combo", num_points=10):
    """
    Create JSON files for beam element responses (forces, stresses, strains)
    """
    JSON_FOLDER, _ = create_output_directories()
    
    # Get all beam elements
    beam_elements = [e for e in ops.getEleTags() if 'ForceBeamColumn3d' in ops.eleType(e)]
    
    # Initialize data storage
    beam_forces_data = {}
    beam_stresses_data = {}
    beam_strains_data = {}
    beam_deflection_data = {}
    
    for ele_tag in beam_elements:
        # Get basic element info
        nodes = ops.eleNodes(ele_tag)
        node_coords = [ops.nodeCoord(n) for n in nodes]
        L = np.linalg.norm(np.array(node_coords[1]) - np.array(node_coords[0]))
        
        # Initialize element data
        ele_force_data = {
            "element_tag": ele_tag,
            "nodes": nodes,
            "node_coords": node_coords,
            "load_combination": load_combination,
            "forces": [],
            "locations": []
        }
        
        ele_stress_data = {
            "element_tag": ele_tag,
            "nodes": nodes,
            "node_coords": node_coords,
            "load_combination": load_combination,
            "stresses": [],
            "locations": []
        }
        
        ele_strain_data = {
            "element_tag": ele_tag,
            "nodes": nodes,
            "node_coords": node_coords,
            "load_combination": load_combination,
            "strains": [],
            "locations": []
        }
        
        ele_deflection_data = {
            "element_tag": ele_tag,
            "nodes": nodes,
            "node_coords": node_coords,
            "load_combination": load_combination,
            "deflections": [],
            "relative_deflections": [],
            "locations": []
        }
        
        # Sample points along element
        for loc in np.linspace(0, 1, num_points):
            x = loc * L
            
            # Get forces
            forces = ops.eleResponse(ele_tag, 'forces', loc)
            ele_force_data["forces"].append(forces)
            ele_force_data["locations"].append(loc)
            
            # Get stresses
            stresses = ops.eleResponse(ele_tag, 'stresses', loc)
            ele_stress_data["stresses"].append(stresses)
            ele_stress_data["locations"].append(loc)
            
            # Get strains
            strains = ops.eleResponse(ele_tag, 'strains', loc)
            ele_strain_data["strains"].append(strains)
            ele_strain_data["locations"].append(loc)
            
            # Calculate deflections
            disp_i = ops.nodeDisp(nodes[0])
            disp_j = ops.nodeDisp(nodes[1])
            
            dx = deflection(disp_i, disp_j, 'dx', x, L)
            dy = deflection(disp_i, disp_j, 'dy', x, L)
            dz = deflection(disp_i, disp_j, 'dz', x, L)
            
            rdy = rel_deflection(disp_i, disp_j, 'dy', x, L)
            rdz = rel_deflection(disp_i, disp_j, 'dz', x, L)
            
            ele_deflection_data["deflections"].append([dx, dy, dz])
            ele_deflection_data["relative_deflections"].append([rdy, rdz])
            ele_deflection_data["locations"].append(loc)
        
        # Store element data
        beam_forces_data[f"element_{ele_tag}"] = ele_force_data
        beam_stresses_data[f"element_{ele_tag}"] = ele_stress_data
        beam_strains_data[f"element_{ele_tag}"] = ele_strain_data
        beam_deflection_data[f"element_{ele_tag}"] = ele_deflection_data
    
    # Save to JSON files
    save_to_json(beam_forces_data, os.path.join(JSON_FOLDER, "member", f"beam_forces_{load_combination}.json"))
    save_to_json(beam_stresses_data, os.path.join(JSON_FOLDER, "member", f"beam_stresses_{load_combination}.json"))
    save_to_json(beam_strains_data, os.path.join(JSON_FOLDER, "member", f"beam_strains_{load_combination}.json"))
    save_to_json(beam_deflection_data, os.path.join(JSON_FOLDER, "member", f"beam_deflections_{load_combination}.json"))

def plot_beam_results(load_combination="combo"):
    """
    Create plots for beam element results from JSON files
    """
    JSON_FOLDER, IMAGE_FOLDER = create_output_directories()
    beam_json_path = os.path.join(JSON_FOLDER, "member")
    
    # Plot forces
    with open(os.path.join(beam_json_path, f"beam_forces_{load_combination}.json"), 'r') as f:
        force_data = json.load(f)
    
    # Create force plots for each component (N, Vy, Vz, T, My, Mz)
    force_components = ['N', 'Vy', 'Vz', 'T', 'My', 'Mz']
    for comp in force_components:
        plt.figure(figsize=(10, 6))
        for ele_id, ele_data in force_data.items():
            locs = ele_data['locations']
            forces = [f[force_components.index(comp)] for f in ele_data['forces']]
            plt.plot(locs, forces, label=f"Element {ele_data['element_tag']}")
        
        plt.xlabel('Position along element')
        plt.ylabel(f'{comp} Force')
        plt.title(f'Beam {comp} Forces - {load_combination}')
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(IMAGE_FOLDER, "member", f"beam_{comp}_forces_{load_combination}.png"))
        plt.close()
    
    # Plot deflections
    with open(os.path.join(beam_json_path, f"beam_deflections_{load_combination}.json"), 'r') as f:
        deflection_data = json.load(f)
    
    directions = ['dx', 'dy', 'dz']
    for dir in directions:
        plt.figure(figsize=(10, 6))
        for ele_id, ele_data in deflection_data.items():
            locs = ele_data['locations']
            deflections = [d[directions.index(dir)] for d in ele_data['deflections']]
            plt.plot(locs, deflections, label=f"Element {ele_data['element_tag']}")
        
        plt.xlabel('Position along element')
        plt.ylabel(f'{dir} Deflection')
        plt.title(f'Beam {dir} Deflections - {load_combination}')
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(IMAGE_FOLDER, "member", f"beam_{dir}_deflections_{load_combination}.png"))
        plt.close()

# =============================================
# Shell Element Functions
# =============================================
def create_shell_response_json(load_combination="combo", samples_per_side=5):
    """
    Create JSON files for shell element responses (forces, stresses, strains, deflections)
    """
    JSON_FOLDER, _ = create_output_directories()
    
    # Get all shell elements
    shell_elements = [e for e in ops.getEleTags() if 'Shell' in ops.eleType(e)]
    
    # Initialize data storage
    shell_forces_data = {}
    shell_stresses_data = {}
    shell_strains_data = {}
    shell_deflection_data = {}
    
    for ele_tag in shell_elements:
        # Get element nodes and coordinates
        nodes = ops.eleNodes(ele_tag)
        node_coords = [ops.nodeCoord(n) for n in nodes]
        
        # Initialize element data
        ele_force_data = {
            "element_tag": ele_tag,
            "nodes": nodes,
            "node_coords": node_coords,
            "load_combination": load_combination,
            "forces": ops.eleResponse(ele_tag, 'forces')
        }
        
        ele_stress_data = {
            "element_tag": ele_tag,
            "nodes": nodes,
            "node_coords": node_coords,
            "load_combination": load_combination,
            "stresses": ops.eleResponse(ele_tag, 'stresses')
        }
        
        ele_strain_data = {
            "element_tag": ele_tag,
            "nodes": nodes,
            "node_coords": node_coords,
            "load_combination": load_combination,
            "strains": ops.eleResponse(ele_tag, 'strains')
        }
        
        # Calculate deflections at nodes
        node_deflections = []
        for node in nodes:
            disp = ops.nodeDisp(node)
            node_deflections.append(disp)
        
        ele_deflection_data = {
            "element_tag": ele_tag,
            "nodes": nodes,
            "node_coords": node_coords,
            "load_combination": load_combination,
            "deflections": node_deflections
        }
        
        # Store element data
        shell_forces_data[f"element_{ele_tag}"] = ele_force_data
        shell_stresses_data[f"element_{ele_tag}"] = ele_stress_data
        shell_strains_data[f"element_{ele_tag}"] = ele_strain_data
        shell_deflection_data[f"element_{ele_tag}"] = ele_deflection_data
    
    # Save to JSON files
    save_to_json(shell_forces_data, os.path.join(JSON_FOLDER, "shell_element", f"shell_forces_{load_combination}.json"))
    save_to_json(shell_stresses_data, os.path.join(JSON_FOLDER, "shell_element", f"shell_stresses_{load_combination}.json"))
    save_to_json(shell_strains_data, os.path.join(JSON_FOLDER, "shell_element", f"shell_strains_{load_combination}.json"))
    save_to_json(shell_deflection_data, os.path.join(JSON_FOLDER, "shell_element", f"shell_deflections_{load_combination}.json"))

def plot_shell_results(load_combination="combo"):
    """
    Create plots for shell element results from JSON files
    """
    JSON_FOLDER, IMAGE_FOLDER = create_output_directories()
    shell_json_path = os.path.join(JSON_FOLDER, "shell_element")
    
    # Plot stresses
    with open(os.path.join(shell_json_path, f"shell_stresses_{load_combination}.json"), 'r') as f:
        stress_data = json.load(f)
    
    stress_components = ['sxx', 'syy', 'sxy', 'vmis']
    for comp in stress_components:
        plt.figure(figsize=(10, 8))
        
        # Collect data for contour plot
        coords = []
        values = []
        for ele_id, ele_data in stress_data.items():
            centroid = np.mean(ele_data['node_coords'], axis=0)
            coords.append(centroid[:2])  # Use only x,y for 2D plot
            if comp == 'vmis':
                # Calculate von Mises stress
                sxx = ele_data['stresses'][0]
                syy = ele_data['stresses'][1]
                sxy = ele_data['stresses'][2]
                vmis = np.sqrt(sxx**2 + syy**2 - sxx*syy + 3*sxy**2)
                values.append(vmis)
            else:
                idx = stress_components.index(comp)
                values.append(ele_data['stresses'][idx])
        
        # Create triangulation and contour plot
        coords = np.array(coords)
        triangulation = mtri.Triangulation(coords[:, 0], coords[:, 1])
        
        plt.tricontourf(triangulation, values, levels=20, cmap='jet')
        plt.colorbar()
        plt.title(f'Shell {comp} Stress - {load_combination}')
        plt.xlabel('X [m]')
        plt.ylabel('Y [m]')
        plt.axis('equal')
        plt.savefig(os.path.join(IMAGE_FOLDER, "shell_element", f"shell_{comp}_stress_{load_combination}.png"))
        plt.close()
    
    # Plot deflections
    with open(os.path.join(shell_json_path, f"shell_deflections_{load_combination}.json"), 'r') as f:
        deflection_data = json.load(f)
    
    directions = ['dx', 'dy', 'dz']
    for dir in directions:
        plt.figure(figsize=(10, 8))
        
        # Collect data for contour plot
        coords = []
        values = []
        for ele_id, ele_data in deflection_data.items():
            centroid = np.mean(ele_data['node_coords'], axis=0)
            coords.append(centroid[:2])  # Use only x,y for 2D plot
            node_deflections = [d[directions.index(dir)] for d in ele_data['deflections']]
            values.append(np.mean(node_deflections))  # Average deflection at element
        
        # Create triangulation and contour plot
        coords = np.array(coords)
        triangulation = mtri.Triangulation(coords[:, 0], coords[:, 1])
        
        plt.tricontourf(triangulation, values, levels=20, cmap='jet')
        plt.colorbar()
        plt.title(f'Shell {dir} Deflection - {load_combination}')
        plt.xlabel('X [m]')
        plt.ylabel('Y [m]')
        plt.axis('equal')
        plt.savefig(os.path.join(IMAGE_FOLDER, "shell_element", f"shell_{dir}_deflection_{load_combination}.png"))
        plt.close()

# =============================================
# Deflection Calculation Functions
# =============================================
def deflection(disp_i, disp_j, direction, x, L):
    """
    Calculate deflection at position x along beam using linear interpolation.
    """
    direction_map = {'dx': 0, 'dy': 1, 'dz': 2}
    idx = direction_map[direction]
    
    # Linear interpolation between nodes
    return disp_i[idx] + (disp_j[idx] - disp_i[idx]) * (x / L)

def rel_deflection(disp_i, disp_j, direction, x, L):
    """
    Calculate relative deflection (deflection relative to chord between end points).
    """
    direction_map = {'dy': 1, 'dz': 2}
    idx = direction_map[direction]
    
    # Deflection at position x
    deflection_x = disp_i[idx] + (disp_j[idx] - disp_i[idx]) * (x / L)
    
    # Chord deflection (linear between end points)
    chord_deflection = disp_i[idx] + (disp_j[idx] - disp_i[idx]) * (x / L)
    
    # Relative deflection is difference from chord
    return deflection_x - chord_deflection

# =============================================
# Main Analysis Workflow
# =============================================
# def run_analysis(load_combination="combo"):
#     """
#     Main function to run the entire analysis workflow
#     """
#     # Create JSON files first
#     create_beam_response_json(load_combination)
#     create_shell_response_json(load_combination)
    
#     # Then create plots
#     plot_beam_results(load_combination)
#     plot_shell_results(load_combination)
    
#     # Create structural model plots
#     structural_model_plot(load_combination=load_combination)

def structural_model_plot(OUTPUT_FOLDER="postprocessing_folder", load_combination="combo"):
    """Create plots of the structural model"""
    JSON_FOLDER, IMAGE_FOLDER = create_output_directories()
    STRUCTURAL_MODEL_FOLDER = os.path.join(IMAGE_FOLDER, "structural_model")
    
    # Model plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Get all elements
    elements = ops.getEleTags()
    
    for ele_tag in elements:
        # Get node coordinates of the element
        ele_nodes = ops.eleNodes(ele_tag)
        node_coords = np.array([ops.nodeCoord(node) for node in ele_nodes])
        
        # Create a filled polygon
        poly = Poly3DCollection([node_coords], alpha=0.5, linewidth=1, edgecolor='k')
        poly.set_facecolor('yellow')
        ax.add_collection3d(poly)

    # Overlay the original model edges
    opsv.plot_model(element_labels=0, node_labels=0, ax=ax, fmt_model={'color': 'k', 'linewidth': 1})
    plt.title(f"Model - {load_combination}")
    filepath = os.path.join(STRUCTURAL_MODEL_FOLDER, f"model_{load_combination}.png")
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Deformation plot
    plt.figure(figsize=(10, 8))
    opsv.plot_model()
    plt.title(f"Deformation - {load_combination}")
    filepath = os.path.join(STRUCTURAL_MODEL_FOLDER, f"deformation_{load_combination}.png")
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Load plot
    plt.figure(figsize=(10, 8))
    opsv.plot_load()
    plt.title(f"Load - {load_combination}")
    filepath = os.path.join(STRUCTURAL_MODEL_FOLDER, f"load_{load_combination}.png")
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Extruded shapes plot
    element_data_path = os.path.join("output_folder", "json_files", "element_data.json")
    with open(element_data_path, 'r') as f:
        data = json.load(f)
    
    ele_shapes = {}
    for ele in data["elements"]:
        tag = ele["eleTag"]
        sec_tag = ele["secTag"]
        
        for section_type in section_definitions:
            sections = section_definitions[section_type]
            for section_name, section_data in sections.items():
                if section_data["section_tag"] != sec_tag:
                    continue
                
                if section_data["type"] == "rectangular":
                    B = section_data["B"]
                    H = section_data["H"]
                    ele_shapes[tag] = ["rect", [B, H]]
                elif section_data["type"] == "circular":
                    D = section_data["D_Sec"]
                    ele_shapes[tag] = ["circ", [D]]
                elif section_data["type"] == "L":
                    H = section_data["H"]
                    B = section_data["B"]
                    t = section_data["t"]
                    ele_shapes[tag] = ["L", [H, B, t]]
                elif section_data["type"] == "I":
                    B = section_data["B"]
                    H = section_data["H"]
                    tf = section_data.get("tf", B / 10.)
                    tw = section_data.get("tw", H / 6.)
                    ele_shapes[tag] = ["I", [B, H, tf, tw]]
                    
    plt.figure(figsize=(10, 8))
    opsv.plot_extruded_shapes_3d(ele_shapes)
    plt.title(f"Extruded Shapes - {load_combination}")
    filepath = os.path.join(STRUCTURAL_MODEL_FOLDER, f"extruded_shapes_{load_combination}.png")
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"All structural model plots saved to {STRUCTURAL_MODEL_FOLDER}")
# import os
# import json
# from matplotlib import pyplot as plt
# import numpy as np
# import pandas as pd
# from test5 import CQC, calculate_floor_masses, calculate_floor_stiffness, extract_and_combine_forces, extract_node_data, model, plot, response_spectrum_analysis
# import openseespy.opensees as ops
# import sys
# from math import pi, sqrt
# import opsvis as opsv



# Tn = [0.0, 0.06, 0.1, 0.12, 0.18, 0.24, 0.3, 0.36, 0.4, 0.42, 
#     0.48, 0.54, 0.6, 0.66, 0.72, 0.78, 0.84, 0.9, 0.96, 1.02, 
#     1.08, 1.14, 1.2, 1.26, 1.32, 1.38, 1.44, 1.5, 1.56, 1.62, 
#     1.68, 1.74, 1.8, 1.86, 1.92, 1.98, 2.04, 2.1, 2.16, 2.22, 
#     2.28, 2.34, 2.4, 2.46, 2.52, 2.58, 2.64, 2.7, 2.76, 2.82, 
#     2.88, 2.94, 3.0, 3.06, 3.12, 3.18, 3.24, 3.3, 3.36, 3.42, 
#     3.48, 3.54, 3.6, 3.66, 3.72, 3.78, 3.84, 3.9, 3.96, 4.02, 
#     4.08, 4.14, 4.2, 4.26, 4.32, 4.38, 4.44, 4.5, 4.56, 4.62, 
#     4.68, 4.74, 4.8, 4.86, 4.92, 4.98, 5.04, 5.1, 5.16, 5.22, 
#     5.28, 5.34, 5.4, 5.46, 5.52, 5.58, 5.64, 5.7, 5.76, 5.82, 
#     5.88, 5.94, 6.0]

# Sa = [1.9612, 3.72628, 4.903, 4.903, 4.903, 4.903, 4.903, 4.903, 4.903, 4.6696172, 
#     4.0861602, 3.6321424, 3.2683398, 2.971218, 2.7241068, 2.5142584, 2.3348086, 2.1788932, 2.0425898, 1.9229566, 
#     1.8160712, 1.7199724, 1.6346602, 1.5562122, 1.485609, 1.4208894, 1.3620534, 1.3071398, 1.2571292, 1.211041, 
#     1.166914, 1.1267094, 1.0894466, 1.054145, 1.0217852, 0.990406, 0.960988, 0.9335312, 0.9080356, 0.8835206, 
#     0.8599862, 0.838413, 0.8168398, 0.7972278, 0.7785964, 0.759965, 0.7432948, 0.7266246, 0.710935, 0.6952454, 
#     0.6805364, 0.666808, 0.6540602, 0.6285646, 0.6040496, 0.5814958, 0.5609032, 0.5403106, 0.5206986, 0.5030478, 
#     0.485397, 0.4697074, 0.4540178, 0.4393088, 0.4255804, 0.411852, 0.3991042, 0.3863564, 0.3755698, 0.3638026, 
#     0.353016, 0.34321, 0.333404, 0.3245786, 0.3157532, 0.3069278, 0.2981024, 0.2902576, 0.2833934, 0.2755486, 
#     0.2686844, 0.2618202, 0.254956, 0.2490724, 0.2431888, 0.2373052, 0.2314216, 0.2265186, 0.220635, 0.215732, 
#     0.210829, 0.205926, 0.2020036, 0.1971006, 0.1931782, 0.1892558, 0.1853334, 0.181411, 0.1774886, 0.1735662, 
#     0.1706244, 0.166702, 0.1637602]


# output_dir = "modal_analysis"
# os.makedirs(output_dir, exist_ok=True)
# model()
# response_spectrum_analysis(output_dir, Tn, Sa, 5)











# # done
# ops.wipe()

# print(f"\nAll output files and images have been saved to the '{output_dir}' directory.")












import csv
import json
from math import sqrt
import os
from post_processing import extract_beam_results
import openseespy.opensees as ops
import opsvis as opsv
from matplotlib import pyplot as plt

def calculate_Cs(S, T, TB, TC, TD, xi, Z, I, R):
    # Calculate the damping correction factor mu
    mu = (10 / (5 + xi)) ** 0.5
    # Ensure mu is not smaller than 0.55
    mu = max(mu, 0.55)

    # Initialize Cs
    Cs = 0

    # Calculate Cs based on the given conditions
    if 0 <= T <= TB:
        Cs = S * (1 + (T / TB) * (2.5 * mu - 1))
    elif TB < T <= TC:  # Fixed: changed from TB <= T to TB < T
        Cs = 2.5 * S * mu
    elif TC < T <= TD:  # Fixed: changed from TC <= T to TC < T
        Cs = 2.5 * S * mu * (TC / T)
    elif TD < T <= 4:   # Fixed: changed from TD <= T to TD < T
        Cs = 2.5 * S * mu * (TC * TD / T ** 2)

    Sa = (2 / 3) * (Z * I / R) * Cs
    return Cs, Sa

def CQC(mu, lambdas, dmp, scalf):
    u = 0.0
    ne = len(lambdas)
    for i in range(ne):
        for j in range(ne):
            di = dmp[i]
            dj = dmp[j]
            bij = lambdas[i]/lambdas[j]
            rho = ((8.0*sqrt(di*dj)*(di+bij*dj)*(bij**(3.0/2.0))) /
                ((1.0-bij**2.0)**2.0 + 4.0*di*dj*bij*(1.0+bij**2.0) + 
                4.0*(di**2.0 + dj**2.0)*bij**2.0))
            u += scalf[i]*mu[i] * scalf[j]*mu[j] * rho;
    return sqrt(u)


def plot(output_dir, JSON_FOLDER, modal_props, Tn, Sa, eigs, nmode, story_heights):
    # Create the modal_analysis directory if it doesn't exist
    # Plot model and save to file
    model_plot_path = os.path.join(output_dir, "model_plot.png")
    opsv.plot_model()
    plt.savefig(model_plot_path)
    plt.title("Structure Model")
    plt.close()


    # Gráficas de modos
    for i in range(1, nmode + 1):
        mode_plot_path = os.path.join(output_dir, f'mode{i}.png')
        
        # Create the mode shape plot
        opsv.plot_mode_shape(i, 
                            endDispFlag=0, 
                            fig_wi_he=(18, 18), 
                            node_supports=False)
        
        # Improve the title and add mode frequency if available
        plt.title(f"Mode Shape {i}", fontsize=16)
        
        # Add grid and adjust layout
        plt.grid(True)
        plt.tight_layout()
        
        # Save the figure with high DPI
        plt.savefig(mode_plot_path, dpi=300, bbox_inches='tight')
        plt.close()

    # Plot Response Spectrum
    response_spectrum_path = os.path.join(output_dir, 'response_spectrum.png')
    plt.figure(figsize=(10, 6))
    plt.plot(Tn, Sa, 'b-', linewidth=2, label='Design Spectrum')
    plt.xlabel('Period (s)')
    plt.ylabel('Spectral Acceleration (g)')
    plt.title('Response Spectrum')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(response_spectrum_path)
    plt.close()


    # ================== NODE DATA EXTRACTION ==================

    # Extract and save node data
    node_coords, node_displacements, node_reactions = extract_node_data(output_dir, eigs)



    # Main calculation sequence
    print("\nCalculating floor properties...")

    # First calculate floor masses
    floor_masses = calculate_floor_masses(node_coords)

    # Then calculate stiffness properties using the floor masses
    floor_stiffness, floor_stiffness_values = calculate_floor_stiffness(
        node_coords, modal_props, eigs, floor_masses)

    

    # Prepare data in JSON-friendly format
    floor_data = []

    for z in sorted(floor_masses.keys()):
        com_x, com_y, _, mass = floor_masses[z]
        if z in floor_stiffness:
            cos_x, cos_y, _ = floor_stiffness[z]
            Kx, Ky, Kr = floor_stiffness_values[z]
        else:
            cos_x, cos_y = com_x, com_y
            Kx, Ky, Kr = 0, 0, 0

        floor_data.append({
            "Floor_Z": round(z, 2),
            "COM_X": round(com_x, 3),
            "COM_Y": round(com_y, 3),
            "Mass": round(mass, 3),
            "COS_X": round(cos_x, 3),
            "COS_Y": round(cos_y, 3),
            "Kx": round(Kx, 3),
            "Ky": round(Ky, 3),
            "Kr": round(Kr, 3)
        })

    # Save to JSON file
    with open(os.path.join(output_dir, "floor_properties.json"), "w") as f:
        json.dump(floor_data, f, indent=4)
    # compute the modal properties
    modal_report_path = os.path.join(output_dir, "ModalReport1.json")
    ops.modalProperties("-print", "-file", modal_report_path, "-unorm")

    # Get modal properties as a dictionary using -return flag
    modal_data = ops.modalProperties("-return")

    # Save to JSON file
    output_path = os.path.join(output_dir, "modal_properties.json")
    with open(output_path, "w") as f:
        json.dump(modal_data, f, indent=4)

    # some settings for the response spectrum analysis
    tsTag = 1 # use the timeSeries 1 as response spectrum function
    direction = 1 # excited DOF = Ux

    # currently we use same damping for each mode
    dmp = [0.05]*len(eigs)
    # we don't want to scale some modes...
    scalf = [1.0]*len(eigs)

    # Run RSA for all modes (your existing code)
    for i in range(len(eigs)):
        ops.responseSpectrumAnalysis(direction, '-Tn', *Tn, '-Sa', *Sa, '-mode', i+1)
        ops.analyze(1)
    

    extract_and_combine_forces_multiple_sections(output_dir, JSON_FOLDER, Tn, Sa, direction, eigs, dmp, scalf, 
                                               num_sections=5, json_filepath='RSA_Forces_MultiSection.json')

    
    extract_and_combine_nodal_responses(output_dir, Tn, Sa, direction, eigs, dmp, scalf,
                                      story_heights,
                                      json_filepath='RSA_Nodal_Responses.json',
                                      csv_filepath='RSA_Nodal_Responses.csv',
                                      drift_json_filepath='RSA_Story_Drifts.json')


def model():
    ops.wipe()

    # define a 3D model
    ops.model("basic","-ndm",3,"-ndf",6)



    # a uniaxial material for transverse shear
    ops.uniaxialMaterial("Elastic",2,938000000.0)

    # the elastic beam section and aggregator
    ops.section("Elastic",1,30000000000.0,0.09,0.0006749999999999999,0.0006749999999999999,12500000000.0,0.0011407499999999994)
    ops.section("Aggregator",3,2,"Vy",2,"Vz","-section",1)

    # nodes and masses
    ops.node(1,0,0,0)
    ops.node(2,0,0,3,"-mass",200,200,200,0,0,0)
    ops.node(3,4,0,3,"-mass",200,200,200,0,0,0)
    ops.node(4,4,0,0)
    ops.node(5,0,0,6,"-mass",200,200,200,0,0,0)
    ops.node(6,4,0,6,"-mass",200,200,200,0,0,0)
    ops.node(7,4,3,6,"-mass",200,200,200,0,0,0)
    ops.node(8,0,3,6,"-mass",200,200,200,0,0,0)
    ops.node(9,0,3,3,"-mass",200,200,200,0,0,0)
    ops.node(10,0,3,0)
    ops.node(11,4,3,3,"-mass",200,200,200,0,0,0)
    ops.node(12,4,3,0)
    ops.node(13,2,1.5,6)
    ops.node(14,2,1.5,3)

    # beam elements
    ops.beamIntegration("Lobatto", 1, 3, 5)
    # beam_column_elements forceBeamColumn
    # Geometric transformation command
    ops.geomTransf("Linear", 1, 1.0, 0.0, -0.0)
    ops.element("forceBeamColumn", 1, 1, 2, 1, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 2, 0.0, 0.0, 1.0)
    ops.element("forceBeamColumn", 2, 2, 3, 2, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 3, 1.0, 0.0, -0.0)
    ops.element("forceBeamColumn", 3, 4, 3, 3, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 4, 1.0, 0.0, -0.0)
    ops.element("forceBeamColumn", 4, 2, 5, 4, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 5, 0.0, 0.0, 1.0)
    ops.element("forceBeamColumn", 5, 5, 6, 5, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 6, 0.0, 0.0, 1.0)
    ops.element("forceBeamColumn", 6, 7, 6, 6, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 7, 0.0, 0.0, 1.0)
    ops.element("forceBeamColumn", 7, 8, 7, 7, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 8, 0.0, 0.0, 1.0)
    ops.element("forceBeamColumn", 8, 9, 2, 8, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 9, 0.0, 0.0, 1.0)
    ops.element("forceBeamColumn", 9, 8, 5, 9, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 10, 1.0, 0.0, -0.0)
    ops.element("forceBeamColumn", 10, 10, 9, 10, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 11, 1.0, 0.0, -0.0)
    ops.element("forceBeamColumn", 11, 3, 6, 11, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 12, 1.0, 0.0, -0.0)
    ops.element("forceBeamColumn", 12, 11, 7, 12, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 13, 0.0, 0.0, 1.0)
    ops.element("forceBeamColumn", 13, 11, 3, 13, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 14, 0.0, 0.0, 1.0)
    ops.element("forceBeamColumn", 14, 9, 11, 14, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 15, 1.0, 0.0, -0.0)
    ops.element("forceBeamColumn", 15, 12, 11, 15, 1)
    # Geometric transformation command
    ops.geomTransf("Linear", 16, 1.0, 0.0, -0.0)
    ops.element("forceBeamColumn", 16, 9, 8, 16, 1)

    # Constraints.sp fix
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(10, 1, 1, 1, 1, 1, 1)
    ops.fix(4, 1, 1, 1, 1, 1, 1)
    ops.fix(12, 1, 1, 1, 1, 1, 1)
    ops.fix(13, 0, 0, 1, 1, 1, 0)
    ops.fix(14, 0, 0, 1, 1, 1, 0)

    # Constraints.mp rigidDiaphragm
    ops.rigidDiaphragm(3, 14, 2, 3, 9, 11)
    ops.rigidDiaphragm(3, 13, 5, 6, 7, 8)

def response_spectrum_analysis(output_dir, JSON_FOLDER, Tn, Sa, nmode):
    # the response spectrum function



    ops.timeSeries("Path",1,"-time", *Tn, "-values", *Sa)
    ops.mass(4, 20.217, 20.217, 20.217, 0, 0, 0)
    # define some analysis settings
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", 0.0001, 10)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 0.0)
    ops.analysis("Static")
    ops.reactions()

    # run the eigenvalue analysis with 7 modes
    # and obtain the eigenvalues
    eigs = ops.eigen("-genBandArpack", nmode)

    # After running ops.eigen(), get modal properties with return option
    modal_props = ops.modalProperties("-return")

    # Define your story heights
    story_heights = {
        1: 0.0,      # Ground level
        2: 10.0,      # First floor at 4m
        # 3: 6.0,      # Second floor at 8m
        # 4: 12.0,     # Third floor at 12m
        # ... add more stories as needed
    }
    plot(output_dir, JSON_FOLDER, modal_props, Tn, Sa, eigs, nmode, story_heights)
    


def extract_and_combine_forces_multiple_sections(output_dir, JSON_FOLDER, Tn, Sa, direction, eigs, dmp, scalf, 
                                               num_sections=10, json_filepath='RSA_Forces_MultiSection.json'):
    """
    Extract all member forces from RSA at multiple integration points and perform CQC combination
    
    Args:
        output_dir (str): Output directory path
        Tn (list): List of periods for response spectrum
        Sa (list): List of spectral accelerations
        direction (int): Excitation direction (1=X, 2=Y, 3=Z)
        eigs (list): Eigenvalues from modal analysis
        dmp (list): Damping ratios for each mode
        scalf (list): Scaling factors for each mode
        JSON_FOLDER (str): Path to folder containing element_data.json
        num_sections (int): Number of evaluation points along elements (default=10)
        json_filepath (str): Path to save JSON results
        csv_filepath (str): Path to save CSV results
    """
    # Initialize dictionaries to store all forces at all sections
    modal_forces = {
        'P': {}, 'Vy': {}, 'Vz': {}, 
        'T': {}, 'My': {}, 'Mz': {}
    }
    element_tags = ops.getEleTags()
    
    # Extract forces for each mode
    for mode in range(1, len(eigs)+1):
        ops.responseSpectrumAnalysis(direction, '-Tn', *Tn, '-Sa', *Sa, '-mode', mode)
        
        # Initialize mode entries
        for force_type in modal_forces:
            modal_forces[force_type][mode] = {}
        
        # Get forces for all elements at all sections
        for ele_tag in element_tags:
            # Initialize element entry
            for force_type in modal_forces:
                modal_forces[force_type][mode][ele_tag] = {}
            
            try:
                # Use extract_beam_results to get forces at multiple points
                beam_results = extract_beam_results(ele_tag, nep=num_sections, JSON_FOLDER=JSON_FOLDER)
                
                # Distribute forces to sections
                for section_idx, force_data in enumerate(beam_results['forces'], 1):
                    modal_forces['P'][mode][ele_tag][section_idx] = force_data['N']
                    modal_forces['Vy'][mode][ele_tag][section_idx] = force_data['Vy']
                    modal_forces['Vz'][mode][ele_tag][section_idx] = force_data['Vz']
                    modal_forces['T'][mode][ele_tag][section_idx] = force_data['T']
                    modal_forces['My'][mode][ele_tag][section_idx] = force_data['My']
                    modal_forces['Mz'][mode][ele_tag][section_idx] = force_data['Mz']
            except Exception as e:
                print(f"Warning: Could not extract forces for element {ele_tag} in mode {mode}: {str(e)}")
                # Initialize with zeros if extraction fails
                for section in range(1, num_sections+1):
                    modal_forces['P'][mode][ele_tag][section] = 0.0
                    modal_forces['Vy'][mode][ele_tag][section] = 0.0
                    modal_forces['Vz'][mode][ele_tag][section] = 0.0
                    modal_forces['T'][mode][ele_tag][section] = 0.0
                    modal_forces['My'][mode][ele_tag][section] = 0.0
                    modal_forces['Mz'][mode][ele_tag][section] = 0.0
    
    # Perform CQC combination for each section
    cqc_forces = {
        'P': {}, 'Vy': {}, 'Vz': {}, 
        'T': {}, 'My': {}, 'Mz': {}
    }
    
    for ele_tag in element_tags:
        # Initialize element entry in CQC results
        for force_type in cqc_forces:
            cqc_forces[force_type][ele_tag] = {}
        
        # CQC combination for each section
        for section in range(1, num_sections+1):
            # Extract modal forces for this element and section
            P_modes = [modal_forces['P'][m][ele_tag][section] for m in range(1, len(eigs)+1)]
            Vy_modes = [modal_forces['Vy'][m][ele_tag][section] for m in range(1, len(eigs)+1)]
            Vz_modes = [modal_forces['Vz'][m][ele_tag][section] for m in range(1, len(eigs)+1)]
            T_modes = [modal_forces['T'][m][ele_tag][section] for m in range(1, len(eigs)+1)]
            My_modes = [modal_forces['My'][m][ele_tag][section] for m in range(1, len(eigs)+1)]
            Mz_modes = [modal_forces['Mz'][m][ele_tag][section] for m in range(1, len(eigs)+1)]
            
            # Perform CQC combination
            cqc_forces['P'][ele_tag][section] = CQC(P_modes, eigs, dmp, scalf)
            cqc_forces['Vy'][ele_tag][section] = CQC(Vy_modes, eigs, dmp, scalf)
            cqc_forces['Vz'][ele_tag][section] = CQC(Vz_modes, eigs, dmp, scalf)
            cqc_forces['T'][ele_tag][section] = CQC(T_modes, eigs, dmp, scalf)
            cqc_forces['My'][ele_tag][section] = CQC(My_modes, eigs, dmp, scalf)
            cqc_forces['Mz'][ele_tag][section] = CQC(Mz_modes, eigs, dmp, scalf)

    # Find critical forces (max absolute values across sections)
    critical_forces = {}
    for ele_tag in element_tags:
        critical_forces[ele_tag] = {
            'P': max(abs(f) for f in cqc_forces['P'][ele_tag].values()),
            'Vy': max(abs(f) for f in cqc_forces['Vy'][ele_tag].values()),
            'Vz': max(abs(f) for f in cqc_forces['Vz'][ele_tag].values()),
            'T': max(abs(f) for f in cqc_forces['T'][ele_tag].values()),
            'My': max(abs(f) for f in cqc_forces['My'][ele_tag].values()),
            'Mz': max(abs(f) for f in cqc_forces['Mz'][ele_tag].values())
        }

    # Save to JSON
    json_filepath = os.path.join(output_dir, json_filepath)
    with open(json_filepath, 'w') as f:
        json.dump({
            'modal_forces': modal_forces,
            'cqc_forces': cqc_forces,
            'critical_forces': critical_forces,
            'eigenvalues': eigs,
            'damping': dmp,
            'scaling_factors': scalf,
            'num_sections': num_sections
        }, f, indent=4)
    
    
    return modal_forces, cqc_forces, critical_forces


def extract_node_data(output_dir, eigs):
    """Extract and save node coordinates, displacements, and reactions"""
    # Get all node tags
    node_tags = ops.getNodeTags()
    
    # Initialize dictionaries to store data
    node_coords = {}
    node_displacements = {}
    node_reactions = {}
    
    # Extract data for each node
    for node_tag in node_tags:
        # Get coordinates
        node_coords[node_tag] = ops.nodeCoord(node_tag)
        
        # Get displacements (modal displacements for each mode)
        node_displacements[node_tag] = {}
        for mode in range(1, len(eigs)+1):
            node_displacements[node_tag][mode] = ops.nodeEigenvector(node_tag, mode)
        
        # Get reactions (only for fixed nodes)
        # if ops.nodeReaction(node_tag) != [0.0]*6:  # Only save if non-zero
        node_reactions[node_tag] = ops.nodeReaction(node_tag)
    
    
    return node_coords, node_displacements, node_reactions

# Function to calculate floor masses (COM and mass)
def calculate_floor_masses(node_coords):
    # Group nodes by floor (Z coordinate)
    floors = {}
    for node_tag, coords in node_coords.items():
        z = coords[2]
        if z not in floors:
            floors[z] = []
        floors[z].append(node_tag)
    
    # Sort floors by elevation
    sorted_floors = sorted(floors.keys())
    
    # Calculate mass properties for each floor
    floor_masses = {}
    for z in sorted_floors:
        total_mass_x = 0
        total_mass_y = 0
        total_mass = 0
        for node_tag in floors[z]:
            mass = ops.nodeMass(node_tag)
            if sum(mass) > 0:  # Only nodes with mass
                x, y, _ = node_coords[node_tag]
                total_mass_x += x * mass[0]  # Using x-mass (DOF 1)
                total_mass_y += y * mass[1]  # Using y-mass (DOF 2)
                total_mass += mass[0]  # Assuming mass is same in x and y
        
        if total_mass > 0:
            com_x = total_mass_x / total_mass
            com_y = total_mass_y / total_mass
            floor_masses[z] = (com_x, com_y, z, total_mass)
        else:
            # For floors with no mass, use geometric center
            x_coords = [node_coords[tag][0] for tag in floors[z]]
            y_coords = [node_coords[tag][1] for tag in floors[z]]
            com_x = sum(x_coords)/len(x_coords) if x_coords else 0
            com_y = sum(y_coords)/len(y_coords) if y_coords else 0
            floor_masses[z] = (com_x, com_y, z, 0)
    
    return floor_masses

# Function to calculate floor stiffness properties
def calculate_floor_stiffness(node_coords, modal_props, eigs, floor_masses):
    # Get modal properties
    parti_factor_MX = modal_props["partiFactorMX"]
    parti_factor_MY = modal_props["partiFactorMY"]
    parti_factor_RMZ = modal_props["partiFactorRMZ"]
    
    floor_stiffness = {}
    floor_stiffness_values = {}
    
    # Find dominant translational modes
    x_modes = [i for i in range(len(parti_factor_MX)) if abs(parti_factor_MX[i]) > 0.3]
    y_modes = [i for i in range(len(parti_factor_MY)) if abs(parti_factor_MY[i]) > 0.3]
    
    if not x_modes or not y_modes:
        print("Warning: Could not identify clear translational modes for stiffness calculation")
        return floor_stiffness, floor_stiffness_values
    
    # We'll use the first significant mode in each direction
    x_mode = x_modes[0]
    y_mode = y_modes[0]
    
    # Calculate global stiffness eccentricities
    ey = parti_factor_RMZ[x_mode]/parti_factor_MX[x_mode] if parti_factor_MX[x_mode] != 0 else 0
    ex = -parti_factor_RMZ[y_mode]/parti_factor_MY[y_mode] if parti_factor_MY[y_mode] != 0 else 0
    
    # Calculate approximate stiffness values using modal frequencies
    omega_x = sqrt(eigs[x_mode])
    omega_y = sqrt(eigs[y_mode])
    
    # Total mass (sum of all floor masses)
    total_mass = sum(mass for _, _, _, mass in floor_masses.values())
    
    # Calculate floor stiffness properties
    for z, (com_x, com_y, _, mass) in floor_masses.items():
        if total_mass > 0:
            # Calculate stiffness values
            Kx = mass * omega_x**2 if mass > 0 else 0
            Ky = mass * omega_y**2 if mass > 0 else 0
            Kr = max(Kx, Ky) * 100  # Simplified rotational stiffness
            
            # Center of stiffness for this floor
            cos_x = com_x + ex
            cos_y = com_y + ey
            
            floor_stiffness[z] = (cos_x, cos_y, z)
            floor_stiffness_values[z] = (Kx, Ky, Kr)
    
    return floor_stiffness, floor_stiffness_values


def extract_and_combine_nodal_responses(output_dir, Tn, Sa, direction, eigs, dmp, scalf,
                                      story_heights,
                                      json_filepath='RSA_Nodal_Responses.json',
                                      csv_filepath='RSA_Nodal_Responses.csv',
                                      drift_json_filepath='RSA_Story_Drifts.json'):
    """
    Extract all nodal displacements and reactions from RSA and perform CQC combination
    
    Args:
        output_dir (str): Directory to save output files
        Tn (list): List of periods for response spectrum
        Sa (list): List of spectral accelerations
        direction (int): Excitation direction (1=X, 2=Y, 3=Z)
        eigs (list): Eigenvalues from modal analysis
        dmp (list): Damping ratios for each mode
        scalf (list): Scaling factors for each mode
        story_heights (dict): Dictionary mapping story level to height {story: height}
        json_filepath (str): Path to save JSON results
        csv_filepath (str): Path to save CSV results
        drift_json_filepath (str): Path to save story drift JSON results
    
    Returns:
        tuple: (modal_displacements, modal_reactions, cqc_displacements, cqc_reactions, story_drifts)
    """
    # Initialize dictionaries to store all nodal responses
    modal_displacements = {
        'Ux': {}, 'Uy': {}, 'Uz': {},
        'Rx': {}, 'Ry': {}, 'Rz': {}
    }
    modal_reactions = {
        'Fx': {}, 'Fy': {}, 'Fz': {},
        'Mx': {}, 'My': {}, 'Mz': {}
    }
    
    # Get all node tags
    node_tags = ops.getNodeTags()
    
    # Extract responses for each mode
    for mode in range(1, len(eigs)+1):
        ops.responseSpectrumAnalysis(direction, '-Tn', *Tn, '-Sa', *Sa, '-mode', mode)
        ops.reactions()  # <-- Add this line to compute reactions

        # Initialize mode entries
        for disp_type in modal_displacements:
            modal_displacements[disp_type][mode] = {}
        for react_type in modal_reactions:
            modal_reactions[react_type][mode] = {}
        
        # Get displacements and reactions for all nodes
        for node_tag in node_tags:
            try:
                # Get nodal displacements
                displacements = ops.nodeDisp(node_tag)
                
                # Store displacements (handle 2D/3D cases)
                if len(displacements) >= 6:  # 3D case with rotations
                    modal_displacements['Ux'][mode][node_tag] = displacements[0]
                    modal_displacements['Uy'][mode][node_tag] = displacements[1]
                    modal_displacements['Uz'][mode][node_tag] = displacements[2]
                    modal_displacements['Rx'][mode][node_tag] = displacements[3]
                    modal_displacements['Ry'][mode][node_tag] = displacements[4]
                    modal_displacements['Rz'][mode][node_tag] = displacements[5]
                elif len(displacements) >= 3:  # 3D case without rotations or 2D with rotation
                    modal_displacements['Ux'][mode][node_tag] = displacements[0]
                    modal_displacements['Uy'][mode][node_tag] = displacements[1]
                    modal_displacements['Uz'][mode][node_tag] = displacements[2] if len(displacements) > 2 else 0.0
                    modal_displacements['Rx'][mode][node_tag] = 0.0
                    modal_displacements['Ry'][mode][node_tag] = 0.0
                    modal_displacements['Rz'][mode][node_tag] = displacements[2] if len(displacements) == 3 else 0.0
                else:  # 2D case
                    modal_displacements['Ux'][mode][node_tag] = displacements[0] if len(displacements) > 0 else 0.0
                    modal_displacements['Uy'][mode][node_tag] = displacements[1] if len(displacements) > 1 else 0.0
                    modal_displacements['Uz'][mode][node_tag] = 0.0
                    modal_displacements['Rx'][mode][node_tag] = 0.0
                    modal_displacements['Ry'][mode][node_tag] = 0.0
                    modal_displacements['Rz'][mode][node_tag] = 0.0
                    
            except:
                print(f"Warning: Could not get displacements for node {node_tag}")
                for disp_type in modal_displacements:
                    modal_displacements[disp_type][mode][node_tag] = 0.0
            
            try:
                # Get nodal reactions
                reactions = ops.nodeReaction(node_tag)
                
                # Store reactions (handle 2D/3D cases)
                if len(reactions) >= 6:  # 3D case with moments
                    modal_reactions['Fx'][mode][node_tag] = reactions[0]
                    modal_reactions['Fy'][mode][node_tag] = reactions[1]
                    modal_reactions['Fz'][mode][node_tag] = reactions[2]
                    modal_reactions['Mx'][mode][node_tag] = reactions[3]
                    modal_reactions['My'][mode][node_tag] = reactions[4]
                    modal_reactions['Mz'][mode][node_tag] = reactions[5]
                elif len(reactions) >= 3:  # 3D case without moments or 2D with moment
                    modal_reactions['Fx'][mode][node_tag] = reactions[0]
                    modal_reactions['Fy'][mode][node_tag] = reactions[1]
                    modal_reactions['Fz'][mode][node_tag] = reactions[2] if len(reactions) > 2 else 0.0
                    modal_reactions['Mx'][mode][node_tag] = 0.0
                    modal_reactions['My'][mode][node_tag] = 0.0
                    modal_reactions['Mz'][mode][node_tag] = reactions[2] if len(reactions) == 3 else 0.0
                else:  # 2D case
                    modal_reactions['Fx'][mode][node_tag] = reactions[0] if len(reactions) > 0 else 0.0
                    modal_reactions['Fy'][mode][node_tag] = reactions[1] if len(reactions) > 1 else 0.0
                    modal_reactions['Fz'][mode][node_tag] = 0.0
                    modal_reactions['Mx'][mode][node_tag] = 0.0
                    modal_reactions['My'][mode][node_tag] = 0.0
                    modal_reactions['Mz'][mode][node_tag] = 0.0
                    
            except:
                print(f"Warning: Could not get reactions for node {node_tag}")
                for react_type in modal_reactions:
                    modal_reactions[react_type][mode][node_tag] = 0.0
    
    # Perform CQC combination
    cqc_displacements = {
        'Ux': {}, 'Uy': {}, 'Uz': {},
        'Rx': {}, 'Ry': {}, 'Rz': {}
    }
    cqc_reactions = {
        'Fx': {}, 'Fy': {}, 'Fz': {},
        'Mx': {}, 'My': {}, 'Mz': {}
    }
    
    for node_tag in node_tags:
        # CQC combination for displacements
        for disp_type in cqc_displacements:
            disp_modes = [modal_displacements[disp_type][m][node_tag] for m in range(1, len(eigs)+1)]
            cqc_displacements[disp_type][node_tag] = CQC(disp_modes, eigs, dmp, scalf)
        
        # CQC combination for reactions
        for react_type in cqc_reactions:
            react_modes = [modal_reactions[react_type][m][node_tag] for m in range(1, len(eigs)+1)]
            cqc_reactions[react_type][node_tag] = CQC(react_modes, eigs, dmp, scalf)
    
    # Calculate story drifts if story heights are provided
    story_drifts = None
    print("story_heights")
    print(story_heights)
    if story_heights is not None:
        story_drifts = extract_story_drifts(cqc_displacements, node_tags, story_heights)
        
        # Save story drifts to separate JSON file
        drift_json_filepath = os.path.join(output_dir, drift_json_filepath)
        with open(drift_json_filepath, 'w') as f:
            json.dump({
                'story_drifts': story_drifts,
                'story_heights': story_heights,
                'direction': direction,
                'eigenvalues': eigs,
                'damping': dmp,
                'scaling_factors': scalf
            }, f, indent=4)
        
        print(f"Story drifts saved to: {drift_json_filepath}")
    
    # Save to JSON
    json_filepath = os.path.join(output_dir, json_filepath)
    with open(json_filepath, 'w') as f:
        json.dump({
            'modal_displacements': modal_displacements,
            'modal_reactions': modal_reactions,
            'cqc_displacements': cqc_displacements,
            'cqc_reactions': cqc_reactions,
            'story_drifts': story_drifts,
            'eigenvalues': eigs,
            'damping': dmp,
            'scaling_factors': scalf,
            'node_tags': node_tags
        }, f, indent=4)
    
    # Calculate base and story shears (now with mode-wise results)
    shear_results = calculate_base_and_story_shears(
        output_dir=output_dir,
        modal_reactions=modal_reactions,
        cqc_reactions=cqc_reactions,
        story_heights=story_heights,
        eigs=eigs,
        json_filepath='RSA_Base_Story_Shears.json'  # You can customize this filename if needed
    )
    
    return modal_displacements, modal_reactions, cqc_displacements, cqc_reactions, story_drifts


def extract_story_drifts(cqc_displacements, node_tags, story_heights):
    """
    Calculate story drifts from nodal displacements
    
    Args:
        cqc_displacements (dict): CQC combined displacements
        node_tags (list): List of node tags
        story_heights (dict): Dictionary mapping story level to height {story: height}
    
    Returns:
        dict: Story drifts for each direction
    """
    story_drifts = {'X': {}, 'Y': {}}
    
    # Group nodes by story level (assuming node tags or coordinates indicate story)
    # This is a simplified approach - you may need to modify based on your node numbering
    story_nodes = {}
    for node_tag in node_tags:
        # Get node coordinates to determine story level
        coords = ops.nodeCoord(node_tag)
        z_coord = coords[2] if len(coords) > 2 else 0.0
        
        # Find which story this node belongs to
        story = None
        for story_level, height in story_heights.items():
            if abs(z_coord - height) < 0.01:  # tolerance for floating point comparison
                story = story_level
                break
        
        if story is not None:
            if story not in story_nodes:
                story_nodes[story] = []
            story_nodes[story].append(node_tag)
    
    # Calculate drifts between stories
    sorted_stories = sorted(story_nodes.keys())
    
    for i in range(1, len(sorted_stories)):
        upper_story = sorted_stories[i]
        lower_story = sorted_stories[i-1]
        
        # Get maximum displacement for each story
        upper_disp_x = max([abs(cqc_displacements['Ux'][node]) for node in story_nodes[upper_story]])
        lower_disp_x = max([abs(cqc_displacements['Ux'][node]) for node in story_nodes[lower_story]])
        
        upper_disp_y = max([abs(cqc_displacements['Uy'][node]) for node in story_nodes[upper_story]])
        lower_disp_y = max([abs(cqc_displacements['Uy'][node]) for node in story_nodes[lower_story]])
        
        # Calculate drift
        story_height = story_heights[upper_story] - story_heights[lower_story]
        drift_x = abs(upper_disp_x - lower_disp_x) / story_height
        drift_y = abs(upper_disp_y - lower_disp_y) / story_height
        
        story_drifts['X'][f"Story_{upper_story}"] = drift_x
        story_drifts['Y'][f"Story_{upper_story}"] = drift_y
    
    return story_drifts



def calculate_base_and_story_shears(output_dir, modal_reactions, cqc_reactions, story_heights, eigs,
                                  json_filepath='RSA_Base_Story_Shears.json'):
    """
    Calculate base shear and story shear forces from nodal reactions (both mode-wise and CQC-combined)
    
    Args:
        output_dir (str): Output directory path
        modal_reactions (dict): Modal reactions from extract_and_combine_nodal_responses
        cqc_reactions (dict): CQC combined reactions from extract_and_combine_nodal_responses
        story_heights (dict): Dictionary mapping story level to height {story: height}
        eigs (list): List of eigenvalues from modal analysis
        json_filepath (str): Path to save JSON results
    
    Returns:
        dict: Dictionary containing base shear and story shear results (both mode-wise and CQC-combined)
    """
    # Initialize results dictionary
    shear_results = {
        'base_shear': {
            'CQC': {'Fx': 0.0, 'Fy': 0.0, 'Fz': 0.0},
            'modal': {}
        },
        'story_shears': {
            'CQC': {},
            'modal': {}
        },
        'story_overturning_moments': {
            'CQC': {},
            'modal': {}
        },
        'modes': [f"Mode {i+1}" for i in range(len(eigs))]
    }
    
    # Get all node tags with reactions
    node_tags = [tag for tag in cqc_reactions['Fx'].keys()]
    
    # Group nodes by story level
    story_nodes = {}
    for node_tag in node_tags:
        # Get node coordinates to determine story level
        coords = ops.nodeCoord(node_tag)
        z_coord = coords[2] if len(coords) > 2 else 0.0
        
        # Find which story this node belongs to
        story = None
        for story_level, height in story_heights.items():
            if abs(z_coord - height) < 0.01:  # tolerance for floating point comparison
                story = story_level
                break
        
        if story is not None:
            if story not in story_nodes:
                story_nodes[story] = []
            story_nodes[story].append(node_tag)
    
    # =============================================
    # CQC-combined results (existing functionality)
    # =============================================
    
    # Calculate base shear (sum of reactions at base nodes)
    base_nodes = story_nodes.get(1, [])  # Assuming story 1 is the base
    for node_tag in base_nodes:
        shear_results['base_shear']['CQC']['Fx'] += cqc_reactions['Fx'][node_tag]
        shear_results['base_shear']['CQC']['Fy'] += cqc_reactions['Fy'][node_tag]
        shear_results['base_shear']['CQC']['Fz'] += cqc_reactions['Fz'][node_tag]
    
    # Calculate story shears and moments (CQC-combined)
    for story, nodes in story_nodes.items():
        story_shear_x = 0.0
        story_shear_y = 0.0
        story_moment_x = 0.0
        story_moment_y = 0.0
        
        for node_tag in nodes:
            # Get node coordinates
            coords = ops.nodeCoord(node_tag)
            z = coords[2] if len(coords) > 2 else 0.0
            
            # Sum shear forces
            story_shear_x += cqc_reactions['Fx'][node_tag]
            story_shear_y += cqc_reactions['Fy'][node_tag]
            
            # Sum moments about base (z=0)
            story_moment_x += cqc_reactions['Fy'][node_tag] * z
            story_moment_y += cqc_reactions['Fx'][node_tag] * z
        
        shear_results['story_shears']['CQC'][story] = {
            'Fx': story_shear_x,
            'Fy': story_shear_y
        }
        
        shear_results['story_overturning_moments']['CQC'][story] = {
            'Mx': story_moment_x,
            'My': story_moment_y
        }
    
    # Calculate cumulative story shears (from top down) - CQC-combined
    sorted_stories = sorted(story_nodes.keys(), reverse=True)
    cumulative_shear_x = 0.0
    cumulative_shear_y = 0.0
    
    for story in sorted_stories:
        cumulative_shear_x += shear_results['story_shears']['CQC'][story]['Fx']
        cumulative_shear_y += shear_results['story_shears']['CQC'][story]['Fy']
        
        shear_results['story_shears']['CQC'][story]['cumulative_Fx'] = cumulative_shear_x
        shear_results['story_shears']['CQC'][story]['cumulative_Fy'] = cumulative_shear_y
    
    # =============================================
    # Mode-wise results (new functionality)
    # =============================================
    
    # Initialize modal results
    for mode in range(1, len(eigs)+1):
        shear_results['base_shear']['modal'][mode] = {'Fx': 0.0, 'Fy': 0.0, 'Fz': 0.0}
        shear_results['story_shears']['modal'][mode] = {}
        shear_results['story_overturning_moments']['modal'][mode] = {}
    
    # Calculate mode-wise base shear
    for mode in range(1, len(eigs)+1):
        for node_tag in base_nodes:
            shear_results['base_shear']['modal'][mode]['Fx'] += modal_reactions['Fx'][mode][node_tag]
            shear_results['base_shear']['modal'][mode]['Fy'] += modal_reactions['Fy'][mode][node_tag]
            shear_results['base_shear']['modal'][mode]['Fz'] += modal_reactions['Fz'][mode][node_tag]
    
    # Calculate mode-wise story shears and moments
    for mode in range(1, len(eigs)+1):
        for story, nodes in story_nodes.items():
            story_shear_x = 0.0
            story_shear_y = 0.0
            story_moment_x = 0.0
            story_moment_y = 0.0
            
            for node_tag in nodes:
                # Get node coordinates
                coords = ops.nodeCoord(node_tag)
                z = coords[2] if len(coords) > 2 else 0.0
                
                # Sum shear forces
                story_shear_x += modal_reactions['Fx'][mode][node_tag]
                story_shear_y += modal_reactions['Fy'][mode][node_tag]
                
                # Sum moments about base (z=0)
                story_moment_x += modal_reactions['Fy'][mode][node_tag] * z
                story_moment_y += modal_reactions['Fx'][mode][node_tag] * z
            
            shear_results['story_shears']['modal'][mode][story] = {
                'Fx': story_shear_x,
                'Fy': story_shear_y
            }
            
            shear_results['story_overturning_moments']['modal'][mode][story] = {
                'Mx': story_moment_x,
                'My': story_moment_y
            }
        
        # Calculate cumulative story shears for this mode (from top down)
        cumulative_shear_x = 0.0
        cumulative_shear_y = 0.0
        
        for story in sorted_stories:
            cumulative_shear_x += shear_results['story_shears']['modal'][mode][story]['Fx']
            cumulative_shear_y += shear_results['story_shears']['modal'][mode][story]['Fy']
            
            shear_results['story_shears']['modal'][mode][story]['cumulative_Fx'] = cumulative_shear_x
            shear_results['story_shears']['modal'][mode][story]['cumulative_Fy'] = cumulative_shear_y
    
    # Save to JSON
    json_filepath = os.path.join(output_dir, json_filepath)
    with open(json_filepath, 'w') as f:
        json.dump(shear_results, f, indent=4)
    
    return shear_results



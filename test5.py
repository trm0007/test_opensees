import csv
import json
from math import sqrt
import os
import openseespy.opensees as ops
import opsvis as opsv
from matplotlib import pyplot as plt


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


def extract_and_combine_forces_multiple_sections(output_dir, Tn, Sa, direction, eigs, dmp, scalf, 
                                               num_sections=10, json_filepath='RSA_Forces_MultiSection.json', 
                                               csv_filepath='RSA_Forces_MultiSection.csv'):
    """
    Extract all member forces from RSA at multiple integration points and perform CQC combination
    
    Args:
        Tn (list): List of periods for response spectrum
        Sa (list): List of spectral accelerations
        direction (int): Excitation direction (1=X, 2=Y, 3=Z)
        eigs (list): Eigenvalues from modal analysis
        dmp (list): Damping ratios for each mode
        scalf (list): Scaling factors for each mode
        num_sections (int): Number of integration points to extract (default=10)
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
            
            # Extract forces at multiple integration points
            for section in range(1, num_sections+1):
                try:
                    forces = ops.eleResponse(ele_tag, 'section', str(section), 'force')
                    modal_forces['P'][mode][ele_tag][section] = forces[0]
                    modal_forces['Vy'][mode][ele_tag][section] = forces[1]
                    modal_forces['Vz'][mode][ele_tag][section] = forces[2]
                    modal_forces['T'][mode][ele_tag][section] = forces[3]
                    modal_forces['My'][mode][ele_tag][section] = forces[4]
                    modal_forces['Mz'][mode][ele_tag][section] = forces[5]
                except:
                    # If section doesn't exist, store zeros or skip
                    print(f"Warning: Section {section} not available for element {ele_tag}")
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

    # Find critical forces
    critical_forces = find_critical_forces(cqc_forces, element_tags)
    
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
    
    # Save to CSV for easier analysis
    csv_filepath = os.path.join(output_dir, csv_filepath)
    with open(csv_filepath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        writer.writerow(['Element', 'Section', 'P', 'Vy', 'Vz', 'T', 'My', 'Mz'])
        
        # Write CQC results
        for ele_tag in element_tags:
            for section in range(1, num_sections+1):
                writer.writerow([
                    ele_tag, section,
                    cqc_forces['P'][ele_tag][section],
                    cqc_forces['Vy'][ele_tag][section],
                    cqc_forces['Vz'][ele_tag][section],
                    cqc_forces['T'][ele_tag][section],
                    cqc_forces['My'][ele_tag][section],
                    cqc_forces['Mz'][ele_tag][section]
                ])
    
    return modal_forces, cqc_forces, critical_forces


def find_critical_forces(cqc_forces, element_tags):
    """
    Find critical (maximum) forces across all sections for each element
    
    Args:
        cqc_forces (dict): CQC combined forces from extract_and_combine_forces_multiple_sections
        element_tags (list): List of element tags
    
    Returns:
        dict: Dictionary with maximum forces and their locations
    """
    critical_forces = {}
    
    for ele_tag in element_tags:
        critical_forces[ele_tag] = {}
        
        for force_type in ['P', 'Vy', 'Vz', 'T', 'My', 'Mz']:
            # Find maximum force and its location
            max_force = 0.0
            max_section = 1
            
            for section in cqc_forces[force_type][ele_tag]:
                if abs(cqc_forces[force_type][ele_tag][section]) > abs(max_force):
                    max_force = cqc_forces[force_type][ele_tag][section]
                    max_section = section
            
            critical_forces[ele_tag][force_type] = {
                'value': max_force,
                'section': max_section
            }
    
    return critical_forces


def get_section_coordinates(num_sections):
    """
    Get normalized coordinates (0 to 1) for integration points along element length
    This assumes Gauss-Lobatto integration points distribution
    
    Args:
        num_sections (int): Number of sections/integration points
    
    Returns:
        list: Normalized coordinates from 0 (start) to 1 (end)
    """
    if num_sections == 1:
        return [0.5]
    elif num_sections == 2:
        return [0.0, 1.0]
    else:
        # For more sections, use equally spaced points including ends
        return [i/(num_sections-1) for i in range(num_sections)]


def plot(output_dir, modal_props, Tn, Sa, eigs, nmode, story_heights):
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


    # Center of mass is directly available
    center_of_mass = modal_props["centerOfMass"]
    print(f"Center of Mass: X={center_of_mass[0]}, Y={center_of_mass[1]}, Z={center_of_mass[2]}")


    # To estimate center of stiffness, we'll use the first two translational modes
    # Get participation factors
    parti_factor_MX = modal_props["partiFactorMX"]
    parti_factor_MY = modal_props["partiFactorMY"]
    parti_factor_RMZ = modal_props["partiFactorRMZ"]

    # Find first strong X and Y translational modes
    x_mode = None
    y_mode = None
    for i in range(len(parti_factor_MX)):
        if abs(parti_factor_MX[i]) > 0.5:  # Threshold for significant X participation
            x_mode = i
            break
    for i in range(len(parti_factor_MY)):
        if abs(parti_factor_MY[i]) > 0.5:  # Threshold for significant Y participation
            y_mode = i
            break

    # Calculate stiffness eccentricities
    if x_mode is not None and y_mode is not None:
        # X-direction stiffness eccentricity (Y coordinate)
        ey = parti_factor_RMZ[x_mode]/parti_factor_MX[x_mode]
        
        # Y-direction stiffness eccentricity (X coordinate)
        ex = -parti_factor_RMZ[y_mode]/parti_factor_MY[y_mode]
        
        # Center of stiffness coordinates
        center_of_stiffness_x = center_of_mass[0] + ex
        center_of_stiffness_y = center_of_mass[1] + ey
        center_of_stiffness_z = center_of_mass[2]  # Typically same as COM in Z
        
        print(f"Center of Stiffness: X={center_of_stiffness_x}, Y={center_of_stiffness_y}, Z={center_of_stiffness_z}")
    else:
        print("Could not identify clear translational modes for stiffness center calculation")
    
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

    # Print results
    print("\nFloor Properties:")
    print("{:<10} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}".format(
        "Floor Z", "COM X", "COM Y", "Mass", "COS X", "COS Y", "Kx", "Ky", "Kr"))
    for z in sorted(floor_masses.keys()):
        com_x, com_y, _, mass = floor_masses[z]
        if z in floor_stiffness:
            cos_x, cos_y, _ = floor_stiffness[z]
            Kx, Ky, Kr = floor_stiffness_values[z]
        else:
            cos_x, cos_y = com_x, com_y  # Default to COM if no stiffness calc
            Kx, Ky, Kr = 0, 0, 0
        
        print("{:<10.2f} {:<15.3f} {:<15.3f} {:<15.3f} {:<15.3f} {:<15.3f} {:<15.3g} {:<15.3g} {:<15.3g}".format(
            z, com_x, com_y, mass, cos_x, cos_y, Kx, Ky, Kr))

    # # Save results to file
    # with open(os.path.join(output_dir, "floor_properties.txt"), "w") as f:
    #     f.write("Floor Properties:\n")
    #     f.write("{:<10} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}\n".format(
    #         "Floor Z", "COM X", "COM Y", "Mass", "COS X", "COS Y", "Kx", "Ky", "Kr"))
    #     for z in sorted(floor_masses.keys()):
    #         com_x, com_y, _, mass = floor_masses[z]
    #         if z in floor_stiffness:
    #             cos_x, cos_y, _ = floor_stiffness[z]
    #             Kx, Ky, Kr = floor_stiffness_values[z]
    #         else:
    #             cos_x, cos_y = com_x, com_y
    #             Kx, Ky, Kr = 0, 0, 0
            
    #         f.write("{:<10.2f} {:<15.3f} {:<15.3f} {:<15.3f} {:<15.3f} {:<15.3f} {:<15.3g} {:<15.3g} {:<15.3g}\n".format(
    #             z, com_x, com_y, mass, cos_x, cos_y, Kx, Ky, Kr))

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

    # define a recorder for the (use a higher precision otherwise the results
    # won't match with those obtained from eleResponse)
    filename = os.path.join(output_dir, 'ele_1_sec_1.txt')
    ops.recorder('Element', '-file', filename, '-closeOnWrite', '-precision', 16, '-ele', 1, 'section', '1', 'force')

    # some settings for the response spectrum analysis
    tsTag = 1 # use the timeSeries 1 as response spectrum function
    direction = 1 # excited DOF = Ux

    # currently we use same damping for each mode
    dmp = [0.05]*len(eigs)
    # we don't want to scale some modes...
    scalf = [1.0]*len(eigs)
    # CQC function

    # ========================================================================
    # TEST 00
    # run a response spectrum analysis for each mode.
    # then do modal combination in post-processing.
    # Use Tn and Sa lists
    # ========================================================================
    ops.responseSpectrumAnalysis(direction, '-Tn', *Tn, '-Sa', *Sa)

    # read the My values [3rd column] for each step 
    # (1 for each mode, they are section forces associated to each modal displacement)
    My = []
    with open(filename, 'r') as f:
        lines = f.read().split('\n')
        for line in lines:
            if len(line) > 0:
                tokens = line.split(' ')
                My.append(float(tokens[2]))

    # post process the results doing the CQC modal combination for the My response (3rd column in section forces)
    MyCQC = CQC(My, eigs, dmp, scalf)

    print('\n\nTEST 00:\nRun a Response Spectrum Analysis for all modes.')
    print('Do CQC combination in post processing.')
    print('Use Tn and Sa lists.\n')
    print('{0: >10}{1: >15}'.format('Mode', 'My'))
    for i in range(len(eigs)):
        print('{0: >10}{1: >15f}'.format(i+1, My[i]))
    print('{0: >10}{1: >15f}'.format('CQC', MyCQC))

    # ========================================================================
    # TEST 01
    # run a response spectrum analysis for each mode.
    # then do modal combination in post-processing.
    # Use a Path timeSeries to store the Tn-Sa pairs
    # ========================================================================
    ops.responseSpectrumAnalysis(tsTag, direction)

    # read the My values [3rd column] for each step 
    # (1 for each mode, they are section forces associated to each modal displacement)
    My = []
    with open(filename, 'r') as f:
        lines = f.read().split('\n')
        for line in lines:
            if len(line) > 0:
                tokens = line.split(' ')
                My.append(float(tokens[2]))

    # post process the results doing the CQC modal combination for the My response (3rd column in section forces)
    MyCQC = CQC(My, eigs, dmp, scalf)

    print('\n\nTEST 01:\nRun a Response Spectrum Analysis for all modes.')
    print('Do CQC combination in post processing.')
    print('Use a Path timeSeries to store Tn-Sa pairs.\n')
    print('{0: >10}{1: >15}'.format('Mode', 'My'))
    for i in range(len(eigs)):
        print('{0: >10}{1: >15f}'.format(i+1, My[i]))
    print('{0: >10}{1: >15f}'.format('CQC', MyCQC))


    ops.recorder('Element', '-file', os.path.join(output_dir, 'ele_1_sec_1.txt'), '-precision', 16, '-ele', 1, 'section', '1', 'force')

    # Create recorders for all elements
    ops.recorder('Element', '-file', os.path.join(output_dir, 'ElementForces.txt'), '-closeOnWrite', 
                '-precision', 16, '-eleRange', 1, 17, 
                'section', '1', 'force')
    
    # # Create recorders for all elements in JSON format
    # element_forces_json = os.path.join(output_dir, 'ElementForces.json')
    # ops.recorder('Element', '-json', '-file', element_forces_json, '-closeOnWrite', 
    #             '-precision', 16, '-eleRange', 1, 17, 
    #             'section', '1', 'force')

    # Run RSA for all modes (your existing code)
    for i in range(len(eigs)):
        ops.responseSpectrumAnalysis(direction, '-Tn', *Tn, '-Sa', *Sa, '-mode', i+1)
        ops.analyze(1)
    


    modal_forces, cqc_forces = extract_and_combine_forces(
        output_dir, Tn, Sa, direction, eigs, dmp, scalf,
        json_filepath='results2100.json'
    )

    extract_and_combine_forces_multiple_sections(output_dir, Tn, Sa, direction, eigs, dmp, scalf, 
                                               num_sections=5, json_filepath='RSA_Forces_MultiSection.json', 
                                               csv_filepath='RSA_Forces_MultiSection.csv')

    
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

def response_spectrum_analysis(output_dir, Tn, Sa, nmode):
    # the response spectrum function



    ops.timeSeries("Path",1,"-time", *Tn, "-values", *Sa)

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
        2: 3.0,      # First floor at 4m
        3: 6.0,      # Second floor at 8m
        # 4: 12.0,     # Third floor at 12m
        # ... add more stories as needed
    }
    plot(output_dir, modal_props, Tn, Sa, eigs, nmode, story_heights)
    

def extract_and_combine_forces(output_dir, Tn, Sa, direction, eigs, dmp, scalf, json_filepath='RSA_Forces.json', csv_filepath='RSA_Forces.csv'):
    """
    Extract all member forces from RSA and perform CQC combination
    
    Args:
        Tn (list): List of periods for response spectrum
        Sa (list): List of spectral accelerations
        direction (int): Excitation direction (1=X, 2=Y, 3=Z)
        eigs (list): Eigenvalues from modal analysis
        dmp (list): Damping ratios for each mode
        scalf (list): Scaling factors for each mode
        json_filepath (str): Path to save JSON results
        csv_filepath (str): Path to save CSV results
    """
    # Initialize dictionaries to store all forces
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
        
        # Get forces for all elements
        for ele_tag in element_tags:
            forces = ops.eleResponse(ele_tag, 'section', '1', 'force')
            modal_forces['P'][mode][ele_tag] = forces[0]
            modal_forces['Vy'][mode][ele_tag] = forces[1]
            modal_forces['Vz'][mode][ele_tag] = forces[2]
            modal_forces['T'][mode][ele_tag] = forces[3]
            modal_forces['My'][mode][ele_tag] = forces[4]
            modal_forces['Mz'][mode][ele_tag] = forces[5]
    
    # Perform CQC combination
    cqc_forces = {
        'P': {}, 'Vy': {}, 'Vz': {}, 
        'T': {}, 'My': {}, 'Mz': {}
    }
    
    for ele_tag in element_tags:
        # Extract modal forces for this element
        P_modes = [modal_forces['P'][m][ele_tag] for m in range(1, len(eigs)+1)]
        Vy_modes = [modal_forces['Vy'][m][ele_tag] for m in range(1, len(eigs)+1)]
        Vz_modes = [modal_forces['Vz'][m][ele_tag] for m in range(1, len(eigs)+1)]
        T_modes = [modal_forces['T'][m][ele_tag] for m in range(1, len(eigs)+1)]
        My_modes = [modal_forces['My'][m][ele_tag] for m in range(1, len(eigs)+1)]
        Mz_modes = [modal_forces['Mz'][m][ele_tag] for m in range(1, len(eigs)+1)]
        
        # Perform CQC combination
        cqc_forces['P'][ele_tag] = CQC(P_modes, eigs, dmp, scalf)
        cqc_forces['Vy'][ele_tag] = CQC(Vy_modes, eigs, dmp, scalf)
        cqc_forces['Vz'][ele_tag] = CQC(Vz_modes, eigs, dmp, scalf)
        cqc_forces['T'][ele_tag] = CQC(T_modes, eigs, dmp, scalf)
        cqc_forces['My'][ele_tag] = CQC(My_modes, eigs, dmp, scalf)
        cqc_forces['Mz'][ele_tag] = CQC(Mz_modes, eigs, dmp, scalf)


    # Save to JSON
    json_filepath = os.path.join(output_dir, json_filepath)
    with open(json_filepath, 'w') as f:
        json.dump({
            'modal_forces': modal_forces,
            'cqc_forces': cqc_forces,
            'eigenvalues': eigs,
            'damping': dmp,
            'scaling_factors': scalf
        }, f, indent=4)
    
    return modal_forces, cqc_forces

def extract_and_combine_forces_multiple_sections(output_dir, Tn, Sa, direction, eigs, dmp, scalf, 
                                               num_sections=10, json_filepath='RSA_Forces_MultiSection.json', 
                                               csv_filepath='RSA_Forces_MultiSection.csv'):
    """
    Extract all member forces from RSA at multiple integration points and perform CQC combination
    
    Args:
        Tn (list): List of periods for response spectrum
        Sa (list): List of spectral accelerations
        direction (int): Excitation direction (1=X, 2=Y, 3=Z)
        eigs (list): Eigenvalues from modal analysis
        dmp (list): Damping ratios for each mode
        scalf (list): Scaling factors for each mode
        num_sections (int): Number of integration points to extract (default=10)
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
            
            # Extract forces at multiple integration points
            for section in range(1, num_sections+1):
                try:
                    forces = ops.eleResponse(ele_tag, 'section', str(section), 'force')
                    modal_forces['P'][mode][ele_tag][section] = forces[0]
                    modal_forces['Vy'][mode][ele_tag][section] = forces[1]
                    modal_forces['Vz'][mode][ele_tag][section] = forces[2]
                    modal_forces['T'][mode][ele_tag][section] = forces[3]
                    modal_forces['My'][mode][ele_tag][section] = forces[4]
                    modal_forces['Mz'][mode][ele_tag][section] = forces[5]
                except:
                    # If section doesn't exist, store zeros or skip
                    print(f"Warning: Section {section} not available for element {ele_tag}")
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

    # Save to JSON
    json_filepath = os.path.join(output_dir, json_filepath)
    with open(json_filepath, 'w') as f:
        json.dump({
            'modal_forces': modal_forces,
            'cqc_forces': cqc_forces,
            'eigenvalues': eigs,
            'damping': dmp,
            'scaling_factors': scalf,
            'num_sections': num_sections
        }, f, indent=4)
    
    # Save to CSV for easier analysis
    csv_filepath = os.path.join(output_dir, csv_filepath)
    with open(csv_filepath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        writer.writerow(['Element', 'Section', 'P', 'Vy', 'Vz', 'T', 'My', 'Mz'])
        
        # Write CQC results
        for ele_tag in element_tags:
            for section in range(1, num_sections+1):
                writer.writerow([
                    ele_tag, section,
                    cqc_forces['P'][ele_tag][section],
                    cqc_forces['Vy'][ele_tag][section],
                    cqc_forces['Vz'][ele_tag][section],
                    cqc_forces['T'][ele_tag][section],
                    cqc_forces['My'][ele_tag][section],
                    cqc_forces['Mz'][ele_tag][section]
                ])
    
    return modal_forces, cqc_forces


def find_critical_forces(cqc_forces, element_tags):
    """
    Find critical (maximum) forces across all sections for each element
    
    Args:
        cqc_forces (dict): CQC combined forces from extract_and_combine_forces_multiple_sections
        element_tags (list): List of element tags
    
    Returns:
        dict: Dictionary with maximum forces and their locations
    """
    critical_forces = {}
    
    for ele_tag in element_tags:
        critical_forces[ele_tag] = {}
        
        for force_type in ['P', 'Vy', 'Vz', 'T', 'My', 'Mz']:
            # Find maximum force and its location
            max_force = 0.0
            max_section = 1
            
            for section in cqc_forces[force_type][ele_tag]:
                if abs(cqc_forces[force_type][ele_tag][section]) > abs(max_force):
                    max_force = cqc_forces[force_type][ele_tag][section]
                    max_section = section
            
            critical_forces[ele_tag][force_type] = {
                'value': max_force,
                'section': max_section
            }
    
    return critical_forces

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
    
    # Save to files
    # 1. Node coordinates
    with open(os.path.join(output_dir, 'node_coordinates.txt'), 'w') as f:
        f.write("Node\tX\tY\tZ\n")
        for node_tag, coords in node_coords.items():
            f.write(f"{node_tag}\t{coords[0]:.4f}\t{coords[1]:.4f}\t{coords[2]:.4f}\n")
    
    # 2. Node displacements (for each mode)
    with open(os.path.join(output_dir, 'node_displacements.txt'), 'w') as f:
        f.write("Node\tMode\tUX\tUY\tUZ\tRX\tRY\tRZ\n")
        for node_tag, modes in node_displacements.items():
            for mode, disp in modes.items():
                f.write(f"{node_tag}\t{mode}\t")
                f.write("\t".join([f"{d:.6e}" for d in disp]))
                f.write("\n")
    
    # 3. Node reactions
    with open(os.path.join(output_dir, 'node_reactions.txt'), 'w') as f:
        f.write("Node\tFX\tFY\tFZ\tMX\tMY\tMZ\n")
        for node_tag, react in node_reactions.items():
            f.write(f"{node_tag}\t")
            f.write("\t".join([f"{r:.6e}" for r in react]))
            f.write("\n")
    
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
    
    # Save to CSV for easier analysis
    csv_filepath = os.path.join(output_dir, csv_filepath)
    with open(csv_filepath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write displacement results
        writer.writerow(['Node', 'Ux', 'Uy', 'Uz', 'Rx', 'Ry', 'Rz'])
        for node_tag in node_tags:
            writer.writerow([
                node_tag,
                cqc_displacements['Ux'][node_tag],
                cqc_displacements['Uy'][node_tag],
                cqc_displacements['Uz'][node_tag],
                cqc_displacements['Rx'][node_tag],
                cqc_displacements['Ry'][node_tag],
                cqc_displacements['Rz'][node_tag]
            ])
        
        # Add separator
        writer.writerow([])
        writer.writerow(['REACTIONS'])
        writer.writerow(['Node', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'])
        
        # Write reaction results
        for node_tag in node_tags:
            writer.writerow([
                node_tag,
                cqc_reactions['Fx'][node_tag],
                cqc_reactions['Fy'][node_tag],
                cqc_reactions['Fz'][node_tag],
                cqc_reactions['Mx'][node_tag],
                cqc_reactions['My'][node_tag],
                cqc_reactions['Mz'][node_tag]
            ])
    
    return modal_displacements, modal_reactions, cqc_displacements, cqc_reactions, story_drifts


def find_critical_nodal_responses(cqc_displacements, cqc_reactions, node_tags):
    """
    Find critical (maximum) nodal displacements and reactions
    
    Args:
        cqc_displacements (dict): CQC combined displacements
        cqc_reactions (dict): CQC combined reactions
        node_tags (list): List of node tags
    
    Returns:
        dict: Dictionary with maximum responses and their locations
    """
    critical_responses = {
        'displacements': {},
        'reactions': {}
    }
    
    # Find critical displacements
    for disp_type in ['Ux', 'Uy', 'Uz', 'Rx', 'Ry', 'Rz']:
        max_disp = 0.0
        max_node = None
        
        for node_tag in node_tags:
            if abs(cqc_displacements[disp_type][node_tag]) > abs(max_disp):
                max_disp = cqc_displacements[disp_type][node_tag]
                max_node = node_tag
        
        critical_responses['displacements'][disp_type] = {
            'value': max_disp,
            'node': max_node
        }
    
    # Find critical reactions
    for react_type in ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']:
        max_react = 0.0
        max_node = None
        
        for node_tag in node_tags:
            if abs(cqc_reactions[react_type][node_tag]) > abs(max_react):
                max_react = cqc_reactions[react_type][node_tag]
                max_node = node_tag
        
        critical_responses['reactions'][react_type] = {
            'value': max_react,
            'node': max_node
        }
    
    return critical_responses


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
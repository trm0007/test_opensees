from import_ import *     
from input import *
from load_combinations import *
from Grid_and_structure_creation import *
from test3 import ModalAnalysis
from wall_meshing import *
from sections_function import *
from units import *
import matplotlib.pyplot as plt
from ipywidgets import widgets
from run_function import *
from post_processing import *





import openseespy.opensees as ops
import numpy as np


def modal_analysis_with_participation():
    """
    Performs modal analysis and calculates mass participation factors.
    Only considers X and Y directions.
    CORRECTED VERSION - fixes the double counting issue in mass participation calculation.
    """
    print("=== MODAL ANALYSIS WITH MASS PARTICIPATION FACTORS ===")

    # Set up the system for eigenvalue analysis
    ops.wipeAnalysis()
    ops.system('BandGeneral')  # or 'FullGeneral' for smaller models
    ops.numberer('RCM')
    ops.constraints('Transformation')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')

    # Perform eigenvalue analysis
    Nmodos = 3
    eigenValues = np.array(ops.eigen(Nmodos))
    ops.modalProperties('-print', '-file', 'ModalReport.txt', '-unorm')
    # Check if eigenvalue analysis was successful
    if len(eigenValues) == 0:
        print("ERROR: Eigenvalue analysis failed!")
        return

    # Calculate frequencies and periods
    w_i = np.sqrt(eigenValues)  # Use abs to handle potential numerical errors
    T_i = 2 * np.pi / w_i
    f_i = w_i / (2 * np.pi)

    print("Modal Analysis Results:")
    for i in range(len(eigenValues)):
        print(f"Mode {i+1}:")
        print(f"  Eigenvalue: {eigenValues[i]:.6e}")
        print(f"  Frequency: {f_i[i]:.4f} Hz")
        print(f"  Period: {T_i[i]:.4f} sec")

    print("\n=== MODAL MASS PARTICIPATION FACTORS ===")

    # Get all node tags
    node_tags = ops.getNodeTags()
    print(f"Total nodes: {len(node_tags)}")

    if len(node_tags) == 0:
        print("ERROR: No nodes found in the model!")
        return

    # Initialize arrays for mass participation calculation
    total_mass_x = 0.0
    total_mass_y = 0.0

    # Arrays to store modal participation factors
    L_x = np.zeros(Nmodos)  # Modal participation factor in X direction
    L_y = np.zeros(Nmodos)  # Modal participation factor in Y direction

    M_eff_x = np.zeros(Nmodos)  # Effective modal mass in X direction
    M_eff_y = np.zeros(Nmodos)  # Effective modal mass in Y direction

    # Get nodal masses and calculate total mass
    nodal_masses = {}
    valid_nodes = []
    
    for node in node_tags:
        mass = ops.nodeMass(node)
        print(f"Node {node}: mass = {mass}")

        
        if mass and len(mass) >= 2:  # At least 2D
            nodal_masses[node] = mass[:2]  # Take only X and Y components
            valid_nodes.append(node)
            
            total_mass_x += mass[0]
            total_mass_y += mass[1]
        else:
            print(f"Warning: Node {node} has no mass or invalid mass definition")

    print(f"Valid nodes with mass: {len(valid_nodes)}")
    print(f"Total mass in X direction: {total_mass_x:.6f}")
    print(f"Total mass in Y direction: {total_mass_y:.6f}")

    if len(valid_nodes) == 0:
        print("ERROR: No nodes with valid mass found!")
        return

    if total_mass_x <= 0 or total_mass_y <= 0:
        print("ERROR: Total mass must be positive in both directions!")
        return

    # Calculate modal participation factors for each mode
    for mode in range(Nmodos):
        mode_num = mode + 1
        
        # Initialize sums for this mode
        sum_m_phi_x = 0.0  # Σ(m_i * φ_i,x)
        sum_m_phi_y = 0.0  # Σ(m_i * φ_i,y)
        
        # CORRECTED: Separate modal mass calculations for each direction
        modal_mass_x = 0.0  # Σ(m_i * φ_i,x^2) - for X direction normalization
        modal_mass_y = 0.0  # Σ(m_i * φ_i,y^2) - for Y direction normalization
        
        # Get mode shape for all nodes
        mode_shapes = {}
        for node in valid_nodes:
            # Get the mode shape (eigenvector) for this node and mode
            phi = ops.nodeEigenvector(node, mode_num)
            
            if phi and len(phi) >= 2:
                mode_shapes[node] = phi[:2]  # Take only X and Y components
            else:
                print(f"Warning: Invalid mode shape for node {node}, mode {mode_num}")
                continue
        
        # Calculate participation factors and modal masses
        for node in valid_nodes:
            if node not in mode_shapes:
                continue
                
            mass = nodal_masses[node]
            phi = mode_shapes[node]
            
            # Modal participation factor components: L = Σ(m_i * φ_i)
            sum_m_phi_x += mass[0] * phi[0]
            sum_m_phi_y += mass[1] * phi[1]
            
            # CORRECTED: Separate modal mass for each direction
            modal_mass_x += mass[0] * phi[0]**2
            modal_mass_y += mass[1] * phi[1]**2
        
        # Store modal participation factors
        L_x[mode] = sum_m_phi_x
        L_y[mode] = sum_m_phi_y
        
        # CORRECTED: Calculate effective modal masses separately for each direction
        # M_eff = L^2 / (φ^T * M * φ) where φ^T * M * φ is calculated per direction
        if modal_mass_x > 1e-12:  # Avoid division by very small numbers
            M_eff_x[mode] = (sum_m_phi_x**2) / modal_mass_x
        else:
            print(f"Warning: Very small modal mass in X direction for mode {mode_num}")
            M_eff_x[mode] = 0.0
            
        if modal_mass_y > 1e-12:  # Avoid division by very small numbers
            M_eff_y[mode] = (sum_m_phi_y**2) / modal_mass_y
        else:
            print(f"Warning: Very small modal mass in Y direction for mode {mode_num}")
            M_eff_y[mode] = 0.0

    # Calculate mass participation ratios (percentages)
    mass_participation_x = (M_eff_x / total_mass_x) * 100 
    mass_participation_y = (M_eff_y / total_mass_y) * 100 

    # Print results
    print("\nModal Mass Participation Factors:")
    print("=" * 60)
    print(f"{'Mode':<6} {'Freq(Hz)':<10} {'Period(s)':<10} {'UX(%)':<10} {'UY(%)':<10}")
    print("=" * 60)

    cumulative_x = 0.0
    cumulative_y = 0.0

    for i in range(Nmodos):
        cumulative_x += mass_participation_x[i]
        cumulative_y += mass_participation_y[i]
        
        print(f"{i+1:<6} {f_i[i]:<10.4f} {T_i[i]:<10.4f} {mass_participation_x[i]:<10.2f} "
              f"{mass_participation_y[i]:<10.2f}")

    print("=" * 60)
    print(f"{'Total':<26} {cumulative_x:<10.2f} {cumulative_y:<10.2f}")

    # Additional detailed output
    print("\nDetailed Modal Information:")
    print("=" * 60)
    for i in range(Nmodos):
        print(f"\nMode {i+1}:")
        print(f"  Eigenvalue: {eigenValues[i]:.6e}")
        print(f"  Frequency: {f_i[i]:.4f} Hz")
        print(f"  Period: {T_i[i]:.4f} sec")
        print(f"  Modal participation factor Lx: {L_x[i]:.6f}")
        print(f"  Modal participation factor Ly: {L_y[i]:.6f}")
        print(f"  Effective modal mass Mx: {M_eff_x[i]:.6f}")
        print(f"  Effective modal mass My: {M_eff_y[i]:.6f}")
        print(f"  Mass participation UX: {mass_participation_x[i]:.2f}%")
        print(f"  Mass participation UY: {mass_participation_y[i]:.2f}%")

    # Verification: Check that effective modal masses sum correctly
    print(f"\nVerification:")
    print(f"Sum of effective modal masses in X: {np.sum(M_eff_x):.6f}")
    print(f"Total mass in X: {total_mass_x:.6f}")
    print(f"Sum of effective modal masses in Y: {np.sum(M_eff_y):.6f}")
    print(f"Total mass in Y: {total_mass_y:.6f}")

    # Check if sufficient mass participation is achieved
    print("\n" + "=" * 60)
    print("MASS PARTICIPATION SUMMARY:")
    print("=" * 60)
    print(f"Cumulative mass participation in X direction: {cumulative_x:.2f}%")
    print(f"Cumulative mass participation in Y direction: {cumulative_y:.2f}%")

    # Check adequacy (typically 90% is required)
    adequacy_threshold = 90.0
    
    if cumulative_x >= adequacy_threshold:
        print(f"✓ X direction: Sufficient mass participation (≥{adequacy_threshold}%)")
    else:
        print(f"⚠ X direction: Insufficient mass participation (<{adequacy_threshold}%)")
        print(f"  Consider increasing the number of modes or checking model stiffness in X direction")
        
    if cumulative_y >= adequacy_threshold:
        print(f"✓ Y direction: Sufficient mass participation (≥{adequacy_threshold}%)")
    else:
        print(f"⚠ Y direction: Insufficient mass participation (<{adequacy_threshold}%)")
        print(f"  Consider increasing the number of modes or checking model stiffness in Y direction")

    # Additional warnings for potential issues
    print("\n" + "=" * 60)
    print("DIAGNOSTIC INFORMATION:")
    print("=" * 60)
    
    # Check for very similar frequencies (potential numerical issues)
    freq_tolerance = 1e-6
    for i in range(len(f_i)-1):
        if abs(f_i[i] - f_i[i+1]) < freq_tolerance:
            print(f"⚠ Warning: Modes {i+1} and {i+2} have very similar frequencies")
            print(f"  This might indicate numerical issues or repeated eigenvalues")
    
    # Check for very low participation in any direction
    if cumulative_x < 50.0:
        print(f"⚠ Warning: Very low mass participation in X direction ({cumulative_x:.2f}%)")
        print(f"  This might indicate excessive stiffness or modeling issues in X direction")
    
    if cumulative_y < 50.0:
        print(f"⚠ Warning: Very low mass participation in Y direction ({cumulative_y:.2f}%)")
        print(f"  This might indicate excessive stiffness or modeling issues in Y direction")

    # Return results for further analysis if needed
    results = {
        'eigenvalues': eigenValues,
        'frequencies': f_i,
        'periods': T_i,
        'modal_participation_factors': {'Lx': L_x, 'Ly': L_y},
        'effective_modal_masses': {'Mx': M_eff_x, 'My': M_eff_y},
        'mass_participation_percentages': {'UX': mass_participation_x, 'UY': mass_participation_y},
        'cumulative_participation': {'X': cumulative_x, 'Y': cumulative_y},
        'total_masses': {'X': total_mass_x, 'Y': total_mass_y}
    }
    
    return results

# # Alternative function for debugging mode shapes and masses
# def debug_modal_analysis(mode_num=1):
#     """
#     Debug function to examine mode shapes and masses in detail for a specific mode.
#     """
#     print(f"=== DEBUGGING MODE {mode_num} ===")
    
#     node_tags = ops.getNodeTags()
#     print(f"Total nodes: {len(node_tags)}")
    
#     # Check a few nodes in detail
#     sample_nodes = node_tags[:min(5, len(node_tags))]
    
#     for node in sample_nodes:
#         mass = ops.nodeMass(node)
#         phi = ops.nodeEigenvector(node, mode_num)
        
#         print(f"\nNode {node}:")
#         print(f"  Mass: {mass}")
#         print(f"  Mode shape: {phi}")
        
#         if mass and phi and len(mass) >= 2 and len(phi) >= 2:
#             print(f"  Mass*phi_x: {mass[0] * phi[0]:.6f}")
#             print(f"  Mass*phi_y: {mass[1] * phi[1]:.6f}")
#             print(f"  Mass*phi_x^2: {mass[0] * phi[0]**2:.6f}")
#             print(f"  Mass*phi_y^2: {mass[1] * phi[1]**2:.6f}")





def final_run(JSON_FOLDER, materials_list, section_definitions, materials_config, sections_config, opensees_element_json_file, combo_name):
    # Initialize OpenSees model
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    
    # =============================================
    # Node Creation
    # =============================================
    node_name_to_id = assign_node_into_opensees(JSON_FOLDER)
    
    # =============================================
    # Support Creation
    # =============================================
    apply_fixity_conditions(fixity_data, JSON_FOLDER)


    
    # =============================================
    # Material Definition
    # =============================================
    for mat in materials_list:
        if mat['ID'] == 'Concrete02':
            ops.uniaxialMaterial('Concrete02', mat['matTag'], 
                                mat['fpc'], mat['epsc0'], 
                                mat['fpcu'], mat['epsU'],
                                mat['lamda'], mat['ft'], mat['Ets'])
        elif mat['ID'] == 'Steel02':
            ops.uniaxialMaterial('Steel02', mat['matTag'],
                                mat['Fy'], mat['E0'],
                                mat['b'], mat['R0'],
                                mat['cR1'], mat['cR2'])
    
    # =============================================
    # Section Assignment
    # =============================================
    section_name_to_id = assign_member_section_into_opensees(section_definitions, materials)
    
    # # =============================================
    # # Shell Material and Section Assignment
    # # =============================================
    material_name_to_id, shell_section_name_to_id = assign_shell_material_and_section_into_opensees(materials_config, sections_config)
    
    # =============================================
    # Beam & Column Element Assignment
    # =============================================
    assign_beam_column_into_opensees(opensees_element_json_file)
    
    # # =============================================
    # # Shell Element Assignment
    # # =============================================
    # assign_shell_element_into_opensees(JSON_FOLDER, shell_section_name_to_id)
    assign_shell_element_into_opensees(JSON_FOLDER, sections_config)
    
    print("Model setup completed successfully")
    # # =============================================
    # # Create rigid diaphragms
    # # =============================================
     # # Step 2: Create rigid diaphragms
    # create_rigid_diaphragms(JSON_FOLDER)

    # # =============================================
    # # Loading and Analysis
    # # =============================================
    # 1. First create the time series and load pattern
    ops.timeSeries('Linear', 1, '-factor', 1.0)
    ops.pattern('Plain', 1, 1)  # Create load pattern with tag 1

    assign_opensees_loads(JSON_FOLDER, combo_name=combo_name)
    calculate_and_apply_nodal_masses(JSON_FOLDER, g)


    # run_gravity(JSON_FOLDER, steps=10, comb_name=combo_name)
    # run_gravity(JSON_FOLDER, steps=10, nep=5, OUTPUT_FOLDER = "postprocessing_folder", load_combination=combo_name)
    
    # # extract_all_element_data(JSON_FOLDER, nep=5, output_folder="postprocessing_folder", load_combination=combo_name)
    # get_node_results(OUTPUT_FOLDER = "postprocessing_folder", load_combination=combo_name)
    # calculate_slab_reinforcement_from_shell_forces(JSON_FOLDER, load_combination=combo_name)
    # structural_model_plot(OUTPUT_FOLDER="postprocessing_folder", load_combination=combo_name)

    
    # modal_analysis_with_participation()
    # debug_modal_analysis(mode_num=1)

        # Run your existing modal analysis first
    results = modal_analysis_with_participation()
    

    # results = analyze_modal_properties(json_folder=JSON_FOLDER, mode_number=1, z_spacing=z_spacing,)
    # for story_num in results['stiffnesses']:
    #     print(f"Story {story_num}:")
    #     print(f"  Stiffness: {results['stiffnesses'][story_num]:.2e}")
    #     print(f"  Mass: {results['masses'][story_num]:.2f}")
    #     print(f"  CoM: {results['centers_of_mass'][story_num]}")
    #     print(f"  CoS: {results['centers_of_stiffness'][story_num]}")
    #     print(f"  Eccentricity: {results['eccentricities'][story_num]}")
 


combo_name = "unfactored_load"
OUTPUT_FOLDER = "output_folder"  # Main output directory
os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist
# Subdirectories
JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
os.makedirs(JSON_FOLDER, exist_ok=True)
os.makedirs(IMAGE_FOLDER, exist_ok=True)
opensees_element_json_file = os.path.join(JSON_FOLDER, "element_data.json")
final_run(JSON_FOLDER, materials, section_definitions, materials_config, sections_config, opensees_element_json_file, combo_name)






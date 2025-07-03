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
from test4 import CQC, calculate_floor_masses, calculate_floor_stiffness, extract_node_data, model, plot, response_spectrum_analysis



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
    # # 1. First create the time series and load pattern
    # ops.timeSeries('Linear', 1, '-factor', 1.0)
    # ops.pattern('Plain', 1, 1)  # Create load pattern with tag 1

    # assign_opensees_loads(JSON_FOLDER, combo_name=combo_name)
    # calculate_and_apply_nodal_masses(JSON_FOLDER, g)


    # run_gravity(JSON_FOLDER, steps=10, comb_name=combo_name)
    # run_gravity(JSON_FOLDER, steps=10, nep=5, OUTPUT_FOLDER = "postprocessing_folder", load_combination=combo_name)

    Tn = [0.0, 0.06, 0.1, 0.12, 0.18, 0.24, 0.3, 0.36, 0.4, 0.42, 
    0.48, 0.54, 0.6, 0.66, 0.72, 0.78, 0.84, 0.9, 0.96, 1.02, 
    1.08, 1.14, 1.2, 1.26, 1.32, 1.38, 1.44, 1.5, 1.56, 1.62, 
    1.68, 1.74, 1.8, 1.86, 1.92, 1.98, 2.04, 2.1, 2.16, 2.22, 
    2.28, 2.34, 2.4, 2.46, 2.52, 2.58, 2.64, 2.7, 2.76, 2.82, 
    2.88, 2.94, 3.0, 3.06, 3.12, 3.18, 3.24, 3.3, 3.36, 3.42, 
    3.48, 3.54, 3.6, 3.66, 3.72, 3.78, 3.84, 3.9, 3.96, 4.02, 
    4.08, 4.14, 4.2, 4.26, 4.32, 4.38, 4.44, 4.5, 4.56, 4.62, 
    4.68, 4.74, 4.8, 4.86, 4.92, 4.98, 5.04, 5.1, 5.16, 5.22, 
    5.28, 5.34, 5.4, 5.46, 5.52, 5.58, 5.64, 5.7, 5.76, 5.82, 
    5.88, 5.94, 6.0]

    Sa = [1.9612, 3.72628, 4.903, 4.903, 4.903, 4.903, 4.903, 4.903, 4.903, 4.6696172, 
        4.0861602, 3.6321424, 3.2683398, 2.971218, 2.7241068, 2.5142584, 2.3348086, 2.1788932, 2.0425898, 1.9229566, 
        1.8160712, 1.7199724, 1.6346602, 1.5562122, 1.485609, 1.4208894, 1.3620534, 1.3071398, 1.2571292, 1.211041, 
        1.166914, 1.1267094, 1.0894466, 1.054145, 1.0217852, 0.990406, 0.960988, 0.9335312, 0.9080356, 0.8835206, 
        0.8599862, 0.838413, 0.8168398, 0.7972278, 0.7785964, 0.759965, 0.7432948, 0.7266246, 0.710935, 0.6952454, 
        0.6805364, 0.666808, 0.6540602, 0.6285646, 0.6040496, 0.5814958, 0.5609032, 0.5403106, 0.5206986, 0.5030478, 
        0.485397, 0.4697074, 0.4540178, 0.4393088, 0.4255804, 0.411852, 0.3991042, 0.3863564, 0.3755698, 0.3638026, 
        0.353016, 0.34321, 0.333404, 0.3245786, 0.3157532, 0.3069278, 0.2981024, 0.2902576, 0.2833934, 0.2755486, 
        0.2686844, 0.2618202, 0.254956, 0.2490724, 0.2431888, 0.2373052, 0.2314216, 0.2265186, 0.220635, 0.215732, 
        0.210829, 0.205926, 0.2020036, 0.1971006, 0.1931782, 0.1892558, 0.1853334, 0.181411, 0.1774886, 0.1735662, 
        0.1706244, 0.166702, 0.1637602]


    output_dir = "modal_analysis"
    os.makedirs(output_dir, exist_ok=True)
    # model()
    response_spectrum_analysis(output_dir, JSON_FOLDER, Tn, Sa, 5)
    
    # # extract_all_element_data(JSON_FOLDER, nep=5, output_folder="postprocessing_folder", load_combination=combo_name)
    # get_node_results(OUTPUT_FOLDER = "postprocessing_folder", load_combination=combo_name)
    # calculate_slab_reinforcement_from_shell_forces(JSON_FOLDER, load_combination=combo_name)
    # structural_model_plot(OUTPUT_FOLDER="postprocessing_folder", load_combination=combo_name)

    



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






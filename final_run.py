from import_ import *     
from input import *
from load_combinations import *
from Grid_and_structure_creation import *
from wall_meshing import *
from sections_function import *
from units import *
import matplotlib.pyplot as plt
from ipywidgets import widgets
from run_function import *
from post_processing import *



def final_run(JSON_FOLDER, materials_list, section_definitions, materials_config, sections_config, opensees_element_json_file):
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
    # # Loading and Analysis
    # # =============================================
    # 1. First create the time series and load pattern
    ops.timeSeries('Linear', 1, '-factor', 1.0)
    ops.pattern('Plain', 1, 1)  # Create load pattern with tag 1

    assign_opensees_loads(JSON_FOLDER, combo_name="Comb2")
    eleTags = [27, 28, 32]  # Replace with actual element tags

    # Define a function to apply loads on multiple elements

    # Apply loads
    # apply_3d_beam_loads(eleTags, Wy=-2.0, Wz=-3.0, Wx=0.0, Py=-5.0, Pz=-10.0, Px=0.0, xL=0.4, pattern_tag=101)




    
    # fig = plt.figure(figsize=(10, 8))
    # ax = fig.add_subplot(111, projection='3d')

    # # Get all shell elements (in this case, just element 1)
    # shell_elements = ops.getEleTags()  # Returns [1]

    # for ele_tag in shell_elements:
    #     # Get node coordinates of the element
    #     ele_nodes = ops.eleNodes(ele_tag)
    #     node_coords = np.array([ops.nodeCoord(node) for node in ele_nodes])
        
    #     # Create a filled polygon
    #     poly = Poly3DCollection([node_coords], alpha=0.5, linewidth=1, edgecolor='k')
        
    #     # Assign color (modify logic as needed)
    #     poly.set_facecolor('yellow')  # Single color for all elements
        
    #     ax.add_collection3d(poly)

    # # Overlay the original model edges (optional)
    # opsv.plot_model(element_labels=0, node_labels=0, ax=ax, fmt_model={'color': 'k', 'linewidth': 1})
    # plt.show()
    # Example: define multiple elements


    run_gravity(JSON_FOLDER, steps=10, comb_name="Comb2")
    
    extract_all_element_data(JSON_FOLDER, nep=5, output_folder="postprocessing_folder", load_combination="combo")
    get_node_results(OUTPUT_FOLDER = "postprocessing_folder", load_combination="combo")
    calculate_slab_reinforcement_from_shell_forces(JSON_FOLDER, load_combination="combo")
    structural_model_plot(OUTPUT_FOLDER="postprocessing_folder", load_combination="combo")
       


OUTPUT_FOLDER = "output_folder"  # Main output directory
os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist
# Subdirectories
JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
os.makedirs(JSON_FOLDER, exist_ok=True)
os.makedirs(IMAGE_FOLDER, exist_ok=True)
opensees_element_json_file = os.path.join(JSON_FOLDER, "element_data.json")
final_run(JSON_FOLDER, materials, section_definitions, materials_config, sections_config, opensees_element_json_file)






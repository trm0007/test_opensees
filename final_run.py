from import_ import *     
from input import *
from load_combinations import merge_structures
from Grid_and_structure_creation import load_structure
from wall_meshing import create_combined_structure_json
from sections_function import *
from units import *
import matplotlib.pyplot as plt
from ipywidgets import widgets
from run_function import *
from post_processing import *

def final_run(JSON_FOLDER, materials_list, section_definitions, materials_config, sections_config, opensees_element_json_file):

    
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    
    # =============================================
    # Node Creation
    # =============================================
    assign_node_into_opensees(JSON_FOLDER)

    # =============================================
    # Support Creation
    # =============================================
    for group in fixity_data.values():
        for nodeTag in group["nodes"]:
            ops.fix(nodeTag, *group["vals"])

    # =============================================
    # Material Definition
    # =============================================
    define_materials(materials_list)

    # # =============================================
    # # Section Assignment
    # # =============================================
    assign_member_section_into_opensees(section_definitions)

    # # =============================================
    # # Beam & Column Element Assignment
    # # =============================================
    assign_beam_column_into_opensees(opensees_element_json_file)
    
    # # =============================================
    # # Shell Material and Section and Shell Element Assignment
    # # =============================================
    try:
        assign_shell_material_and_section_into_opensees(materials_config, sections_config)
    except Exception as e:
        print(f"Error assigning shell material and section: {e}")

    # # =============================================
    # # Shell Element Assignment
    # # =============================================

    try:
        assign_shell_element_into_opensees(JSON_FOLDER, section_config=sections_config)
    except Exception as e:
        print(f"Error assigning shell material and section: {e}")
    # assign_shell_element_into_opensees(JSON_FOLDER, section_config=sections_config)

    # # =============================================
    # # Loading and Analysis
    # # =============================================
    ops.timeSeries('Linear', 1, '-factor', 1.0)
    # Remove the ops.pattern() call here
    # apply_nodal_loads(JSON_FOLDER, load_combination="Comb2", load_pattern_tag=1, time_series_tag=1)
    # apply_member_loads(JSON_FOLDER, load_case_name="Comb2")
    # process_structure_loads(
    #     JSON_FOLDER, 
    #     operation_type="nodal_loads",
    #     load_cases=load_cases,
    #     load_combinations=load_combinations,
    #     numbering=1,
    #     load_combination="Comb2"
    # )

    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Get all shell elements (in this case, just element 1)
    shell_elements = ops.getEleTags()  # Returns [1]

    for ele_tag in shell_elements:
        # Get node coordinates of the element
        ele_nodes = ops.eleNodes(ele_tag)
        node_coords = np.array([ops.nodeCoord(node) for node in ele_nodes])
        
        # Create a filled polygon
        poly = Poly3DCollection([node_coords], alpha=0.5, linewidth=1, edgecolor='k')
        
        # Assign color (modify logic as needed)
        poly.set_facecolor('yellow')  # Single color for all elements
        
        ax.add_collection3d(poly)

    # Overlay the original model edges (optional)
    opsv.plot_model(element_labels=0, node_labels=0, ax=ax, fmt_model={'color': 'k', 'linewidth': 1})
    plt.show()


    # run_gravity(JSON_FOLDER, steps=10, comb_name="Comb2")

    # # # =============================================
    # # # Visualization
    # # # =============================================
    
    # # Extruded shapes visualization
    # # Initialize the ele_shapes dictionary
    # ele_shapes = {}
    # # ele_shapes = {
    # # 1: ['rect', [200, 300]],           # Rectangular
    # # 2: ['I', [200, 400, 12, 20]],      # I-beam  
    # # 3: ['L', [300, 200, 25]],          # L-shape: H=300, B=200, t=25mm
    # # 4: ['L', [400, 400, 30]],          # Equal leg L-shape
    # # }
    # # Extract section types dynamically
    # section_types = section_definitions.keys()

    # # Process each section type
    # for section_type in section_types:
    #     sections = section_definitions[section_type]
        
    #     for section_name, section_data in sections.items():
    #         section_tag = section_data['section_tag']
            
    #         if section_data['type'] == 'rectangular':
    #             # Convert dimensions to feet (assuming inch is defined as 1/12)
    #             B = section_data['B'] / 12.0  # width in feet
    #             H = section_data['H'] / 12.0  # height in feet
    #             ele_shapes[section_tag] = ['rect', [B, H]]
                
    #         elif section_data['type'] == 'circular':
    #             D = section_data['D_Sec'] / 12.0  # diameter in feet
    #             ele_shapes[section_tag] = ['circ', [D]]
            
    #         elif section_data['type'] == 'L':
    #             H = section_data['H'] / 12.0      # height in feet
    #             B = section_data['B'] / 12.0      # width in feet
    #             t = section_data['t'] / 12.0      # thickness in feet
    #             ele_shapes[section_tag] = ['L', [H, B, t]]

    # print("ele_shapes = {")
    # for key, value in sorted(ele_shapes.items()):
    #     print(f"    {key}: {value},")
    # print("}")

    # structural_model_plot(ele_shapes, OUTPUT_FOLDER="postprocessing_folder")

    # # Save node results to JSON
    # get_node_results()

    # print("\nElement Types in Model:")
    # # 3. Plot axial forces

    # plot_beam_forces(force_type='N', sfac=0.5)  # Axial force diagram
    # plot_beam_forces(force_type='Vy', sfac=1.0)  # Shear force diagram in Y direction 
    # plot_beam_forces(force_type='Vz', sfac=0.5)  # Shear force diagram in Z direction
    # plot_beam_forces(force_type='T', sfac=0.1)   # Torsional moment diagram
    # plot_beam_forces(force_type='My', sfac=0.1 , view_angle=(30, 45), figsize=(15, 12))  # Bending moment diagram about Y axis
    # plot_beam_forces(force_type='Mz', sfac=0.5)  # Bending moment diagram about Z axis
    # calculate_beam_stress(load_combination="combo", num_points=10, view_angle=None, figsize=(12, 10))
    # shell_analysis(JSON_FOLDER)
    # deflection_data = calculate_beam_deflection(
    #     load_combination="combo2",
    #     num_points=20,
    #     view_angle=(30, 45),
    #     figsize=(15, 12)
    # )
    # # shell_deflection_data = calculate_shell_deflection(
    # #     load_combination="combo2",
    # #     num_points_per_direction=5,
    # #     figsize=(18, 12)
    # # )
    # analyze_shell_deflections(load_case="combo2", samples_per_side=6)
    # calculate_slab_reinforcement_from_json()
       


OUTPUT_FOLDER = "output_folder"  # Main output directory
os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist
# Subdirectories
JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
os.makedirs(JSON_FOLDER, exist_ok=True)
os.makedirs(IMAGE_FOLDER, exist_ok=True)
opensees_element_json_file = os.path.join(JSON_FOLDER, "element_data.json")
final_run(JSON_FOLDER, materials, section_definitions, materials_config, sections_config, opensees_element_json_file)






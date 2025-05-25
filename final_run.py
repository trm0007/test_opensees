from import_ import *     
from input import *
from load_combinations import merge_structures, process_structure_loads
from Grid_and_structure_creation import load_structure
from wall_meshing import create_combined_structure_json
from sections_function import *
from units import *
import matplotlib.pyplot as plt
import opsvis as opsv
from ipywidgets import widgets
from run_function import *
from post_processing import *

def final_run(JSON_FOLDER, materials_list, section_definitions, materials_config, sections_config, opensees_element_json_file):
    import matplotlib.pyplot as plt
    
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
    assign_shell_material_and_section_into_opensees(materials_config, sections_config)

    # # =============================================
    # # Shell Element Assignment
    # # =============================================
    assign_shell_element_into_opensees(JSON_FOLDER, section_config=sections_config)

    # # =============================================
    # # Loading and Analysis
    # # =============================================
    ops.timeSeries('Linear', 1, '-factor', 1.0)
    ops.pattern('Plain', 1, 1)

    process_structure_loads(
        JSON_FOLDER, 
        operation_type="nodal_loads",
        load_cases=load_cases,
        load_combinations=load_combinations,
        numbering=1,
        load_combination="Comb2"
    )

    run_gravity(JSON_FOLDER, steps=10, comb_name="Comb2")

    # # =============================================
    # # Visualization
    # # =============================================
    
    # Extruded shapes visualization
    ele_shapes = {
        1: ['rect', [1.0, 1.3333333333333333]],  # Section1: B=1.0, H=1.333
        2: ['rect', [1.0, 1.3333333333333333]],  # Section1
        3: ['rect', [1.0, 1.3333333333333333]],  # Section1
        4: ['rect', [1.0, 1.3333333333333333]],  # Section1
        7: ['rect', [1.0, 2.0]],                 # Section2: B=1.0, H=2.0
        8: ['rect', [1.0, 2.0]],                 # Section2
        11: ['rect', [1.0, 2.0]],                # Section2
        12: ['rect', [1.0, 2.0]]                 # Section2
    }

    structural_model_plot(ele_shapes, OUTPUT_FOLDER="postprocessing_folder")

    # Save node results to JSON
    get_node_results()

    print("\nElement Types in Model:")
    # 3. Plot axial forces

    plot_beam_forces(force_type='N', sfac=0.5)  # Axial force diagram
    plot_beam_forces(force_type='Vy', sfac=1.0)  # Shear force diagram in Y direction 
    plot_beam_forces(force_type='Vz', sfac=0.5)  # Shear force diagram in Z direction
    plot_beam_forces(force_type='T', sfac=0.1)   # Torsional moment diagram
    plot_beam_forces(force_type='My', sfac=0.1 , view_angle=(30, 45), figsize=(15, 12))  # Bending moment diagram about Y axis
    plot_beam_forces(force_type='Mz', sfac=0.5)  # Bending moment diagram about Z axis
    calculate_beam_stress(load_combination="combo", num_points=10, view_angle=None, figsize=(12, 10))
    shell_analysis(JSON_FOLDER)
    deflection_data = calculate_beam_deflection(
        load_combination="combo2",
        num_points=20,
        view_angle=(30, 45),
        figsize=(15, 12)
    )
    # shell_deflection_data = calculate_shell_deflection(
    #     load_combination="combo2",
    #     num_points_per_direction=5,
    #     figsize=(18, 12)
    # )
    analyze_shell_deflections(load_case="combo2", samples_per_side=6)
    calculate_slab_reinforcement_from_json()
       


OUTPUT_FOLDER = "output_folder"  # Main output directory
os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist

# Subdirectories
JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
os.makedirs(JSON_FOLDER, exist_ok=True)
os.makedirs(IMAGE_FOLDER, exist_ok=True)
opensees_element_json_file = os.path.join(JSON_FOLDER, "element_data.json")
final_run(JSON_FOLDER, materials, section_definitions, materials_config, sections_config, opensees_element_json_file)






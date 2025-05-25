from import_ import *    
from load_combinations import *
from units import *
from create_structural_model import *
from input import *
from sections_function import *
from wall_meshing import *
from Grid_and_structure_creation import *


def build_model(x_spacing, y_spacing, z_spacing, materials, section_definitions, n_integration_pts=8):
    """
    Builds the 3D RC Frame model with proper section assignments and transformations
    
    Parameters:
    -----------
    x_spacing : list
        Spacing between nodes along X-axis
    y_spacing : list
        Spacing between nodes along Y-axis
    z_spacing : list
        Spacing between nodes along Z-axis (story heights)
    materials : dict
        Material properties for concrete and steel
    section_definitions : dict
        Section definitions for different structural elements
    n_integration_pts : int, optional
        Number of integration points for fiber sections (default=8)
    """
    
    # =============================================
    # Create Output Directory Structure
    # =============================================
    OUTPUT_FOLDER = "output_folder"  # Main output directory
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist

    # Subdirectories for different output types
    JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")  # For JSON data files
    IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")      # For visualization images
    os.makedirs(JSON_FOLDER, exist_ok=True)
    os.makedirs(IMAGE_FOLDER, exist_ok=True)

    # =============================================
    # Create Basic Structural Grid
    # =============================================
    # Generate nodes and members based on input spacing
    # nodes, members, x_points, y_points, z_points = create_structure(x_spacing, y_spacing, z_spacing)
    
    # # Save initial structure data to JSON
    # save_structure_to_json(nodes, members, JSON_FOLDER, "structure_data_json.json")

    
    # =============================================
    # Structure Modifications
    # =============================================
   
    
    # Update tolerance for geometric comparisons
     
    # Apply modifications to structure
    # updated_structure = update_structure(
    #     JSON_FOLDER,
    #     "structure_data_json.json",
    #     create_new_nodes=create_new_nodes,
    #     create_new_members=create_new_members,
    #     delete_nodes=delete_nodes,
    #     delete_members=delete_members,
    #   )
    update_structure(JSON_FOLDER, x_spacing, y_spacing, z_spacing, create_new_nodes={}, create_new_members={}, delete_nodes=[], delete_members=[], tolerance=1e-6)

    print(f"Structure updated successfully with {len(create_new_nodes)} new nodes and {len(create_new_members)} new members.")
    
    # # =============================================
    # # Load and Filter Structure Data
    # # =============================================
    # # Reload structure from JSON file
    # # nodes, members = load_structure_from_json(JSON_FOLDER, "structure_data_json.json")

    # nodes, members, structure_data = load_structure(JSON_FOLDER)
        
    # # Get unique Z levels (story heights)
    # z_points = get_z_levels(nodes)
    # print(f'z_points={z_points}')
    
    # # Filter structural elements by Z level
    # filtered_nodes = filter_nodes_by_z(nodes, z_points)
    # filtered_columns = filter_column_members_by_z(members, nodes, z_points)
    # filtered_beams_x = filter_beam_x_by_z(members, nodes, z_points)
    # filtered_beams_y = filter_beam_y_by_z(members, nodes, z_points)
        
    # # Save filtered data to separate JSON files
    # save_filtered_data(filtered_nodes, JSON_FOLDER, "filtered_nodes.json")
    # save_filtered_data(filtered_columns, JSON_FOLDER, "filtered_columns.json")
    # save_filtered_data(filtered_beams_x, JSON_FOLDER, "beamX.json")
    # save_filtered_data(filtered_beams_y, JSON_FOLDER, "beamY.json")

    # # =============================================
    # # Visualization of Structural Elements
    # # =============================================
    # plot_structure(nodes, members, IMAGE_FOLDER, title="Full Structure", save_path="full_structure.png")
    # # Generate plots for each Z level showing different element types
    # for z_level in z_points:
    #     # Plot nodes at current Z level
    #     plot_filtered_elements(
    #         {z_level: filtered_nodes.get(z_level, {})},
    #         nodes,
    #         IMAGE_FOLDER,
    #         title_prefix=f"Nodes at Z Level",
    #         save_path=f"grid/nodes_z{z_level}.png"
    #     )
        
    #     # Plot columns (only for levels below the top)
    #     if z_level < max(z_points):
    #         plot_filtered_elements(
    #             {z_level: filtered_columns.get(z_level, {})},
    #             nodes,
    #             IMAGE_FOLDER,
    #             title_prefix=f"Columns starting at Z Level",
    #             save_path=f"grid/columns_z{z_level}.png"
    #         )
        
    #     # Plot beams in X direction
    #     plot_filtered_elements(
    #         {z_level: filtered_beams_x.get(z_level, {})},
    #         nodes,
    #         IMAGE_FOLDER,
    #         title_prefix=f"X Beams at Z Level",
    #         save_path=f"grid/beams_x_z{z_level}.png"
    #     )
        
    #     # Plot beams in Y direction
    #     plot_filtered_elements(
    #         {z_level: filtered_beams_y.get(z_level, {})},
    #         nodes,
    #         IMAGE_FOLDER,
    #         title_prefix=f"Y Beams at Z Level",
    #         save_path=f"grid/beams_y_z{z_level}.png"
    #     )

    # =============================================
    # Section Definitions and Visualization
    # =============================================
    # Process all section definitions and generate visualizations
    for category in section_definitions.values():
        for section_name, section_info in category.items():
            if section_info["type"] == "rectangular":
                # Draw rectangular RC sections
                draw_rc_section10(
                    section_info["section_tag"],
                    section_info["core_tag"],
                    section_info["cover_tag"],
                    section_info["steel_tag"],
                    section_info["H"],
                    section_info["B"],
                    section_info["cover_H"],
                    section_info["cover_B"],
                    section_info["offset"],
                    section_info["n_bars_top"],
                    section_info["dia_top"],
                    section_info["n_bars_bot"],
                    section_info["dia_bot"],
                    section_info["n_bars_secondary_top"],
                    section_info["dia_sec_top"],
                    section_info["n_bars_secondary_bot"],
                    section_info["dia_sec_bot"],
                    section_info["n_bars_int"],
                    section_info["dia_int"],
                    IMAGE_FOLDER
                )
            elif section_info["type"] == "circular":
                # Draw circular RC sections
                draw_circular_rc_section10(
                    section_info["section_tag"],
                    section_info["core_tag"],
                    section_info["cover_tag"],
                    section_info["steel_tag"],
                    section_info["D_Sec"],
                    section_info["cover_Sec"],
                    section_info["num_Bars_Sec"],
                    section_info["bar_dia_Sec"],
                    IMAGE_FOLDER,
                    section_info["ri"],
                    section_info["nf_Core_R"],
                    section_info["nf_Core_T"],
                    section_info["nf_Cover_R"],
                    section_info["nf_Cover_T"]
                )

    # =============================================
    # Create Member-Section Mapping
    # =============================================
    mapping = create_member_section_mapping(
        JSON_FOLDER,
        'member_section_mapping.json'
    )

    # =============================================
    # Generate OpenSees Structural Model
    # =============================================
    opensees_element_json_file = os.path.join(JSON_FOLDER, "element_data.json")
    # Create the actual OpenSees model elements
    shapes, memberLengths = create_structural_model(section_definitions, JSON_FOLDER, opensees_element_json_file, numIntgrPts=8)

    # =============================================
    # Shell Mesh Generation (if needed)
    # =============================================
    # Process each surface configuration for shell elements
    if surface_configurations is not None and surface_configurations:
        for config_name, config_data in surface_configurations.items():
            print(f"\nProcessing {config_name} configuration ({config_data['section_name']})...")
            
            # Create shell mesh for this configuration
            mesh_data = create_proper_mesh_for_closed_area_3d1(
                points=config_data["points"],
                predefined_points=config_data["predefined_points"],
                shell_section_name=config_data["section_name"],
                JSON_FOLDER=JSON_FOLDER,
                IMAGE_FOLDER=IMAGE_FOLDER,
                num_x_div=config_data["num_x_div"],
                num_y_div=config_data["num_y_div"],
                numbering=list(surface_configurations.keys()).index(config_name) + 1,
                add_shell=config_data["add_shell"],
                remove_shell=config_data["remove_shell"]
            )

    # =============================================
    # Final Data Consolidation
    # =============================================
    # Combine all structural data into single JSON file
    # create_combined_structure_json(JSON_FOLDER, output_file="combined_structure_data.json")
    shell_elements_data = create_combined_structure_json(JSON_FOLDER)
    print(f'shell_elements_data={shell_elements_data}')
    # =============================================
    # Final Data Consolidation
    # =============================================

    # apply_shell_load(JSON_FOLDER, load_case_name="DL", numbering=1, pressure=-100.0)

    # Process each surface configuration for shell elements
    if surface_configurations is not None and surface_configurations:
        for config_name, config_data in surface_configurations.items():           
            # Apply shell loads for each load case
            numbering = list(surface_configurations.keys()).index(config_name) + 1
            apply_shell_load(
                JSON_FOLDER=JSON_FOLDER,
                load_case_names=config_data["load_case_names"],  # Plural
                pressures=config_data["pressures"],              # Plural
                numbering=numbering
            )
    apply_self_weight(JSON_FOLDER, load_case_names="self_weight")
    
    print("\nAll configurations processed successfully!")
    print("Model Built Successfully!")

build_model(x_spacing, y_spacing, z_spacing, materials, section_definitions, n_integration_pts=8) 


OUTPUT_FOLDER = "output_folder"  # Main output directory
os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist

# Subdirectories
JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
os.makedirs(JSON_FOLDER, exist_ok=True)
os.makedirs(IMAGE_FOLDER, exist_ok=True)
# opensees_element_json_file = os.path.join(JSON_FOLDER, "element_data.json")

# final_run(JSON_FOLDER, IMAGE_FOLDER, materials, section_definitions, opensees_element_json_file, json_filename="structure_data_json.json")

# OUTPUT_FOLDER = "postprocessing_folder"  # Main output directory
# os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist

# # Subdirectories
# JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
# IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
# os.makedirs(JSON_FOLDER, exist_ok=True)
# os.makedirs(IMAGE_FOLDER, exist_ok=True)


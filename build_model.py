from import_ import *    
from load_combinations import *
from units import *
from create_beam_column_member import *
from input import *
from sections_function import *
from wall_meshing import *
from Grid_and_structure_creation import *
from lateral_load_distribution import *


# def build_model(x_spacing, y_spacing, z_spacing, materials, section_definitions, OUTPUT_FOLDER = "output_folder", n_integration_pts=8):
#     """
#     Builds the 3D RC Frame model with proper section assignments and transformations
    
#     Parameters:
#     -----------
#     x_spacing : list
#         Spacing between nodes along X-axis
#     y_spacing : list
#         Spacing between nodes along Y-axis
#     z_spacing : list
#         Spacing between nodes along Z-axis (story heights)
#     materials : dict
#         Material properties for concrete and steel
#     section_definitions : dict
#         Section definitions for different structural elements
#     n_integration_pts : int, optional
#         Number of integration points for fiber sections (default=8)
#     """
    
#     # =============================================
#     # Create Output Directory Structure
#     # =============================================
#     # Main output directory
#     os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist

#     # Subdirectories for different output types
#     JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")  # For JSON data files
#     IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")      # For visualization images
#     os.makedirs(JSON_FOLDER, exist_ok=True)
#     os.makedirs(IMAGE_FOLDER, exist_ok=True)

#     # =============================================
#     # Create Basic Structural Grid
#     # =============================================
#     update_structure(JSON_FOLDER, x_spacing, y_spacing, z_spacing, create_new_nodes={}, create_new_members={}, delete_nodes=[], delete_members=[], tolerance=1e-6)

#     print(f"Structure updated successfully with {len(create_new_nodes)} new nodes and {len(create_new_members)} new members.")
    
#     # =============================================
#     # Load and Filter Structure Data
#     # =============================================
#     # Reload structure from JSON file
#     # nodes, members = load_structure_from_json(JSON_FOLDER, "structure_data_json.json")

#     nodes, members, structure_data = load_structure(JSON_FOLDER)
        
#     # Get unique Z levels (story heights)
#     z_points = get_z_levels(nodes)
#     print(f'z_points={z_points}')
    
#     # Filter structural elements by Z level
#     filtered_nodes = filter_nodes_by_z(nodes, z_points)
#     filtered_columns = filter_column_members_by_z(members, nodes, z_points)
#     filtered_beams_x = filter_beam_x_by_z(members, nodes, z_points)
#     filtered_beams_y = filter_beam_y_by_z(members, nodes, z_points)
        
#     # Save filtered data to separate JSON files
#     save_filtered_data(filtered_nodes, JSON_FOLDER, "filtered_nodes.json")
#     save_filtered_data(filtered_columns, JSON_FOLDER, "filtered_columns.json")
#     save_filtered_data(filtered_beams_x, JSON_FOLDER, "beamX.json")
#     save_filtered_data(filtered_beams_y, JSON_FOLDER, "beamY.json")

#     # =============================================
#     # Visualization of Structural Elements
#     # =============================================
#     plot_structure(nodes, members, IMAGE_FOLDER, title="Full Structure", save_path="full_structure.png")
#     # Generate plots for each Z level showing different element types
#     for z_level in z_points:
#         # Plot nodes at current Z level
#         plot_filtered_elements(
#             {z_level: filtered_nodes.get(z_level, {})},
#             nodes,
#             IMAGE_FOLDER,
#             title_prefix=f"Nodes at Z Level",
#             save_path=f"grid/nodes_z{z_level}.png"
#         )
        
#         # Plot columns (only for levels below the top)
#         if z_level < max(z_points):
#             plot_filtered_elements(
#                 {z_level: filtered_columns.get(z_level, {})},
#                 nodes,
#                 IMAGE_FOLDER,
#                 title_prefix=f"Columns starting at Z Level",
#                 save_path=f"grid/columns_z{z_level}.png"
#             )
        
#         # Plot beams in X direction
#         plot_filtered_elements(
#             {z_level: filtered_beams_x.get(z_level, {})},
#             nodes,
#             IMAGE_FOLDER,
#             title_prefix=f"X Beams at Z Level",
#             save_path=f"grid/beams_x_z{z_level}.png"
#         )
        
#         # Plot beams in Y direction
#         plot_filtered_elements(
#             {z_level: filtered_beams_y.get(z_level, {})},
#             nodes,
#             IMAGE_FOLDER,
#             title_prefix=f"Y Beams at Z Level",
#             save_path=f"grid/beams_y_z{z_level}.png"
#         )

#     # =============================================
#     # Section Definitions and Visualization
#     # =============================================
#     # Process all section definitions and generate visualizations
#     for category in section_definitions.values():
#         for section_name, section_info in category.items():
#             if section_info["type"] == "rectangular":
#                 # Draw rectangular RC sections
#                 draw_rc_section10(
#                     section_info["section_tag"],
#                     section_info["core_tag"],
#                     section_info["cover_tag"],
#                     section_info["steel_tag"],
#                     section_info["H"],
#                     section_info["B"],
#                     section_info["cover_H"],
#                     section_info["cover_B"],
#                     section_info["offset"],
#                     section_info["n_bars_top"],
#                     section_info["dia_top"],
#                     section_info["n_bars_bot"],
#                     section_info["dia_bot"],
#                     section_info["n_bars_secondary_top"],
#                     section_info["dia_sec_top"],
#                     section_info["n_bars_secondary_bot"],
#                     section_info["dia_sec_bot"],
#                     section_info["n_bars_int"],
#                     section_info["dia_int"],
#                     IMAGE_FOLDER
#                 )
#             elif section_info["type"] == "circular":
#                 # Draw circular RC sections
#                 draw_circular_rc_section10(
#                     section_info["section_tag"],
#                     section_info["core_tag"],
#                     section_info["cover_tag"],
#                     section_info["steel_tag"],
#                     section_info["D_Sec"],
#                     section_info["cover_Sec"],
#                     section_info["num_Bars_Sec"],
#                     section_info["bar_dia_Sec"],
#                     IMAGE_FOLDER,
#                     section_info["ri"],
#                     section_info["nf_Core_R"],
#                     section_info["nf_Core_T"],
#                     section_info["nf_Cover_R"],
#                     section_info["nf_Cover_T"]
#                 )

#     # =============================================
#     # Create Member-Section Mapping
#     # =============================================
#     mapping = create_member_section_mapping(
#         JSON_FOLDER,
#         'member_section_mapping.json'
#     )

#     # =============================================
#     # Generate OpenSees Structural Model
#     # =============================================
#     opensees_element_json_file = os.path.join(JSON_FOLDER, "element_data.json")
#     # Create the actual OpenSees model elements
#     shapes, memberLengths = create_structural_model(section_definitions, JSON_FOLDER, opensees_element_json_file, numIntgrPts=8)

#     # =============================================
#     # Shell Mesh Generation (if needed)
#     # =============================================
#     # Process each surface configuration for shell elements
#     if surface_configurations is not None and surface_configurations:
#         for config_name, config_data in surface_configurations.items():
#             print(f"\nProcessing {config_name} configuration ({config_data['section_name']})...")
            
#             # Create shell mesh for this configuration
#             mesh_data = create_proper_mesh_for_closed_area_3d1(
#                 points=config_data["points"],
#                 predefined_points=config_data["predefined_points"],
#                 shell_section_name=config_data["section_name"],
#                 JSON_FOLDER=JSON_FOLDER,
#                 IMAGE_FOLDER=IMAGE_FOLDER,
#                 num_x_div=config_data["num_x_div"],
#                 num_y_div=config_data["num_y_div"],
#                 numbering=list(surface_configurations.keys()).index(config_name) + 1,
#                 add_shell=config_data["add_shell"],
#                 remove_shell=config_data["remove_shell"]
#             )

#     # =============================================
#     # Final Data Consolidation
#     # =============================================
#     # Combine all structural data into single JSON file
#     # create_combined_structure_json(JSON_FOLDER, output_file="combined_structure_data.json")
#     shell_elements_data = create_combined_structure_json(JSON_FOLDER)
#     print(f'shell_elements_data={shell_elements_data}')
#     # =============================================
#     # Final Data Consolidation
#     # =============================================

#     # apply_shell_load(JSON_FOLDER, load_case_name="DL", numbering=1, pressure=-100.0)
#     apply_shell_load(JSON_FOLDER, load_case_names=["DL"], pressures=[-100.0], numbering=1)


#     # Process each surface configuration for shell elements
#     if surface_configurations is not None and surface_configurations:
#         for config_name, config_data in surface_configurations.items():           
#             # Apply shell loads for each load case
#             numbering = list(surface_configurations.keys()).index(config_name) + 1
#             apply_shell_load(
#                 JSON_FOLDER=JSON_FOLDER,
#                 load_case_names=config_data["load_case_names"],  # Plural
#                 pressures=config_data["pressures"],              # Plural
#                 numbering=numbering
#             )
#     apply_self_weight(JSON_FOLDER, load_case_names="self_weight")
    
#     print("\nAll configurations processed successfully!")
#     print("Model Built Successfully!")





                        
# Combine all element loads (without self-weight since it's not defined)

# Generate combined loads


def reset_model(OUTPUT_FOLDER="output_folder"):
    """
    Resets the model by deleting all generated files in the output folder without confirmation.
    
    Parameters:
    -----------
    OUTPUT_FOLDER : str, optional
        The main output directory to be cleared (default="output_folder")
        
    Returns:
    --------
    bool
        True if reset was successful, False if an error occurred
    """
    try:
        # Check if the output folder exists
        if os.path.exists(OUTPUT_FOLDER):
            print(f"Resetting model by clearing contents of {OUTPUT_FOLDER}...")
            
            # Remove all files and subdirectories in the output folder
            for filename in os.listdir(OUTPUT_FOLDER):
                file_path = os.path.join(OUTPUT_FOLDER, filename)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    print(f"Failed to delete {file_path}. Reason: {e}")
            
            print("Model reset successfully!")
            return True
        else:
            print(f"Output folder {OUTPUT_FOLDER} does not exist. Nothing to reset.")
            return True
            
    except Exception as e:
        print(f"Error resetting model: {e}")
        return False
    
OUTPUT_FOLDER = "output_folder"  # Main output directory
os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist

# build_model(x_spacing, y_spacing, z_spacing, materials, section_definitions, OUTPUT_FOLDER = "output_folder", n_integration_pts=8) 

def build_model(x_spacing, y_spacing, z_spacing, materials, section_definitions, 
                OUTPUT_FOLDER="output_folder", n_integration_pts=8, this_run="all"):
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
    this_run : str or list, optional
        Controls which parts of the code to run (default="all")
        Options: 
        - "all" - Run everything
        - "create_structure" - Create basic structural grid
        - "filter_data" - Filter structure data
        - "visualization" - Generate visualizations
        - "sections" - Process section definitions
        - "mapping" - Create member-section mapping
        - "structural_model" - Generate OpenSees structural model
        - "shell_mesh" - Generate shell meshes
        - "loads" - Apply loads
        - "final_consolidation" - Perform final data consolidation
        Can also pass a list of options to run multiple specific parts
    """
    
    # Convert this_run to list if it's a string
    if isinstance(this_run, str):
        run_parts = [this_run]
    else:
        run_parts = this_run
    
    # =============================================
    # Create Output Directory Structure (always run)
    # =============================================
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
    IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
    os.makedirs(JSON_FOLDER, exist_ok=True)
    os.makedirs(IMAGE_FOLDER, exist_ok=True)

    # =============================================
    # Create Basic Structural Grid
    # =============================================
    if "all" in run_parts or "create_structure" in run_parts:
        # reset_model(OUTPUT_FOLDER=OUTPUT_FOLDER)
        update_structure(JSON_FOLDER, x_spacing, y_spacing, z_spacing, 
                        create_new_nodes={}, create_new_members={}, 
                        delete_nodes=[], delete_members=[], tolerance=1e-6)
        print(f"Structure updated successfully.")

    # # =============================================
    # # Load and Filter Structure Data
    # # =============================================
    # if "all" in run_parts or "filter_data" in run_parts:
    #     nodes, members, structure_data = load_structure(JSON_FOLDER)
    #     z_points = get_z_levels(nodes)
    #     print(f'z_points={z_points}')
        
    #     filtered_nodes = filter_nodes_by_z(nodes, z_points)
    #     filtered_columns = filter_column_members_by_z(members, nodes, z_points)
    #     filtered_beams_x = filter_beam_x_by_z(members, nodes, z_points)
    #     filtered_beams_y = filter_beam_y_by_z(members, nodes, z_points)
            
    #     save_filtered_data(filtered_nodes, JSON_FOLDER, "filtered_nodes.json")
    #     save_filtered_data(filtered_columns, JSON_FOLDER, "filtered_columns.json")
    #     save_filtered_data(filtered_beams_x, JSON_FOLDER, "beamX.json")
    #     save_filtered_data(filtered_beams_y, JSON_FOLDER, "beamY.json")

    # # =============================================
    # # Visualization of Structural Elements
    # # =============================================
    # if "all" in run_parts or "visualization" in run_parts:
    #     if 'nodes' not in locals():
    #         nodes, members, structure_data = load_structure(JSON_FOLDER)
        
    #     plot_structure(nodes, members, IMAGE_FOLDER, title="Full Structure", save_path="full_structure.png")
        
    #     for z_level in z_points:
    #         plot_filtered_elements(
    #             {z_level: filtered_nodes.get(z_level, {})},
    #             nodes,
    #             IMAGE_FOLDER,
    #             title_prefix=f"Nodes at Z Level",
    #             save_path=f"grid/nodes_z{z_level}.png"
    #         )
            
    #         if z_level < max(z_points):
    #             plot_filtered_elements(
    #                 {z_level: filtered_columns.get(z_level, {})},
    #                 nodes,
    #                 IMAGE_FOLDER,
    #                 title_prefix=f"Columns starting at Z Level",
    #                 save_path=f"grid/columns_z{z_level}.png"
    #             )
            
    #         plot_filtered_elements(
    #             {z_level: filtered_beams_x.get(z_level, {})},
    #             nodes,
    #             IMAGE_FOLDER,
    #             title_prefix=f"X Beams at Z Level",
    #             save_path=f"grid/beams_x_z{z_level}.png"
    #         )
            
    #         plot_filtered_elements(
    #             {z_level: filtered_beams_y.get(z_level, {})},
    #             nodes,
    #             IMAGE_FOLDER,
    #             title_prefix=f"Y Beams at Z Level",
    #             save_path=f"grid/beams_y_z{z_level}.png"
    #         )

    # =============================================
    # Section Definitions and Visualization
    # =============================================
    if "all" in run_parts or "sections" in run_parts:
        for category in section_definitions.values():
            for section_name, section_info in category.items():
                if section_info["type"] == "rectangular":
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
                elif section_info["type"] == "L":
                    draw_L_rc_section10({
                        'sec_tag': section_info["section_tag"],
                        'core_tag': section_info["core_tag"],
                        'cover_tag': section_info["cover_tag"],
                        'steel_tag': section_info["steel_tag"],
                        'H': section_info["H"],
                        'B': section_info["B"],
                        't': section_info["t"],
                        'cover_H': section_info["cover_H"],
                        'cover_B': section_info["cover_B"],
                        'n_bars_vertical_outer': section_info["n_bars_vertical_outer"],
                        'dia_vertical_outer': section_info["dia_vertical_outer"],
                        'n_bars_horizontal_outer': section_info["n_bars_horizontal_outer"],
                        'dia_horizontal_outer': section_info["dia_horizontal_outer"],
                        'n_bars_vertical_inner': section_info["n_bars_vertical_inner"],
                        'dia_vertical_inner': section_info["dia_vertical_inner"],
                        'n_bars_horizontal_inner': section_info["n_bars_horizontal_inner"],
                        'dia_horizontal_inner': section_info["dia_horizontal_inner"],
                        'corner_bar_dia': section_info["corner_bar_dia"],
                        'IMAGE_FOLDER': IMAGE_FOLDER
                    })


    # =============================================
    # Create Member-Section Mapping
    # =============================================
    if "all" in run_parts or "member_section_mapping" in run_parts:
        mapping = create_member_section_mapping(
            JSON_FOLDER,
            'member_section_mapping.json'
        )

    # =============================================
    # Generate OpenSees Structural Model
    # =============================================
    if "all" in run_parts or "create_beam_column_member" in run_parts:
        opensees_element_json_file = os.path.join(JSON_FOLDER, "element_data.json")
        shapes, memberLengths = create_beam_column_member(section_definitions, JSON_FOLDER, 
                                                       opensees_element_json_file, numIntgrPts=8)

    # =============================================
    # Shell Mesh Generation (if needed)
    # =============================================
    if ("all" in run_parts or "shell_mesh" in run_parts) and surface_configurations is not None and surface_configurations:
        for config_name, config_data in surface_configurations.items():
            print(f"\nProcessing {config_name} configuration ({config_data['section_name']})...")
            
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
    if "all" in run_parts or "final_consolidation" in run_parts:
        shell_elements_data = create_combined_structure_json(JSON_FOLDER)
        # print(f'shell_elements_data={shell_elements_data}')

    # =============================================
    # Apply Loads
    # =============================================
    if "all" in run_parts or "apply_loads" in run_parts:
        
        
    #     if surface_configurations is not None and surface_configurations:
    #         for config_name, config_data in surface_configurations.items():           
    #             numbering = list(surface_configurations.keys()).index(config_name) + 1
    #             apply_shell_load(
    #                 JSON_FOLDER=JSON_FOLDER,
    #                 load_case_names=config_data["load_case_names"],
    #                 pressures=config_data["pressures"],
    #                 numbering=numbering
    #             )
    #     apply_self_weight(JSON_FOLDER, load_case_names="self_weight")
    # save_nodal_load_cases(JSON_FOLDER, nodal_load_entries)
    # create_load_combinations(JSON_FOLDER, load_combinations)
    # combined_loads = create_combined_loads(load_combinations, *all_element_loads, JSON_FOLDER=JSON_FOLDER)
    # convert_member_loads_to_nodal(JSON_FOLDER)
    # nodal_mass_combo(JSON_FOLDER, load_combinations)
        process_loads(JSON_FOLDER, surface_configurations, nodal_load_entries, 
                    load_combinations, all_element_loads, run_parts)
        print("\nOperation completed successfully!")
        print(f"Ran the following parts: {run_parts}")

# Example usage:
# Run only the structure creation and filtering
# build_model(x_spacing, y_spacing, z_spacing, materials, section_definitions, 
#             OUTPUT_FOLDER="output_folder", this_run="create_structure")

# Run only the visualization
# build_model(x_spacing, y_spacing, z_spacing, materials, section_definitions, 
#             OUTPUT_FOLDER="output_folder", this_run="visualization")

# Run everything (default)
# build_model(x_spacing, y_spacing, z_spacing, materials, section_definitions, 
#             OUTPUT_FOLDER="output_folder")


# # Run everything (default)
# build_model(x_spacing, y_spacing, z_spacing, materials, section_definitions)

# # Run just structure creation and filtering
# build_model(x_spacing, y_spacing, z_spacing, materials, section_definitions,
#             this_run=["create_structure", "filter_data"])

# # Run just visualization and sections
# build_model(x_spacing, y_spacing, z_spacing, materials, section_definitions,
#             this_run=["visualization", "sections"])

# Run everything except loads
build_model(x_spacing, y_spacing, z_spacing, materials, section_definitions,
            this_run=["create_structure", "filter_data", "visualization", 
                      "sections", "member_section_mapping", "create_beam_column_member", "shell_mesh",
                      "final_consolidation", "apply_loads"])



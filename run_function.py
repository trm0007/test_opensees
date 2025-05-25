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



# Extruded shapes visualization - Fixed version

# Updated visualization function with better error handling

# =============================================
# Section Visualization Setup
# =============================================

def extract_materials(materials_config):
    all_materials = []
    for mat in materials_config["materials"]:
        mat_data = {
            "name": mat["name"],
            "type": mat["type"],
            "id": mat["id"],
            **mat["properties"]
        }
        all_materials.append(mat_data)
    return all_materials

def extract_sections(sections_config):
    all_sections = []
    for sec in sections_config["sections"]:
        sec_data = {
            "name": sec["name"],
            "type": sec["type"],
            "id": sec["id"],
            "material_id": sec["material_id"],
            **sec["properties"]
        }
        all_sections.append(sec_data)
    return all_sections




def assign_node_into_opensees(JSON_FOLDER):
    """
    Creates nodes in OpenSees directly from the combined JSON file.

    Args:
        JSON_FOLDER (str): Path to folder containing the combined JSON file
        json_filename (str): Name of the JSON file (default matches create_combined_structure_json output)

    Returns:
        dict: Mapping from node names to OpenSees node IDs
    """
    # Load the combined JSON file
    # Call the first function to load structure data
    nodes_dict, members_dict, structure = load_structure(JSON_FOLDER)

    # Call the second function to create combined structure JSON
    combined_data = create_combined_structure_json(JSON_FOLDER)
    # 2. Load the combined structure data to get all node coordinates
    # Merge structural data
    data = merge_structures(structure, combined_data)
    # print(f'Processed data structure: {data}')

    # Process nodes and create OpenSees nodes
    nodes = data['nodes']
    all_nodes = {}  # Dictionary to store node information

    for node in nodes:
        try:
            # Validate required fields exist
            node_id = node['id']
            node_name = node['name']
            
            # Require explicit coordinate values - no defaults
            if 'x' not in node or 'y' not in node or 'z' not in node:
                raise ValueError(f"Node {node_id} ({node_name}) is missing coordinate data")
                
            node_coordinates = {
                'x': node['x'],
                'y': node['y'], 
                'z': node['z']
            }
            
            # Store node information
            all_nodes[node_name] = {
                'id': node_id,
                'coordinates': node_coordinates
            }
            
            # Create node in OpenSees
            ops.node(node_id, node_coordinates['x'], node_coordinates['y'], node_coordinates['z'])
            
            print(f"Created node - ID: {node_id}, Name: {node_name}, "
                f"Coordinates: ({node_coordinates['x']}, {node_coordinates['y']}, {node_coordinates['z']})")
        
        except KeyError as e:
            print(f"Error: Node {node.get('id', 'UNKNOWN')} is missing required field {str(e)}")
            continue
        except ValueError as e:
            print(f"Error: {str(e)}")
            continue


def define_materials(materials_list):
    for mat in materials_list:
        ops.uniaxialMaterial(*mat.values())


def assign_member_section_into_opensees(section_definitions):
    # =============================================
    # Section Definitions (using build_RC_rect_section)
    # =============================================
    # Create section_tags dynamically
    section_tags = {}
    tag_counter = 1

    # Iterate through each category in section_definitions
    for category in section_definitions.values():
        # Iterate through each section in the category
        for section_name in category:
            section_tags[section_name] = tag_counter
            tag_counter += 1
    
    for section_name, tag in section_tags.items():
        if section_name in section_definitions["rectangular_sections"]:
            section_info = section_definitions["rectangular_sections"][section_name]
            define_rc_section10_opensees(
                tag,
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
                section_info["dia_int"]
            )
        elif section_name in section_definitions["circular_sections"]:
            section_info = section_definitions["circular_sections"][section_name]
            define_circular_rc_section10_opensees(
                tag,
                section_info["core_tag"],
                section_info["cover_tag"],
                section_info["steel_tag"],
                section_info["D_Sec"],
                section_info["cover_Sec"],
                section_info["num_Bars_Sec"],
                section_info["bar_dia_Sec"],
                section_info["ri"],
                section_info["nf_Core_R"],
                section_info["nf_Core_T"],
                section_info["nf_Cover_R"],
                section_info["nf_Cover_T"]
            )


def assign_beam_column_into_opensees(opensees_element_json_file):
    with open(opensees_element_json_file, "r") as f:
        d = json.load(f)

    for ele in d["elements"]:
        ops.geomTransf(ele["transType"], ele["transTag"],
                       float(ele["vecxz"][0]), float(ele["vecxz"][1]), float(ele["vecxz"][2]))

        ops.element('nonlinearBeamColumn', ele["eleTag"], ele["node_i_id"], ele["node_j_id"],
                    ele["numIntgrPts"], ele["secTag"], ele["transTag"], '-integration', 'Lobatto')



def assign_shell_material_and_section_into_opensees(materials_config, sections_config):
    materials = extract_materials(materials_config)
    sections = extract_sections(sections_config)

    # print("Extracted Materials:")
    # for material in materials:
    #     print("  Material:")
    #     for key, value in material.items():
    #         print(f"    {key}: {value}")
    #     print()

    # print("Extracted Sections:")
    # for section in sections:
    #     print("  Section:")
    #     for key, value in section.items():
    #         print(f"    {key}: {value}")
    #     print()

    # Create materials
    for mat in materials_config['materials']:
        mat_id = mat['id']
        mat_type = mat['type']
        props = mat['properties']

        if mat_type == "ENT":
            ops.uniaxialMaterial('ENT', mat_id, props['E'])
        elif mat_type == "ElasticIsotropic":
            ops.nDMaterial('ElasticIsotropic', mat_id, props['E'], props['nu'], props['rho'])

    # Create sections
    for sec in sections_config['sections']:
        sec_id = sec['id']
        mat_id = sec['material_id']
        props = sec['properties']

        if sec['type'] == "PlateFiber":
            ops.section('PlateFiber', sec_id, mat_id, props['thickness'])

    return True


def assign_shell_element_into_opensees(JSON_FOLDER, section_config=None):
    """
    Assigns shell elements into OpenSees using section definitions from section_config.
    
    Args:
        JSON_FOLDER (str): Path to folder containing the combined structure JSON
        output_file (str): Name of the combined JSON file
        section_config (dict): Section configuration dictionary containing section IDs
    """
    # Call the first function to load structure data
    nodes_dict, members_dict, structure = load_structure(JSON_FOLDER)

    # Call the second function to create combined structure JSON
    combined_data = create_combined_structure_json(JSON_FOLDER)
    # 2. Load the combined structure data to get all node coordinates
    combined_data = merge_structures(structure, combined_data)
    
    # Create node name to OpenSees tag mapping
    node_name_to_tag = {node['name']: node['id'] for node in combined_data['nodes']}
    
    # Create section name to ID mapping from section_config
    section_name_to_id = {}
    if section_config:
        for section in section_config.get('sections', []):
            section_name_to_id[section['name']] = section['id']
    
    # Initialize counters
    created_elements = 0
    skipped_elements = 0
    missing_sections = set()
    
    # Assign shell elements
    for shell_elem in combined_data['shell_elements']:
        elem_id = shell_elem['id']
        elem_name = shell_elem.get('name', f"element_{elem_id}")
        node_names = shell_elem['node_names']
        section_name = shell_elem.get('shell_section', 'default')
        
        # Get section ID from config
        section_id = section_name_to_id.get(section_name)
        if section_id is None:
            print(f"Warning: Section '{section_name}' not found in section_config for element {elem_name}")
            missing_sections.add(section_name)
            skipped_elements += 1
            continue
        
        # Get node tags and check for missing nodes
        node_tags = []
        missing = False
        
        for node_name in node_names:
            if node_name not in node_name_to_tag:
                print(f"Warning: Node {node_name} not found for element {elem_name}")
                missing = True
                break
            node_tags.append(node_name_to_tag[node_name])
        
        if missing:
            skipped_elements += 1
            continue
        
        # Create element if all nodes were resolved
        try:
            if len(node_tags) == 4:
                ops.element("ShellMITC4", elem_id, *node_tags, section_id)
                # print(f"Created ShellMITC4 element {elem_name} (ID:{elem_id}) with section '{section_name}' (ID:{section_id})")
                created_elements += 1
            elif len(node_tags) == 3:
                ops.element("ShellDKGT", elem_id, *node_tags, section_id)
                # print(f"Created ShellDKGT element {elem_name} (ID:{elem_id}) with section '{section_name}' (ID:{section_id})")
                created_elements += 1
            else:
                print(f"Warning: Element {elem_name} has unsupported number of nodes ({len(node_tags)})")
                skipped_elements += 1
        except Exception as e:
            print(f"Error creating element {elem_name}: {str(e)}")
            skipped_elements += 1
    
    # Print summary
    print(f"\nElement creation summary:")
    print(f"  Successfully created: {created_elements}")
    print(f"  Skipped: {skipped_elements}")
    print(f"  Total processed: {len(combined_data['shell_elements'])}")
    
    if missing_sections:
        print("\nWarning: Missing section definitions for:")
        for section in missing_sections:
            print(f"  - {section}")
        print("Please ensure these sections are defined in your section_config")
    
    return created_elements > 0


# Analysis Functions - Improved with better recorder handling
def ensure_output_dir():
    """Ensure output directory exists."""
    os.makedirs('FGU_RC3DF_files', exist_ok=True)
# =============================================
# 6. run_gravity
# =============================================
def run_gravity(JSON_FOLDER, steps=10, comb_name=None):
    """Run gravity analysis with dynamic recorder naming."""
    ensure_output_dir()
    reaction_file = os.path.join(JSON_FOLDER, "post_processing", f"Gravity_Reactions_{comb_name}.out" if comb_name else "Gravity_Reactions.out")
    # reaction_file = f"Gravity_Reactions_{comb_name}.out" if comb_name else "Gravity_Reactions.out"
    reaction_path = os.path.join('FGU_RC3DF_files', reaction_file)
    
    # ops.recorder('Node', '-file', reaction_path,
    #             '-time', '-node', '-dof', 'reaction')

    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-6, 100, 0, 2)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1/steps)
    ops.analysis('Static')
    # ops.analyze(steps)  # e.g., steps=10
    
    ops.record()
    
    ok = ops.analyze(steps)
    ops.reactions()     # Must call to update reactions!
    
    if ok == 0:
        print(f"Gravity analysis {'for ' + comb_name if comb_name else ''} completed successfully")
        return True
    else:
        print(f"Gravity analysis {'for ' + comb_name if comb_name else ''} failed")
        return False


# =============================================
# 7. run_modal
# =============================================
def run_modal(n_evs=3, comb_name=None):
    """
    Runs Modal analysis with support for load combinations.
    
    Args:
        n_evs (int): Number of eigenvalues to compute
        comb_name (str): Name of the load combination (for output files)
        
    Returns:
        np.array: Array of eigenvalues
    """
    ensure_output_dir()
    
    # Create unique recorder names if combination is specified
    if comb_name:
        eigen_files = [
            os.path.join('FGU_RC3DF_files', f'ModalAnalysis_EigenVec1_{comb_name}.out'),
            os.path.join('FGU_RC3DF_files', f'ModalAnalysis_EigenVec2_{comb_name}.out'),
            os.path.join('FGU_RC3DF_files', f'ModalAnalysis_EigenVec3_{comb_name}.out')
        ]
        eigenval_file = os.path.join('FGU_RC3DF_files', f'ModalAnalysis_EigenVal_{comb_name}.out')
    else:
        eigen_files = [
            os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVec1.out'),
            os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVec2.out'),
            os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVec3.out')
        ]
        eigenval_file = os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVal.out')

    # Set up recorders for each mode
    for i, file in enumerate(eigen_files[:n_evs], 1):
        ops.recorder('Node', '-file', file,
                    '-node', *list(range(5,9)), '-dof', 1, 2, f'eigen {i}')

    # Analysis configuration
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('BandGen')
    ops.test('NormDispIncr', 1.0e-12, 25, 0, 2)
    ops.algorithm('Newton')
    ops.analysis('Transient')  # Needed for eigen commands
    
    # Compute eigenvalues
    lamda = np.array(ops.eigen(n_evs))
    
    # Write eigenvalues to file
    with open(eigenval_file, "w") as eig_file:
        eig_file.write("lambda omega period frequency\n")
        for l in lamda:
            omega = l**0.5
            period = 2*np.pi/omega
            freq = omega/(2*np.pi)
            eig_file.write(f"{l:2.6e} {omega:2.6e} {period:2.6e} {freq:2.6e}\n")

    ops.record()
    print(f"Modal analysis {'for ' + comb_name if comb_name else ''} completed")
    return lamda
# =============================================
# 8. run_pushover
# =============================================
def run_pushover(steps=10000, direction='X', comb_name=None):
    """Run pushover analysis with dynamic recorder naming."""
    ensure_output_dir()
    
    reaction_file = f"Pushover_Horizontal_Reactions{direction}_{comb_name}.out" if comb_name else f"Pushover_Horizontal_Reactions{direction}.out"
    disp_file = f"Pushover_Story_Displacement{direction}_{comb_name}.out" if comb_name else f"Pushover_Story_Displacement{direction}.out"
    
    reaction_path = os.path.join('FGU_RC3DF_files', reaction_file)
    disp_path = os.path.join('FGU_RC3DF_files', disp_file)

    d_o_f = 1 if direction == 'X' else 2
    phi = 1.0
    
    ops.recorder('Node', '-file', reaction_path,
                '-time', '-node', *list(range(1,5)), '-dof', d_o_f, 'reaction')
    ops.recorder('Node', '-file', disp_path,
                '-time', '-node', *list(range(5,9)), '-dof', d_o_f, 'disp')

    # Create lateral load pattern
    pattern_tag = 100 if comb_name else 2  # Use high tag for combination cases
    ops.pattern('Plain', pattern_tag, 1)
    step = 1.0e-05
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGen')
    ops.test('NormDispIncr', 0.000001, 100)
    ops.algorithm('NewtonLineSearch', True, 0.8, 1000, 0.1, 10.0)
    ops.integrator('DisplacementControl', 5, d_o_f, step)
    ops.analysis('Static')
    ops.record()

    ok = ops.analyze(steps)
    
    if ok == 0:
        print(f'Pushover Analysis in {direction} {"for " + comb_name if comb_name else ""} completed successfully')
        return True
    else:
        print(f'Pushover Analysis in {direction} {"for " + comb_name if comb_name else ""} failed')
        return False
# =============================================
# 9. run_time_history
# =============================================
def run_time_history(direction='X', g_motion_id=1, scaling_id=1,
                    lamda=1.0, acc_file='FGU_RC3DF_files/acc_1.txt',
                    comb_name=None, analysis_duration_ratio=0.29):
    """
    Runs Time history analysis with support for load combinations.
    
    Args:
        direction (str): Direction of excitation ('X' or 'Y')
        g_motion_id (int): Ground motion identifier
        scaling_id (int): Scaling factor identifier
        lamda (float): Scaling factor for ground motion
        acc_file (str): Path to acceleration file
        comb_name (str): Name of the load combination
        analysis_duration_ratio (float): Fraction of full duration to analyze (0-1)
        
    Returns:
        bool: True if analysis succeeded, False otherwise
    """
    ensure_output_dir()
    
    # Create output file names
    if comb_name:
        reaction_file = os.path.join('FGU_RC3DF_files', 
                                   f'TimeHistory_Reactions_{direction}_{g_motion_id}_{scaling_id}_{comb_name}.out')
        disp_file = os.path.join('FGU_RC3DF_files', 
                               f'TimeHistory_Displacement_{direction}_{g_motion_id}_{scaling_id}_{comb_name}.out')
        accel_file = os.path.join('FGU_RC3DF_files', 
                                 f'TimeHistory_Acceleration_{direction}_{g_motion_id}_{scaling_id}_{comb_name}.out')
    else:
        reaction_file = os.path.join('FGU_RC3DF_files', 
                                    f'TimeHistory_Reactions_{direction}_{g_motion_id}.out')
        disp_file = os.path.join('FGU_RC3DF_files', 
                               f'TimeHistory_Displacement_{direction}_{g_motion_id}.out')
        accel_file = os.path.join('FGU_RC3DF_files', 
                                f'TimeHistory_Acceleration_{direction}_{g_motion_id}.out')

    # Determine DOF and damping parameters
    dof = 1 if direction == 'X' else 2
    
    # Get modal properties (assuming first mode dominates)
    try:
        eigenvals = np.loadtxt(os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVal.out'), 
                             skiprows=1)
        omega = eigenvals[0,1] if direction == 'X' else eigenvals[2,1]
    except:
        print("Warning: Could not read eigenvalue file, using default omega=2π")
        omega = 2*np.pi  # Fallback value
    
    xi = 0.05  # Damping ratio (5%)
    alpha_M = 0.0       # Mass proportional damping
    beta_K = 2*xi/omega # Stiffness proportional damping
    
    # Set up recorders
    ops.recorder('Node', '-file', reaction_file,
                '-time', '-node', *list(range(1,5)), '-dof', dof, 'reaction')
    ops.recorder('Node', '-file', disp_file,
                '-time', '-node', *list(range(5,9)), '-dof', dof, 'disp')
    ops.recorder('Node', '-file', accel_file,
                '-time', '-node', *list(range(5,9)), '-dof', dof, 'accel')

    # Load acceleration time history
    accelerogram = np.loadtxt(acc_file)
    dt = 0.02  # Time step in acceleration file
    n_steps = len(accelerogram)
    
    # Analysis parameters
    tol = 1.0e-6
    max_iter = 500
    analysis_dt = 0.01  # Analysis time step (should be ≤ dt/2 for accuracy)

    # Define time series and pattern
    ops.timeSeries('Path', 2, '-dt', dt, '-values', *accelerogram, '-factor', lamda)
    ops.pattern('UniformExcitation', 3, dof, '-accel', 2)
    
    # Analysis configuration
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', tol, max_iter, 0, 2)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.rayleigh(alpha_M, beta_K, 0.0, 0.0)
    ops.analysis('Transient')

    # Run analysis
    print(f"Running Time-History analysis (λ={lamda}) {'for ' + comb_name if comb_name else ''}")
    start_time = time.time()
    
    ok = 0
    current_time = ops.getTime()
    final_time = n_steps * dt
    target_time = analysis_duration_ratio * final_time
    
    while ok == 0 and current_time < target_time:
        ok = ops.analyze(1, analysis_dt)
        current_time = ops.getTime()
        
        # Optional: Print progress
        if int(current_time/dt) % 100 == 0:
            print(f"Time: {current_time:.2f}s ({current_time/target_time:.1%})")
    
    elapsed_time = time.time() - start_time
    
    if ok == 0:
        print(f"Time-History completed in {elapsed_time:.2f}s")
        return True
    else:
        print(f"Time-History failed at {current_time:.2f}s")
        return False
# =============================================
# 10. reset_analysis
# =============================================
def reset_analysis():
    """Reset the analysis state."""
    ops.setTime(0.0)
    ops.loadConst()
    ops.remove('recorders')
    ops.wipeAnalysis()

     



# OUTPUT_FOLDER = "output_folder"  # Main output directory
# os.makedirs(OUTPUT_FOLDER, exist_ok=True)  # Create if doesn't exist

# Subdirectories
# JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "json_files")
# IMAGE_FOLDER = os.path.join(OUTPUT_FOLDER, "images")
# os.makedirs(JSON_FOLDER, exist_ok=True)
# os.makedirs(IMAGE_FOLDER, exist_ok=True)
# opensees_element_json_file = os.path.join(JSON_FOLDER, "element_data.json")
# final_run(JSON_FOLDER, materials, section_definitions, materials_config, sections_config, opensees_element_json_file)






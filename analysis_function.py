from import_ import *
from units import *

# Analysis Functions - Improved with better recorder handling
def ensure_output_dir():
    """Ensure output directory exists."""
    os.makedirs('FGU_RC3DF_files', exist_ok=True)

def run_gravity(steps=10, comb_name=None):
    """Run gravity analysis with dynamic recorder naming."""
    ensure_output_dir()
    
    reaction_file = f"Gravity_Reactions_{comb_name}.out" if comb_name else "Gravity_Reactions.out"
    reaction_path = os.path.join('FGU_RC3DF_files', reaction_file)
    
    ops.recorder('Node', '-file', reaction_path,
                '-time', '-node', *list(range(1,5)), '-dof', *list(range(1,7)), 'reaction')

    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-6, 100, 0, 2)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1/steps)
    ops.analysis('Static')
    ops.record()
    
    ok = ops.analyze(steps)
    
    if ok == 0:
        print(f"Gravity analysis {'for ' + comb_name if comb_name else ''} completed successfully")
        return True
    else:
        print(f"Gravity analysis {'for ' + comb_name if comb_name else ''} failed")
        return False

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

def run_pushover(m_1, steps=10000, direction='X', comb_name=None):
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
    for node in range(5, 9):
        if direction == 'X':
            ops.load(node, m_1*phi, 0, 0, 0, 0, 0)
        else:
            ops.load(node, 0, m_1*phi, 0, 0, 0, 0)

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
    
def reset_analysis():
    """Reset the analysis state."""
    ops.setTime(0.0)
    ops.loadConst()
    ops.remove('recorders')
    ops.wipeAnalysis()

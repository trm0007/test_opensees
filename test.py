# from import_ import *
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# # Clear any existing model
# ops.wipe()

# # Create model
# ops.model('basic', '-ndm', 3, '-ndf', 6)

# # Define material properties
# E = 200000.0      # Young's modulus (MPa)
# nu = 0.3          # Poisson's ratio
# h = 0.2           # Shell thickness (m)
# rho = 2500.0      # Density (kg/m³)

# # Create shell section
# ops.section('ElasticMembranePlateSection', 1, E, nu, h, rho)

# # Define nodes for a simple rectangular shell element
# ops.node(1, 0.0, 0.0, 0.0)
# ops.node(2, 2.0, 0.0, 0.0)
# ops.node(3, 2.0, 2.0, 0.0)
# ops.node(4, 0.0, 2.0, 0.0)

# # Create shell element (ShellMITC4)
# ops.element('ShellMITC4', 1, 1, 2, 3, 4, 1)

# # Apply boundary conditions (fix one edge)
# ops.fix(1, 1, 1, 1, 1, 1, 1)
# ops.fix(4, 1, 1, 1, 1, 1, 1)

# # Optional: Add some loading for visualization
# ops.timeSeries('Linear', 1)
# ops.pattern('Plain', 1, 1)
# ops.load(2, 0, 0, -1000, 0, 0, 0)  # Downward load at node 2
# ops.load(3, 0, 0, -1000, 0, 0, 0)  # Downward load at node 3

# # Create analysis
# ops.constraints('Plain')
# ops.numberer('RCM')
# ops.system('BandGeneral')
# ops.algorithm('Linear')
# ops.integrator('LoadControl', 1.0)
# ops.analysis('Static')

# # Run analysis
# ops.analyze(1)

# # ====================== VISUALIZATION ====================== #
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

# # Set axis labels and title
# ax.set_xlabel('X (m)')
# ax.set_ylabel('Y (m)')
# ax.set_zlabel('Z (m)')
# ax.set_title('3D Shell Element with Fill Color')

# plt.tight_layout()
# plt.show()




# import openseespy.opensees as ops
# import opsvis as ops_vis
# import numpy as np
# import matplotlib.pyplot as plt




# # Clear any existing model
# ops.wipe()

# # Create model builder
# ops.model('basic', '-ndm', 3, '-ndf', 6)

# # Define materials
# ops.uniaxialMaterial('Elastic', 1, 200000.0)  # Steel E = 200 GPa

# # Define nodes (simple 3-story frame)
# # Node tag, x, y, z coordinates
# ops.node(1, 0.0, 0.0, 0.0)    # Base
# ops.node(2, 0.0, 3.0, 0.0)    # 1st floor
# ops.node(3, 0.0, 6.0, 0.0)    # 2nd floor  
# ops.node(4, 0.0, 9.0, 0.0)    # 3rd floor

# # Define boundary conditions (fix base)
# ops.fix(1, 1, 1, 1, 1, 1, 1)

# # Define masses (concentrated at floors)
# mass = 1000.0  # kg
# ops.mass(2, mass, mass, 0.0, 0.0, 0.0, 0.0)
# ops.mass(3, mass, mass, 0.0, 0.0, 0.0, 0.0)
# ops.mass(4, mass, mass, 0.0, 0.0, 0.0, 0.0)

# # Define geometric transformation
# ops.geomTransf('Linear', 1, 0.0, 0.0, 1.0)

# # Define elements (columns)
# # Element tag, node i, node j, area, E, G, J, Iy, Iz, transformation
# A = 0.01      # Cross-sectional area (m²)
# E = 200e9     # Elastic modulus (Pa)
# G = 80e9      # Shear modulus (Pa)
# J = 8.33e-6   # Torsional constant (m⁴)
# Iy = 8.33e-6  # Moment of inertia about y-axis (m⁴)
# Iz = 8.33e-6  # Moment of inertia about z-axis (m⁴)

# ops.element('elasticBeamColumn', 1, 1, 2, A, E, G, J, Iy, Iz, 1)
# ops.element('elasticBeamColumn', 2, 2, 3, A, E, G, J, Iy, Iz, 1)
# ops.element('elasticBeamColumn', 3, 3, 4, A, E, G, J, Iy, Iz, 1)

# # MODAL ANALYSIS
# print("=== MODAL ANALYSIS ===")

# # Set up the system for eigenvalue analysis
# ops.system('BandGeneral')
# ops.numberer('RCM')
# ops.constraints('Plain')
# ops.integrator('LoadControl', 1.0)
# ops.algorithm('Linear')
# ops.analysis('Static')

# # Perform eigenvalue analysis
# Nmodos = 3
# eigenValues = ops.eigen(Nmodos)

# # Calculate frequencies and periods
# Nmodos = 3
# w_i = np.sqrt(ops.eigen(Nmodos))
# T_i = 2 * np.pi / w_i

# print("Modal Analysis Results:")
# for i in range(Nmodos):
#     print(f"Mode {i+1}:")
#     print(f"  Frequency: {w_i[i]/(2*np.pi):.4f} Hz")
#     print(f"  Period: {T_i[i]:.4f} sec")
#     print()

# # Get mode shapes (simplified - just for reference)
# modeShapes = []
# for mode in range(1, Nmodos + 1):
#     nodeDisps = []
#     for node in [2, 3, 4]:  # Only free nodes
#         disp = ops.nodeEigenvector(node, mode, 1)  # X-direction displacement
#         nodeDisps.append(disp)
#     modeShapes.append(nodeDisps)

# # RESPONSE SPECTRUM ANALYSIS
# print("\n=== RESPONSE SPECTRUM ANALYSIS ===")

# # Define response spectrum (simplified)
# # Periods and corresponding spectral accelerations
# spectrum_periods = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0])
# spectrum_sa = np.array([0.4, 0.4, 0.6, 1.0, 0.8, 0.4, 0.2])  # in g

# # Interpolate spectrum for calculated periods
# spectral_accelerations = np.interp(T_i, spectrum_periods, spectrum_sa)

# print("Response Spectrum Results:")
# print("Mode | Period (s) | Frequency (Hz) | Sa (g) | Modal Response")
# print("-" * 60)

# modal_responses = []
# for i in range(Nmodos):
#     sa_g = spectral_accelerations[i]
#     # Modal response = Sa * participation factor (simplified as 1.0)
#     participation_factor = 1.0  # Simplified assumption
#     modal_response = sa_g * 9.81 * participation_factor  # Convert to m/s²
#     modal_responses.append(modal_response)
    
#     print(f"{i+1:4d} | {T_i[i]:8.4f} | {w_i[i]/(2*np.pi):10.4f} | {sa_g:6.2f} | {modal_response:8.2f} m/s²")

# # Calculate total response using SRSS (Square Root of Sum of Squares)
# total_response = np.sqrt(sum([resp**2 for resp in modal_responses]))
# print(f"\nTotal Response (SRSS): {total_response:.2f} m/s²")

# # Gráficas de modos
# for i in range(1, Nmodos + 1):
#     ops_vis.plot_mode_shape(i, endDispFlag=0, fig_wi_he=(18, 18), node_supports=False)
#     plt.show()

# # Plot Response Spectrum
# plt.figure(figsize=(10, 6))
# plt.plot(spectrum_periods, spectrum_sa, 'b-', linewidth=2, label='Design Spectrum')
# plt.scatter(T_i, spectral_accelerations, color='red', s=100, 
#            label='Structure Periods', zorder=5)
# plt.xlabel('Period (s)')
# plt.ylabel('Spectral Acceleration (g)')
# plt.title('Response Spectrum')
# plt.legend()
# plt.grid(True, alpha=0.3)
# plt.show()

# # Clean up
# ops.wipe()

# print("\nAnalysis completed successfully!")
# print(f"The structure has {Nmodos} modes analyzed.")
# print(f"Fundamental period: {T_i[0]:.4f} seconds")
# print(f"Total response (SRSS): {total_response:.2f} m/s²")



# # import numpy as np
# # import matplotlib.pyplot as plt

# def calculate_response_spectrum(S, TB, TC, TD, xi, Z, I, R, file_name='response_spec.txt',
#                                 plot_name='response_spectrum.png'):
#     # Define the function to calculate Cs based on given parameters
#     def calculate_Cs(S, T, TB, TC, TD, xi):
#         # Calculate the damping correction factor mu
#         mu = (10 / (5 + xi)) ** 0.5
#         # Ensure mu is not smaller than 0.55
#         mu = max(mu, 0.55)

#         # Initialize Cs
#         Cs = 0

#         # Calculate Cs based on the given conditions
#         if 0 <= T <= TB:
#             Cs = S * (1 + (T / TB) * (2.5 * mu - 1))
#         elif TB < T <= TC:  # Fixed: changed from TB <= T to TB < T
#             Cs = 2.5 * S * mu
#         elif TC < T <= TD:  # Fixed: changed from TC <= T to TC < T
#             Cs = 2.5 * S * mu * (TC / T)
#         elif TD < T <= 4:   # Fixed: changed from TD <= T to TD < T
#             Cs = 2.5 * S * mu * (TC * TD / T ** 2)

#         return Cs

#     # T values from 0 to 4 seconds with interval 0.005
#     T_values = np.arange(0, 4.005, 0.005)  # Fixed: changed to 4.005 to include T=4

#     # Calculate Cs for each T value
#     Cs_values = []
#     for T in T_values:
#         Cs = calculate_Cs(S, T, TB, TC, TD, xi)
#         Cs_values.append(Cs)

#     # Calculate Sa values
#     Cs_values = np.array(Cs_values)
#     Sa_values = (2 / 3) * (Z * I / R) * Cs_values

#     # Save results to file
#     with open(file_name, 'w') as f:
#         f.write("T (sec)\tSa\n")
#         for i in range(len(T_values)):
#             f.write(f"{T_values[i]:.3f}\t{Sa_values[i]:.6f}\n")  # Fixed: changed to 3 decimal places for T

#     # Additional data points to plot in red
#     red_points_T = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050]
#     red_points_Sa = [0.060000, 0.060900, 0.061800, 0.062700, 0.063600, 0.064500, 0.065400, 0.066300, 0.067200, 0.068100, 0.069000]
    
#     # Plot T vs Sa (Fixed: plotting Sa instead of Cs)
#     plt.figure(figsize=(10, 6))
#     plt.plot(T_values, Sa_values, 'b-', label='Sa (Spectral Acceleration)', linewidth=2)
#     plt.plot(red_points_T, red_points_Sa, 'ro', label='Specific Points', markersize=6, markerfacecolor='red', markeredgecolor='darkred')
#     plt.xlabel('Period T (sec)')
#     plt.ylabel('Spectral Acceleration Sa')
#     plt.title('Response Spectrum')
#     plt.grid(True, alpha=0.3)
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(plot_name, dpi=300, bbox_inches='tight')
#     plt.show()

#     return T_values, Sa_values, Cs_values  # Added return values for further analysis

# # Parameters
# S = 1.5   # Soil factor
# TB = 0.5  # Lower limit of the period
# TC = 1.5  # Upper limit of the period
# TD = 2.0  # Lower limit of the period for displacement
# xi = 5    # Viscous damping ratio (%)
# Z = 0.2   # Seismic zone factor
# I = 1.5   # Importance factor
# R = 5.0   # Response reduction factor

# # Call the function
# T_vals, Sa_vals, Cs_vals = calculate_response_spectrum(S, TB, TC, TD, xi, Z, I, R)







def generate_seismic_load_combinations(Z, S):
    Ev = 0.5 * (2/3) * Z * S
    DL_12 = 1.2 + Ev
    DL_09 = 0.9 + Ev

    combos = {}

    eq_signs = [1.0, -1.0]
    eqx_factors = [1.0, 0.3]
    eqy_factors = [1.0, 0.3]

    idx = 1
    # First group: 1.2DL + EQx ± 0.3EQy + LL
    for sx in eq_signs:
        for sy in [0.3, -0.3]:
            combos[f"LC{idx}"] = [("DL", DL_12), ("EQx", sx), ("EQy", sy), ("LL", 1.0)]
            idx += 1

    # Second group: 1.2DL + 0.3EQx ± EQy + LL
    for sx in [0.3, -0.3]:
        for sy in eq_signs:
            combos[f"LC{idx}"] = [("DL", DL_12), ("EQx", sx), ("EQy", sy), ("LL", 1.0)]
            idx += 1

    # Third group: 0.9DL + EQx ± 0.3EQy + LL
    for sx in eq_signs:
        for sy in [0.3, -0.3]:
            combos[f"LC{idx}"] = [("DL", DL_09), ("EQx", sx), ("EQy", sy), ("LL", 1.0)]
            idx += 1

    # Fourth group: 0.9DL + 0.3EQx ± EQy + LL
    for sx in [0.3, -0.3]:
        for sy in eq_signs:
            combos[f"LC{idx}"] = [("DL", DL_09), ("EQx", sx), ("EQy", sy), ("LL", 1.0)]
            idx += 1

    return combos

# Example usage:
Z = 0.28
S = 1.15
load_combinationsiiii = generate_seismic_load_combinations(Z, S)
print(load_combinationsiiii)

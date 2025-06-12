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











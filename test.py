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


# Define original list of nodes
original_nodes = [
    {"id": 1, "name": "n1", "coord": [0.0, 0.0, 0.0]},
    {"id": 2, "name": "n2", "coord": [1.0, 0.0, 0.0]},
    {"id": 3, "name": "n3", "coord": [1.0, 1.0, 0.0]},
    {"id": 4, "name": "n4", "coord": [0.0, 1.0, 0.0]},
]

# Create duplicate nodes
duplicate_nodes = []
start_id = 101
for i, node in enumerate(original_nodes):
    duplicate_nodes.append({
        "id": start_id + i,
        "name": f"{node['name']}_dup",
        "coord": node["coord"]
    })

# Create zeroLength element definitions
zero_length_elements = []
start_ele_id = 1001
for i in range(len(original_nodes)):
    ele = {
        "eleTag": start_ele_id + i,
        "eleNodes": [original_nodes[i]["id"], duplicate_nodes[i]["id"]],
        "matTags": [1, 2],
        "dirs": [1, 2],
        "rFlag": 1,
        "vecx": [1.0, 0.0, 0.0],
        "vecyp": [0.0, 1.0, 0.0]
    }
    zero_length_elements.append(ele)

# Function to create zeroLength elements
def create_zero_length_elements(elements_data):
    for data in elements_data:
        cmd = ["element", "zeroLength", data["eleTag"]]
        cmd += data["eleNodes"]
        cmd += ["-mat"] + data["matTags"]
        cmd += ["-dir"] + data["dirs"]
        if "rFlag" in data:
            cmd += ["-doRayleigh", data["rFlag"]]
        if "vecx" in data and "vecyp" in data:
            cmd += ["-orient"] + data["vecx"] + data["vecyp"]
        print(" ".join(map(str, cmd)))

# Call and print
create_zero_length_elements(zero_length_elements)

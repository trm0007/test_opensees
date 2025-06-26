# import math

# def calculate_reinforcement_spacing(moment_kipft, width_inch, d_inch, fy_psi, bar_dia_inch, phi=0.9):
#     Mu_inlb = moment_kipft * 12 * 1000  # kip-ft to in-lb
#     d = d_inch
#     fy = fy_psi
#     Ast_required = Mu_inlb / (phi * 0.87 * fy * d * (1 - (0.42 * d / d)))  # simplified rectangular
#     area_single_bar = math.pi * (bar_dia_inch ** 2) / 4  # in²
#     spacing_inch = (area_single_bar * width_inch) / Ast_required
#     return spacing_inch

# # Input values
# moments = [3.5, 4.0, 12.0]  # in kip-ft
# width = 12  # inch (1 ft strip)
# d = 6  # inch effective depth
# fy = 60000  # psi
# bar_dia = 5/8  # inch (e.g., #4 bar = 1/2 inch)

# for m in moments:
#     spacing = calculate_reinforcement_spacing(m, width, d, fy, bar_dia)
#     spacing_rounded = max(3, min(round(spacing), 12))  # round to nearest inch, limit 3" to 12"
#     print(f"Provide {bar_dia*25.4:.0f} mm bar @ {spacing_rounded} inch c/c for moment {m} kip-ft")


# import numpy as np
# from rich.pretty import pprint
# from concreteproperties.material import Concrete, SteelBar
# from concreteproperties.stress_strain_profile import (
#     EurocodeNonLinear,
#     RectangularStressBlock,
#     SteelElasticPlastic,
# )
# from sectionproperties.pre.library.concrete_sections import concrete_circular_section
# from concreteproperties.concrete_section import ConcreteSection

# concrete = Concrete(
#     name="40 MPa Concrete",
#     density=2.4e-6,
#     stress_strain_profile=EurocodeNonLinear(
#         elastic_modulus=32.8e3,
#         ultimate_strain=0.0035,
#         compressive_strength=40,
#         compressive_strain=0.0023,
#         tensile_strength=3.8,
#         tension_softening_stiffness=10e3,
#     ),
#     ultimate_stress_strain_profile=RectangularStressBlock(
#         compressive_strength=40,
#         alpha=0.79,
#         gamma=0.87,
#         ultimate_strain=0.003,
#     ),
#     flexural_tensile_strength=3.8,
#     colour="lightgrey",
# )

# steel = SteelBar(
#     name="500 MPa Steel",
#     density=7.85e-6,
#     stress_strain_profile=SteelElasticPlastic(
#         yield_strength=500,
#         elastic_modulus=200e3,
#         fracture_strain=0.05,
#     ),
#     colour="grey",)

# geom = concrete_circular_section(
#     d=600,
#     area_conc=np.pi * 600 * 600 / 4,
#     n_conc=32,
#     dia_bar=20,
#     area_bar=310,
#     n_bar=10,
#     cover=45,
#     n_circle=4,
#     conc_mat=concrete,
#     steel_mat=steel,
# )


# conc_sec = ConcreteSection(geom)
# conc_sec.plot_section()

# uncr_stress_res_1 = conc_sec.calculate_uncracked_stress(m_x=50e6)
# uncr_stress_res_2 = conc_sec.calculate_uncracked_stress(m_x=25e6, m_y=35e6, n=200e3)


# # stresses from first analysis section
# stresses = uncr_stress_res_1.concrete_stresses[0]
# print("Concrete Stresses (first analysis section):")
# print(stresses)

# # total concrete force in X direction from first analysis section
# N, Mx, My = uncr_stress_res_1.concrete_forces[0]
# print("\nConcrete Forces (first analysis section):")
# print(f"Axial Force (N): {N}")
# print(f"Moment about X (Mx): {Mx}")
# print(f"Moment about Y (My): {My}")

# # reinforcement stress and strain at index 0
# stress_0 = uncr_stress_res_1.lumped_reinforcement_stresses[0]
# strain_0 = uncr_stress_res_1.lumped_reinforcement_strains[0]
# print("\nReinforcement Stress at Index 0:")
# print(stress_0)

# print("Reinforcement Strain at Index 0:")
# print(strain_0)

# cracked_res = conc_sec.calculate_cracked_properties(theta=0)
# print(f"M_cr = {cracked_res.m_cr / 1e6:.2f} kN.m")
# print(f"d_n,c = {cracked_res.d_nc:.2f} mm")


def create_slabs(node_data):
    slabs = {}
    
    # Extract all nodes and sort them by y then x coordinates for proper ordering
    nodes = sorted(node_data.keys(), key=lambda n: (node_data[n][1], node_data[n][0]))
    
    # Group nodes by their y-coordinate to identify rows
    y_coords = sorted({node_data[n][1] for n in nodes})
    rows = {}
    for y in y_coords:
        rows[y] = sorted([n for n in nodes if node_data[n][1] == y], 
                         key=lambda n: node_data[n][0])
    
    # Create rectangular slabs (4-noded) where possible
    rect_slab_count = 1
    tri_slab_count = 1
    
    for i in range(len(y_coords)-1):
        y1 = y_coords[i]
        y2 = y_coords[i+1]
        row1 = rows[y1]
        row2 = rows[y2]
        
        # Try to match nodes between rows to form rectangles
        j = 0
        while j < len(row1)-1 and j < len(row2)-1:
            n1 = row1[j]
            n2 = row1[j+1]
            n3 = row2[j+1]
            n4 = row2[j]
            
            # Check if these nodes form a proper quadrilateral
            if (node_data[n1][0] == node_data[n4][0] and  # x-coord matches vertically
                node_data[n2][0] == node_data[n3][0] and  # x-coord matches vertically
                node_data[n1][1] == node_data[n2][1] and  # y-coord matches horizontally
                node_data[n4][1] == node_data[n3][1]):   # y-coord matches horizontally
                
                # Create rectangular slab with nodes ordered anticlockwise
                slab_name = f"rect_slab_{rect_slab_count}"
                slabs[slab_name] = {
                    "nodes": [n1, n2, n3, n4],  # Already ordered anticlockwise
                    "coordinates": {
                        n1: node_data[n1],
                        n2: node_data[n2],
                        n3: node_data[n3],
                        n4: node_data[n4]
                    },
                    "type": "rectangular",
                    # "area": calculate_area([n1, n2, n3, n4], node_data)
                }
                rect_slab_count += 1
                j += 1  # Move to next potential rectangle
            else:
                # Can't form rectangle, try triangles instead
                # Create first triangle (n1, n2, n4)
                slab_name = f"tri_slab_{tri_slab_count}"
                nodes_tri = order_nodes_anticlockwise([n1, n2, n4], node_data)
                slabs[slab_name] = {
                    "nodes": nodes_tri,
                    "coordinates": {n: node_data[n] for n in nodes_tri},
                    "type": "triangular",
                    # "area": calculate_area(nodes_tri, node_data)
                }
                tri_slab_count += 1
                
                # Create second triangle (n2, n3, n4) if possible
                if j+1 < len(row2):
                    slab_name = f"tri_slab_{tri_slab_count}"
                    nodes_tri = order_nodes_anticlockwise([n2, n3, n4], node_data)
                    slabs[slab_name] = {
                        "nodes": nodes_tri,
                        "coordinates": {n: node_data[n] for n in nodes_tri},
                        "type": "triangular",
                        # "area": calculate_area(nodes_tri, node_data)
                    }
                    tri_slab_count += 1
                j += 1
    
    return slabs

def order_nodes_anticlockwise(nodes, node_data):
    """Order nodes in anticlockwise direction"""
    if len(nodes) == 3:
        # For triangles, calculate centroid
        x = sum(node_data[n][0] for n in nodes) / 3
        y = sum(node_data[n][1] for n in nodes) / 3
        
        # Calculate angles from centroid
        def angle_from_centroid(n):
            dx = node_data[n][0] - x
            dy = node_data[n][1] - y
            return math.atan2(dy, dx)
        
        return sorted(nodes, key=angle_from_centroid, reverse=True)
    elif len(nodes) == 4:
        # For quadrilaterals, order is already anticlockwise in our creation method
        return nodes
    return nodes



# Example usage with your provided data
import math

node_data = {
    "n36": [0.0, 0.0, 10.0], "n37": [10.0, 0.0, 10.0], "n38": [22.0, 0.0, 10.0], 
    "n39": [32.0, 0.0, 10.0], "n40": [44.0, 0.0, 10.0], "n41": [0.0, 11.0, 10.0], 
    "n42": [10.0, 11.0, 10.0], "n43": [22.0, 11.0, 10.0], "n44": [32.0, 11.0, 10.0], 
    "n45": [44.0, 11.0, 10.0], "n46": [0.0, 24.0, 10.0], "n47": [10.0, 24.0, 10.0], 
    "n48": [22.0, 24.0, 10.0], "n49": [32.0, 24.0, 10.0], "n50": [44.0, 24.0, 10.0], 
    "n51": [0.0, 35.0, 10.0], "n52": [10.0, 35.0, 10.0], "n53": [22.0, 35.0, 10.0], 
    "n54": [32.0, 35.0, 10.0], "n55": [44.0, 35.0, 10.0], "n56": [0.0, 48.0, 10.0], 
    "n57": [10.0, 48.0, 10.0], "n58": [22.0, 48.0, 10.0], "n59": [32.0, 48.0, 10.0], 
    "n60": [44.0, 48.0, 10.0], "n61": [0.0, 59.0, 10.0], "n62": [10.0, 59.0, 10.0], 
    "n63": [22.0, 59.0, 10.0], "n64": [32.0, 59.0, 10.0], "n65": [44.0, 59.0, 10.0], 
    "n66": [0.0, 72.0, 10.0], "n67": [10.0, 72.0, 10.0], "n68": [22.0, 72.0, 10.0], 
    "n69": [32.0, 72.0, 10.0], "n70": [44.0, 72.0, 10.0]
}

slabs = create_slabs(node_data)

# Print the results
for slab_name, slab_data in slabs.items():
    print(f"{slab_name}:")
    print(f"  Type: {slab_data['type']}")
    print(f"  Nodes (anticlockwise): {slab_data['nodes']}")
    print("  Coordinates:")
    for node, coord in slab_data['coordinates'].items():
        print(f"    {node}: {coord}")
    print()
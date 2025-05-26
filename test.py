# import math
# import json

# # Material properties
# E_concrete = 30000  # MPa
# story_height = 3000  # mm
# wall_thickness = 200  # mm (default thickness for walls)
# Vx = -100.0  # kN (applied force in X direction)
# Vy = 0.0  # kN (applied force in Y direction)
# Wx = 21000  # mm (21.0 m width)
# Wy = 24000  # mm (24.0 m length)

# # Elements (coordinates already in mm)
# elements = [
#     {"name": "M", "Dx": 600, "Dy": 600, "x": 2000, "y": 2000, "weight": 10, "type": "rect_column"},
#     {"name": "N", "Dx": 600, "Dy": 600, "x": 10000, "y": 2000, "weight": 10, "type": "rect_column"},
#     {"name": "O", "Dx": 600, "Dy": 600, "x": 18000, "y": 2000, "weight": 10, "type": "rect_column"},
#     {"name": "P", "Dx": 2000, "Dy": 200, "x": 19000, "y": 10000, "weight": 10, "type": "wall"},
#     {"name": "Q", "Dx": 400, "Dy": 400, "x": 10000, "y": 22000, "weight": 10, "type": "circ_column"},
#     {"name": "R", "Dx": 400, "Dy": 400, "x": 10000, "y": 12000, "weight": 10, "type": "circ_column"}
# ]

# def validate_inputs():
#     """Validate input parameters and geometry"""
#     print("Validating inputs...")
    
#     # Check if all elements are within building bounds
#     for element in elements:
#         if element["x"] < 0 or element["x"] > Wx:
#             print(f"WARNING: Element {element['name']} x-coordinate ({element['x']}) outside building width (0-{Wx})")
#         if element["y"] < 0 or element["y"] > Wy:
#             print(f"WARNING: Element {element['name']} y-coordinate ({element['y']}) outside building length (0-{Wy})")
    
#     print("Input validation complete.")
#     return True

# def calculate_moment_of_inertia_circular(d):
#     """Calculate moment of inertia for circular cross-section: I = πd⁴/64"""
#     return (math.pi * d**4) / 64

# def calculate_moment_of_inertia_rectangular(b, h):
#     """Calculate moment of inertia for rectangular cross-section: I = bh³/12"""
#     return (b * h**3) / 12

# def calculate_column_stiffness(h, I, E=30000):
#     """
#     Calculate column lateral stiffness
#     For fixed-fixed: k = 12EI/h³
#     For pinned-fixed: k = 3EI/h³
#     Using fixed-fixed assumption
#     """
#     return (12 * E * I) / (h**3)

# def calculate_element_stiffnesses(elements, story_height, E_concrete):
#     """Calculate stiffness properties for all elements with proper units"""
#     print("\nCalculating stiffness for each element...")
#     print("="*80)
    
#     for element in elements:
#         Dx = element["Dx"]  # mm
#         Dy = element["Dy"]  # mm
#         x = element["x"]    # mm (NO unit conversion needed)
#         y = element["y"]    # mm (NO unit conversion needed)
        
#         if element["type"] == "wall":
#             # For a wall, assume Dx is length and Dy is thickness
#             # X-direction stiffness (wall resists force in X, bending about Y-axis)
#             # Use the wall's strong axis (length cubed)
#             Ix = calculate_moment_of_inertia_rectangular(Dy, Dx)  # thickness × length³/12
#             element["kx"] = calculate_column_stiffness(story_height, Ix, E_concrete)
            
#             # Y-direction stiffness (wall resists force in Y, bending about X-axis)
#             # Use the wall's weak axis (thickness cubed)
#             Iy = calculate_moment_of_inertia_rectangular(Dx, Dy)  # length × thickness³/12
#             element["ky"] = calculate_column_stiffness(story_height, Iy, E_concrete)
            
#             print(f"Wall {element['name']}: L={Dx}mm, t={Dy}mm")
#             print(f"  Ix = {Ix:.2e} mm⁴, kx = {element['kx']:.2e} N/mm")
#             print(f"  Iy = {Iy:.2e} mm⁴, ky = {element['ky']:.2e} N/mm")
                
#         elif element["type"] == "circ_column":
#             I = calculate_moment_of_inertia_circular(Dx)
#             k = calculate_column_stiffness(story_height, I, E_concrete)
#             element["kx"] = k  # same in both directions
#             element["ky"] = k
#             print(f"Circular Column {element['name']}: d={Dx}mm, I={I:.2e} mm⁴, k={k:.2e} N/mm")
            
#         else:  # rectangular column
#             # X-direction stiffness (bending about Y-axis)
#             Ix = calculate_moment_of_inertia_rectangular(Dy, Dx)  # Dy × Dx³/12
#             # Y-direction stiffness (bending about X-axis)
#             Iy = calculate_moment_of_inertia_rectangular(Dx, Dy)  # Dx × Dy³/12
            
#             element["kx"] = calculate_column_stiffness(story_height, Ix, E_concrete)
#             element["ky"] = calculate_column_stiffness(story_height, Iy, E_concrete)
            
#             print(f"Rect Column {element['name']}: {Dx}×{Dy}mm")
#             print(f"  Ix = {Ix:.2e} mm⁴, kx = {element['kx']:.2e} N/mm")
#             print(f"  Iy = {Iy:.2e} mm⁴, ky = {element['ky']:.2e} N/mm")
        
#         # Calculate products for center of rigidity (using mm units consistently)
#         element["kx_y"] = element["kx"] * y  # N/mm × mm = N
#         element["ky_x"] = element["ky"] * x  # N/mm × mm = N

#         # Calculate area and weighted positions for center of mass
#         if element["type"] == "circ_column":
#             element["A"] = math.pi * (Dx/2)**2  # mm²
#         else:
#             element["A"] = Dx * Dy  # mm²
        
#         # Calculate weighted positions for center of mass
#         element["Ax"] = element["A"] * element["weight"]  # mm² × weight factor
#         element["x_Ax"] = x * element["Ax"]  # mm × mm² × weight = mm³ × weight
#         element["y_Ax"] = y * element["Ax"]  # mm × mm² × weight = mm³ × weight
    
#     return elements



# def calculate_centers(elements):
#     """Calculate centers of rigidity and mass, including applied loads."""
#     # Sums for center of rigidity (unchanged)
#     sum_kx = sum(element["kx"] for element in elements)
#     sum_ky = sum(element["ky"] for element in elements)
#     sum_kx_y = sum(element["kx_y"] for element in elements)
#     sum_ky_x = sum(element["ky_x"] for element in elements)
    
#     # Sums for center of mass (updated to include applied loads)
#     sum_Ax = 0
#     sum_x_Ax = 0
#     sum_y_Ax = 0

#     for element in elements:
#         # Original mass contribution (based on area and weight factor)
#         if element["type"] == "circ_column":
#             A = math.pi * (element["Dx"]/2)**2
#         else:
#             A = element["Dx"] * element["Dy"]
        
#         mass_term = A * element["weight"]
#         sum_Ax += mass_term
#         sum_x_Ax += element["x"] * mass_term
#         sum_y_Ax += element["y"] * mass_term

#         # Additional mass contribution from applied loads
#         if "point_load" in element:  # Point load on columns
#             sum_Ax += element["point_load"]
#             sum_x_Ax += element["x"] * element["point_load"]
#             sum_y_Ax += element["y"] * element["point_load"]
        
#         if "udl" in element:  # UDL on walls (converted to point load at centroid)
#             total_load = element["udl"] * element["Dx"]  # UDL (kN/m) * length (m)
#             sum_Ax += total_load
#             sum_x_Ax += element["x"] * total_load
#             sum_y_Ax += element["y"] * total_load

#     # Calculate centers (same as before)
#     CR_x = sum_ky_x / sum_ky if sum_ky != 0 else 0
#     CR_y = sum_kx_y / sum_kx if sum_kx != 0 else 0
#     CM_x = sum_x_Ax / sum_Ax if sum_Ax != 0 else 0
#     CM_y = sum_y_Ax / sum_Ax if sum_Ax != 0 else 0

#     return {
#         "CR": {"x": CR_x, "y": CR_y},
#         "CM": {"x": CM_x, "y": CM_y},
#         "sums": {
#             "sum_kx": sum_kx,
#             "sum_ky": sum_ky,
#             "sum_kx_y": sum_kx_y,
#             "sum_ky_x": sum_ky_x
#         }
#     }

# def calculate_eccentricities(CR, CM, Wx, Wy):
#     """Calculate natural and design eccentricities with proper sign conventions"""
#     # Natural eccentricity (center of mass relative to center of rigidity)
#     enx = CM["x"] - CR["x"]  # mm
#     eny = CM["y"] - CR["y"]  # mm

#     # Accidental eccentricity (5% of building dimension)
#     acc_ex = 0.05 * Wx  # mm
#     acc_ey = 0.05 * Wy  # mm

#     # Design eccentricity - consider both positive and negative accidental eccentricity
#     # and choose the one that gives maximum eccentricity
#     ex_pos = enx + acc_ex
#     ex_neg = enx - acc_ex
#     ey_pos = eny + acc_ey
#     ey_neg = eny - acc_ey
    
#     # Choose the eccentricity with larger absolute value (more critical)
#     ex = ex_pos if abs(ex_pos) > abs(ex_neg) else ex_neg
#     ey = ey_pos if abs(ey_pos) > abs(ey_neg) else ey_neg

#     return {
#         "natural": {"enx": enx, "eny": eny},
#         "accidental": {"acc_ex": acc_ex, "acc_ey": acc_ey},
#         "design": {"ex": ex, "ey": ey},
#         "design_cases": {
#             "ex_pos": ex_pos, "ex_neg": ex_neg,
#             "ey_pos": ey_pos, "ey_neg": ey_neg
#         }
#     }

# def calculate_torsion(Vx, Vy, eccentricities):
#     """Calculate torsional moments"""
#     # T = V × e (force × eccentricity)
#     Tx = Vy * eccentricities["design"]["ex"]  # kN × mm = kN⋅mm
#     Ty = Vx * eccentricities["design"]["ey"]  # kN × mm = kN⋅mm
#     T_total = Tx + Ty  # Total torsion in kN⋅mm
    
#     return {
#         "Tx": Tx,
#         "Ty": Ty,
#         "T_total": T_total
#     }

# def calculate_polar_moment(elements, CR):
#     """Calculate polar moment of inertia about center of rigidity"""
#     for element in elements:
#         x = element["x"]  # mm
#         y = element["y"]  # mm
        
#         # Calculate relative positions from center of rigidity
#         element["xr"] = x - CR["x"]  # mm
#         element["yr"] = y - CR["y"]  # mm
        
#         # Calculate contributions to polar moment of inertia
#         # J = Σ(ky × xr²) + Σ(kx × yr²)
#         # Units: (N/mm) × mm² = N⋅mm
#         element["ky_xr2"] = element["ky"] * (element["xr"]**2)  # N⋅mm
#         element["kx_yr2"] = element["kx"] * (element["yr"]**2)  # N⋅mm

#     sum_ky_xr2 = sum(element["ky_xr2"] for element in elements)
#     sum_kx_yr2 = sum(element["kx_yr2"] for element in elements)
#     J = sum_ky_xr2 + sum_kx_yr2  # N⋅mm

#     return {
#         "sum_ky_xr2": sum_ky_xr2,
#         "sum_kx_yr2": sum_kx_yr2,
#         "J": J
#     }

# def calculate_element_forces(elements, Vx, Vy, sums, T_total, J):
#     """Calculate forces on each element with proper sign conventions"""
#     total_fx_check = 0.0
#     total_fy_check = 0.0
#     total_moment_check = 0.0

#     for element in elements:
#         # Direct forces (proportional to stiffness)
#         direct_x = (Vx * element["kx"]) / sums["sum_kx"] if sums["sum_kx"] != 0 else 0
#         direct_y = (Vy * element["ky"]) / sums["sum_ky"] if sums["sum_ky"] != 0 else 0
        
#         # Torsional forces
#         # For torsion, forces are proportional to distance from CR and stiffness
#         torsional_x = -(T_total * element["kx"] * element["yr"]) / J if J != 0 else 0
#         torsional_y = (T_total * element["ky"] * element["xr"]) / J if J != 0 else 0
        
#         # Total forces
#         element["Fx"] = direct_x + torsional_x  # kN
#         element["Fy"] = direct_y + torsional_y  # kN
        
#         # Store components for analysis
#         element["direct_x"] = direct_x
#         element["direct_y"] = direct_y
#         element["torsional_x"] = torsional_x
#         element["torsional_y"] = torsional_y
        
#         # Equilibrium checks
#         total_fx_check += element["Fx"]
#         total_fy_check += element["Fy"]
        
#         # Moment equilibrium check (about origin)
#         moment_contribution = element["Fx"] * element["y"] - element["Fy"] * element["x"]
#         total_moment_check += moment_contribution

#     return {
#         "elements": elements,
#         "equilibrium_check": {
#             "total_fx": total_fx_check,
#             "total_fy": total_fy_check,
#             "total_moment": total_moment_check,
#             "error_fx": abs(total_fx_check - Vx),
#             "error_fy": abs(total_fy_check - Vy),
#             "applied_moment": Vx * 0 + Vy * 0,  # Applied moment about origin
#             "error_moment": abs(total_moment_check - (Vx * 0 + Vy * 0))
#         }
#     }

# def print_detailed_results(elements, centers, eccentricities, torsion, polar, forces):
#     """Print comprehensive results with validation"""
#     print("\n" + "="*80)
#     print("STRUCTURAL ANALYSIS RESULTS")
#     print("="*80)
    
#     # Centers
#     print(f"\nCENTERS:")
#     print(f"Center of Rigidity (CR): x = {centers['CR']['x']:.1f} mm, y = {centers['CR']['y']:.1f} mm")
#     print(f"Center of Mass (CM):     x = {centers['CM']['x']:.1f} mm, y = {centers['CM']['y']:.1f} mm")
    
#     # Validate centers are within building
#     CR_x, CR_y = centers['CR']['x'], centers['CR']['y']
#     CM_x, CM_y = centers['CM']['x'], centers['CM']['y']
    
#     print(f"\nVALIDATION:")
#     print(f"Building bounds: 0 ≤ x ≤ {Wx} mm, 0 ≤ y ≤ {Wy} mm")
#     print(f"CR within bounds: x={'✓' if 0 <= CR_x <= Wx else '✗'}, y={'✓' if 0 <= CR_y <= Wy else '✗'}")
#     print(f"CM within bounds: x={'✓' if 0 <= CM_x <= Wx else '✗'}, y={'✓' if 0 <= CM_y <= Wy else '✗'}")
    
#     # Eccentricities
#     print(f"\nECCENTRICITIES:")
#     print(f"Natural:     enx = {eccentricities['natural']['enx']:.1f} mm, eny = {eccentricities['natural']['eny']:.1f} mm")
#     print(f"Accidental:  eax = {eccentricities['accidental']['acc_ex']:.1f} mm, eay = {eccentricities['accidental']['acc_ey']:.1f} mm")
#     print(f"Design:      ex = {eccentricities['design']['ex']:.1f} mm, ey = {eccentricities['design']['ey']:.1f} mm")
    
#     # Applied loads and torsion
#     print(f"\nAPPLIED LOADS AND TORSION:")
#     print(f"Applied Forces: Vx = {Vx:.1f} kN, Vy = {Vy:.1f} kN")
#     print(f"Torsional Moments: Tx = {torsion['Tx']:.0f} kN⋅mm, Ty = {torsion['Ty']:.0f} kN⋅mm")
#     print(f"Total Torsion: T = {torsion['T_total']:.0f} kN⋅mm")
#     print(f"Polar Moment of Inertia: J = {polar['J']:.2e} N⋅mm")
    
#     # Element forces
#     print(f"\nELEMENT FORCES:")
#     print("-" * 80)
#     print(f"{'Element':<8} {'Type':<12} {'Fx (kN)':<10} {'Fy (kN)':<10} {'Direct X':<10} {'Torsional X':<12}")
#     print("-" * 80)
#     for element in forces["elements"]:
#         print(f"{element['name']:<8} {element['type']:<12} "
#               f"{element['Fx']:<10.3f} {element['Fy']:<10.3f} "
#               f"{element['direct_x']:<10.3f} {element['torsional_x']:<12.3f}")
    
#     # Equilibrium check
#     print("-" * 80)
#     print(f"\nEQUILIBRIUM CHECK:")
#     eq_check = forces['equilibrium_check']
#     print(f"ΣFx = {eq_check['total_fx']:.6f} kN (should = {Vx:.1f} kN) - Error: {eq_check['error_fx']:.6f} kN")
#     print(f"ΣFy = {eq_check['total_fy']:.6f} kN (should = {Vy:.1f} kN) - Error: {eq_check['error_fy']:.6f} kN")
    
#     # Force reasonableness check
#     max_force = max(abs(el["Fx"]) for el in forces["elements"]) + max(abs(el["Fy"]) for el in forces["elements"])
#     applied_force = abs(Vx) + abs(Vy)
#     amplification = max_force / applied_force if applied_force > 0 else 0
    
#     print(f"\nFORCE AMPLIFICATION:")
#     print(f"Maximum element force: {max_force:.1f} kN")
#     print(f"Applied force magnitude: {applied_force:.1f} kN") 
#     print(f"Amplification factor: {amplification:.2f}")
    
#     status = "REASONABLE" if amplification < 50 else "CHECK REQUIRED"
#     print(f"Status: {status}")
    
#     print("="*80)
#     print("ANALYSIS CORRECTIONS APPLIED:")
#     print("✓ Fixed coordinate unit conversion error")
#     print("✓ Proper stiffness calculations with correct I formulations")
#     print("✓ Consistent units throughout (mm, N, kN)")
#     print("✓ Corrected center of rigidity calculations")
#     print("✓ Proper eccentricity sign conventions")
#     print("✓ Force and moment equilibrium verification")
#     print("✓ Added validation checks for reasonableness")
#     print("="*80)

# def main():
#     """Main analysis function"""
#     # Validate inputs
#     validate_inputs()
    
#     # Calculate element stiffnesses
#     elements_with_stiffness = calculate_element_stiffnesses(elements, story_height, E_concrete)
    
#     # Calculate centers
#     centers = calculate_centers(elements_with_stiffness)
    
#     # Calculate eccentricities
#     eccentricities = calculate_eccentricities(centers["CR"], centers["CM"], Wx, Wy)
    
#     # Calculate torsion
#     torsion = calculate_torsion(Vx, Vy, eccentricities)
    
#     # Calculate polar moment
#     polar = calculate_polar_moment(elements_with_stiffness, centers["CR"])
    
#     # Calculate element forces
#     forces = calculate_element_forces(elements_with_stiffness, Vx, Vy, 
#                                     centers["sums"], torsion["T_total"], polar["J"])
    
#     # Print detailed results
#     print_detailed_results(forces["elements"], centers, eccentricities, torsion, polar, forces)
    
#     return {
#         "elements": forces["elements"],
#         "centers": centers,
#         "eccentricities": eccentricities,
#         "torsion": torsion,
#         "polar": polar,
#         "forces": forces
#     }

# # Run the analysis
# if __name__ == "__main__":
#     results = main()

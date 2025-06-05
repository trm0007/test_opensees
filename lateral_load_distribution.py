import math

def lateral_load_distribution(Lx, Ly,E, poisson_ratio, fc, Vx, Vy, components, elements):
    """Perform all structural calculations in one function."""
    
    # 1. Calculate center of mass
    # print("\n1. CENTER OF MASS CALCULATION:")
    total_wt = sum(comp["wt"] for comp in components)
    sum_wt_xi = sum(comp["wt"] * comp["xi"] for comp in components)
    sum_wt_yi = sum(comp["wt"] * comp["yi"] for comp in components)
    Xm = sum_wt_xi / total_wt
    Ym = sum_wt_yi / total_wt
    COM =(Xm, Ym)
    total_Kx = 0
    total_Ky = 0
    sum_Kx_x = 0
    sum_Ky_y = 0
    element_details = []
    
    for element in elements:
        Dx = element["Dx"]
        Dy = element["Dy"]
        x = element["x"]
        y = element["y"]
        height = element["height"]
        element_type = element.get("type")
        
        # Calculate stiffness in both directions
        Ix = (Dy * Dx**3) / 12
        Iy = (Dx * Dy**3) / 12
        
        
        if element_type == "column":
            Kx = (12 * E * Ix) / (height**3)
            Ky = (12 * E * Iy) / (height**3)
        else:
            denom_x = 1 + 0.6 * (1 + poisson_ratio) * (Dx**2 / height**2)
            denom_y = 1 + 0.6 * (1 + poisson_ratio) * (Dy**2 / height**2)
            Kx = (3 * E * Ix) / (height**3 * denom_x)
            Ky = (3 * E * Iy) / (height**3 * denom_y)
        
        element_details.append({
            "name": element.get("name", "Unknown"),
            "Kx": Kx,
            "Ky": Ky,
            "x": x,
            "y": y
        })
        
        total_Kx += Kx
        total_Ky += Ky
        sum_Kx_x += Kx * x
        sum_Ky_y += Ky * y
    
    Xr = sum_Kx_x / total_Kx if total_Kx > 0 else 0
    Yr = sum_Ky_y / total_Ky if total_Ky > 0 else 0
    COR = (Xr, Yr)


    # 3. Calculate eccentricity
    # print("\n3. ECCENTRICITY CALCULATION:")
    ex = Xm - Xr
    ey = Ym - Yr
    min_ex = 0.05 * Lx
    min_ey = 0.05 * Ly
    ex_total = ex + (-min_ex if ex < 0 else min_ex)
    ey_total = ey + (-min_ey if ey < 0 else min_ey)
    eccentricity = (ex, ey)
    eccentricity_total = (ex_total, ey_total)

    # 4. Calculate torsional moment
    T = Vy * ex_total + Vx * ey_total

    # 5. Calculate torsional rigidity
    # print("\n4. TORSIONAL RIGIDITY CALCULATION:")
    Jr_total = 0
    Jr_X = 0
    Jr_Y = 0
    

    for element in element_details:
        dx = element["x"] - Xr
        dy = element["y"] - Yr
        Kx = element["Kx"]
        Ky = element["Ky"]

        Jr_contribution_x = Kx * dy**2
        Jr_contribution_y = Ky * dx**2
        
        Jr_X += Jr_contribution_x
        Jr_Y += Jr_contribution_y

        # print(f"{element['name']:<8} {Kx/1e6:<12.2f} {Ky/1e6:<12.2f} {dx:<8.2f} {dy:<8.2f} {Jr_contribution_x/1e6:<12.2f} {Jr_contribution_y/1e6:<12.2f}")

    Jr_total = Jr_X + Jr_Y
    torsional_rigidity = (Jr_X, Jr_Y, Jr_total)

    # 6. Calculate floor shear distribution
    print("\n5. FLOOR SHEAR DISTRIBUTION:")
    results = []
    total_Fx = 0
    total_Fy = 0
    
    for element in element_details:
        dx = element['x'] - Xr
        dy = element['y'] - Yr
        Kx = element['Kx']
        Ky = element['Ky']
        
        Fy = (Vy * Ky / total_Ky) + (T * Kx * dx / Jr_total)
        Fx = (Vx * Kx / total_Kx) - (T * Ky * dy / Jr_total)
        
        results.append({
            'name': element['name'],
            'Fx': Fx,
            'Fy': Fy,
            'dx': dx,
            'dy': dy
        })
        total_Fx += Fx
        total_Fy += Fy

    return {
        "center_of_mass": COM,
        "center_of_rigidity": COR,
        "total_stiffness": (total_Kx, total_Ky),
        "element_details": element_details,
        "eccentricity": eccentricity,
        "eccentricity_total": eccentricity_total,
        "torsional_moment": T,
        "torsional_rigidity": torsional_rigidity,
        "shear_distribution": results,
        "total_shear": (total_Fx, total_Fy)
    }




# Example usage with converted units (inches and kips)
components = [
    {"wt": 629.4652, "xi": 314.9608, "yi": 246.063},  # 2800 kg = 629.4652 kips, 8m = 314.9608in, 6.25m = 246.063in
    {"wt": -18.88396, "xi": 551.1814, "yi": 462.598},  # -84 kg = -18.88396 kips, 14m = 551.1814in, 11.75m = 462.598in
    {"wt": 25.898, "xi": 5.905515, "yi": 236.2206},    # 115.2 kg = 25.898 kips, 0.15m = 5.905515in, 6m = 236.2206in
    {"wt": 25.898, "xi": 314.9608, "yi": 236.2206},
    {"wt": 25.898, "xi": 623.9961, "yi": 236.2206},    # 15.85m = 623.9961in
    {"wt": 25.898, "xi": 78.7402, "yi": 5.905515},      # 2m = 78.7402in, 0.15m = 5.905515in
    {"wt": 25.898, "xi": 551.1814, "yi": 5.905515},
    {"wt": 25.898, "xi": 78.7402, "yi": 486.2206},      # 12.35m = 486.2206in
    {"wt": 25.898, "xi": 551.1814, "yi": 427.1655},     # 10.85m = 427.1655in
]

elements = [
    {"name": "W2-1", "type": "wall", "Dx": 11.811, "Dy": 157.4804, "height": 157.4804, "x": 78.7402, "y": 5.905515},
    {"name": "W2-2", "type": "wall", "Dx": 11.811, "Dy": 157.4804, "height": 157.4804, "x": 551.1814, "y": 5.905515},
    {"name": "W2-3", "type": "wall", "Dx": 11.811, "Dy": 157.4804, "height": 157.4804, "x": 78.7402, "y": 486.2206},
    {"name": "W2-4", "type": "wall", "Dx": 11.811, "Dy": 157.4804, "height": 157.4804, "x": 551.1814, "y": 427.1655},
    {"name": "W1-1", "type": "wall", "Dx": 157.4804, "Dy": 11.811, "height": 157.4804, "x": 5.905515, "y": 236.2206},
    {"name": "W1-2", "type": "wall", "Dx": 157.4804, "Dy": 11.811, "height": 157.4804, "x": 314.9608, "y": 236.2206},
    {"name": "W1-3", "type": "wall", "Dx": 157.4804, "Dy": 11.811, "height": 157.4804, "x": 623.9961, "y": 236.2206},
]

Lx = 629.9216  # 16m = 629.9216in
Ly = 492.126   # 12.5m = 492.126in
fc = 4.35  # in ksi
E = 57000 * math.sqrt(fc)  # E in ksi
poisson_ratio = 0.25
Vx = 0.0
Vy = -34.163   # -34.163 kN = -7.6806 kips
results = lateral_load_distribution(Lx, Ly, E, poisson_ratio, fc, Vx, Vy, components, elements)
# The rest of the structural analysis code would remain the same, 
# just interpreting all values as inches and kips instead of meters and kN

# Print all results in US Customary Units (inches and kips)
print("\n" + "="*50)
print("STRUCTURAL ANALYSIS RESULTS (US Customary Units)".center(50))
print("="*50)

# 1. Center of Mass and Rigidity
print("\n1. CENTER OF MASS & RIGIDITY")
print("-"*50)
print(f"{'Center of Mass (COM)':<30}: X = {results['center_of_mass'][0]:.3f} in, Y = {results['center_of_mass'][1]:.3f} in")
print(f"{'Center of Rigidity (COR)':<30}: X = {results['center_of_rigidity'][0]:.3f} in, Y = {results['center_of_rigidity'][1]:.3f} in")

# 2. Stiffness Properties
print("\n2. STIFFNESS PROPERTIES")
print("-"*50)
print(f"{'Total Stiffness in X-direction':<30}: {results['total_stiffness'][0]:.3f} kips/in")
print(f"{'Total Stiffness in Y-direction':<30}: {results['total_stiffness'][1]:.3f} kips/in")

# 3. Eccentricity
print("\n3. ECCENTRICITY")
print("-"*50)
print(f"{'Natural Eccentricity (ex, ey)':<30}: ex = {results['eccentricity'][0]:.3f} in, ey = {results['eccentricity'][1]:.3f} in")
print(f"{'Design Eccentricity (ex_total, ey_total)':<30}: ex = {results['eccentricity_total'][0]:.3f} in, ey = {results['eccentricity_total'][1]:.3f} in")

# 4. Torsional Properties
print("\n4. TORSIONAL PROPERTIES")
print("-"*50)
print(f"{'Torsional Moment (T)':<30}: {results['torsional_moment']:.3f} kip·in")
print(f"{'Torsional Rigidity (Jr_X, Jr_Y)':<30}: {results['torsional_rigidity'][0]:.3f}, {results['torsional_rigidity'][1]:.3f} kip·in²")
print(f"{'Total Torsional Rigidity (Jr_total)':<30}: {results['torsional_rigidity'][2]:.3f} kip·in²")

# 5. Element Stiffness Details
print("\n5. ELEMENT STIFFNESS DETAILS")
print("-"*50)
print(f"{'Element':<10}{'Kx (kips/in)':<15}{'Ky (kips/in)':<15}{'X (in)':<10}{'Y (in)':<10}")
print("-"*50)
for element in results['element_details']:
    print(f"{element['name']:<10}{element['Kx']:<15.3f}{element['Ky']:<15.3f}{element['x']:<10.3f}{element['y']:<10.3f}")

# 6. Shear Distribution
print("\n6. SHEAR DISTRIBUTION")
print("-"*50)
print(f"{'Element':<10}{'Fx (kips)':<15}{'Fy (kips)':<15}")
print("-"*50)
for result in results['shear_distribution']:
    print(f"{result['name']:<10}{result['Fx']:<15.3f}{result['Fy']:<15.3f}")

# Print totals
print("-"*50)
print(f"{'TOTAL':<10}{results['total_shear'][0]:<15.3f}{results['total_shear'][1]:<15.3f}")
print("="*50)
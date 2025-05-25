# =============================================
# Unit Definitions (from your input)
# =============================================

# SI (MKS) Base Units
meter = 1.0          # Length (m)
kilogram = 1.0       # Mass (kg)
second = 1.0         # Time (s)
newton = 1.0         # Force (N)
pascal = 1.0         # Pressure (Pa)

# SI Derived Units
kN = 1e3 * newton
MPa = 1e6 * pascal
GPa = 1e9 * pascal
kNm = kN * meter
Nm = newton * meter
milimeter = 1e-3 * meter
# FPS Base Units
inch = 0.0254 * meter        # 1 inch = 0.0254 m
foot = 12.0 * inch           # 1 ft = 12 inches = 0.3048 m
pound = 0.45359237 * kilogram   # 1 lb = 0.45359237 kg
second_fps = second             # Same time unit
lbf = 4.44822 * newton          # 1 lbf = 4.44822 N

# FPS Derived Units
ksi = 1e3 * (lbf / (inch**2))
kip = 1e3 * lbf
kipft = kip * foot
kipin = kip * inch
psf = lbf / (foot**2)
psi = lbf / (inch**2)


# # =============================================
# # Unit Comparison (Differentiation)
# # =============================================

# def compare_units(unit1, unit2, unit1_name, unit2_name, quantity):
#     ratio = unit1 / unit2
#     print(f"1 {unit1_name} = {ratio:.6g} {unit2_name} ({quantity})")

# print("\n===== Length Units =====")
# compare_units(meter, foot, "meter", "foot", "Length")
# compare_units(foot, meter, "foot", "meter", "Length")
# compare_units(inch, meter, "inch", "meter", "Length")
# compare_units(meter, milimeter, "meter", "milimeter", "Length")

# print("\n===== Mass Units =====")
# compare_units(kilogram, pound, "kg", "lb", "Mass")
# compare_units(pound, kilogram, "lb", "kg", "Mass")

# print("\n===== Force Units =====")
# compare_units(newton, lbf, "N", "lbf", "Force")
# compare_units(lbf, newton, "lbf", "N", "Force")
# compare_units(kN, kip, "kN", "kip", "Force")
# compare_units(kip, kN, "kip", "kN", "Force")

# print("\n===== Stress/Pressure Units =====")
# compare_units(pascal, psi, "Pa", "psi", "Stress")
# compare_units(psi, pascal, "psi", "Pa", "Stress")
# compare_units(MPa, ksi, "MPa", "ksi", "Stress")
# compare_units(ksi, MPa, "ksi", "MPa", "Stress")

# print("\n===== Moment/Torque Units =====")
# compare_units(Nm, lbf*foot, "Nm", "lbf·ft", "Moment")
# compare_units(kipft, kNm, "kip·ft", "kN·m", "Moment")
# compare_units(kipin, Nm, "kip·in", "Nm", "Moment")

# print("\n===== Summary of Key Conversions =====")
# print("1 meter = 3.28084 feet")
# print("1 kg = 2.20462 lb")
# print("1 N = 0.224809 lbf")
# print("1 kN = 0.224809 kip")
# print("1 MPa = 0.145038 ksi")
# print("1 Nm = 0.737562 lbf·ft")
# print("1 kip·ft = 1.35582 kN·m")
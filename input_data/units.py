# =============================================
# Unit Definitions (Consistent FPS System with Kip as Base Force Unit)
# =============================================

# FPS Base Units - Modified to use kip instead of pound_force
kip = 1.0                     # Force (k)
foot = 1.0                    # Length (ft)
second = 1.0                  # Time (s)
g = 32.174 * foot/(second**2) # Gravitational acceleration (ft/s²)

# Derived units
# pound_force = kip / 1000.0    # Force (lbf)


# FPS Derived Units
inch = foot / 12.0
ksi = kip / (inch**2)
ksf = kip / (foot**2)
kcf = kip / (foot**3)

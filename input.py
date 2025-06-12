# =============================================
# Unit Definitions (Consistent FPS System with Kip as Base Force Unit)
# =============================================

# FPS Base Units - Modified to use kip instead of pound_force
kip = 1.0                     # Force (k)
foot = 1.0                    # Length (ft)
second = 1.0                  # Time (s)
g = 32.174 * foot/(second**2) # Gravitational acceleration (ft/s²)

# Derived units
pound_force = kip / 1000.0    # Force (lbf)
pound_mass = pound_force / g  # Mass (lbm = lbf·s²/ft)
slug = pound_force * second**2 / foot  # Consistent mass unit in FPS
kip_mass = kip / g            # Mass equivalent of 1 kip (k·s²/ft)

# FPS Derived Units
inch = foot / 12.0
ksi = kip / (inch**2)
psi = pound_force / (inch**2)
psf = pound_force / (foot**2)
ksf = kip / (foot**2)
kipin = kip * inch
kipft = kip * foot

# =============================================
# Module Imports and Initial Setup
# =============================================
import numpy as np        # Numerical computing library

# =============================================
# Material Property Definitions (FPS Units with Kip)
# =============================================
# Concrete properties:
f_c_1 = -3.6 * ksi              # Compressive strength of unconfined concrete (-3.6 ksi = -3600 psi)
f_c_2 = -4.0 * ksi              # Compressive strength of confined concrete (-4.0 ksi = -4000 psi)
eps_c = -0.002                  # Strain at peak compressive stress (dimensionless)
eps_u = -0.02                   # Ultimate compressive strain (dimensionless)

# Steel reinforcement properties:
f_y = 60.0 * ksi                # Yield strength of reinforcing steel (60 ksi)
E_s = 29000.0 * ksi             # Elastic modulus of steel (29000 ksi)
rebar_dia = 1.0 * inch          # #8 bar diameter (1 inch)
rebar = 0.25 * np.pi * (rebar_dia)**2  # Cross-sectional area of rebar

# =============================================
# Material (BEAM COLUMN MEMBER MATERIALS)
# =============================================
materials = [
    # Unconfined concrete (Concrete02 model)
    {
        'ID': 'Concrete02',
        'matTag': 1,            # Material tag for reference
        'fpc': f_c_1,           # Compressive strength
        'epsc0': eps_c,         # Strain at peak strength (dimensionless)
        'fpcu': 0.2 * f_c_1,    # Crushing strength
        'epsU': eps_u,          # Strain at crushing strength (dimensionless)
        'lamda': 0.1,           # Ratio between unloading/loading stiffness (dimensionless)
        'ft': -0.1 * f_c_1,     # Tensile strength (10% of compressive)
        'Ets': (-0.1 * f_c_1) / 0.002  # Tension softening stiffness
    },
    
    # Confined concrete (Concrete02 model)
    {
        'ID': 'Concrete02',
        'matTag': 2,            # Different tag for confined concrete
        'fpc': f_c_2,           # Higher compressive strength
        'epsc0': eps_c,         # Dimensionless
        'fpcu': 0.2 * f_c_2,
        'epsU': eps_u,          # Dimensionless
        'lamda': 0.1,           # Dimensionless
        'ft': -0.1 * f_c_2,
        'Ets': (-0.1 * f_c_2) / 0.002
    },
    
    # Reinforcing steel (Steel02 model - Giuffre-Menegotto-Pinto)
    {
        'ID': 'Steel02',
        'matTag': 3,
        'Fy': f_y,              # Yield strength
        'E0': E_s,              # Initial elastic modulus
        'b': 0.005,             # Strain-hardening ratio (dimensionless)
        'R0': 20.0,             # Initial transition radius (dimensionless)
        'cR1': 0.925,           # Curve fitting parameter (dimensionless)
        'cR2': 0.15             # Curve fitting parameter (dimensionless)
    }
]

# =============================================
# Section Definitions (BEAM COLUMN MEMBER SECTIONS)
# =============================================
# section_definitions = {
#     # Rectangular RC sections
#     "rectangular_sections": {
#         "section1": {  # Example: 16×12 inch beam
#             "type": "rectangular",
#             "section_tag": 1,          # Unique section identifier
#             "H": 16.0 * inch,          # Section height
#             "B": 12.0 * inch,          # Section width
#             "cover_H": 1.5 * inch,     # Cover thickness (vertical)
#             "cover_B": 1.5 * inch,     # Cover thickness (horizontal)
#             "core_tag": 2,             # Confined concrete material
#             "cover_tag": 1,            # Unconfined concrete material
#             "steel_tag": 3,            # Steel material
#             # Main reinforcement (top and bottom)
#             "n_bars_top": 3,           # Number of top bars
#             "dia_top": 1.0 * inch,     # Diameter of top bars
#             "n_bars_bot": 3,           # Number of bottom bars
#             "dia_bot": 1.0 * inch,     # Diameter of bottom bars
#             # Secondary reinforcement (sides)
#             "n_bars_secondary_top": 2,
#             "dia_sec_top": 0.625 * inch,  # 5/8 inch (#5 bar)
#             "n_bars_secondary_bot": 2,
#             "dia_sec_bot": 0.625 * inch,  # 5/8 inch (#5 bar)
#             # Internal ties/stirrups
#             "n_bars_int": 4,
#             "dia_int": 0.5 * inch,     # 1/2 inch (#4 bar)
#             "offset": 2.0 * inch       # Offset for secondary bars from edge
#         },
#         "section2": {  # Larger beam: 24×12 inch
#             "type": "rectangular",
#             "section_tag": 2,
#             "H": 24.0 * inch,
#             "B": 12.0 * inch,
#             "cover_H": 1.5 * inch,
#             "cover_B": 1.5 * inch,
#             "core_tag": 2,
#             "cover_tag": 1,
#             "steel_tag": 3,
#             "n_bars_top": 4,
#             "dia_top": 1.0 * inch,
#             "n_bars_bot": 4,
#             "dia_bot": 1.0 * inch,
#             "n_bars_secondary_top": 2,
#             "dia_sec_top": 0.625 * inch,
#             "n_bars_secondary_bot": 2,
#             "dia_sec_bot": 0.625 * inch,
#             "n_bars_int": 6,
#             "dia_int": 0.5 * inch,
#             "offset": 2.0 * inch
#         }
#     },
    
#     # Circular RC sections (columns)
#     "circular_sections": {
#         "section3": {  # 20-inch diameter solid column
#             "type": "circular",
#             "section_tag": 3,
#             "D_Sec": 20.0 * inch,      # Diameter
#             "cover_Sec": 1.5 * inch,   # Cover thickness
#             "num_Bars_Sec": 8,         # Number of longitudinal bars
#             "bar_dia_Sec": 1.0 * inch, # Bar diameter
#             "core_tag": 2,             # Confined concrete
#             "cover_tag": 1,            # Unconfined concrete
#             "steel_tag": 3,            # Steel material
#             "ri": 0.0 * inch,          # Inner radius (0 for solid)
#             # Fiber discretization parameters:
#             "nf_Core_R": 8,            # Radial fibers in core (dimensionless)
#             "nf_Core_T": 8,            # Theta fibers in core (dimensionless)
#             "nf_Cover_R": 4,           # Radial fibers in cover (dimensionless)
#             "nf_Cover_T": 8            # Theta fibers in cover (dimensionless)
#         },
#         "section4": {  # 32-inch diameter hollow column
#             "type": "circular",
#             "section_tag": 4,
#             "D_Sec": 32.0 * inch,
#             "cover_Sec": 2.0 * inch,
#             "num_Bars_Sec": 12,
#             "bar_dia_Sec": 1.25 * inch,  # 1-1/4 inch (#10 bar)
#             "core_tag": 2,
#             "cover_tag": 1,
#             "steel_tag": 3,
#             "ri": 8.0 * inch,          # Inner radius (hollow section)
#             "nf_Core_R": 8,            # Dimensionless
#             "nf_Core_T": 12,           # Dimensionless
#             "nf_Cover_R": 4,           # Dimensionless
#             "nf_Cover_T": 12           # Dimensionless
#         }
#     }
# }
section_definitions = {
    # Rectangular RC sections
    "rectangular_sections": {
        "col1": {  # Example: 16×12 inch beam
            "type": "rectangular",
            "section_tag": 1,          # Unique section identifier
            "H": 10.0 * inch,          # Section height
            "B": 10.0 * inch,          # Section width
            "cover_H": 1.5 * inch,     # Cover thickness (vertical)
            "cover_B": 1.5 * inch,     # Cover thickness (horizontal)
            "core_tag": 2,             # Confined concrete material
            "cover_tag": 1,            # Unconfined concrete material
            "steel_tag": 3,            # Steel material
            # Main reinforcement (top and bottom)
            "n_bars_top": 2,           # Number of top bars
            "dia_top": 0.625 * inch,     # Diameter of top bars
            "n_bars_bot": 2,           # Number of bottom bars
            "dia_bot": 0.625 * inch,     # Diameter of bottom bars
            # Secondary reinforcement (sides)
            "n_bars_secondary_top": 0,
            "dia_sec_top": 0.625 * inch,  # 5/8 inch (#5 bar)
            "n_bars_secondary_bot": 0,
            "dia_sec_bot": 0.625 * inch,  # 5/8 inch (#5 bar)
            # Internal ties/stirrups
            "n_bars_int": 4,
            "dia_int": 0.5 * inch,     # 1/2 inch (#4 bar)
            "offset": 2.0 * inch,      # Offset for secondary bars from edge
            "area": 192.0 * inch * inch,  # Cross-sectional area (H × B)
            "unit_weight": 0.150,       # Unit weight in pcf (pounds per cubic foot)
            "rotation": 0.0,
        },
        "beam1": {  # Larger beam: 24×12 inch
            "type": "rectangular",
            "section_tag": 2,
            "H": 15.0 * inch,
            "B": 10.0 * inch,
            "cover_H": 1.5 * inch,
            "cover_B": 1.5 * inch,
            "core_tag": 2,
            "cover_tag": 1,
            "steel_tag": 3,
            # Main reinforcement (top and bottom)
            "n_bars_top": 2,           # Number of top bars
            "dia_top": 0.625 * inch,     # Diameter of top bars
            "n_bars_bot": 2,           # Number of bottom bars
            "dia_bot": 0.625 * inch,     # Diameter of bottom bars
            # Secondary reinforcement (sides)
            "n_bars_secondary_top": 0,
            "dia_sec_top": 0.625 * inch,  # 5/8 inch (#5 bar)
            "n_bars_secondary_bot": 0,
            "dia_sec_bot": 0.625 * inch,  # 5/8 inch (#5 bar)
            # Internal ties/stirrups
            "n_bars_int": 4,
            "dia_int": 0.5 * inch,     # 1/2 inch (#4 bar)
            "offset": 2.0 * inch,      # Offset for secondary bars from edge
            "area": 192.0 * inch * inch,  # Cross-sectional area (H × B)
            "unit_weight": 0.150,       # Unit weight in pcf (pounds per cubic foot)
            "rotation": 0.0,
        }
    },
    
    # Circular RC sections (columns)
    "circular_sections": {
        "col2": {  # 20-inch diameter solid column
            "type": "circular",
            "section_tag": 3,
            "D_Sec": 10.0 * inch,      # Diameter
            "cover_Sec": 1.5 * inch,   # Cover thickness
            "num_Bars_Sec": 8,         # Number of longitudinal bars
            "bar_dia_Sec": 1.0 * inch, # Bar diameter
            "core_tag": 2,             # Confined concrete
            "cover_tag": 1,            # Unconfined concrete
            "steel_tag": 3,            # Steel material
            "ri": 0.0 * inch,          # Inner radius (0 for solid)
            # Fiber discretization parameters:
            "nf_Core_R": 8,            # Radial fibers in core (dimensionless)
            "nf_Core_T": 8,            # Theta fibers in core (dimensionless)
            "nf_Cover_R": 4,           # Radial fibers in cover (dimensionless)
            "nf_Cover_T": 8,           # Theta fibers in cover (dimensionless)
            "area": 3.14159 * (10.0/2 * inch) * (10.0/2 * inch),  # π × r²
            "unit_weight": 0.150,       # Unit weight in pcf (pounds per cubic foot)
            "rotation": 0.0,
        },
        "section4": {  # 32-inch diameter hollow column
            "type": "circular",
            "section_tag": 4,
            "D_Sec": 12.0 * inch,
            "cover_Sec": 2.0 * inch,
            "num_Bars_Sec": 12,
            "bar_dia_Sec": 1.25 * inch,  # 1-1/4 inch (#10 bar)
            "core_tag": 2,
            "cover_tag": 1,
            "steel_tag": 3,
            "ri": 8.0 * inch,          # Inner radius (hollow section)
            "nf_Core_R": 8,            # Dimensionless
            "nf_Core_T": 12,           # Dimensionless
            "nf_Cover_R": 4,           # Dimensionless
            "nf_Cover_T": 12,          # Dimensionless
            "area": 3.14159 * ((32.0/2 * inch) * (32.0/2 * inch) - (8.0 * inch) * (8.0 * inch)),  # π × (R² - r²)
            "unit_weight": 0.150,       # Unit weight in pcf (pounds per cubic foot)
            "rotation": 0.0,
        }
    }
}

member_section_mapping = {
    'col1': 
    ["cz4","cz7","cz13","cz16","cz19","cz22","cz25","cz28","cz31","cz34","cz40","cz43",
             
    "cz5","cz8","cz14","cz17","cz20","cz23","cz26","cz29","cz32","cz35","cz41","cz44",
    
    "cz6","cz9","cz15","cz18","cz21","cz24","cz27","cz30","cz33","cz36","cz42", "cz45"
    ],

    'col2': 
    ["cz1", "cz10", "cz37", "cz46",
             
    "cz2", "cz11", "cz38", "cz47",
    
    "cz3", "cz12", "cz39", "cz48"
    ],

    'beam1':
    [

    "bx61","bx62","bx63","bx64","bx65","bx66","bx67","bx68","bx69","bx70","bx71","bx72"
    ,
    "bx73","bx74","bx75","bx76","bx77","bx78","bx79","bx80","bx81","bx82","bx83","bx84"
    ,
    "bx85","bx86","bx87","bx88","bx89","bx90","bx91","bx92","bx93","bx94","bx95","bx96",
    

    "by109","by110","by111","by112","by113","by114","by115","by116","by117","by118","by119","by120",

    "by121","by122","by123","by124","by125","by126","by127","by128","by129","by130","by131","by132",

    "by133","by134","by135","by136","by137","by138","by139","by140","by141","by142","by143","by144"
    
]
}

# # Parameters (in mm)
# params = {
#     'sec_tag': 1,
#     'core_tag': 1,
#     'cover_tag': 2,
#     'steel_tag': 3,
#     'H': 800,
#     'B': 800,
#     't': 200,
#     'cover_H': 40,
#     'cover_B': 40,
#     'n_bars_vertical_outer': 5,
#     'dia_vertical_outer': 16,
#     'n_bars_horizontal_outer': 6,
#     'dia_horizontal_outer': 16,
#     'n_bars_vertical_inner': 4,
#     'dia_vertical_inner': 12,
#     'n_bars_horizontal_inner': 5,
#     'dia_horizontal_inner': 12,
#     'corner_bar_dia': 16,
#     'IMAGE_FOLDER': "output"
# }

# # Generate section
# draw_L_rc_section10(**params)

# =============================================
# Structural Grid Definition
# =============================================
# Node spacing in X, Y, and Z directions
x_spacing = [10.0 * foot, 12.0 * foot, 11.0 * foot]     # X-direction spacing
y_spacing = [10.0 * foot, 12.0 * foot, 11.0 * foot]     # Y-direction spacing 
z_spacing = [10.0 * foot, 12.0 * foot, 11.0 * foot]     # Story height

# =============================================
# Node and Member Management
# =============================================
# Define any custom nodes/members to add (empty dictionary for now)
create_new_nodes = {}
create_new_members = {}

# Define nodes/members to remove from model
delete_nodes = []
delete_members = ['bx49', 'bx50', 'bx51', 'bx52', 'bx53', 'bx54', 'bx55', 'bx56', 'bx57', 'bx58', 'bx59', 'bx60', 'by97', 'by98', 'by99', 'by100', 'by101', 'by102', 'by103', 'by104', 'by105', 'by106', 'by107', 'by108']  # Members to be removed

# =============================================
# FIXITY DATA
# =============================================
fixity_data = {
    "base_nodes1": {"nodes": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], "vals": [1, 1, 1, 1, 1, 1]},
    # "base_nodes2": {"nodes": [], "vals": [0, 1, 1, 0, 1, 0]},
    # "base_nodes3": {"nodes": [], "vals": [1, 0, 0, 1, 0, 1]},
}

# =============================================
# Material (SHELL MATERIALS)
# =============================================
materials_config = {
    "materials": [
        {
            "name": "soil_spring",
            "type": "ENT",
            "id": 100000,
            "properties": {
                "E": 3.0 * ksi  # Elastic modulus for soil spring (3000 psi = 3.0 ksi)
            }
        },
        {
            "name": "concrete_material",
            "type": "ElasticIsotropic",
            "id": 100001,
            "properties": {
                "E": 3.0 * ksi,           # Elastic modulus for concrete (3000 psi = 3.0 ksi)
                "nu": 0.2,                # Poisson's ratio (dimensionless)
                "rho": 0.15 * (kip / foot**3) / g  # Density (150 pcf = 0.15 kcf converted to consistent mass units)
            }
        }
    ]
}
# K = 1e8          # Modulus of subgrade (soil spring stiffness) reaction (Pa/m)
# ops.uniaxialMaterial("ENT", 1, K)            # No Tension Soil
# =============================================
# Section Definitions (Shell Elements)
# =============================================
sections_config = {
    "sections": [
        {
            "name": "foundation_section",
            "type": "PlateFiber",
            "id": 100000,
            "material_id": 100001,
            "properties": {
                "thickness": 20.0 * inch  # 20-inch thick foundation
            }
        },
        {
            "name": "floor_slab1",
            "type": "PlateFiber",
            "id": 100001,
            "material_id": 100001,
            "properties": {
                "thickness": 6.0 * inch  # 12-inch thick floor slab
            }
        },
        # {
        #     "name": "floor_slab2",
        #     "type": "PlateFiber",
        #     "id": 100002,
        #     "material_id": 100001,
        #     "properties": {
        #         "thickness": 6.0 * inch  # 12-inch thick floor slab
        #     }
        # },
        # {
        #     "name": "floor_slab3",
        #     "type": "PlateFiber",
        #     "id": 100003,
        #     "material_id": 100001,
        #     "properties": {
        #         "thickness": 6.0 * inch  # 12-inch thick floor slab
        #     }
        # },
    ]
}

# =============================================
# Surface Mesh Configuration
# =============================================
# Define shell element surfaces
surface_configurations = {
    "first_floor_slab1": {
        "points": [ "n17", "n18", "n22", "n21"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "first_floor_slab2": {
        "points": [ "n18", "n19", "n23", "n22"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "first_floor_slab3": {
        "points": [ "n19", "n20", "n24", "n23"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "first_floor_slab4": {
        "points": [ "n21", "n22", "n26", "n25"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "first_floor_slab5": {
        "points": [ "n22", "n23", "n27", "n26"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "first_floor_slab6": {
        "points": [ "n25", "n26", "n30", "n29"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },

    "first_floor_slab7": {
        "points": [ "n26", "n27", "n31", "n30"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "first_floor_slab8": {
        "points": [ "n27", "n28", "n32", "n31"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },

    "second_floor_slab1": {
        "points": [ "n33", "n34", "n38", "n37"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "second_floor_slab2": {
        "points": [ "n34", "n35", "n39", "n38"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "second_floor_slab3": {
        "points": [ "n35", "n36", "n40", "n39"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "second_floor_slab4": {
        "points": [ "n37", "n38", "n42", "n41"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "second_floor_slab5": {
        "points": [ "n38", "n39", "n43", "n42"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "second_floor_slab6": {
        "points": [ "n41", "n42", "n46", "n45"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },

    "second_floor_slab7": {
        "points": [ "n42", "n43", "n47", "n46"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "second_floor_slab8": {
        "points": [ "n43", "n44", "n48", "n47"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },

    
"third_floor_slab1": {
        "points": [ "n49", "n50", "n54", "n53"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "third_floor_slab2": {
        "points": [ "n50", "n51", "n55", "n54"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "third_floor_slab3": {
        "points": [ "n51", "n52", "n56", "n55"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "third_floor_slab4": {
        "points": [ "n53", "n54", "n58", "n57"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "third_floor_slab5": {
        "points": [ "n54", "n55", "n59", "n58"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "third_floor_slab6": {
        "points": [ "n57", "n58", "n62", "n61"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },

    "third_floor_slab7": {
        "points": [ "n58", "n59", "n63", "n62"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    "third_floor_slab8": {
        "points": [ "n59", "n60", "n64", "n63"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
    
    "shear_wall1": {
        "points": [ "n1", "n2", "n18", "n17"],
        "section_name": "floor_slab1",  # Reference to section in sections_config
        "add_shell": {},     # Can specify predefined shell elements here
        "remove_shell": [],  # Shells to remove
        "predefined_points": {},
        "num_x_div": 3,     # Number of mesh divisions in X (dimensionless)
        "num_y_div": 3,     # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL', 'self_weight'],
        "pressures": [-0.02 * ksf, -0.045 * ksf, -0.075 * ksf],  # Dead load (20 psf = 0.02 ksf) and live load (40 psf = 0.04 ksf)
        "thickness": 6*inch
    },
}

# surface_configurations = {
#     "horizontal": {
#         "points": [
#             [0, 0, 0],
#             [2, 0, 0],
#             [2.5, 1.5, 0],
#             [1.5, 2.5, 0],
#             [0.5, 2, 0],
#             [-0.5, 1, 0]
#         ],
#         "section_name": "floor_slab",
#         "add_shell": {
#             # "Shell1": ["P1", "P2", "P3", "P4"]
#         },
#         "remove_shell": [],
#         "predefined_points": {
#             # 'P1': np.array([0, 0, 0]),
#             # 'P2': np.array([2, 0, 1.5])
#         },
#         "num_x_div": 4,
#         "num_y_div": 4
#     },
#     "vertical": {
#         "points": [
#             [0, 0, 0],    # Bottom-front
#             [3, 0, 0],    # Bottom-back
#             [3, 0, 2],    # Top-back
#             [0, 0, 2]     # Top-front
#         ],
#         "section_name": "exterior_wall",
#         "add_shell": {
#             # "Shell2": ["P5", "P6", "P7"]
#         },
#         "remove_shell": ["Shell4"],
#         "predefined_points": {
#             # 'P3': np.array([2, 2, 1.5]),
#             # 'P4': np.array([0, 2, 0])
#         },
#         "num_x_div": 4,
#         "num_y_div": 4
#     },
#     "inclined": {
#         "points": [
#             [0, 0, 0],    # Base point 1
#             [3, 0, 0],    # Base point 2
#             [2, 2, 2],    # Top point 1
#             [0, 2, 2]     # Top point 2
#         ],
#         "section_name": "roof_slab",
#         "add_shell": {
#             # "Shell3": ["P8", "P9", "P10", "P11"]
#         },
#         "remove_shell": [],
#         "predefined_points": {
#             # 'P5': np.array([1, 1, 1]),
#             # 'P6': np.array([1.5, 0.5, 1.2])
#         },
#         "num_x_div": 6,  # More divisions for better approximation of slope
#         "num_y_div": 6
#     },
#     "complex_inclined": {
#         "points": [
#             [0, 0, 0],     # Base point
#             [4, 0, 1],     # Right point (slightly elevated)
#             [3, 3, 3],     # Top point
#             [1, 3, 2],     # Left point
#             [0, 2, 1.5]    # Front point
#         ],
#         "section_name": "curved_wall",
#         "add_shell": {
#             # "Shell4": ["P1", "P2", "P3"],
#             # "Shell5": ["P4", "P5", "P6"]
#         },
#         "remove_shell": ["Shell5"],
#         "predefined_points": {
#             # 'P7': np.array([0.5, 1.5, 0.8]),
#             # 'P8': np.array([0.5, 0, 0.5])
#         },
#         "num_x_div": 8,  # More divisions for complex geometry
#         "num_y_div": 8
#     }
# }

# =============================================
# 13. Exampples
# =============================================
# Define load cases with mass information
load_cases={
  "load_cases": [
    {"name": "DL","type": "gravity","factor": 1.0,"direction": None },
    {"name": "LL","type": "gravity","factor": 1.0,"direction": None},
    {"name": "self_weight","type": "gravity","factor": 1.0,"direction": None},
    # {"name": "mass","type": "gravity","factor": 1.0,"direction": None},
    # {"name": "WindX","type": "lateral","factor": 1.0,"direction": "X"},
    # {"name": "QuakeY","type": "lateral","factor": 1.0,"direction": "Y"}
  ]
}

nodal_load_entries = [
    ["n10002", "DL", 1.2, -0.5, -1.15, 0.3, -0.4, 0.1],
    ["n10003", "DL", -2.1, 0.7, -2.3, -0.2, 0.5, -0.3],
    ["n10004", "DL", 0.8, 1.0, -1.15, 0.0, -0.6, 0.2],
    ["n6", "DL", -1.4, 0.9, -0.575, 0.1, 0.3, -0.1],
    ["n5", "DL", 1.1, -0.6, -0.575, -0.5, 0.2, 0.3],
    ["n10006", "DL", -0.9, 1.2, -1.15, 0.4, -0.1, -0.2],
    ["n8", "DL", -1.7, 0.3, -0.575, -0.2, -0.3, 0.1],
    ["n7", "DL", 0.5, -0.4, -0.575, 0.2, 0.1, -0.1],
    ["n10002", "LL", -1.3, 0.6, -2.3, 0.0, -0.4, 0.2],
    ["n10003", "LL", 0.9, -0.8, -4.6, -0.1, 0.5, -0.2],
    ["n10004", "LL", 1.5, 0.4, -2.3, 0.3, -0.3, 0.1],
    ["n6", "LL", -1.0, 1.1, -1.15, -0.4, 0.2, 0.0],
    ["n5", "LL", 0.6, -0.2, -1.15, 0.1, -0.5, 0.4],
    ["n10006", "LL", -0.7, 0.8, -2.3, 0.2, 0.3, -0.1],
    ["n8", "LL", 0.2, 0.5, -1.15, 0.0, 0.0, -0.2],
    ["n7", "LL", -0.6, -0.3, -1.15, 0.4, -0.1, 0.2]
]

nodal_load_entries1 = [
    ["n10007", "DL", 0.7, -0.3, -1.15, 0.2, -0.5, 0.0],
    ["n10008", "DL", -1.6, 0.8, -0.575, -0.3, 0.1, -0.2],
    ["n10007", "LL", 1.0, -0.5, -2.3, 0.1, -0.2, 0.3],
    ["n10008", "LL", -0.8, 0.6, -1.15, -0.2, 0.4, -0.1]
]

nodal_load_entries2 = [
    ["n10009", "DL", -1.2, 1.0, -0.575, 0.0, 0.2, 0.1],
    ["n10010", "DL", 0.6, -0.9, -2.3, 0.3, -0.1, -0.3],
    ["n10009", "LL", -1.1, 0.4, -1.15, -0.1, 0.0, 0.2],
    ["n10010", "LL", 1.3, -0.7, -2.3, 0.5, -0.3, -0.2]
]

nodal_load_entries3 = [
    ["n10011", "DL", 0.9, -0.2, -0.575, -0.4, 0.3, -0.1],
    ["n10012", "DL", -1.4, 1.1, -1.15, 0.2, -0.2, 0.0],
    ["n10011", "LL", 0.4, 0.7, -2.3, -0.3, 0.1, -0.2],
    ["n10012", "LL", -0.5, -0.6, -1.15, 0.0, 0.5, 0.1]
]

nodal_load_entries4 = [
    ["n10013", "DL", -0.8, 0.5, -2.3, 0.1, 0.0, -0.4],
    ["n10014", "DL", 1.2, -0.4, -0.575, -0.2, 0.3, 0.2],
    ["n10013", "LL", 0.6, 0.3, -1.15, 0.4, -0.1, -0.3],
    ["n10014", "LL", -1.0, -0.9, -2.3, -0.1, 0.2, 0.0]
]

# Combine all load entries into a single list
# all_nodal_load_entries = nodal_load_entries + nodal_load_entries1 + nodal_load_entries2 + nodal_load_entries3 + nodal_load_entries4

all_nodal_load_entries = []

    
# Process and save the load cases

# Element load definitions (unchanged from your input)
element_loads1 = {
    "[5,6,7,8]": [
        {
            "LoadCase": ["DL"],    
            "uniform": {"x": -0.0, "y": -0.0, "z": -10.0},
            "point": {"x": -0.0, "y": -0.0, "z": -0.0, "location": 0.5},
            "temperature_points": [
                {"temp": 0.0, "y": 0.0},
                {"temp": 0.0, "y": 0.0},
                {"temp": 0.0, "y": 0.0}
            ]
        }
    ]
}

element_loads2 = {
    "[2,4,7,8]": [
        {
            "LoadCase": ["LL"],    
            "uniform": {"x": -0.0, "y": -0.0, "z": -10.0},
            "point": {"x": -0.0, "y": -0.0, "z": -10.0, "location": 0.5},
            "temperature_points": [
                {"temp": 0.0, "y": 0.0},
                {"temp": 0.0, "y": 0.0},
                {"temp": 0.0, "y": 0.0}
            ]
        }
    ]
}

element_loads2 = {
    "[15,26,7,8,1,22]": [
        {
            "LoadCase": ["LL"],    
            "uniform": {"x": -0.0, "y": -0.0, "z": -10.0},
            "point": {"x": -0.0, "y": -0.0, "z": -10.0, "location": 0.5},
            "temperature_points": [
                {"temp": 0.0, "y": 0.0},
                {"temp": 0.0, "y": 0.0},
                {"temp": 0.0, "y": 0.0}
            ]
        }
    ]
}

all_element_loads = [element_loads1, element_loads2, element_loads2]

simple_list = [
    [5, ["DL"], [0.0, 0.0, -10.0], [0.0, 0.0, 0.0, 0.5], [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]],
    [6, ["DL"], [0.0, 0.0, -10.0], [0.0, 0.0, 0.0, 0.5], [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]],
    [7, ["DL"], [0.0, 0.0, -10.0], [0.0, 0.0, 0.0, 0.5], [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]],
    [8, ["DL"], [0.0, 0.0, -10.0], [0.0, 0.0, 0.0, 0.5], [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]]
]



all_element_loads = [element_loads1, element_loads2]

element_loads1 = {
    "['bx62', 'bx63', 'by109', 'by110']": [
        {
            "LoadCase": ["DL"],    
            "uniform": {"x": 0.0, "y": 0.0, "z": -10.0},
            "point": {"x": 0.0, "y": 0.0, "z": 0.0, "location": 0.5},
            "temperature_points": [
                {"temp": 0.0, "y": 0.0},
                {"temp": 0.0, "y": 0.0},
                {"temp": 0.0, "y": 0.0}
            ]
        }
    ]
}

element_loads2 = {
    "['cz2', 'bx61', 'by109', 'by110']": [
        {
            "LoadCase": ["LL"],    
            "uniform": {"x": 0.0, "y": 0.0, "z": -10.0},
            "point": {"x": 0.0, "y": 0.0, "z": -10.0, "location": 0.5},
            "temperature_points": [
                {"temp": 0.0, "y": 0.0},
                {"temp": 0.0, "y": 0.0},
                {"temp": 0.0, "y": 0.0}
            ]
        }
    ]
}

element_loads3 = {
    "['by117', 'by126', 'by109', 'by110', 'cz1', 'by122']": [
        {
            "LoadCase": ["LL"],    
            "uniform": {"x": 0.0, "y": 0.0, "z": -10.0},
            "point": {"x": 0.0, "y": 0.0, "z": -10.0, "location": 0.5},
            "temperature_points": [
                {"temp": 0.0, "y": 0.0},
                {"temp": 0.0, "y": 0.0},
                {"temp": 0.0, "y": 0.0}
            ]
        }
    ]
}

all_element_loads = [element_loads1, element_loads2, element_loads3]

load_combinations = {
    "mass": [("DL", 1.0), ("LL", 0.5), ("self_weight", 0.0)],
    "Comb2": [("DL", 1.2), ("LL", 1.6), ("self_weight", 1.2)],
    "Comb1": [("DL", 1.4)]
}

# Example 1: Apply nodal loads
# applied_nodal_loads = process_structure_loads(
#     JSON_FOLDER, 
#     operation_type="nodal_loads",
#     load_combinations=load_combinations,
#     numbering=1,
#     load_combination="Comb2"
# )
# print("Applied nodal loads:", applied_nodal_loads)

# Example 2: Apply member loads (uniform and point loads)
# member_loads = [
#     # Uniform load
#     {
#         'element_names': ["main_beam1", "main_beam2"],
#         'load_type': 'beamUniform',
#         'load_case': 'DL',
#         'wy': -10.0,  # Vertical load (negative = downward)
#         'wz': -5.0   # Lateral load
#     },
#     # Point load
#     {
#         'element_names': ["cantilever_beam"],
#         'load_type': 'beamPoint',
#         'Py': -15.0,  # Vertical point load
#         'x': 0.7       # 70% along the beam length
#     }
# ]
# process_structure_loads(
#     JSON_FOLDER,
#     operation_type="member_loads",
#     load_pattern_tag=1,
#     loads=member_loads,
#     load_combinations=load_combinations,
#     load_combination="Comb2"
# )

# Example 3: Apply masses using node names
# process_structure_loads(
#     JSON_FOLDER,
#     operation_type="masses",
#     loaded_nodes=["base_node1", "base_node2", "roof_node1", "roof_node2"],
#     m_1=1500.0,  # Total mass (kg)
#     mass_distribution=[0.4, 0.4, 0.1, 0.1]  # More mass at base nodes
# )

# Example 4: Thermal and mechanical combined loads
# combined_loads = [
#     # Thermal gradient load
#     {
#         'element_names': ["exposed_beam1", "exposed_beam2"],
#         'load_type': 'beamThermal',
#         'temp_points': [20.0, -0.1, 25.0, 0.0, 30.0, 0.1]  # Temp (C) and y-coord pairs
#     },
#     # Simultaneous uniform load
#     {
#         'element_names': ["exposed_beam1", "exposed_beam2"],
#         'load_type': 'beamUniform',
#         'wy': -3.0  # Snow load
#     }
# ]
# process_structure_loads(
#     JSON_FOLDER,
#     operation_type="member_loads",
#     load_pattern_tag=2,  # Different pattern for thermal effects
#     loads=combined_loads
# )

'''
# =============================================
# Unit Definitions (Consistent MKS System with kN as Base Force Unit)
# =============================================

# MKS Base Units - Modified to use kN instead of newton
kN = 1.0                      # Force (kN)
meter = 1.0                   # Length (m)
second = 1.0                  # Time (s)
g = 9.81 * meter/(second**2)  # Gravitational acceleration (m/s²)

# Derived units
newton = kN / 1000.0          # Force (N)
kilogram = newton * second**2 / meter  # Mass (kg = N·s²/m)
ton_mass = 1000.0 * kilogram  # Metric ton (1000 kg)
kN_mass = kN * second**2 / meter  # Mass equivalent of 1 kN (kg)

# MKS Derived Units
millimeter = 1e-3 * meter     # Length
centimeter = 1e-2 * meter     # Length
MN = 1e3 * kN                 # Meganewton (1000 kN)
kPa = kN / (meter**2) / 1000.0 # Kilopascal (kN/m²/1000)
MPa = kN / (meter**2)         # Megapascal (kN/m²)
GPa = 1e3 * MPa               # Gigapascal (1000 MPa)
kNm = kN * meter              # Kilonewton-meter
Nm = newton * meter           # Newton-meter

# =============================================
# Module Imports and Initial Setup
# =============================================
import numpy as np            # Numerical computing library

# =============================================
# Material Property Definitions (MKS Units with kN)
# =============================================
# Concrete properties:
f_c_1 = -25 * MPa             # Compressive strength of unconfined concrete (negative for compression)
f_c_2 = -28 * MPa             # Compressive strength of confined concrete (negative for compression)
eps_c = -0.002                # Strain at peak compressive stress (dimensionless)
eps_u = -0.02                 # Ultimate compressive strain (dimensionless)

# Steel reinforcement properties:
f_y = 420.0 * MPa             # Yield strength of reinforcing steel
E_s = 210.0 * GPa             # Elastic modulus of steel
rebar_dia = 25 * millimeter   # Bar diameter
rebar = 0.25 * np.pi * (rebar_dia)**2  # Cross-sectional area of 25mm diameter rebar

# =============================================
# Material Model Definitions (OpenSees)
# =============================================
materials = [
    # Unconfined concrete (Concrete02 model)
    {
        'ID': 'Concrete02',
        'matTag': 1,          # Material tag for reference
        'fpc': f_c_1,         # Compressive strength
        'epsc0': eps_c,       # Strain at peak strength (dimensionless)
        'fpcu': 0.2 * f_c_1,  # Crushing strength
        'epsU': eps_u,        # Strain at crushing strength (dimensionless)
        'lamda': 0.1,         # Ratio between unloading/loading stiffness (dimensionless)
        'ft': -0.1 * f_c_1,   # Tensile strength (10% of compressive)
        'Ets': (-0.1 * f_c_1) / 0.002  # Tension softening stiffness
    },
    
    # Confined concrete (Concrete02 model)
    {
        'ID': 'Concrete02',
        'matTag': 2,          # Different tag for confined concrete
        'fpc': f_c_2,         # Higher compressive strength
        'epsc0': eps_c,       # Dimensionless
        'fpcu': 0.2 * f_c_2,
        'epsU': eps_u,        # Dimensionless
        'lamda': 0.1,         # Dimensionless
        'ft': -0.1 * f_c_2,
        'Ets': (-0.1 * f_c_2) / 0.002
    },
    
    # Reinforcing steel (Steel02 model - Giuffre-Menegotto-Pinto)
    {
        'ID': 'Steel02',
        'matTag': 3,
        'Fy': f_y,            # Yield strength
        'E0': E_s,            # Initial elastic modulus
        'b': 0.005,           # Strain-hardening ratio (dimensionless)
        'R0': 20.0,           # Initial transition radius (dimensionless)
        'cR1': 0.925,         # Curve fitting parameter (dimensionless)
        'cR2': 0.15           # Curve fitting parameter (dimensionless)
    }
]

# =============================================
# Section Definitions (Fiber Sections)
# =============================================
section_definitions = {
    # Rectangular RC sections
    "rectangular_sections": {
        "section1": {  # Example: 400x300mm beam
            "type": "rectangular",
            "section_tag": 1,             # Unique section identifier
            "H": 400 * millimeter,        # Section height
            "B": 300 * millimeter,        # Section width
            "cover_H": 40 * millimeter,   # Cover thickness (vertical)
            "cover_B": 40 * millimeter,   # Cover thickness (horizontal)
            "core_tag": 2,                # Confined concrete material
            "cover_tag": 1,               # Unconfined concrete material
            "steel_tag": 3,               # Steel material
            # Main reinforcement (top and bottom)
            "n_bars_top": 3,              # Number of top bars
            "dia_top": 25 * millimeter,   # Diameter of top bars
            "n_bars_bot": 3,              # Number of bottom bars
            "dia_bot": 25 * millimeter,   # Diameter of bottom bars
            # Secondary reinforcement (sides)
            "n_bars_secondary_top": 2,
            "dia_sec_top": 16 * millimeter,
            "n_bars_secondary_bot": 2,
            "dia_sec_bot": 16 * millimeter,
            # Internal ties/stirrups
            "n_bars_int": 4,
            "dia_int": 12 * millimeter,
            "offset": 50 * millimeter     # Offset for secondary bars from edge
        },
        "section2": {  # Larger beam: 600x300mm
            "type": "rectangular",
            "section_tag": 2,
            "H": 600 * millimeter,
            "B": 300 * millimeter,
            "cover_H": 40 * millimeter,
            "cover_B": 40 * millimeter,
            "core_tag": 2,
            "cover_tag": 1,
            "steel_tag": 3,
            "n_bars_top": 4,
            "dia_top": 25 * millimeter,
            "n_bars_bot": 4,
            "dia_bot": 25 * millimeter,
            "n_bars_secondary_top": 2,
            "dia_sec_top": 16 * millimeter,
            "n_bars_secondary_bot": 2,
            "dia_sec_bot": 16 * millimeter,
            "n_bars_int": 6,
            "dia_int": 12 * millimeter,
            "offset": 50 * millimeter
        }
    },
    
    # Circular RC sections (columns)
    "circular_sections": {
        "section3": {  # 500mm diameter solid column
            "type": "circular",
            "section_tag": 3,
            "D_Sec": 500 * millimeter,     # Diameter
            "cover_Sec": 40 * millimeter,  # Cover thickness
            "num_Bars_Sec": 8,             # Number of longitudinal bars
            "bar_dia_Sec": 25 * millimeter,# Bar diameter
            "core_tag": 2,                 # Confined concrete
            "cover_tag": 1,                # Unconfined concrete
            "steel_tag": 3,                # Steel material
            "ri": 0.0 * millimeter,        # Inner radius (0 for solid)
            # Fiber discretization parameters:
            "nf_Core_R": 8,                # Radial fibers in core (dimensionless)
            "nf_Core_T": 8,                # Theta fibers in core (dimensionless)
            "nf_Cover_R": 4,               # Radial fibers in cover (dimensionless)
            "nf_Cover_T": 8                # Theta fibers in cover (dimensionless)
        },
        "section4": {  # 800mm diameter hollow column
            "type": "circular",
            "section_tag": 4,
            "D_Sec": 800 * millimeter,
            "cover_Sec": 50 * millimeter,
            "num_Bars_Sec": 12,
            "bar_dia_Sec": 32 * millimeter,
            "core_tag": 2,
            "cover_tag": 1,
            "steel_tag": 3,
            "ri": 200 * millimeter,        # Inner radius (hollow section)
            "nf_Core_R": 8,                # Dimensionless
            "nf_Core_T": 12,               # Dimensionless
            "nf_Cover_R": 4,               # Dimensionless
            "nf_Cover_T": 12               # Dimensionless
        }
    }
}

# =============================================
# Structural Grid Definition
# =============================================
# Node spacing in X, Y, and Z directions
x_spacing = [6.0 * meter]      # X-direction spacing
y_spacing = [7.0 * meter]      # Y-direction spacing 
z_spacing = [3.5 * meter]      # Story height

# =============================================
# Node and Member Management
# =============================================
# Define any custom nodes/members to add (empty dictionary for now)
create_new_nodes = {}
create_new_members = {}

# Define nodes/members to remove from model
delete_nodes = []
delete_members = ['bx5', 'bx6', 'by9', 'by10']  # Members to be removed

# =============================================
# FIXITY DATA
# =============================================
fixity_data = {
    "base_nodes1": {"nodes": [1, 2, 3, 4], "vals": [1, 1, 1, 1, 1, 1]},
    # "base_nodes2": {"nodes": [], "vals": [0, 1, 1, 0, 1, 0]},
    # "base_nodes3": {"nodes": [], "vals": [1, 0, 0, 1, 0, 1]},
}

# =============================================
# Material Configurations for Shell Elements
# =============================================
materials_config = {
    "materials": [
        {
            "name": "soil_spring",
            "type": "ENT",
            "id": 100000,
            "properties": {
                "E": 21 * GPa  # Elastic modulus for soil spring
            }
        },
        {
            "name": "concrete_material",
            "type": "ElasticIsotropic",
            "id": 100001,
            "properties": {
                "E": 21 * GPa,                           # Elastic modulus for concrete
                "nu": 0.2,                               # Poisson's ratio (dimensionless)
                "rho": 2.4 * (kN_mass / meter**3)        # Density (2400 kg/m³ = 2.4 t/m³)
            }
        }
    ]
}

# =============================================
# Section Definitions for Shell Elements
# =============================================
sections_config = {
    "sections": [
        {
            "name": "foundation_section",
            "type": "PlateFiber",
            "id": 100000,
            "material_id": 100001,
            "properties": {
                "thickness": 500 * millimeter  # 500mm thick foundation
            }
        },
        {
            "name": "floor_slab",
            "type": "PlateFiber",
            "id": 100001,
            "material_id": 100001,
            "properties": {
                "thickness": 300 * millimeter  # 300mm thick floor slab
            }
        }
    ]
}

# =============================================
# Surface Mesh Configuration
# =============================================
# Define shell element surfaces
surface_configurations = {
    "horizontal": {
        "points": [
            # Define a rectangular floor slab at z=3.5m
            [0.0 * meter, 0.0 * meter, 3.5 * meter],
            [6.0 * meter, 0.0 * meter, 3.5 * meter],
            [6.0 * meter, 7.0 * meter, 3.5 * meter],
            [0.0 * meter, 7.0 * meter, 3.5 * meter],
        ],
        "section_name": "floor_slab",  # Reference to section in sections_config
        "add_shell": {},               # Can specify predefined shell elements here
        "remove_shell": [],            # Shells to remove
        "predefined_points": {
            'n5': np.array([0.0 * meter, 0.0 * meter, 3.5 * meter]),
            'n6': np.array([6.0 * meter, 0.0 * meter, 3.5 * meter]),
            'n7': np.array([0.0 * meter, 7.0 * meter, 3.5 * meter]),
            'n8': np.array([6.0 * meter, 7.0 * meter, 3.5 * meter]),
        },
        "num_x_div": 2,                # Number of mesh divisions in X (dimensionless)
        "num_y_div": 2,                # Number of mesh divisions in Y (dimensionless)
        "load_case_names": ['DL', 'LL'],
        "pressures": [5.0 * kPa, 5.0 * kPa]  # Dead load and live load pressures
    }
}

'''
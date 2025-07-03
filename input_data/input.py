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
# =============================================
# Module Imports and Initial Setup
# =============================================
import numpy as np        # Numerical computing library

# =============================================
# Structural Grid Definition
# =============================================
# Node spacing in X, Y, and Z directions
x_spacing = [10.0 * foot]     # X-direction spacing
y_spacing = [10.0 * foot]     # Y-direction spacing 
z_spacing = [10.0 * foot]     # Story height


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
        'name': 'Concrete1',
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
        'name': 'Concrete2',
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
        'name': 'Steel1',
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
# Material (SHELL MATERIALS)
# =============================================
materials_config = {
    "materials": [
        {
            "name": "soil_spring",
            "type": "ENT",
            "id": 100,
            "properties": {
                "E": 3.0 * ksi  # Elastic modulus for soil spring (3000 psi = 3.0 ksi)
            }
        },
        {
            "name": "slab_concrete_material",
            "type": "ElasticIsotropic",
            "id": 101,
            "properties": {
                "E": 3.0 * ksi,           # Elastic modulus for concrete (3000 psi = 3.0 ksi)
                "nu": 0.2,                # Poisson's ratio (dimensionless)
                "rho": 0.15 * kcf  # Density (150 pcf = 0.15 kcf converted to consistent mass units)
            }
        }
    ]
}

# K = 1e8          # Modulus of subgrade (soil spring stiffness) reaction (Pa/m)
# ops.uniaxialMaterial("ENT", 1, K)            # No Tension Soil

# =============================================
# Section Definitions (BEAM COLUMN MEMBER SECTIONS)
# =============================================


section_definitions = {
    # Rectangular RC Columns
    "rectangular_columns": {
        "col1": {  # 10×10 inch column
            "type": "rectangular",
            "section_tag": 1,
            "H": 10.0 * inch,
            "B": 10.0 * inch,
            "cover_H": 1.5 * inch,
            "cover_B": 1.5 * inch,
            "core_tag": "Concrete2",
            "cover_tag": "Concrete1",
            "steel_tag": "Steel1",
            "n_bars_top": 4,  # Columns typically have equal reinforcement on all sides
            "dia_top": 0.625 * inch,
            "n_bars_bot": 4,
            "dia_bot": 0.625 * inch,
            "n_bars_secondary_top": 2,
            "dia_sec_top": 0.625 * inch,
            "n_bars_secondary_bot": 2,
            "dia_sec_bot": 0.625 * inch,
            "n_bars_int": 4,
            "dia_int": 0.5 * inch,
            "offset": 2.0 * inch,
            "area": 100.0 * inch * inch,
            "unit_weight": 0.150 * kcf,
            "rotation": 0.0,
        },
        "col2": {  # 12×12 inch column
            "type": "rectangular",
            "section_tag": 2,
            "H": 24.0 * inch,
            "B": 10.0 * inch,
            "cover_H": 1.5 * inch,
            "cover_B": 1.5 * inch,
            "core_tag": "Concrete2",
            "cover_tag": "Concrete1",
            "steel_tag": "Steel1",
            "n_bars_top": 4,
            "dia_top": 0.75 * inch,
            "n_bars_bot": 4,
            "dia_bot": 0.75 * inch,
            "n_bars_secondary_top": 2,
            "dia_sec_top": 0.625 * inch,
            "n_bars_secondary_bot": 2,
            "dia_sec_bot": 0.625 * inch,
            "n_bars_int": 4,
            "dia_int": 0.5 * inch,
            "offset": 2.0 * inch,
            "area": 144.0 * inch * inch,
            "unit_weight": 0.150 * kcf,
            "rotation": 0.0,
        }
    },
    
    # Rectangular RC Beams
    "rectangular_beams": {
        "beam1": {  # 10×15 inch beam (original)
            "type": "rectangular",
            "section_tag": 3,
            "H": 15.0 * inch,  # Height (depth)
            "B": 10.0 * inch,   # Width
            "cover_H": 1.5 * inch,
            "cover_B": 1.5 * inch,
            "core_tag": "Concrete2",
            "cover_tag": "Concrete1",
            "steel_tag": "Steel1",
            # Main reinforcement (top and bottom)
            "n_bars_top": 2,
            "dia_top": 0.625 * inch,
            "n_bars_bot": 3,  # Typically more at bottom for positive moment
            "dia_bot": 0.625 * inch,
            # Secondary reinforcement (sides)
            "n_bars_secondary_top": 0,
            "dia_sec_top": 0.625 * inch,
            "n_bars_secondary_bot": 0,
            "dia_sec_bot": 0.625 * inch,
            # Internal ties/stirrups
            "n_bars_int": 2,  # Beams typically have 2-leg stirrups
            "dia_int": 0.5 * inch,
            "offset": 2.0 * inch,
            "area": 150.0 * inch * inch,
            "unit_weight": 0.150 * kcf,
            "rotation": 0.0,
        },
        "beam2": {  # 12×18 inch beam
            "type": "rectangular",
            "section_tag": 4,
            "H": 15.0 * inch,
            "B": 12.0 * inch,
            "cover_H": 1.5 * inch,
            "cover_B": 1.5 * inch,
            "core_tag": "Concrete2",
            "cover_tag": "Concrete1",
            "steel_tag": "Steel1",
            "n_bars_top": 3,
            "dia_top": 0.75 * inch,
            "n_bars_bot": 4,
            "dia_bot": 0.75 * inch,
            "n_bars_secondary_top": 0,
            "dia_sec_top": 0.625 * inch,
            "n_bars_secondary_bot": 0,
            "dia_sec_bot": 0.625 * inch,
            "n_bars_int": 2,
            "dia_int": 0.5 * inch,
            "offset": 2.0 * inch,
            "area": 216.0 * inch * inch,
            "unit_weight": 0.150 * kcf,
            "rotation": 0.0,
        }
    },
    
    # Circular RC Columns
    "circular_columns": {
        "col3": {  # 12-inch diameter solid column
            "type": "circular",
            "section_tag": 5,
            "D_Sec": 10.0 * inch,
            "cover_Sec": 1.5 * inch,
            "num_Bars_Sec": 6,
            "bar_dia_Sec": 0.75 * inch,
            "core_tag": "Concrete2",
            "cover_tag": "Concrete1",
            "steel_tag": "Steel1",
            "ri": 0.0 * inch,
            "nf_Core_R": 8,
            "nf_Core_T": 8,
            "nf_Cover_R": 4,
            "nf_Cover_T": 8,
            "area": 113.1 * inch * inch,
            "unit_weight": 0.150 * kcf,
            "rotation": 0.0,
        },
        "col4": {  # 18-inch diameter solid column
            "type": "circular",
            "section_tag": 6,
            "D_Sec": 12.0 * inch,
            "cover_Sec": 1.5 * inch,
            "num_Bars_Sec": 8,
            "bar_dia_Sec": 0.875 * inch,
            "core_tag": "Concrete2",
            "cover_tag": "Concrete1",
            "steel_tag": "Steel1",
            "ri": 0.0 * inch,
            "nf_Core_R": 10,
            "nf_Core_T": 12,
            "nf_Cover_R": 4,
            "nf_Cover_T": 12,
            "area": 254.5 * inch * inch,
            "unit_weight": 0.150 * kcf,
            "rotation": 0.0,
        }
    },
    "L_columns": {
    "col5": {
        "type": "L",
        "section_tag": 7,
        "core_tag": "Concrete2",
        "cover_tag": "Concrete1",
        "steel_tag": "Steel1",
        "H": 24.0 * inch,
        "B": 24.0 * inch,
        "t": 6.0 * inch,  # leg thickness
        "cover_H": 1.5 * inch,
        "cover_B": 1.5 * inch,
        "n_bars_vertical_outer": 4,
        "dia_vertical_outer": 0.75 * inch,
        "n_bars_horizontal_outer": 4,
        "dia_horizontal_outer": 0.75 * inch,
        "n_bars_vertical_inner": 2,
        "dia_vertical_inner": 0.625 * inch,
        "n_bars_horizontal_inner": 2,
        "dia_horizontal_inner": 0.625 * inch,
        "corner_bar_dia": 0.75 * inch,
        "area": (24.0 * 6.0 + 18.0 * 6.0) * inch * inch,  # L-shape area
        "unit_weight": 0.150 * kcf,
        "rotation": 0.0,
    }
}

}

member_section_mapping = {

    'col1': 
    [
        {
            "name": "cz1",
            "rotation": 0
        },
        {
            "name": "cz2",
            "rotation": 0
        },
        {
            "name": "cz3",
            "rotation": 0
        },
        {
            "name": "cz4",
            "rotation": 0
        }
    ],
    
 

    'beam1':
    [
        {
            "name": "bx7",
            "rotation": 0
        },
        {
            "name": "bx8",
            "rotation": 0
        }
    ],
       
       'beam2': 
[
        {
            "name": "by11",
            "rotation": 0
        },
        {
            "name": "by12",
            "rotation": 0
        }
    ],
    }



# =============================================
# Node and Member Management
# =============================================
# Define any custom nodes/members to add (empty dictionary for now)
create_new_nodes = {}
create_new_members = {}

# Define nodes/members to remove from model
delete_nodes = []
delete_members = [
        "bx5",
        "bx6",
        "by9",
        "by10",
        ]  # Members to be removed

# =============================================
# FIXITY DATA
# =============================================
fixity_data = {
    "base_nodes1": {"nodes": ["n1", "n2", "n3", "n4"], "vals": [1, 1, 1, 1, 1, 1]},
    # "base_nodes2": {"nodes": [], "vals": [0, 1, 1, 0, 1, 0]},
    # "base_nodes3": {"nodes": [], "vals": [1, 0, 0, 1, 0, 1]},
}


# =============================================
# Section Definitions (Shell Elements)
# =============================================
sections_config = {
    "sections": [
        # {
        #     "name": "foundation_section",
        #     "type": "PlateFiber",
        #     "material_id": "soil_spring",
        #     "section_id": 100,
        #     "properties": {
        #         "thickness": 20.0 * inch  # 20-inch thick foundation
        #     }
        # },
        {
            "name": "floor_slab1",
            "type": "PlateFiber",
            "material_id": "slab_concrete_material",
            "section_id": 101,
            "properties": {
                "thickness": 6.0 * inch  # 12-inch thick floor slab
            }
        },
        {
            "name": "floor_slab2",
            "type": "PlateFiber",
            "material_id": "slab_concrete_material",
            "section_id": 102,
            "properties": {
                "thickness": 6.0 * inch  # 12-inch thick floor slab
            }
        },

    ]
}

# =============================================
# Surface Mesh Configuration
# =============================================
# Define shell element surfaces
surface_configurations = {
    
    "10.0_slab1": {
        "points": [
            "n5",
            "n6",
            "n8",
            "n7"
        ],
        "section_name": "floor_slab1",
        "add_shell": {},
        "remove_shell": [],
        "predefined_points": {},
        "num_x_div": 3,
        "num_y_div": 3,
        "load_case_names": [
            "DL",
            "LL",
            "self_weight"
        ],
        "pressures": [
            -0.02 * ksf,
            -0.045 * ksf,
            -0.075 * ksf
        ],
        "thickness": 6
    }
}

# =============================================
# 13. Exampples
# =============================================


nodal_load_entries = [
    # ["n10002", "DL", 1.2, -0.5, -1.15, 0.3, -0.4, 0.1],
    # ["n10003", "DL", -2.1, 0.7, -2.3, -0.2, 0.5, -0.3],

]

nodal_load_entries1 = [
    # ["n10007", "DL", 0.7, -0.3, -1.15, 0.2, -0.5, 0.0],
    # ["n10008", "DL", -1.6, 0.8, -0.575, -0.3, 0.1, -0.2],

]

nodal_load_entries2 = [
    # ["n10009", "DL", -1.2, 1.0, -0.575, 0.0, 0.2, 0.1],
    # ["n10010", "DL", 0.6, -0.9, -2.3, 0.3, -0.1, -0.3],

]

nodal_load_entries3 = [
    # ["n10011", "FL", 0.9, -0.2, -0.575, -0.4, 0.3, -0.1],
    # ["n10012", "DL", -1.4, 1.1, -1.15, 0.2, -0.2, 0.0],

]

nodal_load_entries4 = [
    # ["n10013", "LL", -0.8, 0.5, -2.3, 0.1, 0.0, -0.4],
    # ["n10014", "LL", 1.2, -0.4, -0.575, -0.2, 0.3, 0.2],

]

# Combine all load entries into a single list
# all_nodal_load_entries = nodal_load_entries + nodal_load_entries1 + nodal_load_entries2 + nodal_load_entries3 + nodal_load_entries4

all_nodal_load_entries = nodal_load_entries1 + nodal_load_entries2 + nodal_load_entries3 + nodal_load_entries4

    
# Process and save the load cases


beamx= [
        "bx7",
        "bx8",
    ]
beamy= [
        "by11",
        "by12",
    ]

loading1 = [{
    "LoadCase": ["DL"],
    "uniform": {"x": 0.0, "y": 0.0, "z": -1.0},
    "point": {"x": 0.0, "y": 0.0, "z": 0.0, "location": 0.5},
    "temperature_points": [
        {"temp": 0.0, "y": 0.0},
        {"temp": 0.0, "y": 0.0},
        {"temp": 0.0, "y": 0.0}
    ]
}]

loading2 = [{
    "LoadCase": ["LL"],
    "uniform": {"x": 0.0, "y": 0.0, "z": -1.0},
    "point": {"x": 0.0, "y": 0.0, "z": 0.0, "location": 0.5},
    "temperature_points": [
        {"temp": 0.0, "y": 0.0},
        {"temp": 0.0, "y": 0.0},
        {"temp": 0.0, "y": 0.0}
    ]
}]

loading_mapping = [(beamx, loading1), (beamy, loading2)]


# all_element_loads = []

load_combinations = {
    "mass": [("DL", 1.0), ("LL", 0.25), ("self_weight", 1.0)],
    "Comb2": [("DL", 1.2), ("LL", 1.6), ("self_weight", 1.2)],
    "Comb1": [("DL", 1.4), ("self_weight", 1.4)],
    "unfactored_load": [("DL", 1.0), ("LL", 1.0),  ("self_weight", 1.0)],
}

load_combinationsiiii = {
    "LC1": [("DL", 1.4)],
    "LC2": [("DL", 1.2), ("LL", 1.6), ("Lr", 0.5)],
    "LC3": [("DL", 1.2), ("Lr", 1.6), ("LL", 1.0)],
    "LC4": [("DL", 1.2), ("Lr", 1.6), ("Wx", 0.8)],
    "LC5": [("DL", 1.2), ("Lr", 1.6), ("Wx", -0.8)],
    "LC6": [("DL", 1.2), ("Lr", 1.6), ("Wy", 0.8)],
    "LC7": [("DL", 1.2), ("Lr", 1.6), ("Wy", -0.8)],
    "LC8": [("DL", 1.2), ("W", 1.6), ("LL", 1.0), ("Lr", 0.5)],
    "LC9": [("DL", 1.2), ("Wx", -1.6), ("LL", 1.0), ("Lr", 0.5)],
    "LC10": [("DL", 1.2), ("Wy", 1.6), ("LL", 1.0), ("Lr", 0.5)],
    "LC11": [("DL", 1.2), ("Wy", -1.6), ("LL", 1.0), ("Lr", 0.5)],
    "LC12": [("DL", 1.2), ("EQx", 1.0), ("LL", 1.0)],
    "LC13": [("DL", 1.2), ("EQx", -1.0), ("LL", 1.0)],
    "LC14": [("DL", 1.2), ("EQy", 1.0), ("LL", 1.0)],
    "LC15": [("DL", 1.2), ("EQy", -1.0), ("LL", 1.0)],
    "LC16": [("DL", 0.9), ("Wx", 1.6)],
    "LC17": [("DL", 0.9), ("Wx", -1.6)],
    "LC18": [("DL", 0.9), ("Wy", 1.6)],
    "LC19": [("DL", 0.9), ("Wy", -1.6)],
    "LC20": [("DL", 0.9), ("EQx", 1.0)],
    "LC21": [("DL", 0.9), ("EQx", -1.0)],
    "LC22": [("DL", 0.9), ("EQy", 1.0)],
    "LC23": [("DL", 0.9), ("EQy", -1.0)],

    # Seismic Load Combinations for Category C(v) and D (Ev = 0.11*DL)
    "LC24": [("DL", 1.31), ("EQx", 1.0), ("EQy", 0.3), ("LL", 1.0)],
    "LC25": [("DL", 1.31), ("EQx", -1.0), ("EQy", 0.3), ("LL", 1.0)],
    "LC26": [("DL", 1.31), ("EQx", 1.0), ("EQy", -0.3), ("LL", 1.0)],
    "LC27": [("DL", 1.31), ("EQx", -1.0), ("EQy", -0.3), ("LL", 1.0)],
    "LC28": [("DL", 1.31), ("EQx", 0.3), ("EQy", 1.0), ("LL", 1.0)],
    "LC29": [("DL", 1.31), ("EQx", -0.3), ("EQy", 1.0), ("LL", 1.0)],
    "LC30": [("DL", 1.31), ("EQx", 0.3), ("EQy", -1.0), ("LL", 1.0)],
    "LC31": [("DL", 1.31), ("EQx", -0.3), ("EQy", -1.0), ("LL", 1.0)],
    "LC32": [("DL", 1.01), ("EQx", 1.0), ("EQy", 0.3), ("LL", 1.0)],
    "LC33": [("DL", 1.01), ("EQx", -1.0), ("EQy", 0.3), ("LL", 1.0)],
    "LC34": [("DL", 1.01), ("EQx", 1.0), ("EQy", -0.3), ("LL", 1.0)],
    "LC35": [("DL", 1.01), ("EQx", -1.0), ("EQy", -0.3), ("LL", 1.0)],
    "LC36": [("DL", 1.01), ("EQx", 0.3), ("EQy", 1.0), ("LL", 1.0)],
    "LC37": [("DL", 1.01), ("EQx", -0.3), ("EQy", 1.0), ("LL", 1.0)],
    "LC38": [("DL", 1.01), ("EQx", 0.3), ("EQy", -1.0), ("LL", 1.0)],
    "LC39": [("DL", 1.01), ("EQx", -0.3), ("EQy", -1.0), ("LL", 1.0)]
}


zero_length_nodes = [
    {"name": "n1",  "Kx": "1e8", "Ky": "1e8", "Kz": "1e8"},
    {"name": "n2",  "Kx": "1e8", "Ky": "1e8", "Kz": "1e8"},
    {"name": "n3",  "Kx": "1e8", "Ky": "1e8", "Kz": "1e8"},
    {"name": "n4",  "Kx": "1e8", "Ky": "1e8", "Kz": "1e8"},
]

# create_zero_length_elements(JSON_FOLDER, zero_length_nodes)

def load_comb(type, Z, S):
    Ev = 0.5 * (2 / 3) * Z * S

    if type == "B":
        load_combinationsiiii = {
            "LC1": [("DL", 1.4)],
            "LC2": [("DL", 1.2), ("LL", 1.6), ("Lr", 0.5)],
            "LC3": [("DL", 1.2), ("Lr", 1.6), ("LL", 1.0)],
            "LC4": [("DL", 1.2), ("Lr", 1.6), ("Wx", 0.8)],
            "LC5": [("DL", 1.2), ("Lr", 1.6), ("Wx", -0.8)],
            "LC6": [("DL", 1.2), ("Lr", 1.6), ("Wy", 0.8)],
            "LC7": [("DL", 1.2), ("Lr", 1.6), ("Wy", -0.8)],
            "LC8": [("DL", 1.2), ("W", 1.6), ("LL", 1.0), ("Lr", 0.5)],
            "LC9": [("DL", 1.2), ("Wx", -1.6), ("LL", 1.0), ("Lr", 0.5)],
            "LC10": [("DL", 1.2), ("Wy", 1.6), ("LL", 1.0), ("Lr", 0.5)],
            "LC11": [("DL", 1.2), ("Wy", -1.6), ("LL", 1.0), ("Lr", 0.5)],
            "LC12": [("DL", 1.2), ("EQx", 1.0), ("LL", 1.0)],
            "LC13": [("DL", 1.2), ("EQx", -1.0), ("LL", 1.0)],
            "LC14": [("DL", 1.2), ("EQy", 1.0), ("LL", 1.0)],
            "LC15": [("DL", 1.2), ("EQy", -1.0), ("LL", 1.0)],
            "LC16": [("DL", 0.9), ("Wx", 1.6)],
            "LC17": [("DL", 0.9), ("Wx", -1.6)],
            "LC18": [("DL", 0.9), ("Wy", 1.6)],
            "LC19": [("DL", 0.9), ("Wy", -1.6)],
            "LC20": [("DL", 0.9), ("EQx", 1.0)],
            "LC21": [("DL", 0.9), ("EQx", -1.0)],
            "LC22": [("DL", 0.9), ("EQy", 1.0)],
            "LC23": [("DL", 0.9), ("EQy", -1.0)],
        }
    else:
        load_combinationsiiii = {
            "LC1": [("DL", 1.4)],
            "LC2": [("DL", 1.2), ("LL", 1.6), ("Lr", 0.5)],
            "LC3": [("DL", 1.2), ("Lr", 1.6), ("LL", 1.0)],
            "LC4": [("DL", 1.2), ("Lr", 1.6), ("Wx", 0.8)],
            "LC5": [("DL", 1.2), ("Lr", 1.6), ("Wx", -0.8)],
            "LC6": [("DL", 1.2), ("Lr", 1.6), ("Wy", 0.8)],
            "LC7": [("DL", 1.2), ("Lr", 1.6), ("Wy", -0.8)],
            "LC8": [("DL", 1.2), ("W", 1.6), ("LL", 1.0), ("Lr", 0.5)],
            "LC9": [("DL", 1.2), ("Wx", -1.6), ("LL", 1.0), ("Lr", 0.5)],
            "LC10": [("DL", 1.2), ("Wy", 1.6), ("LL", 1.0), ("Lr", 0.5)],
            "LC11": [("DL", 1.2), ("Wy", -1.6), ("LL", 1.0), ("Lr", 0.5)],
            "LC12": [("DL", 0.9), ("Wx", 1.6)],
            "LC13": [("DL", 0.9), ("Wx", -1.6)],
            "LC14": [("DL", 0.9), ("Wy", 1.6)],
            "LC15": [("DL", 0.9), ("Wy", -1.6)],
            "LC16": [("DL", 1.2 + Ev), ("EQx", 1.0), ("EQy", 0.3), ("LL", 1.0)],
            "LC17": [("DL", 1.2 + Ev), ("EQx", -1.0), ("EQy", 0.3), ("LL", 1.0)],
            "LC18": [("DL", 1.2 + Ev), ("EQx", 1.0), ("EQy", -0.3), ("LL", 1.0)],
            "LC19": [("DL", 1.2 + Ev), ("EQx", -1.0), ("EQy", -0.3), ("LL", 1.0)],
            "LC20": [("DL", 1.2 + Ev), ("EQx", 0.3), ("EQy", 1.0), ("LL", 1.0)],
            "LC21": [("DL", 1.2 + Ev), ("EQx", -0.3), ("EQy", 1.0), ("LL", 1.0)],
            "LC22": [("DL", 1.2 + Ev), ("EQx", 0.3), ("EQy", -1.0), ("LL", 1.0)],
            "LC23": [("DL", 1.2 + Ev), ("EQx", -0.3), ("EQy", -1.0), ("LL", 1.0)],
            "LC24": [("DL", 0.9 + Ev), ("EQx", 1.0), ("EQy", 0.3), ("LL", 1.0)],
            "LC25": [("DL", 0.9 + Ev), ("EQx", -1.0), ("EQy", 0.3), ("LL", 1.0)],
            "LC26": [("DL", 0.9 + Ev), ("EQx", 1.0), ("EQy", -0.3), ("LL", 1.0)],
            "LC27": [("DL", 0.9 + Ev), ("EQx", -1.0), ("EQy", -0.3), ("LL", 1.0)],
            "LC28": [("DL", 0.9 + Ev), ("EQx", 0.3), ("EQy", 1.0), ("LL", 1.0)],
            "LC29": [("DL", 0.9 + Ev), ("EQx", -0.3), ("EQy", 1.0), ("LL", 1.0)],
            "LC30": [("DL", 0.9 + Ev), ("EQx", 0.3), ("EQy", -1.0), ("LL", 1.0)],
            "LC31": [("DL", 0.9 + Ev), ("EQx", -0.3), ("EQy", -1.0), ("LL", 1.0)],
        }

    return load_combinationsiiii

combinations = load_comb("C", Z=0.28, S=1.15)
print(combinations)


























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
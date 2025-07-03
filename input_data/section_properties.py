import numpy as np
from input_data.units import *

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
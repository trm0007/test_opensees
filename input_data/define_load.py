import numpy as np
from input_data.units import *



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



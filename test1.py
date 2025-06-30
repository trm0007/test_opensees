import os

def find_warning_in_all_py_files(target_text):
    cwd = os.getcwd()
    
    py_files = []

    for dirpath, _, filenames in os.walk(cwd):
        for file in filenames:
            if file.endswith(".py"):
                py_files.append(os.path.join(dirpath, file))

    for file_path in py_files:
        with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
            for line_num, line in enumerate(f, 1):
                # Only match lines that contain the warning, but not the definition of the string
                if target_text in line and "target_text" not in line:
                    print(f"{file_path} (Line {line_num}): {line.strip()}")

# target_text = "Found X,Y match"
target_text = "save_element_loads"

# find_warning_in_all_py_files(target_text)



import os
print("Current Working Directory:", os.getcwd())


def find_txt_in_py_files(base_folder, target_txt):
    for root, dirs, files in os.walk(base_folder):
        for file in files:
            # print(file)
            if file.endswith('.py') and file != '__init__.py':
                file_path = os.path.join(root, file)
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                    if target_txt in content:
                        print(f'Found in: {file_path}')

# base_path = "sectionproperties"
base_path = "C:\\Users\\User\\Desktop\\abcd\\opensees_final"

search_text = 'loading_mapping'
find_txt_in_py_files(base_path, search_text)



# import os

# def get_all_filenames(folder_path):
#     filenames = []
#     for file in os.listdir(folder_path):
#         full_path = os.path.join(folder_path, file)
#         if os.path.isfile(full_path):
#             filenames.append(file)
#     return filenames

# # Example usage
# folder = 'output_folder\images'
# all_files = get_all_filenames(folder)
# for name in all_files:
#     print(name)

# nodes_data = [
#     {"id": 1, "name": "n1"},
#     {"id": 2, "name": "n2"},
#     {"id": 3, "name": "n3"},
# ]

# nodal_masses = {
#     1: [1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
#     2: [2.0, 2.0, 0.0, 0.0, 0.0, 0.0],
# }

# scaling_factors = {
#     1: 2.5,
#     2: 3.6,
#     # default for others
# }

# mass_commands = []
# for node in nodes_data:
#     node_id = node["id"]
#     original_mass = nodal_masses.get(node_id, [0] * 6)
#     factor = scaling_factors.get(node_id, 1.0)
#     scaled_mass = [round(m * factor, 4) for m in original_mass]
#     mass_commands.append(f"mass {node_id} {' '.join(map(str, scaled_mass))}")

# for cmd in mass_commands:
#     print(cmd)

slab = 150*(0.02+0.045+0.075)
column = 4 * 10 * 10 * 10 * 0.15 /144
beam1 = 2*10*15*10*0.15/144
beam2 = 2*12*15*15*0.15/144
bu1 = 2*1*10*1
bu2 = 2*1*15*1
print( slab + column + beam1 + beam2 + bu1+bu2)


load_data1 = {
  "n10002": [
    0.0,
    0.0,
    -1.166666666666667,
    0.0,
    0.0,
    0.0
  ],
  "n10003": [
    0.0,
    0.0,
    -2.333333333333334,
    0.0,
    0.0,
    0.0
  ],
  "n10004": [
    0.0,
    0.0,
    -1.166666666666667,
    0.0,
    0.0,
    0.0
  ],
  "n6": [
    0.0,
    0.0,
    -3.2916666666666665,
    0.0,
    0.0,
    0.0
  ],
  "n10006": [
    0.0,
    0.0,
    -1.1666666666666665,
    0.0,
    0.0,
    0.0
  ],
  "n10005": [
    0.0,
    0.0,
    -2.333333333333333,
    0.0,
    0.0,
    0.0
  ],
  "n5": [
    0.0,
    0.0,
    -3.2916666666666665,
    0.0,
    0.0,
    0.0
  ],
  "n10008": [
    0.0,
    0.0,
    -1.1666666666666665,
    0.0,
    0.0,
    0.0
  ],
  "n10009": [
    0.0,
    0.0,
    -2.333333333333333,
    0.0,
    0.0,
    0.0
  ],
  "n10010": [
    0.0,
    0.0,
    -1.166666666666667,
    0.0,
    0.0,
    0.0
  ],
  "n10011": [
    0.0,
    0.0,
    -2.333333333333332,
    0.0,
    0.0,
    0.0
  ],
  "n10012": [
    0.0,
    0.0,
    -1.1666666666666663,
    0.0,
    0.0,
    0.0
  ],
  "n10013": [
    0.0,
    0.0,
    -1.1666666666666665,
    0.0,
    0.0,
    0.0
  ],
  "n8": [
    0.0,
    0.0,
    -3.2916666666666665,
    0.0,
    0.0,
    0.0
  ],
  "n10015": [
    0.0,
    0.0,
    -1.166666666666666,
    0.0,
    0.0,
    0.0
  ],
  "n7": [
    0.0,
    0.0,
    -3.291666666666666,
    0.0,
    0.0,
    0.0
  ],
  "n1": [
    0.0,
    0.0,
    -0.5208333333333333,
    0.0,
    0.0,
    0.0
  ],
  "n2": [
    0.0,
    0.0,
    -0.5208333333333333,
    0.0,
    0.0,
    0.0
  ],
  "n3": [
    0.0,
    0.0,
    -0.5208333333333333,
    0.0,
    0.0,
    0.0
  ],
  "n4": [
    0.0,
    0.0,
    -0.5208333333333333,
    0.0,
    0.0,
    0.0
  ]
}

load_data={
  "DL": {
    "n10002": [
      0.0,
      0.0,
      -0.16666666666666669,
      0.0,
      0.0,
      0.0
    ],
    "n10003": [
      0.0,
      0.0,
      -0.33333333333333337,
      0.0,
      0.0,
      0.0
    ],
    "n10004": [
      0.0,
      0.0,
      -0.16666666666666674,
      0.0,
      0.0,
      0.0
    ],
    "n6": [
      0.0,
      0.0,
      -0.08333333333333337,
      0.0,
      0.0,
      0.0
    ],
    "n10006": [
      0.0,
      0.0,
      -0.16666666666666663,
      0.0,
      0.0,
      0.0
    ],
    "n10005": [
      0.0,
      0.0,
      -0.33333333333333326,
      0.0,
      0.0,
      0.0
    ],
    "n5": [
      0.0,
      0.0,
      -0.08333333333333333,
      0.0,
      0.0,
      0.0
    ],
    "n10008": [
      0.0,
      0.0,
      -0.16666666666666666,
      0.0,
      0.0,
      0.0
    ],
    "n10009": [
      0.0,
      0.0,
      -0.33333333333333337,
      0.0,
      0.0,
      0.0
    ],
    "n10010": [
      0.0,
      0.0,
      -0.1666666666666667,
      0.0,
      0.0,
      0.0
    ],
    "n10011": [
      0.0,
      0.0,
      -0.3333333333333332,
      0.0,
      0.0,
      0.0
    ],
    "n10012": [
      0.0,
      0.0,
      -0.16666666666666663,
      0.0,
      0.0,
      0.0
    ],
    "n10013": [
      0.0,
      0.0,
      -0.16666666666666666,
      0.0,
      0.0,
      0.0
    ],
    "n8": [
      0.0,
      0.0,
      -0.08333333333333336,
      0.0,
      0.0,
      0.0
    ],
    "n10015": [
      0.0,
      0.0,
      -0.16666666666666657,
      0.0,
      0.0,
      0.0
    ],
    "n7": [
      0.0,
      0.0,
      -0.08333333333333329,
      0.0,
      0.0,
      0.0
    ]
  },
  "LL": {
    "n10002": [
      0.0,
      0.0,
      -0.375,
      0.0,
      0.0,
      0.0
    ],
    "n10003": [
      0.0,
      0.0,
      -0.75,
      0.0,
      0.0,
      0.0
    ],
    "n10004": [
      0.0,
      0.0,
      -0.3750000000000001,
      0.0,
      0.0,
      0.0
    ],
    "n6": [
      0.0,
      0.0,
      -0.18750000000000006,
      0.0,
      0.0,
      0.0
    ],
    "n10006": [
      0.0,
      0.0,
      -0.3749999999999999,
      0.0,
      0.0,
      0.0
    ],
    "n10005": [
      0.0,
      0.0,
      -0.7499999999999998,
      0.0,
      0.0,
      0.0
    ],
    "n5": [
      0.0,
      0.0,
      -0.18749999999999997,
      0.0,
      0.0,
      0.0
    ],
    "n10008": [
      0.0,
      0.0,
      -0.37499999999999994,
      0.0,
      0.0,
      0.0
    ],
    "n10009": [
      0.0,
      0.0,
      -0.7499999999999999,
      0.0,
      0.0,
      0.0
    ],
    "n10010": [
      0.0,
      0.0,
      -0.3750000000000001,
      0.0,
      0.0,
      0.0
    ],
    "n10011": [
      0.0,
      0.0,
      -0.7499999999999997,
      0.0,
      0.0,
      0.0
    ],
    "n10012": [
      0.0,
      0.0,
      -0.3749999999999999,
      0.0,
      0.0,
      0.0
    ],
    "n10013": [
      0.0,
      0.0,
      -0.37499999999999994,
      0.0,
      0.0,
      0.0
    ],
    "n8": [
      0.0,
      0.0,
      -0.18750000000000006,
      0.0,
      0.0,
      0.0
    ],
    "n10015": [
      0.0,
      0.0,
      -0.3749999999999998,
      0.0,
      0.0,
      0.0
    ],
    "n7": [
      0.0,
      0.0,
      -0.1874999999999999,
      0.0,
      0.0,
      0.0
    ]
  },
  "self_weight": {
    "n10002": [
      0.0,
      0.0,
      -0.6250000000000001,
      0.0,
      0.0,
      0.0
    ],
    "n10003": [
      0.0,
      0.0,
      -1.2500000000000002,
      0.0,
      0.0,
      0.0
    ],
    "n10004": [
      0.0,
      0.0,
      -0.6250000000000002,
      0.0,
      0.0,
      0.0
    ],
    "n6": [
      0.0,
      0.0,
      -3.020833333333333,
      0.0,
      0.0,
      0.0
    ],
    "n10006": [
      0.0,
      0.0,
      -0.6249999999999999,
      0.0,
      0.0,
      0.0
    ],
    "n10005": [
      0.0,
      0.0,
      -1.2499999999999998,
      0.0,
      0.0,
      0.0
    ],
    "n5": [
      0.0,
      0.0,
      -3.020833333333333,
      0.0,
      0.0,
      0.0
    ],
    "n10008": [
      0.0,
      0.0,
      -0.6249999999999999,
      0.0,
      0.0,
      0.0
    ],
    "n10009": [
      0.0,
      0.0,
      -1.25,
      0.0,
      0.0,
      0.0
    ],
    "n10010": [
      0.0,
      0.0,
      -0.6250000000000002,
      0.0,
      0.0,
      0.0
    ],
    "n10011": [
      0.0,
      0.0,
      -1.2499999999999996,
      0.0,
      0.0,
      0.0
    ],
    "n10012": [
      0.0,
      0.0,
      -0.6249999999999998,
      0.0,
      0.0,
      0.0
    ],
    "n10013": [
      0.0,
      0.0,
      -0.6249999999999999,
      0.0,
      0.0,
      0.0
    ],
    "n8": [
      0.0,
      0.0,
      -3.020833333333333,
      0.0,
      0.0,
      0.0
    ],
    "n10015": [
      0.0,
      0.0,
      -0.6249999999999998,
      0.0,
      0.0,
      0.0
    ],
    "n7": [
      0.0,
      0.0,
      -3.020833333333333,
      0.0,
      0.0,
      0.0
    ],
    "n1": [
      0.0,
      0.0,
      -0.5208333333333333,
      0.0,
      0.0,
      0.0
    ],
    "n2": [
      0.0,
      0.0,
      -0.5208333333333333,
      0.0,
      0.0,
      0.0
    ],
    "n3": [
      0.0,
      0.0,
      -0.5208333333333333,
      0.0,
      0.0,
      0.0
    ],
    "n4": [
      0.0,
      0.0,
      -0.5208333333333333,
      0.0,
      0.0,
      0.0
    ]
  }
}

def sum_z_direction_loads(load_data):
    total_load = 0.0
    
    # Iterate through each load category (DL, LL, self_weight)
    for category in load_data.values():
        # Iterate through each node in the category
        for node_load in category.values():
            # The z-direction load is the 3rd element (index 2)
            z_load = node_load[2]
            total_load += z_load
    
    return total_load

total_z_load = sum_z_direction_loads(load_data)
print("Total Z Direction Load:", total_z_load)

# print("load_data:", sum_z_direction_loads(load_data))
total_z_load = sum(values[2] for values in load_data1.values())
print("Total Z Direction Load:", total_z_load)
print(18.14+ 8.44+8.51-1.14)

data = [
    {
        "id": 1,
        "reactions": [5.268703470110725, 8.097712461040633, 18.1400094730076]
    },
    {
        "id": 2,
        "reactions": [4.71610255079849, 7.963348267163205, 8.448467914261016]
    },
    {
        "id": 3,
        "reactions": [5.296525395484837, 7.033546536468834, 8.513513875907012]
    },
    {
        "id": 4,
        "reactions": [4.718668583599228, 6.905392735316315, -1.1853245965089476]
    }
]

sum_reactions = [0, 0, 0]
for item in data:
    for i in range(3):
        sum_reactions[i] += item["reactions"][i]

print(sum_reactions)

total = sum(sum_reactions)
print(total)

ratio = 1000  # kg/m²
L = 10.0/3.28
B = 15.0/3.28

Peso_s = ratio * L * B
Masa_nodal = Peso_s / 8 / 9.81 / 1000  # kN·s²/m

print("Peso_s (kg):", Peso_s)
print("Masa_nodal (kN·s²/m):", Masa_nodal)

# Convert kN·s²/m to kip·s²/ft
# 1 kN = 0.2248089431 kip, 1 m = 3.280839895 ft
Masa_nodal_kip_ft = Masa_nodal * 0.2248089431 / 3.280839895

print("Masa_nodal (kip·s²/ft):", Masa_nodal_kip_ft)


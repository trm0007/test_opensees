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


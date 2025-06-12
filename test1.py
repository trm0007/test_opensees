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
target_text = "member_load_"

find_warning_in_all_py_files(target_text)




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

# nodes, members, structure_data = load_structure(JSON_FOLDER)

import json
import os
import re


# Call the first function to load structure data
nodes_dict, members_dict, structure = load_structure(JSON_FOLDER)

# Call the second function to create combined structure JSON
combined_data = create_combined_structure_json(JSON_FOLDER)
# 2. Load the combined structure data to get all node coordinates
data = merge_structures(structure, combined_data)
print(f'data={data}')
# Extract and sort nodes using 'coordinates'
# Extract node details
nodes = data['nodes']
for node in nodes:
    node_id = node['id']
    node_name = node['name']
    node_coordinates = {
        'x': node.get('x', 0.0),  # Default to 0.0 if not present
        'y': node.get('y', 0.0),
        'z': node.get('z', 0.0)
    }
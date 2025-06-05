data = {
    "elements": {
        "R10001": {
            "id": 10001,
            "type": "rectangle",
            "nodes": [
                "N10002",
                "N10003",
                "N10004",
                "n6"
            ],
            "shell_section": "floor_slab",
            "thickness": 10,
        },
        "R10002": {
            "id": 10002,
            "type": "rectangle",
            "nodes": [
                "n5",
                "N10006",
                "N10003",
                "N10002"
            ],
            "shell_section": "floor_slab",
            "thickness": 10,
        },
        "R10003": {
            "id": 10003,
            "type": "rectangle",
            "nodes": [
                "N10003",
                "N10007",
                "n8",
                "N10004"
            ],
            "shell_section": "floor_slab",
            "thickness": 10,
        },
        "R10004": {
            "id": 10004,
            "type": "rectangle",
            "nodes": [
                "N10006",
                "n7",
                "N10007",
                "N10003"
            ],
            "shell_section": "floor_slab",
            "thickness": 10,
        }
    },
    "nodes": {
        "N10002": {
            "id": 10002,
            "coordinates": [
                10.0,
                0.0,
                11.5
            ],
            "is_predefined": False  # Changed from false to False
        },
        "N10003": {
            "id": 10003,
            "coordinates": [
                10.0,
                11.5,
                11.5
            ],
            "is_predefined": False  # Changed from false to False
        },
        "N10004": {
            "id": 10004,
            "coordinates": [
                20.0,
                11.5,
                11.5
            ],
            "is_predefined": False  # Changed from false to False
        },
        "N10006": {
            "id": 10006,
            "coordinates": [
                0.0,
                11.5,
                11.5
            ],
            "is_predefined": False  # Changed from false to False
        },
        "N10007": {
            "id": 10007,
            "coordinates": [
                10.0,
                23.0,
                11.5
            ],
            "is_predefined": False  # Changed from false to False
        }
    }
}


def check_wall_orientation(data):
    nodes = data['nodes']
    x_values = [node['coordinates'][0] for node in nodes.values()]
    y_values = [node['coordinates'][1] for node in nodes.values()]
    z_values = [node['coordinates'][2] for node in nodes.values()]
    
    if all(z == z_values[0] for z in z_values):
        # Horizontal case (all z equal)
        Dz = next(iter(data['elements'].values()))['thickness']
        Dx = max(x_values) - min(x_values)
        Dy = max(y_values) - min(y_values)
        return {
            "wall_type": "horizontal", 
            "z_value": z_values[0],
            "Dx": Dx,  # Planar x-dimension
            "Dy": Dy,  # Planar y-dimension
            "Dz": Dz   # Thickness (z-direction)
        }
    elif all(x == x_values[0] for x in x_values):
        # Vertical along y-axis (all x equal)
        Dx = next(iter(data['elements'].values()))['thickness']  # Thickness in x-direction
        Dy = max(y_values) - min(y_values)  # Length along y-axis
        Dz = max(z_values) - min(z_values)  # Height
        return {
            "wall_type": "vertical_y", 
            "x_value": x_values[0],
            "Dx": Dx,  # Thickness
            "Dy": Dy,  # Length
            "Dz": Dz   # Height
        }
    elif all(y == y_values[0] for y in y_values):
        # Vertical along x-axis (all y equal)
        Dy = next(iter(data['elements'].values()))['thickness']  # Thickness in y-direction
        Dx = max(x_values) - min(x_values)  # Length along x-axis
        Dz = max(z_values) - min(z_values)  # Height
        return {
            "wall_type": "vertical_x", 
            "y_value": y_values[0],
            "Dx": Dx,  # Length
            "Dy": Dy,  # Thickness
            "Dz": Dz   # Height
        }
    else:
        return {
            "wall_type": "irregular",
            "x_values": x_values,
            "y_values": y_values,
            "z_values": z_values
        }

result = check_wall_orientation(data)
print(f"Wall orientation: {result['wall_type']}")

if result['wall_type'] == "horizontal":
    print(f"All nodes at Z = {result['z_value']}")
    print(f"Dx (length): {result['Dx']}, Dy (width): {result['Dy']}, Dz (thickness): {result['Dz']}")
elif result['wall_type'] == "vertical_y":
    print(f"All nodes at X = {result['x_value']}")
    print(f"Dx (thickness): {result['Dx']}, Dy (length): {result['Dy']}, Dz (height): {result['Dz']}")
elif result['wall_type'] == "vertical_x":
    print(f"All nodes at Y = {result['y_value']}")
    print(f"Dx (length): {result['Dx']}, Dy (thickness): {result['Dy']}, Dz (height): {result['Dz']}")
else:
    print("Irregular wall with varying coordinates")
    print(f"X range: {min(result['x_values'])} to {max(result['x_values'])}")
    print(f"Y range: {min(result['y_values'])} to {max(result['y_values'])}")
    print(f"Z range: {min(result['z_values'])} to {max(result['z_values'])}")
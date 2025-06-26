from import_ import *
from Grid_and_structure_creation import load_structure
from units import *


# def create_structural_model(section_definitions, structure_data_path, member_section_mapping_path, json_file_path, numIntgrPts=8):
#     """Create structural model with robust file handling"""
    
#     try:
#         # 1. Load input data with proper error handling
#         try:
#             with open(structure_data_path, 'r') as f:
#                 structure_data = json.load(f)
            
#             with open(member_section_mapping_path, 'r') as f:
#                 member_section_mapping = json.load(f)
#         except FileNotFoundError as e:
#             print(f"❌ Error: Input file not found: {e.filename}")
#             raise
#         except json.JSONDecodeError as e:
#             print(f"❌ Error: Invalid JSON in input file: {e}")
#             raise
        
#         # 2. Process nodes and members
#         nodes_data = structure_data["nodes"]
#         members_data = structure_data["members"]
        
#         nodes = np.zeros((len(nodes_data), 3))
#         for node in nodes_data:
#             node_id = node["id"]
#             nodes[node_id-1] = [node["x"], node["y"], node["z"]]
        
#         # 3. Create section tags dynamically
#         section_tags = {}
#         tag_counter = 1
#         for category in section_definitions.values():
#             for section_name in category:
#                 section_tags[section_name] = tag_counter
#                 tag_counter += 1
        
#         # 4. Initialize result containers
#         shapes = {}
#         memberLengths = []
#         element_data = []
        
#         # 5. Process each member
#         for mbr_data in members_data:
#             mbr_id = mbr_data["id"]
#             mbr_name = mbr_data["name"]
#             node_i_id = mbr_data["start_node_id"] 
#             node_j_id = mbr_data["end_node_id"]
            
            
#             # Get node coordinates
#             node_i = node_i_id - 1
#             node_j = node_j_id - 1
            
#             # Calculate member length and orientation
#             dx, dy, dz = nodes[node_j] - nodes[node_i]
#             length = math.sqrt(dx**2 + dy**2 + dz**2)
#             memberLengths.append(length)
            
#             # Determine transformation type
#             transType = 'PDelta' if mbr_name.startswith("cz") else 'Linear'
            
#             # Get section tag
#             section_name = member_section_mapping[mbr_name]
#             secTag = section_tags[section_name]
#             # Extract area and unit weight from section definitions
#             area = None
#             unit_weight = None
            
#             # Find the section in section_definitions to get area and unit_weight
#             for category, sections in section_definitions.items():
#                 if section_name in sections:
#                     area = sections[section_name]["area"]
#                     unit_weight = sections[section_name]["unit_weight"]
#                     break
                
            
#             # Calculate local coordinate system
#             local_x_unit = np.array([dx, dy, dz]) / length if length > 0 else np.array([1.0, 0.0, 0.0])
#             reference_vector = np.array([0.0, 1.0, 0.0])
            
#             if abs(np.dot(local_x_unit, reference_vector)) > 0.99:
#                 reference_vector = np.array([1.0, 0.0, 0.0])
            
#             local_z_unit = np.cross(local_x_unit, reference_vector)
#             local_z_unit /= np.linalg.norm(local_z_unit)
            
#             # Store element data
#             element_data.append({
#                 "transType": transType,
#                 "transTag": mbr_id,
#                 "vecxz": local_z_unit.tolist(),
#                 "eleTag": mbr_id,
#                 "node_i_id": node_i_id,
#                 "node_j_id": node_j_id,
#                 "numIntgrPts": numIntgrPts,
#                 "secTag": secTag,
#                 "length": length,
#                 "area": area, 
#                 "unit_weight": unit_weight

#             })
        
#         # 6. Save output data with proper error handling
#         try:
#             # Ensure directory exists
#             os.makedirs(os.path.dirname(json_file_path), exist_ok=True)
            
#             with open(json_file_path, "w") as f:
#                 json.dump({
#                     "elements": element_data,
#                     "nodes": nodes_data,
#                     "sections": section_definitions
#                 }, f, indent=4)
            
#             print(f"✅ Successfully saved element data to {json_file_path}")
            
#         except PermissionError:
#             print(f"❌ Error: Permission denied when writing to {json_file_path}")
#             print("Possible solutions:")
#             print("1. Run the script as administrator")
#             print("2. Choose a different output directory")
#             print("3. Check file/folder permissions")
#             raise
#         except Exception as e:
#             print(f"❌ Error saving output file: {e}")
#             raise
        
#         return shapes, memberLengths
    
#     except Exception as e:
#         print(f"❌ Error in create_structural_model: {e}")
#         raise

# def create_beam_column_member(section_definitions, JSON_FOLDER, json_file_path, numIntgrPts=8):
#     """Create structural model with robust file handling"""
    
#     try:
#         # 1. Load input data with proper error handling
#         try:
#             nodes, members, structure_data = load_structure(JSON_FOLDER)
#             member_section_mapping_path = os.path.join(JSON_FOLDER, 'member_section_mapping.json')
#             with open(member_section_mapping_path, 'r') as f:
#                 member_section_mapping = json.load(f)
#         except FileNotFoundError as e:
#             print(f"❌ Error: Input file not found: {e.filename}")
#             raise
#         except json.JSONDecodeError as e:
#             print(f"❌ Error: Invalid JSON in input file: {e}")
#             raise
        
#         # 2. Process nodes and members
#         nodes_data = structure_data["nodes"]
#         members_data = structure_data["members"]
#         # print(f'members_name1470={members_data}')
#         # members_name = members_data["name"]
        
#         nodes = np.zeros((len(nodes_data), 3))
#         for node in nodes_data:
#             node_id = node["id"]
#             nodes[node_id-1] = [node["x"], node["y"], node["z"]]
        
#         # 3. Calculate moments of inertia for all sections
#         for category, sections in section_definitions.items():
#             for section_name, section_data in sections.items():
#                 if section_data["type"] == "rectangular":
#                     B = section_data["B"]
#                     H = section_data["H"]
#                     Iy = (B * H**3) / 12  # Weak axis
#                     Iz = (H * B**3) / 12  # Strong axis
#                     section_data["Iy"] = Iy
#                     section_data["Iz"] = Iz
#                 elif section_data["type"] == "circular":
#                     # print("avs10", section_data["section_tag"])
#                     D = section_data["D_Sec"]
#                     ri = section_data["ri"]
#                     if ri == 0:  # Solid section
#                         Iy = (math.pi * D**4) / 64
#                         Iz = Iy
#                     else:  # Hollow section
#                         di = 2 * ri  # Inner diameter
#                         Iy = (math.pi * (D**4 - di**4)) / 64
#                         Iz = Iy
#                     section_data["Iy"] = Iy
#                     section_data["Iz"] = Iz
#                     B = D  # For compatibility with existing code
#                     H = D   # For compatibility with existing code
#                 elif section_data["type"] == "L":
#                     B = section_data["B"]
#                     H = section_data["H"]
#                     t = section_data["t"]
#                     # Assume L-section formed by two rectangles: vertical and horizontal legs
#                     A1 = B * t
#                     A2 = (H - t) * t
#                     y1 = t / 2
#                     y2 = (H - t) / 2 + t
#                     z1 = B / 2
#                     z2 = t / 2
#                     A_total = A1 + A2
#                     y_bar = (A1 * y1 + A2 * y2) / A_total
#                     z_bar = (A1 * z1 + A2 * z2) / A_total

#                     Iy = (B * t**3) / 12 + A1 * (y1 - y_bar)**2 + (t * (H - t)**3) / 12 + A2 * (y2 - y_bar)**2
#                     Iz = (t * B**3) / 12 + A1 * (z1 - z_bar)**2 + ((H - t) * t**3) / 12 + A2 * (z2 - z_bar)**2
#                     section_data["Iy"] = Iy
#                     section_data["Iz"] = Iz

        
#         # 4. Create section tags dynamically
#         section_tags = {}
#         tag_counter = 1
#         for category in section_definitions.values():
#             for section_name in category:
#                 section_tags[section_name] = tag_counter
#                 tag_counter += 1
        
#         # 5. Initialize result containers
#         shapes = {}
#         memberLengths = []
#         element_data = []
        
#         # 6. Process each member
#         for mbr_data in members_data:
#             mbr_id = mbr_data["id"]
#             mbr_name = mbr_data["name"]
#             # print(f'mbr_name={mbr_name}')
#             node_i_id = mbr_data["start_node_id"] 
#             node_j_id = mbr_data["end_node_id"]
            
#             # Get node coordinates
#             node_i = node_i_id - 1
#             node_j = node_j_id - 1
            
#             # Calculate member length and orientation
#             dx, dy, dz = nodes[node_j] - nodes[node_i]
#             length = math.sqrt(dx**2 + dy**2 + dz**2)
#             memberLengths.append(length)
            
#             # Determine transformation type
#             transType = 'PDelta' if mbr_name.startswith("cz") else 'Linear'
            
#             # Get section tag
#             # section_name = member_section_mapping[mbr_name]
#             section_name = member_section_mapping[mbr_name]["section"]
#             rotation = member_section_mapping[mbr_name]["rotation"]

#             secTag = section_tags[section_name]
            
#             # Extract section properties
#             area = None
#             unit_weight = None
#             Iy = None
#             Iz = None
#             # rotation = None
#             type = None
#             B = None
#             H = None
#             # Find the section in section_definitions to get properties
#             for category, sections in section_definitions.items():
#                 if section_name in sections:
#                     if sections[section_name]["type"] == "rectangular":
#                         B = sections[section_name]["B"]
#                         H = sections[section_name]["H"]
#                         # rotation = sections[section_name]["rotation"]
#                         rotation = rotation
#                         area = B * H
#                         type = sections[section_name]["type"]
#                     elif sections[section_name]["type"] == "circular":
#                         D = sections[section_name]["D_Sec"]
#                         area = 0.785 * D * D
#                         # rotation = sections[section_name]["rotation"]
#                         rotation = rotation
#                         type = sections[section_name]["type"]
#                     elif section_data["type"] == "L":
#                         B = section_data["B"]
#                         H = section_data["H"]
#                         t = section_data["t"]
#                         # Assume L-section formed by two rectangles: vertical and horizontal legs
#                         A1 = B * t
#                         A2 = (H - t) * t
#                         y1 = t / 2
#                         y2 = (H - t) / 2 + t
#                         z1 = B / 2
#                         z2 = t / 2
#                         area = A1 + A2
#                         y_bar = (A1 * y1 + A2 * y2) / A_total
#                         z_bar = (A1 * z1 + A2 * z2) / A_total

#                         Iy = (B * t**3) / 12 + A1 * (y1 - y_bar)**2 + (t * (H - t)**3) / 12 + A2 * (y2 - y_bar)**2
#                         Iz = (t * B**3) / 12 + A1 * (z1 - z_bar)**2 + ((H - t) * t**3) / 12 + A2 * (z2 - z_bar)**2
#                         section_data["Iy"] = Iy
#                         section_data["Iz"] = Iz
#                         rotation = rotation
#                         type = sections[section_name]["type"]

#                     unit_weight = sections[section_name]["unit_weight"]
#                     Iy = sections[section_name]["Iy"]
#                     Iz = sections[section_name]["Iz"]
#                     break
                
#             local_x_unit = np.array([dx, dy, dz]) / length 
#             if mbr_name.startswith("cz"):
#                 angle = rotation
#                 # User-defined rotation angle (degrees)
#                 angle = math.radians(angle)  # 30° → 0.5236 rad
#                 local_z_unit = np.array([math.sin(angle), math.cos(angle), 1.0])
#             else:
#                 # For non-cz elements, use the more detailed approach
#                 ix = nodes[node_i_id-1, 0]  # X-coord of node i of this member
#                 iy = nodes[node_i_id-1, 1]  # y-coord of node i of this member
#                 iz = nodes[node_i_id-1, 2]  # z-coord of node i of this member
#                 jx = nodes[node_j_id-1, 0]  # X-coord of node j of this member
#                 jy = nodes[node_j_id-1, 1]  # y-coord of node j of this member
#                 jz = nodes[node_j_id-1, 2]  # z-coord of node j of this member

#                 # Local y-vector in global RF using Gram-Schmidt process
#                 i_offset = np.array([ix, iy, iz+1])  # Offset node i by 1m in global z-direction
#                 j_offset = np.array([jx, jy, jz+1])  # Offset node j by 1m in global z-direction
#                 node_k = i_offset + 0.5*(j_offset-i_offset)  # Point in the local x-y plane
#                 node_i_coords = np.array([ix, iy, iz])
#                 vector_in_plane = node_k - node_i_coords  # Vector in the x-y plane
#                 local_y_vector = vector_in_plane - np.dot(vector_in_plane, local_x_unit)*local_x_unit 
#                 local_y_unit = local_y_vector / np.linalg.norm(local_y_vector)  # Local unit vector defining local y-axis

#                 # Local z-vector in global RF using matrix cross product
#                 local_z_unit = np.cross(local_x_unit, local_y_unit)  # Local unit vector defining local z-axis

#             # Store element data
#             element_data.append({
#                 "name": mbr_name,
#                 "type": type,
#                 "transType": transType,
#                 "transTag": mbr_id,
#                 "vecxz": local_z_unit.tolist(),
#                 "eleTag": mbr_id,
#                 "node_i_id": node_i_id,
#                 "node_j_id": node_j_id,
#                 "numIntgrPts": numIntgrPts,
#                 "secTag": secTag,
#                 "length": length,
#                 "B": B,
#                 "H": H,
#                 "area": area, 
#                 "unit_weight": unit_weight,
#                 "Iy": Iy,
#                 "Iz": Iz,
#                 "section_name": section_name,
#                 "rotation": rotation
#             })
#         # 7. Save output data with proper error handling
#         try:
#             # Ensure directory exists
#             os.makedirs(os.path.dirname(json_file_path), exist_ok=True)
            
#             with open(json_file_path, "w") as f:
#                 json.dump({
#                     "elements": element_data,
#                     "nodes": nodes_data,
#                     "sections": section_definitions
#                 }, f, indent=4)
            
#             print(f"✅ Successfully saved element data to {json_file_path}")
            
#         except PermissionError:
#             print(f"❌ Error: Permission denied when writing to {json_file_path}")
#             print("Possible solutions:")
#             print("1. Run the script as administrator")
#             print("2. Choose a different output directory")
#             print("3. Check file/folder permissions")
#             raise
#         except Exception as e:
#             print(f"❌ Error saving output file: {e}")
#             raise
        
#         return shapes, memberLengths
    
#     except Exception as e:
#         print(f"❌ Error in create_structural_model: {e}")
#         raise

def create_beam_column_member(section_definitions, JSON_FOLDER, json_file_path, numIntgrPts=8):
    """Create structural model with robust error handling and debugging"""
    
    try:
        # ====================== DEBUG: INITIAL SECTION DATA ======================
        print("\n" + "="*60)
        print("DEBUG: SECTION DEFINITIONS INPUT")
        print("="*60)
        for category, sections in section_definitions.items():
            print(f"\nCategory: {category}")
            for name, props in sections.items():
                print(f"  Section: {name}")
                print(f"  Type: {props['type']}")
                if props['type'] == 'rectangular':
                    print(f"  Dimensions: B={props['B']}, H={props['H']}")
                elif props['type'] == 'circular':
                    print(f"  Dimensions: D_Sec={props['D_Sec']}, ri={props['ri']}")
                elif props['type'] == 'L':
                    print(f"  Dimensions: B={props['B']}, H={props['H']}, t={props['t']}")
                print(f"  Unit Weight: {props['unit_weight']}")

        # ====================== LOAD INPUT FILES ======================
        try:
            print("\n" + "="*60)
            print("DEBUG: LOADING INPUT FILES")
            print("="*60)
            
            nodes, members, structure_data = load_structure(JSON_FOLDER)
            print(f"Loaded {len(structure_data['nodes'])} nodes and {len(structure_data['members'])} members")
            
            member_section_mapping_path = os.path.join(JSON_FOLDER, 'member_section_mapping.json')
            with open(member_section_mapping_path, 'r') as f:
                member_section_mapping = json.load(f)
            print(f"Loaded member-section mapping with {len(member_section_mapping)} entries")
            
            # Debug print first 5 mappings
            print("\nSample member-section mappings:")
            for i, (mbr, mapping) in enumerate(member_section_mapping.items()):
                if i >= 5: break
                print(f"  {mbr} -> Section: {mapping['section']}, Rotation: {mapping['rotation']}")
                
        except Exception as e:
            print(f"\n❌ ERROR LOADING INPUT FILES: {str(e)}")
            raise

        # ====================== PROCESS NODES ======================
        nodes_data = structure_data["nodes"]
        members_data = structure_data["members"]
        
        nodes = np.zeros((len(nodes_data), 3))
        for node in nodes_data:
            node_id = node["id"]
            nodes[node_id-1] = [node["x"], node["y"], node["z"]]
        print(f"\nProcessed {len(nodes_data)} nodes")

        # ====================== CALCULATE SECTION PROPERTIES ======================
        print("\n" + "="*60)
        print("DEBUG: SECTION PROPERTY CALCULATIONS")
        print("="*60)
        
        for category, sections in section_definitions.items():
            for section_name, section_data in sections.items():
                print(f"\nProcessing section: {section_name} (Type: {section_data['type']})")
                
                try:
                    if section_data["type"] == "rectangular":
                        B = section_data["B"]
                        H = section_data["H"]
                        Iy = (B * H**3) / 12
                        Iz = (H * B**3) / 12
                        section_data["Iy"] = Iy
                        section_data["Iz"] = Iz
                        area = B * H
                        print(f"  Rectangular: B={B}, H={H}")
                        print(f"  Calculated: Iy={Iy:.2e}, Iz={Iz:.2e}, Area={area:.2f}")
                        
                    elif section_data["type"] == "circular":
                        D = section_data["D_Sec"]
                        ri = section_data["ri"]
                        if ri == 0:  # Solid section
                            Iy = (math.pi * D**4) / 64
                            Iz = Iy
                            print(f"  Solid Circular: D={D}")
                        else:  # Hollow section
                            di = 2 * ri
                            Iy = (math.pi * (D**4 - di**4)) / 64
                            Iz = Iy
                            print(f"  Hollow Circular: D={D}, di={di}")
                        section_data["Iy"] = Iy
                        section_data["Iz"] = Iz
                        area = math.pi * (D**2 - (di**2 if ri != 0 else 0)) / 4
                        print(f"  Calculated: Iy=Iz={Iy:.2e}, Area={area:.2f}")
                        
                    elif section_data["type"] == "L":
                        B = section_data["B"]
                        H = section_data["H"]
                        t = section_data["t"]
                        A1 = B * t
                        A2 = (H - t) * t
                        y1 = t / 2
                        y2 = (H - t) / 2 + t
                        z1 = B / 2
                        z2 = t / 2
                        A_total = A1 + A2
                        y_bar = (A1 * y1 + A2 * y2) / A_total
                        z_bar = (A1 * z1 + A2 * z2) / A_total

                        Iy = (B * t**3)/12 + A1*(y1-y_bar)**2 + (t*(H-t)**3)/12 + A2*(y2-y_bar)**2
                        Iz = (t * B**3)/12 + A1*(z1-z_bar)**2 + ((H-t)*t**3)/12 + A2*(z2-z_bar)**2
                        section_data["Iy"] = Iy
                        section_data["Iz"] = Iz
                        print(f"  L-section: B={B}, H={H}, t={t}")
                        print(f"  Calculated: Iy={Iy:.2e}, Iz={Iz:.2e}, Area={A_total:.2f}")
                        
                except KeyError as e:
                    print(f"❌ Missing property in section {section_name}: {str(e)}")
                    raise
                except Exception as e:
                    print(f"❌ Error calculating properties for {section_name}: {str(e)}")
                    raise

        # ====================== CREATE SECTION TAGS ======================
        section_tags = {}
        tag_counter = 1
        for category in section_definitions.values():
            for section_name in category:
                section_tags[section_name] = tag_counter
                tag_counter += 1
        print(f"\nCreated {len(section_tags)} section tags")

        # ====================== PROCESS MEMBERS ======================
        print("\n" + "="*60)
        print("DEBUG: MEMBER PROCESSING")
        print("="*60)
        
        shapes = {}
        memberLengths = []
        element_data = []
        error_count = 0
        
        for mbr_data in members_data:
            try:
                mbr_id = mbr_data["id"]
                mbr_name = mbr_data["name"]
                node_i_id = mbr_data["start_node_id"] 
                node_j_id = mbr_data["end_node_id"]
                
                print(f"\nProcessing member {mbr_id}: {mbr_name}")
                
                # Get section mapping
                section_mapping = member_section_mapping.get(mbr_name)
                if not section_mapping:
                    print(f"❌ No section mapping found for member {mbr_name}")
                    error_count += 1
                    continue
                    
                section_name = section_mapping["section"]
                rotation = section_mapping["rotation"]
                print(f"  Assigned section: {section_name}, Rotation: {rotation}°")

                # Find section in definitions
                section_props = None
                for category, sections in section_definitions.items():
                    if section_name in sections:
                        section_props = sections[section_name]
                        break
                
                if not section_props:
                    print(f"❌ Section '{section_name}' not found in definitions")
                    error_count += 1
                    continue
                
                # Calculate member geometry
                node_i = node_i_id - 1
                node_j = node_j_id - 1
                dx, dy, dz = nodes[node_j] - nodes[node_i]
                length = math.sqrt(dx**2 + dy**2 + dz**2)
                memberLengths.append(length)
                
                # Determine transformation type
                transType = 'PDelta' if mbr_name.startswith("cz") else 'Linear'
                
                # Extract section properties
                secTag = section_tags[section_name]
                section_type = section_props["type"]
                unit_weight = section_props["unit_weight"]
                Iy = section_props["Iy"]
                Iz = section_props["Iz"]
                
                # Section-specific properties
                if section_type == "rectangular":
                    B = section_props["B"]
                    H = section_props["H"]
                    area = B * H
                    print(f"  Rectangular section: B={B}, H={H}, Area={area:.2f}")
                    
                elif section_type == "circular":
                    D = section_props["D_Sec"]
                    area = math.pi * (D**2) / 4
                    B = D  # For compatibility
                    H = D  # For compatibility
                    print(f"  Circular section: D={D}, Area={area:.2f}")
                    
                elif section_type == "L":
                    B = section_props["B"]
                    H = section_props["H"]
                    t = section_props["t"]
                    area = (B * t) + ((H - t) * t)
                    print(f"  L-section: B={B}, H={H}, t={t}, Area={area:.2f}")
                    
                else:
                    print(f"❌ Unknown section type: {section_type}")
                    error_count += 1
                    continue
                
                # Validate dimensions
                if B <= 0 or H <= 0:
                    print(f"❌ Invalid dimensions for member {mbr_name}: B={B}, H={H}")
                    error_count += 1
                    continue
                
                # Transformation calculations
                local_x_unit = np.array([dx, dy, dz]) / length
                if mbr_name.startswith("cz"):
                    angle = math.radians(rotation)
                    local_z_unit = np.array([math.sin(angle), math.cos(angle), 1.0])
                else:
                    ix, iy, iz = nodes[node_i_id-1]
                    jx, jy, jz = nodes[node_j_id-1]
                    i_offset = np.array([ix, iy, iz+1])
                    j_offset = np.array([jx, jy, jz+1])
                    node_k = i_offset + 0.5*(j_offset-i_offset)
                    node_i_coords = np.array([ix, iy, iz])
                    vector_in_plane = node_k - node_i_coords
                    local_y_vector = vector_in_plane - np.dot(vector_in_plane, local_x_unit)*local_x_unit
                    local_y_unit = local_y_vector / np.linalg.norm(local_y_vector)
                    local_z_unit = np.cross(local_x_unit, local_y_unit)

                # Store element data
                element_data.append({
                    "name": mbr_name,
                    "type": section_type,
                    "transType": transType,
                    "transTag": mbr_id,
                    "vecxz": local_z_unit.tolist(),
                    "eleTag": mbr_id,
                    "node_i_id": node_i_id,
                    "node_j_id": node_j_id,
                    "numIntgrPts": numIntgrPts,
                    "secTag": secTag,
                    "length": length,
                    "B": B,
                    "H": H,
                    "t": t,
                    "area": area,
                    "unit_weight": unit_weight,
                    "Iy": Iy,
                    "Iz": Iz,
                    "section_name": section_name,
                    "rotation": rotation
                })
                print(f"  Successfully processed member {mbr_name}")
                
            except Exception as e:
                print(f"❌ Error processing member {mbr_name}: {str(e)}")
                error_count += 1
                continue

        # ====================== SAVE OUTPUT ======================
        print("\n" + "="*60)
        print("DEBUG: SAVING OUTPUT")
        print("="*60)
        
        print(f"\nProcessing complete with {error_count} errors out of {len(members_data)} members")
        
        try:
            os.makedirs(os.path.dirname(json_file_path), exist_ok=True)
            
            output_data = {
                "elements": element_data,
                "nodes": nodes_data,
                "sections": section_definitions
            }
            
            with open(json_file_path, "w") as f:
                json.dump(output_data, f, indent=4)
            
            print(f"\n✅ Successfully saved element data to {json_file_path}")
            print(f"  Contains {len(element_data)} elements")
            
            # Debug print first element
            print("\nSample element data:")
            print(json.dumps(element_data[0], indent=2))
            
        except Exception as e:
            print(f"\n❌ ERROR SAVING OUTPUT: {str(e)}")
            raise

        return shapes, memberLengths
    
    except Exception as e:
        print("\n" + "="*60)
        print("❌ CRITICAL ERROR IN create_beam_column_member")
        print("="*60)
        print(f"Error Type: {type(e).__name__}")
        print(f"Error Message: {str(e)}")
        print("Stack Trace:")
        # traceback.print_exc()
        raise

# # Example usage:
# if __name__ == "__main__":
#     structure_data_path = "structure_data_json.json"
#     member_section_mapping_path = "member_section_mapping.json"
    
#     # Initialize OpenSees model
#     ops.wipe()
#     ops.model('basic', '-ndm', 3, '-ndf', 6)  # 3D model with 6 DOF per node
    
#     # Create nodes from structure data
#     with open(structure_data_path, 'r') as f:
#         structure_data = json.load(f)
    
#     for node in structure_data["nodes"]:
#         node_id = node["id"]
#         x, y, z = node["x"], node["y"], node["z"]
#         ops.node(node_id, x, y, z)
    
#     # Define materials and sections (simplified example)
#     # In a real application, you would define proper sections based on your requirements
#     E = 200000.0  # Young's modulus (MPa)
#     nu = 0.3      # Poisson's ratio
#     G = E / (2 * (1 + nu))  # Shear modulus
    
#     # Example section properties - would be different for columns and beams
#     # Section 1 (for columns)
#     ops.section('Elastic', 1, E, 0.01, 1.0e-4, 1.0e-4, G, 1.0e-6)
    
#     # Section 2 (for beams)
#     ops.section('Elastic', 2, E, 0.008, 8.0e-5, 8.0e-5, G, 8.0e-7)
    
#     # Create the structural model elements
#     shapes, memberLengths = create_structural_model(structure_data_path, member_section_mapping_path)
    
#     print("Model created successfully!")
#     print(f"Number of members: {len(memberLengths)}")
#     print(f"Member lengths: {memberLengths}")
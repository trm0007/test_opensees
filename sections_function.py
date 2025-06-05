from import_ import *
from units import *


def draw_rc_section10(sec_tag, core_tag, cover_tag, steel_tag, H, B, cover_H, cover_B, offset,
                    n_bars_top, dia_top, n_bars_bot, dia_bot,
                    n_bars_secondary_top, dia_sec_top, n_bars_secondary_bot, dia_sec_bot,
                    n_bars_int, dia_int, IMAGE_FOLDER):
    """
    Draw a rectangular RC section with visualization.
    
    Parameters:
    -----------
    sec_tag : int
        Tag for the section
    core_tag : int
        Material tag for core concrete
    cover_tag : int
        Material tag for cover concrete
    steel_tag : int
        Material tag for reinforcement steel
    H : float
        Height of the section
    B : float
        Width of the section
    cover_H : float
        Cover thickness in H direction
    cover_B : float
        Cover thickness in B direction
    offset : float
        Offset for secondary reinforcement from primary reinforcement
    n_bars_top : int
        Number of primary reinforcement bars at top
    dia_top : float
        Diameter of primary reinforcement bars at top
    n_bars_bot : int
        Number of primary reinforcement bars at bottom
    dia_bot : float
        Diameter of primary reinforcement bars at bottom
    n_bars_secondary_top : int
        Number of secondary reinforcement bars at top
    dia_sec_top : float
        Diameter of secondary reinforcement bars at top
    n_bars_secondary_bot : int
        Number of secondary reinforcement bars at bottom
    dia_sec_bot : float
        Diameter of secondary reinforcement bars at bottom
    n_bars_int : int
        Number of intermediate reinforcement bars
    dia_int : float
        Diameter of intermediate reinforcement bars
    """
    # Calculate coordinates for core and cover
    core_y = H / 2 - cover_H
    core_z = B / 2 - cover_B
    cover_y = H / 2
    cover_z = B / 2

    # # Create the OpenSees fiber section
    # ops.section('Fiber', sec_tag, '-GJ', 1.0e9)
    
    # # Generate a quadrilateral shaped patch (core patch)
    # ops.patch('quad', core_tag,
    #       10, 10,  # nfIJ, nfJK
    #       *[-core_y, core_z],
    #       *[-core_y, -core_z],
    #       *[core_y, -core_z],
    #       *[core_y, core_z])
    
    # # Define the four cover patches
    # ops.patch('quad', cover_tag,
    #       2, 10,
    #       *[-cover_y, cover_z],
    #       *[-core_y, core_z],
    #       *[core_y, core_z],
    #       *[cover_y, cover_z])    
    # ops.patch('quad', cover_tag,
    #       2, 10,
    #       *[-core_y, -core_z],
    #       *[-cover_y, -cover_z],
    #       *[cover_y, -cover_z],
    #       *[core_y, -core_z])    
    # ops.patch('quad', cover_tag,
    #       10, 2,
    #       *[-cover_y, cover_z],
    #       *[-cover_y, -cover_z],
    #       *[-core_y, -core_z],
    #       *[-core_y, core_z])     
    # ops.patch('quad', cover_tag,
    #       10, 2,
    #       *[core_y, core_z],
    #       *[core_y, -core_z],
    #       *[cover_y, -cover_z],
    #       *[cover_y, cover_z]) 

    # # Construct straight lines of fibers for reinforcement
    # # Top primary reinforcement
    # ops.layer('straight', steel_tag, n_bars_top, 0.25*np.pi*dia_top**2,
    #       *[core_y, core_z], *[core_y, -core_z])
    
    # # Bottom primary reinforcement
    # ops.layer('straight', steel_tag, n_bars_bot, 0.25*np.pi*dia_bot**2,
    #   *[-core_y, core_z], *[-core_y, -core_z])
    
    # # Secondary top reinforcement
    # ops.layer('straight', steel_tag, n_bars_secondary_top, 0.25*np.pi*dia_sec_top**2,
    #       *[core_y - offset, core_z], *[core_y - offset, -core_z])
    
    # # Secondary bottom reinforcement
    # ops.layer('straight', steel_tag, n_bars_secondary_bot, 0.25*np.pi*dia_sec_bot**2,
    #       *[-core_y + offset, core_z], *[-core_y + offset, -core_z])
    
    # # Intermediate reinforcement (top and bottom)
    # ops.layer('straight', steel_tag, n_bars_int // 2, 0.25*np.pi*dia_int**2,
    #       *[-core_y, core_z], *[core_y, core_z])
    # ops.layer('straight', steel_tag, n_bars_int // 2, 0.25*np.pi*dia_int**2,
    #       *[-core_y, -core_z], *[core_y, -core_z])

    # VISUALIZATION
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Draw the concrete patches
    # Core concrete
    core_rect = plt.Rectangle((-core_y, -core_z), 2*core_y, 2*core_z, 
                            facecolor='darkgray', edgecolor='black', label='Core Concrete')
    ax.add_patch(core_rect)
    
    # Cover concrete - top
    top_cover = plt.Rectangle((-cover_y, core_z), 2*cover_y, cover_H, 
                            facecolor='lightgray', edgecolor='black', label='Cover Concrete')
    ax.add_patch(top_cover)
    
    # Cover concrete - bottom
    bottom_cover = plt.Rectangle((-cover_y, -cover_z), 2*cover_y, cover_B, 
                               facecolor='lightgray', edgecolor='black')
    ax.add_patch(bottom_cover)
    
    # Cover concrete - left
    left_cover = plt.Rectangle((-cover_y, -core_z), cover_B, 2*core_z, 
                             facecolor='lightgray', edgecolor='black')
    ax.add_patch(left_cover)
    
    # Cover concrete - right
    right_cover = plt.Rectangle((core_y, -core_z), cover_B, 2*core_z, 
                              facecolor='lightgray', edgecolor='black')
    ax.add_patch(right_cover)
    
    # Function to distribute rebars evenly along a line
    def distribute_bars(x1, y1, x2, y2, n_bars):
        if n_bars == 1:
            return [(0.5*(x1+x2), 0.5*(y1+y2))]
        
        points = []
        for i in range(n_bars):
            t = i / (n_bars - 1) if n_bars > 1 else 0.5
            x = x1 + t * (x2 - x1)
            y = y1 + t * (y2 - y1)
            points.append((x, y))
        return points
    
    # Draw reinforcement bars
    # Primary top reinforcement
    top_bar_positions = distribute_bars(core_y, core_z, core_y, -core_z, n_bars_top)
    for x, y in top_bar_positions:
        bar = plt.Circle((x, y), dia_top/2, color='black', 
                       label='Primary Reinforcement' if (x, y) == top_bar_positions[0] else None)
        ax.add_patch(bar)
    
    # Primary bottom reinforcement
    bottom_bar_positions = distribute_bars(-core_y, core_z, -core_y, -core_z, n_bars_bot)
    for x, y in bottom_bar_positions:
        bar = plt.Circle((x, y), dia_bot/2, color='black')
        ax.add_patch(bar)
    
    # Secondary top reinforcement
    sec_top_positions = distribute_bars(core_y-offset, core_z, core_y-offset, -core_z, n_bars_secondary_top)
    for x, y in sec_top_positions:
        bar = plt.Circle((x, y), dia_sec_top/2, color='blue', 
                       label='Secondary Reinforcement' if (x, y) == sec_top_positions[0] else None)
        ax.add_patch(bar)
    
    # Secondary bottom reinforcement
    sec_bottom_positions = distribute_bars(-core_y+offset, core_z, -core_y+offset, -core_z, n_bars_secondary_bot)
    for x, y in sec_bottom_positions:
        bar = plt.Circle((x, y), dia_sec_bot/2, color='blue')
        ax.add_patch(bar)
    
    # Intermediate reinforcement (top and bottom)
    int_top_positions = distribute_bars(-core_y, core_z, core_y, core_z, n_bars_int//2)
    for x, y in int_top_positions:
        bar = plt.Circle((x, y), dia_int/2, color='red', 
                       label='Intermediate Reinforcement' if (x, y) == int_top_positions[0] else None)
        ax.add_patch(bar)
    
    int_bottom_positions = distribute_bars(-core_y, -core_z, core_y, -core_z, n_bars_int//2)
    for x, y in int_bottom_positions:
        bar = plt.Circle((x, y), dia_int/2, color='red')
        ax.add_patch(bar)
    
    # Set plot properties
    ax.set_aspect('equal')
    ax.set_xlim(-cover_y*1.2, cover_y*1.2)
    ax.set_ylim(-cover_z*1.2, cover_z*1.2)
    ax.set_title(f'Rectangular RC Section (tag: {sec_tag})')
    ax.set_xlabel('y-axis')
    ax.set_ylabel('z-axis')
    ax.grid(True)
    ax.legend(loc='upper right')
    
    section_folder = os.path.join(IMAGE_FOLDER, "sections")
    os.makedirs(section_folder, exist_ok=True)
    path = os.path.join(section_folder, f"rect_section_{sec_tag}.png")
    plt.savefig(path, bbox_inches='tight', dpi=300)
    plt.close()

    
    # plt.show()
    return fig, ax

def draw_circular_rc_section10(sec_tag, core_tag, cover_tag, steel_tag,
                             D_Sec, cover_Sec, num_Bars_Sec, bar_dia_Sec, IMAGE_FOLDER,
                             ri=0.0, nf_Core_R=8, nf_Core_T=8, nf_Cover_R=4, nf_Cover_T=8):
    """
    Draw a circular RC section with visualization.
    
    Parameters:
    -----------
    sec_tag : int
        Tag for the section
    core_tag : int
        Material tag for core concrete
    cover_tag : int
        Material tag for cover concrete
    steel_tag : int
        Material tag for reinforcement steel
    D_Sec : float
        Outer diameter of the section
    cover_Sec : float
        Cover thickness to reinforcing steel
    num_Bars_Sec : int
        Number of longitudinal reinforcement bars
    bar_dia_Sec : float
        Diameter of longitudinal reinforcement bars
    ri : float, optional
        Inner radius (for hollow sections), default = 0.0
    nf_Core_R : int, optional
        Number of radial divisions in core, default = 8
    nf_Core_T : int, optional
        Number of theta divisions in core, default = 8
    nf_Cover_R : int, optional
        Number of radial divisions in cover, default = 4
    nf_Cover_T : int, optional
        Number of theta divisions in cover, default = 8
    """
    # Calculate key dimensions
    ro = D_Sec / 2  # outer radius
    rc = ro - cover_Sec  # core radius
    bar_area_Sec = 0.25 * np.pi * bar_dia_Sec**2  # bar area
    
    # constructs a FiberSection 
    # ops.section('Fiber', sec_tag, '-GJ', 1.0e9)
    
    # # Define the core patch (confined concrete)
    # ops.patch('circ', core_tag, nf_Core_T, nf_Core_R, 0, 0, ri, rc, 0, 360)
    
    # # Define the cover patch (unconfined concrete)
    # ops.patch('circ', cover_tag, nf_Cover_T, nf_Cover_R, 0, 0, rc, ro, 0, 360)
    
    # # Define the reinforcing layer
    # theta = 360.0 / num_Bars_Sec  # angle increment between bars
    # ops.layer('circ', steel_tag, num_Bars_Sec, bar_area_Sec, 0, 0, rc, 0, 360)
    
    # Create figure for visualization
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Draw concrete patches
    # Core concrete (in gray)
    core_circle = plt.Circle((0, 0), rc, color='darkgray', ec='black', label='Core Concrete')
    ax.add_patch(core_circle)
    
    # Cover concrete (in light gray)
    cover_ring = plt.Circle((0, 0), ro, color='lightgray', ec='black', label='Cover Concrete')
    ax.add_patch(cover_ring)
    ax.add_patch(core_circle)  # Add core again to create ring effect
    
    # If there's a hollow section
    if ri > 0:
        hollow_circle = plt.Circle((0, 0), ri, color='white', ec='black')
        ax.add_patch(hollow_circle)
    
    # Draw reinforcement bars
    angles = np.linspace(0, 360, num_Bars_Sec, endpoint=False)
    rads = np.radians(angles)
    x_coords = rc * np.cos(rads)
    y_coords = rc * np.sin(rads)
    
    # Draw bars as circles
    for i, (x, y) in enumerate(zip(x_coords, y_coords)):
        bar_circle = plt.Circle((x, y), bar_dia_Sec/2, color='black', 
                              label='Reinforcement' if i == 0 else None)
        ax.add_patch(bar_circle)
    
    # Set plot properties
    ax.set_aspect('equal')
    ax.set_xlim(-ro*1.2, ro*1.2)
    ax.set_ylim(-ro*1.2, ro*1.2)
    ax.set_title(f'Circular RC Section (tag: {sec_tag})')
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.grid(True)
    ax.legend(loc='upper right')

    section_folder = os.path.join(IMAGE_FOLDER, "sections")
    os.makedirs(section_folder, exist_ok=True)
    path = os.path.join(section_folder, f"circ_section_{sec_tag}.png")
    plt.savefig(path, bbox_inches='tight', dpi=300)
    plt.close()

    # plt.show()
    return fig, ax

def draw_L_rc_section10(sec_tag, core_tag, cover_tag, steel_tag, H, B, t, cover_H, cover_B, 
                       n_bars_vertical_outer, dia_vertical_outer, n_bars_horizontal_outer, dia_horizontal_outer,
                       n_bars_vertical_inner, dia_vertical_inner, n_bars_horizontal_inner, dia_horizontal_inner,
                       corner_bar_dia, IMAGE_FOLDER):
    """
    Draw an L-shaped RC section with two layers of reinforcement and proper dimension annotations.
    
    Parameters:
    -----------
    sec_tag : int
        Section tag
    core_tag : int
        Core concrete material tag
    cover_tag : int
        Cover concrete material tag
    steel_tag : int
        Steel material tag
    H, B : float
        Total height and width (mm)
    t : float
        Thickness of legs (mm)
    cover_H, cover_B : float
        Concrete cover dimensions (mm)
    n_bars_* : int
        Number of bars in each layer
    dia_* : float
        Bar diameters (mm)
    corner_bar_dia : float
        Corner bar diameter (mm)
    IMAGE_FOLDER : str
        Output directory
    """
    fig, ax = plt.subplots(figsize=(12, 12))
    
    # Core dimensions
    core_vert_width = t - cover_B
    core_horiz_height = t - cover_H
    
    # Draw concrete - Core (dark gray)
    vert_core = plt.Rectangle((cover_B, cover_H), core_vert_width, H-cover_H, 
                            facecolor='darkgray', edgecolor='black')
    horiz_core = plt.Rectangle((cover_B, cover_H), B-cover_B, core_horiz_height, 
                             facecolor='darkgray', edgecolor='black')
    ax.add_patch(vert_core)
    ax.add_patch(horiz_core)
    
    # Draw concrete - Cover (light gray)
    # Vertical leg covers
    left_cover = plt.Rectangle((0, 0), cover_B, H, facecolor='lightgray', edgecolor='black')
    right_cover = plt.Rectangle((t-cover_B, 0), cover_B, H, facecolor='lightgray', edgecolor='black')
    # Horizontal leg covers
    bottom_cover = plt.Rectangle((0, 0), B, cover_H, facecolor='lightgray', edgecolor='black')
    top_cover = plt.Rectangle((0, t-cover_H), B, cover_H, facecolor='lightgray', edgecolor='black')
    
    ax.add_patch(left_cover)
    ax.add_patch(right_cover)
    ax.add_patch(bottom_cover)
    ax.add_patch(top_cover)
    
    # Calculate bar positions (two layers)
    def calc_bar_positions(start, end, n_bars, offset_from_edge):
        if n_bars == 1:
            return [start + (end-start)/2]
        return [start + offset_from_edge + i*(end-start-2*offset_from_edge)/(n_bars-1) for i in range(n_bars)]
    
    # Vertical bars - Outer layer (near outer face)
    outer_vert_x = t - cover_B/2
    outer_vert_y = calc_bar_positions(cover_H, H, n_bars_vertical_outer, dia_vertical_outer)
    
    # Vertical bars - Inner layer (near inner face)
    inner_vert_x = cover_B + dia_horizontal_inner
    inner_vert_y = calc_bar_positions(cover_H, H, n_bars_vertical_inner, dia_vertical_inner)
    
    # Horizontal bars - Outer layer (near top face)
    outer_horiz_y = t - cover_H/2
    outer_horiz_x = calc_bar_positions(cover_B, B, n_bars_horizontal_outer, dia_horizontal_outer)
    
    # Horizontal bars - Inner layer (near bottom face)
    inner_horiz_y = cover_H + dia_vertical_inner
    inner_horiz_x = calc_bar_positions(cover_B, B, n_bars_horizontal_inner, dia_horizontal_inner)
    
    # Plot reinforcement
    # Outer layer (black)
    for y in outer_vert_y:
        ax.add_patch(plt.Circle((outer_vert_x, y), dia_vertical_outer/2, color='black'))
    for x in outer_horiz_x:
        ax.add_patch(plt.Circle((x, outer_horiz_y), dia_horizontal_outer/2, color='black'))
    
    # Inner layer (red)
    for y in inner_vert_y:
        ax.add_patch(plt.Circle((inner_vert_x, y), dia_vertical_inner/2, color='red'))
    for x in inner_horiz_x:
        ax.add_patch(plt.Circle((x, inner_horiz_y), dia_horizontal_inner/2, color='red'))
    
    # Corner bars (blue)
    ax.add_patch(plt.Circle((inner_vert_x, inner_horiz_y), corner_bar_dia/2, color='blue'))
    ax.add_patch(plt.Circle((outer_vert_x, outer_horiz_y), corner_bar_dia/2, color='blue'))
    
    # Add dimension lines
    def add_dimension(ax, start, end, text, offset, orientation='horizontal'):
        if orientation == 'horizontal':
            ax.plot([start[0], end[0]], [start[1], end[1]], 'k-', lw=0.5)
            ax.text((start[0]+end[0])/2, start[1]+offset, text, ha='center')
        else:
            ax.plot([start[0], end[0]], [start[1], end[1]], 'k-', lw=0.5)
            ax.text(start[0]+offset, (start[1]+end[1])/2, text, va='center')
    
    # Add dimensions
    dim_offset = 50
    add_dimension(ax, (0, -dim_offset), (B, -dim_offset), f"{B}mm", -10)
    add_dimension(ax, (-dim_offset, 0), (-dim_offset, H), f"{H}mm", -10)
    add_dimension(ax, (0, t+dim_offset), (t, t+dim_offset), f"{t}mm", 10)
    
    # Finalize plot
    ax.set_aspect('equal')
    ax.set_xlim(-100, B+100)
    ax.set_ylim(-100, H+100)
    ax.set_title(f'L-Shaped RC Section (Tag: {sec_tag})\nDouble Layer Reinforcement')
    ax.legend([plt.Rectangle((0,0),1,1,fc='darkgray'), 
               plt.Rectangle((0,0),1,1,fc='lightgray'),
               plt.Circle((0,0),5,fc='black'), 
               plt.Circle((0,0),5,fc='red'),
               plt.Circle((0,0),5,fc='blue')],
              ['Core Concrete', 'Cover Concrete', 
               'Outer Layer Steel', 'Inner Layer Steel',
               'Corner Bars'], loc='upper right')
    
    # Save figure
    os.makedirs(os.path.join(IMAGE_FOLDER, "sections"), exist_ok=True)
    plt.savefig(os.path.join(IMAGE_FOLDER, "sections", f"L_section_{sec_tag}_double_layer.png"), 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    return fig, ax


  

# # Example call with parameters from the TCL file
# draw_circular_rc_section10(
#     sec_tag=1, 
#     core_tag=1, 
#     cover_tag=2, 
#     steel_tag=3,
#     D_Sec=60,           # 5 ft diameter converted to inches
#     cover_Sec=5,        # 5 inches cover
#     num_Bars_Sec=16,    # 16 bars
#     bar_dia_Sec=1.69,   # equivalent to 2.25 sq in area
#     ri=0.0,             # solid section
#     nf_Core_R=8,
#     nf_Core_T=8,
#     nf_Cover_R=4,
#     nf_Cover_T=8
# )

# # Example call for a hollow circular section
# draw_circular_rc_section10(
#     sec_tag=1, 
#     core_tag=1, 
#     cover_tag=2, 
#     steel_tag=3,
#     D_Sec=60,           # 5 ft diameter converted to inches
#     cover_Sec=5,        # 5 inches cover
#     num_Bars_Sec=16,    # 16 bars
#     bar_dia_Sec=1.69,   # equivalent to 2.25 sq in area
#     ri=15.0,            # hollow section with inner radius of 15 inches
#     nf_Core_R=8,
#     nf_Core_T=8,
#     nf_Cover_R=4,
#     nf_Cover_T=8
# )

# draw_rc_section10(
#     sec_tag=1, core_tag=1, cover_tag=2, steel_tag=3,
#     H=600, B=300, cover_H=40, cover_B=40, offset=25,
#     n_bars_top=3, dia_top=16,
#     n_bars_bot=4, dia_bot=20,
#     n_bars_secondary_top=5, dia_sec_top=16,
#     n_bars_secondary_bot=6, dia_sec_bot=25,
#     n_bars_int=6, dia_int=16
# )

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
# draw_L_rc_section10({
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
# })



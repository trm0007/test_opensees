from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch

# Create PDF document
doc = SimpleDocTemplate("lateral_load_analysis_report.pdf", pagesize=letter)
styles = getSampleStyleSheet()

# Create custom styles
title_style = ParagraphStyle(
    'CustomTitle',
    parent=styles['Title'],
    fontSize=16,
    spaceAfter=30,
    alignment=1  # Center alignment
)

heading_style = ParagraphStyle(
    'CustomHeading',
    parent=styles['Heading1'],
    fontSize=14,
    spaceAfter=12
)

# Content list
content = []

# Title
content.append(Paragraph("Structural Analysis Report", title_style))
content.append(Paragraph("Lateral Load Distribution with Mesh-Based Walls", title_style))
content.append(Spacer(1, 20))

# Theoretical Background
content.append(Paragraph("Theoretical Background", heading_style))
background_text = """
• Mesh-based wall modeling<br/>
• Rotated element handling<br/>
• Complete stiffness coupling<br/>
• Displacement-based distribution
"""
content.append(Paragraph(background_text, styles['Normal']))
content.append(PageBreak())

# Key Equations
content.append(Paragraph("Key Equations", heading_style))
equations_text = """
<b>Shoelace Formula:</b><br/>
A = 1/2 * |Σ(x<sub>i</sub>*y<sub>i+1</sub> - x<sub>i+1</sub>*y<sub>i</sub>)|<br/><br/>

<b>Stiffness Transformation:</b><br/>
K<sub>xx</sub> = K<sub>x_local</sub>*cos²θ + K<sub>y_local</sub>*sin²θ<br/>
K<sub>yy</sub> = K<sub>x_local</sub>*sin²θ + K<sub>y_local</sub>*cos²θ<br/>
K<sub>xy</sub> = (K<sub>x_local</sub> - K<sub>y_local</sub>)*sinθ*cosθ
"""
content.append(Paragraph(equations_text, styles['Normal']))
content.append(PageBreak())

# Example Calculations
content.append(Paragraph("Example Calculations", heading_style))
example_text = """
<b>Example 1: Quadrilateral Mesh Wall</b><br/>
Area = 21600 in² (15 ft²)<br/>
Centroid = 100 in<br/>
I<sub>x</sub> = 5,760,000 in⁴
"""
content.append(Paragraph(example_text, styles['Normal']))

# Build PDF
doc.build(content)
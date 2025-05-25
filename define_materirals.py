from import_ import *
from units import *

def define_materials(materials_list):
    for mat in materials_list:
        ops.uniaxialMaterial(*mat.values())
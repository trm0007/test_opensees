import json
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import openseespy.opensees as ops
import sys
from math import sqrt

ops.wipe()

# define a 3D model
ops.model("basic","-ndm",3,"-ndf",6)

# the response spectrum function
Tn = [0.0, 0.06, 0.1, 0.12, 0.18, 0.24, 0.3, 0.36, 0.4, 0.42, 
	0.48, 0.54, 0.6, 0.66, 0.72, 0.78, 0.84, 0.9, 0.96, 1.02, 
	1.08, 1.14, 1.2, 1.26, 1.32, 1.38, 1.44, 1.5, 1.56, 1.62, 
	1.68, 1.74, 1.8, 1.86, 1.92, 1.98, 2.04, 2.1, 2.16, 2.22, 
	2.28, 2.34, 2.4, 2.46, 2.52, 2.58, 2.64, 2.7, 2.76, 2.82, 
	2.88, 2.94, 3.0, 3.06, 3.12, 3.18, 3.24, 3.3, 3.36, 3.42, 
	3.48, 3.54, 3.6, 3.66, 3.72, 3.78, 3.84, 3.9, 3.96, 4.02, 
	4.08, 4.14, 4.2, 4.26, 4.32, 4.38, 4.44, 4.5, 4.56, 4.62, 
	4.68, 4.74, 4.8, 4.86, 4.92, 4.98, 5.04, 5.1, 5.16, 5.22, 
	5.28, 5.34, 5.4, 5.46, 5.52, 5.58, 5.64, 5.7, 5.76, 5.82, 
	5.88, 5.94, 6.0]
Sa = [1.9612, 3.72628, 4.903, 4.903, 4.903, 4.903, 4.903, 4.903, 4.903, 4.6696172, 
	4.0861602, 3.6321424, 3.2683398, 2.971218, 2.7241068, 2.5142584, 2.3348086, 2.1788932, 2.0425898, 1.9229566, 
	1.8160712, 1.7199724, 1.6346602, 1.5562122, 1.485609, 1.4208894, 1.3620534, 1.3071398, 1.2571292, 1.211041, 
	1.166914, 1.1267094, 1.0894466, 1.054145, 1.0217852, 0.990406, 0.960988, 0.9335312, 0.9080356, 0.8835206, 
	0.8599862, 0.838413, 0.8168398, 0.7972278, 0.7785964, 0.759965, 0.7432948, 0.7266246, 0.710935, 0.6952454, 
	0.6805364, 0.666808, 0.6540602, 0.6285646, 0.6040496, 0.5814958, 0.5609032, 0.5403106, 0.5206986, 0.5030478, 
	0.485397, 0.4697074, 0.4540178, 0.4393088, 0.4255804, 0.411852, 0.3991042, 0.3863564, 0.3755698, 0.3638026, 
	0.353016, 0.34321, 0.333404, 0.3245786, 0.3157532, 0.3069278, 0.2981024, 0.2902576, 0.2833934, 0.2755486, 
	0.2686844, 0.2618202, 0.254956, 0.2490724, 0.2431888, 0.2373052, 0.2314216, 0.2265186, 0.220635, 0.215732, 
	0.210829, 0.205926, 0.2020036, 0.1971006, 0.1931782, 0.1892558, 0.1853334, 0.181411, 0.1774886, 0.1735662, 
	0.1706244, 0.166702, 0.1637602]
ops.timeSeries("Path",1,"-time", *Tn, "-values", *Sa)

# a uniaxial material for transverse shear
ops.uniaxialMaterial("Elastic",2,938000000.0)

# the elastic beam section and aggregator
ops.section("Elastic",1,30000000000.0,0.09,0.0006749999999999999,0.0006749999999999999,12500000000.0,0.0011407499999999994)
ops.section("Aggregator",3,2,"Vy",2,"Vz","-section",1)

# nodes and masses
ops.node(1,0,0,0)
ops.node(2,0,0,3,"-mass",200,200,200,0,0,0)
ops.node(3,4,0,3,"-mass",200,200,200,0,0,0)
ops.node(4,4,0,0)
ops.node(5,0,0,6,"-mass",200,200,200,0,0,0)
ops.node(6,4,0,6,"-mass",200,200,200,0,0,0)
ops.node(7,4,3,6,"-mass",200,200,200,0,0,0)
ops.node(8,0,3,6,"-mass",200,200,200,0,0,0)
ops.node(9,0,3,3,"-mass",200,200,200,0,0,0)
ops.node(10,0,3,0)
ops.node(11,4,3,3,"-mass",200,200,200,0,0,0)
ops.node(12,4,3,0)
ops.node(13,2,1.5,6)
ops.node(14,2,1.5,3)

# beam elements
ops.beamIntegration("Lobatto", 1, 3, 5)
# beam_column_elements forceBeamColumn
# Geometric transformation command
ops.geomTransf("Linear", 1, 1.0, 0.0, -0.0)
ops.element("forceBeamColumn", 1, 1, 2, 1, 1)
# Geometric transformation command
ops.geomTransf("Linear", 2, 0.0, 0.0, 1.0)
ops.element("forceBeamColumn", 2, 2, 3, 2, 1)
# Geometric transformation command
ops.geomTransf("Linear", 3, 1.0, 0.0, -0.0)
ops.element("forceBeamColumn", 3, 4, 3, 3, 1)
# Geometric transformation command
ops.geomTransf("Linear", 4, 1.0, 0.0, -0.0)
ops.element("forceBeamColumn", 4, 2, 5, 4, 1)
# Geometric transformation command
ops.geomTransf("Linear", 5, 0.0, 0.0, 1.0)
ops.element("forceBeamColumn", 5, 5, 6, 5, 1)
# Geometric transformation command
ops.geomTransf("Linear", 6, 0.0, 0.0, 1.0)
ops.element("forceBeamColumn", 6, 7, 6, 6, 1)
# Geometric transformation command
ops.geomTransf("Linear", 7, 0.0, 0.0, 1.0)
ops.element("forceBeamColumn", 7, 8, 7, 7, 1)
# Geometric transformation command
ops.geomTransf("Linear", 8, 0.0, 0.0, 1.0)
ops.element("forceBeamColumn", 8, 9, 2, 8, 1)
# Geometric transformation command
ops.geomTransf("Linear", 9, 0.0, 0.0, 1.0)
ops.element("forceBeamColumn", 9, 8, 5, 9, 1)
# Geometric transformation command
ops.geomTransf("Linear", 10, 1.0, 0.0, -0.0)
ops.element("forceBeamColumn", 10, 10, 9, 10, 1)
# Geometric transformation command
ops.geomTransf("Linear", 11, 1.0, 0.0, -0.0)
ops.element("forceBeamColumn", 11, 3, 6, 11, 1)
# Geometric transformation command
ops.geomTransf("Linear", 12, 1.0, 0.0, -0.0)
ops.element("forceBeamColumn", 12, 11, 7, 12, 1)
# Geometric transformation command
ops.geomTransf("Linear", 13, 0.0, 0.0, 1.0)
ops.element("forceBeamColumn", 13, 11, 3, 13, 1)
# Geometric transformation command
ops.geomTransf("Linear", 14, 0.0, 0.0, 1.0)
ops.element("forceBeamColumn", 14, 9, 11, 14, 1)
# Geometric transformation command
ops.geomTransf("Linear", 15, 1.0, 0.0, -0.0)
ops.element("forceBeamColumn", 15, 12, 11, 15, 1)
# Geometric transformation command
ops.geomTransf("Linear", 16, 1.0, 0.0, -0.0)
ops.element("forceBeamColumn", 16, 9, 8, 16, 1)

# Constraints.sp fix
ops.fix(1, 1, 1, 1, 1, 1, 1)
ops.fix(10, 1, 1, 1, 1, 1, 1)
ops.fix(4, 1, 1, 1, 1, 1, 1)
ops.fix(12, 1, 1, 1, 1, 1, 1)
ops.fix(13, 0, 0, 1, 1, 1, 0)
ops.fix(14, 0, 0, 1, 1, 1, 0)

# Constraints.mp rigidDiaphragm
ops.rigidDiaphragm(3, 14, 2, 3, 9, 11)
ops.rigidDiaphragm(3, 13, 5, 6, 7, 8)

# define some analysis settings
ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("UmfPack")
ops.test("NormUnbalance", 0.0001, 10)
ops.algorithm("Linear")
ops.integrator("LoadControl", 0.0)
ops.analysis("Static")
import opsvis as opsv
opsv.plot_model()
# plt.show()
# run the eigenvalue analysis with 7 modes
# and obtain the eigenvalues
eigs = ops.eigen("-genBandArpack", 7)

# compute the modal properties
ops.modalProperties("-print", "-file", "ModalReport.txt", "-unorm")

# define a recorder for the (use a higher precision otherwise the results
# won't match with those obtained from eleResponse)
filename = 'ele_1_sec_1.txt'
ops.recorder('Element', '-file', filename, '-closeOnWrite', '-precision', 16, '-ele', 1, 'section', '1', 'force')

# some settings for the response spectrum analysis
tsTag = 1 # use the timeSeries 1 as response spectrum function
direction = 1 # excited DOF = Ux

# currently we use same damping for each mode
dmp = [0.05]*len(eigs)
# we don't want to scale some modes...
scalf = [1.0]*len(eigs)
# CQC function
def CQC(mu, lambdas, dmp, scalf):
	u = 0.0
	ne = len(lambdas)
	for i in range(ne):
		for j in range(ne):
			di = dmp[i]
			dj = dmp[j]
			bij = lambdas[i]/lambdas[j]
			rho = ((8.0*sqrt(di*dj)*(di+bij*dj)*(bij**(3.0/2.0))) /
				((1.0-bij**2.0)**2.0 + 4.0*di*dj*bij*(1.0+bij**2.0) + 
				4.0*(di**2.0 + dj**2.0)*bij**2.0))
			u += scalf[i]*mu[i] * scalf[j]*mu[j] * rho;
	return sqrt(u)

# ========================================================================
# TEST 00
# run a response spectrum analysis for each mode.
# then do modal combination in post-processing.
# Use Tn and Sa lists
# ========================================================================
ops.responseSpectrumAnalysis(direction, '-Tn', *Tn, '-Sa', *Sa)

# read the My values [3rd column] for each step 
# (1 for each mode, they are section forces associated to each modal displacement)
My = []
with open(filename, 'r') as f:
	lines = f.read().split('\n')
	for line in lines:
		if len(line) > 0:
			tokens = line.split(' ')
			My.append(float(tokens[2]))

# post process the results doing the CQC modal combination for the My response (3rd column in section forces)
MyCQC = CQC(My, eigs, dmp, scalf)

print('\n\nTEST 00:\nRun a Response Spectrum Analysis for all modes.')
print('Do CQC combination in post processing.')
print('Use Tn and Sa lists.\n')
print('{0: >10}{1: >15}'.format('Mode', 'My'))
for i in range(len(eigs)):
	print('{0: >10}{1: >15f}'.format(i+1, My[i]))
print('{0: >10}{1: >15f}'.format('CQC', MyCQC))

# ========================================================================
# TEST 01
# run a response spectrum analysis for each mode.
# then do modal combination in post-processing.
# Use a Path timeSeries to store the Tn-Sa pairs
# ========================================================================
ops.responseSpectrumAnalysis(tsTag, direction)

# read the My values [3rd column] for each step 
# (1 for each mode, they are section forces associated to each modal displacement)
My = []
with open(filename, 'r') as f:
	lines = f.read().split('\n')
	for line in lines:
		if len(line) > 0:
			tokens = line.split(' ')
			My.append(float(tokens[2]))

# post process the results doing the CQC modal combination for the My response (3rd column in section forces)
MyCQC = CQC(My, eigs, dmp, scalf)

print('\n\nTEST 01:\nRun a Response Spectrum Analysis for all modes.')
print('Do CQC combination in post processing.')
print('Use a Path timeSeries to store Tn-Sa pairs.\n')
print('{0: >10}{1: >15}'.format('Mode', 'My'))
for i in range(len(eigs)):
	print('{0: >10}{1: >15f}'.format(i+1, My[i]))
print('{0: >10}{1: >15f}'.format('CQC', MyCQC))

# ========================================================================
# TEST 02
# run a response spectrum analysis mode-by-mode.
# grab results during the loop, not using the recorder
# then do modal combination in post-processing.
# ========================================================================
ops.remove('recorder', 0)
My = []
for i in range(len(eigs)):
	ops.responseSpectrumAnalysis(direction, '-Tn', *Tn, '-Sa', *Sa, '-mode', i+1)
	force = ops.eleResponse(1, 'section', '1', 'force')
	My.append(force[2])

# post process the results doing the CQC modal combination for the My response (3rd column in section forces)
MyCQC = CQC(My, eigs, dmp, scalf)

print('\n\nTEST 02:\nRun a Response Spectrum Analysis mode-by-mode.')
print('Grab results during the loop and do CQC combination with them.\n')
print('{0: >10}{1: >15}'.format('Mode', 'My'))
for i in range(len(eigs)):
	print('{0: >10}{1: >15f}'.format(i+1, My[i]))
print('{0: >10}{1: >15f}'.format('CQC', MyCQC))

# Gráficas de modos
for i in range(1, 7 + 1):
    opsv.plot_mode_shape(i, endDispFlag=0, fig_wi_he=(18, 18), node_supports=False)
    # plt.show()

# Plot Response Spectrum
plt.figure(figsize=(10, 6))
plt.plot(Tn, Sa, 'b-', linewidth=2, label='Design Spectrum')
plt.xlabel('Period (s)')
plt.ylabel('Spectral Acceleration (g)')
plt.title('Response Spectrum')
plt.legend()
plt.grid(True, alpha=0.3)
# plt.show()



ops.recorder('Element', '-file', 'ele_1_sec_1.txt', '-precision', 16, '-ele', 1, 'section', '1', 'force')

# Create recorders for all elements
ops.recorder('Element', '-file', 'ElementForces.txt', '-closeOnWrite', 
             '-precision', 16, '-eleRange', 1, 1, 
             'section', '1', 'force')

# Run RSA for all modes (your existing code)
for i in range(len(eigs)):
    ops.responseSpectrumAnalysis(direction, '-Tn', *Tn, '-Sa', *Sa, '-mode', i+1)
    ops.analyze(1)

def extract_and_combine_forces(Tn, Sa, direction, eigs, dmp, scalf, json_filepath='RSA_Forces.json', csv_filepath='RSA_Forces.csv'):
    """
    Extract all member forces from RSA and perform CQC combination
    
    Args:
        Tn (list): List of periods for response spectrum
        Sa (list): List of spectral accelerations
        direction (int): Excitation direction (1=X, 2=Y, 3=Z)
        eigs (list): Eigenvalues from modal analysis
        dmp (list): Damping ratios for each mode
        scalf (list): Scaling factors for each mode
        json_filepath (str): Path to save JSON results
        csv_filepath (str): Path to save CSV results
    """
    # Initialize dictionaries to store all forces
    modal_forces = {
        'P': {}, 'Vy': {}, 'Vz': {}, 
        'T': {}, 'My': {}, 'Mz': {}
    }
    element_tags = ops.getEleTags()
    
    # Extract forces for each mode
    for mode in range(1, len(eigs)+1):
        ops.responseSpectrumAnalysis(direction, '-Tn', *Tn, '-Sa', *Sa, '-mode', mode)
        
        # Initialize mode entries
        for force_type in modal_forces:
            modal_forces[force_type][mode] = {}
        
        # Get forces for all elements
        for ele_tag in element_tags:
            forces = ops.eleResponse(ele_tag, 'section', '1', 'force')
            modal_forces['P'][mode][ele_tag] = forces[0]
            modal_forces['Vy'][mode][ele_tag] = forces[1]
            modal_forces['Vz'][mode][ele_tag] = forces[2]
            modal_forces['T'][mode][ele_tag] = forces[3]
            modal_forces['My'][mode][ele_tag] = forces[4]
            modal_forces['Mz'][mode][ele_tag] = forces[5]
    
    # Perform CQC combination
    cqc_forces = {
        'P': {}, 'Vy': {}, 'Vz': {}, 
        'T': {}, 'My': {}, 'Mz': {}
    }
    
    for ele_tag in element_tags:
        # Extract modal forces for this element
        P_modes = [modal_forces['P'][m][ele_tag] for m in range(1, len(eigs)+1)]
        Vy_modes = [modal_forces['Vy'][m][ele_tag] for m in range(1, len(eigs)+1)]
        Vz_modes = [modal_forces['Vz'][m][ele_tag] for m in range(1, len(eigs)+1)]
        T_modes = [modal_forces['T'][m][ele_tag] for m in range(1, len(eigs)+1)]
        My_modes = [modal_forces['My'][m][ele_tag] for m in range(1, len(eigs)+1)]
        Mz_modes = [modal_forces['Mz'][m][ele_tag] for m in range(1, len(eigs)+1)]
        
        # Perform CQC combination
        cqc_forces['P'][ele_tag] = CQC(P_modes, eigs, dmp, scalf)
        cqc_forces['Vy'][ele_tag] = CQC(Vy_modes, eigs, dmp, scalf)
        cqc_forces['Vz'][ele_tag] = CQC(Vz_modes, eigs, dmp, scalf)
        cqc_forces['T'][ele_tag] = CQC(T_modes, eigs, dmp, scalf)
        cqc_forces['My'][ele_tag] = CQC(My_modes, eigs, dmp, scalf)
        cqc_forces['Mz'][ele_tag] = CQC(Mz_modes, eigs, dmp, scalf)
    
    # Save to JSON
    with open(json_filepath, 'w') as f:
        json.dump({
            'modal_forces': modal_forces,
            'cqc_forces': cqc_forces,
            'eigenvalues': eigs,
            'damping': dmp,
            'scaling_factors': scalf
        }, f, indent=4)
    

    
    return modal_forces, cqc_forces
modal_forces, cqc_forces = extract_and_combine_forces(
    Tn, Sa, direction, eigs, dmp, scalf,
    json_filepath='results2100.json',

)

# done
ops.wipe()

# coding: UTF-8
# from isqtplot.plotter import Plotter
import numpy as np
from ShapeFunctionGenerator import ShapeFunctionGenerator
from assemble import assemble

E = 1
nu = 0.25
stiffnessTensor = E/(1-nu)*np.array([
	[1, nu, 0],
	[nu, 1, 0],
	[0, 0, (1-nu)/2],
	])

nodex = np.array([0,1,0,1]) # x coord of nodes
nodey = np.array([0,0,1,1]) # y coord of nodes
connectivityArray = np.array([
	[0,1,2],
	[2,1,3],
])



xil = np.array([1./3, 2./15, 2./15, 11./15]);
etal = np.array([1./3, 11./15, 2./15, 2./15]);
weight = np.array([-27, 25, 25, 25])/96;
nint = len(xil)
sfGen = ShapeFunctionGenerator()
# Raw shape functions for each variable
psi,dxi,deta = sfGen.generateShapeFunctionsTri(xil,etal)

assemble(nodex,nodey,connectivityArray,
	psi,dxi,deta,weight,stiffnessTensor)

# A = np.ones([2,2,2])
# f = np.array([3,2])
# B = A*f
# print B[0,0,:]
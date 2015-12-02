# coding: UTF-8
import numpy as np

class ShapeFunctionGenerator:
	dimensions = None
	order = None

	def generateShapeFunctionsTri(self,xi,eta,order=1):
		'''
		Generates shape functions on triangular elements
		input:
			xi - 1d array - xi coordinate
			eta - 1d array - eta coordinate
		output:
			phi - values of shape functions
			dxi - values of derivatives wrt xi
			deta - values of derivatives wrt eta
		'''
		if len(xi) != len(eta): 
			raise ValueError('input arrays must be the same length')
		if order==1:
			psi,dxi,deta = self.generateTri1(xi,eta)
		elif order==2:
			psi,dxi,deta = self.generateTri2(xi,eta)
		else:
			raise ValueError('This order of shape functions is not supported')
		return psi,dxi,deta

	def generateTri1(self,xi,eta):
		'''
		generates 1st order shape functions on triangles
		'''
		nint = len(xi)
		psi = np.zeros([3,nint])
		dxi = np.zeros([3,nint])
		deta = np.zeros([3,nint])
		# Shape functions
		psi[0,:] = 1 - xi - eta;
		psi[1,:] = xi;
		psi[2,:] = eta;
		# Derivatives wrt xi
		dxi[0,:] = - np.ones([1,nint])
		dxi[1,:] =   np.ones([1,nint])
		dxi[2,:] =   np.zeros([1,nint])
		# Derivatives wrt eta
		deta[0,:] = - np.ones([1,nint])
		deta[1,:] =   np.zeros([1,nint])
		deta[2,:] =   np.ones([1,nint])
		return psi,dxi,deta

	def generateTri2(self,xi,eta):
		'''
		generates 2nd order shape functions on triangles
		'''
		pass

	# def assebleShapeFunctions


if __name__ == '__main__':
	xil = np.array([1./3, 2./15, 2./15, 11./15]);
	etal = np.array([1./3, 11./15, 2./15, 2./15]);
	sfGen = ShapeFunctionGenerator()
	sfGen.generateShapeFunctionsTri(xil,etal)
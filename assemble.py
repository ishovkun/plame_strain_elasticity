# coding: UTF-8
import numpy as np

def dot(A,B,axis=2):
	'''
	dot product of each slice of 3d matrices
	shapes: (A1,A2,A3) & (B1,B2,B3)
	the function also accepts:
		(A1,A2) & (B1,B2,B3)
		and (A1,A2,A3) & (B1,B2)
	A2 = B1
	A3 = B3
	'''
	if axis==2:
		n1 = A.shape[0]
		n2 = B.shape[1]
		n3 = B.shape[2]
		if (A.shape[1]!=B.shape[0]):
			if len(A.shape)!=2 and len(B.shape)!=2:
				if A.shape[2]!=B.shape[2]:
					raise ValueError('shapes not aligned')
		product = np.empty([n1,n2,n3])
		if len(A.shape)==2:
			for i in xrange(B.shape[2]):	
				product[:,:,i] = np.dot(A,B[:,:,i])
		elif len(B.shape)==2:
			for i in xrange(A.shape[2]):	
				product[:,:,i] = np.dot(A[:,:,i],B)
		else:
			for i in xrange(A.shape[2]):
				product[:,:,i] = np.dot(A[:,:,i],B[:,:,i])
	else:
		raise NotImplementedError
	return product

def div(dxi,deta,dxidx,dxidy,detadx,detady):
	nint = dxi.shape[1]
	dphidx = dxi*dxidx + deta*detadx
	dphidy = dxi*dxidy + deta*detady
	divphi = np.empty([3,6,nint])
	divphi[0,0,:] = dphidx[0]
	divphi[0,1,:] = np.zeros(nint)
	divphi[0,2,:] = dphidx[1]
	divphi[0,3,:] = np.zeros(nint)
	divphi[0,4,:] = dphidx[2]
	divphi[0,5,:] = np.zeros(nint)

	divphi[1,0,:] = np.zeros(nint)
	divphi[1,1,:] = dphidy[0]
	divphi[1,2,:] = np.zeros(nint)
	divphi[1,3,:] = dphidy[1]
	divphi[1,4,:] = np.zeros(nint)
	divphi[1,5,:] = dphidy[2]
	
	divphi[2,0,:] = dphidy[0]
	divphi[2,1,:] = dphidx[0]
	divphi[2,2,:] = dphidy[1]
	divphi[2,3,:] = dphidx[1]
	divphi[2,4,:] = dphidy[2]
	divphi[2,5,:] = dphidx[2]

	return divphi


def assemble(x,y,conn,
		psi,dxi,deta,weight,stiffnessTensor
		):
	nint = psi.shape[1]	# number of int points
	nnodes = len(x)			# number of nodes
	nelem = conn.shape[0] 	# number of elements
	stiffnessMatrix = np.zeros([nnodes,nnodes])
	loadVector = np.zeros([nnodes,1])
	for e in xrange(1):
	# for e in xrange(nelem):
		# global coords of the integration points
		# x = sum(x_i*psi(xi))
		xx = np.dot(psi.T,x[conn[e,:]])
		yy = np.dot(psi.T,y[conn[e,:]])
		# Jacobian
		J11 = np.dot(dxi.T,x[conn[e,:]])
		J12 = np.dot(deta.T,x[conn[e,:]])
		J21 = np.dot(dxi.T,y[conn[e,:]])
		J22 = np.dot(deta.T,y[conn[e,:]])
		detJ = J11*J22 - J21*J12
		# Inverse Jacobian
		dxidx = J22/abs(detJ)
		dxidy = - J12/abs(detJ)
		detadx = - J21/abs(detJ)
		detady = J11/abs(detJ)
		# diffirential operator
		divphi = div(dxi,deta,dxidx,dxidy,detadx,detady)
		# print dot(stiffnessTensor,divphi)[:,:,0]/1*(1-0.25)
		EDphi = dot(stiffnessTensor,divphi)
		# EDphiDphi
		integrand = dot(np.transpose(EDphi, (1, 0, 2)),divphi)
		print integrand[:,:,0]/1*(1-0.25)
		integrand = integrand*weight*detJ
		ke = np.sum(integrand,2)
		# print ke



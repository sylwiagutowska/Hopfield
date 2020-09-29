#this function do tetrahedra exactly like in QE (/~/qe-6.0/PW/src/kpoint_grid.f90)
import numpy as np
def tetrahedra(nk1,nk2,nk3,equiv):
	tetra=[[0 for k in range(6*nk1*nk2*nk3)] for i in range(4)]
	for i in range(nk1):
		for j in range(nk2):
			for k in range(nk3):
           #  n1-n8 are the indices of k-point 1-8 forming a cube
				ip1 = (i+1)%nk1 #np.mod(i+1,nk1)
				jp1 = (j+1)%nk2 #np.mod(j+1,nk2)
				kp1 = (k+1)%nk3 #np.mod(k+1,nk3)
				n1 = k   + j  *nk3 + i   *nk2*nk3 
				n2 = k   + j  *nk3 + ip1 *nk2*nk3 
				n3 = k   + jp1*nk3 + i   *nk2*nk3 
				n4 = k   + jp1*nk3 + ip1 *nk2*nk3 
				n5 = kp1 + j  *nk3 + i   *nk2*nk3 
				n6 = kp1 + j  *nk3 + ip1 *nk2*nk3 
				n7 = kp1 + jp1*nk3 + i   *nk2*nk3 
				n8 = kp1 + jp1*nk3 + ip1 *nk2*nk3 
           #  there are 6 tetrahedra per cube (and nk1*nk2*nk3 cubes)
				n  = 6 * ( k + j*nk3 + i*nk3*nk2 )
				tetra [0][n] = equiv[n1]
				tetra [1][n] = equiv[n2]
				tetra [2][n] = equiv[n3]
				tetra [3][n] = equiv[n6]

				tetra [0][n+1] = equiv[n2]
				tetra [1][n+1] = equiv[n3]
				tetra [2][n+1] = equiv[n4]
				tetra [3][n+1] = equiv[n6]

				tetra [0][n+2] = equiv[n1]
				tetra [1][n+2] = equiv[n3]
				tetra [2][n+2] = equiv[n5]
				tetra [3][n+2] = equiv[n6]

				tetra [0][n+3] = equiv[n3]
				tetra [1][n+3] = equiv[n4]
				tetra [2][n+3] = equiv[n6]
				tetra [3][n+3] = equiv[n8]

				tetra [0][n+4] = equiv[n3]
				tetra [1][n+4] = equiv[n6]
				tetra [2][n+4] = equiv[n7]
				tetra [3][n+4] = equiv[n8]

				tetra [0][n+5] = equiv[n3]
				tetra [1][n+5] = equiv[n5]
				tetra [2][n+5] = equiv[n6]
				tetra [3][n+5] = equiv[n7]
	return tetra
############################

#these  2 functions do dos_t exactly like QE (~/qe-6.0/Modules/dost.f90)
def e_tetra(et,nbnd,ntetra,tetra,weights,Noofatoms):
	etetra=[]
	wtetra=[]
	for nt in range(ntetra):
		etetra.append([])
		wtetra.append([])
		for ibnd in range(nbnd):
			m  =np.array([ et[ibnd][tetra[i][nt]]  for i in range(4) ])
			ind=m.argsort()
			m  =[m[i] for i in ind]
			etetra[nt].append(m)

			n  =np.array([ sum([weights[na][ibnd][tetra[i][nt]]  for i in ind])/4. for na in range(Noofatoms) ])
			wtetra[nt].append(n)

	return [np.array(etetra),np.array(wtetra)]

def dos_t(etetra,wtetra,nbnd,ntetra,e,Noofatoms):
	dost=0
	weight=[0 for i in range(Noofatoms)]
	k=0
	for nt in range(ntetra):
		m1    = etetra[nt]
		for ibnd in range(nbnd):
			m2 = wtetra[nt][ibnd]
			[e1,e2,e3,e4]=m1[ibnd] #energie - wierzcholki czworoscianow
			if  (e<e4 and e>=e3):
				k      = 1./ntetra * 3.*(e4-e)**2 /(e4-e1) /(e4-e2) /(e4-e3)
				dost   = dost+k
				weight =[ weight[i]+k*m2[i] for i in range(Noofatoms)]
			elif(e<e3 and e>=e2):	
				k      = 1./ntetra /(e3-e1) /(e4-e1) * ( 3.*(e2-e1) + 6.*(e-e2) - 3.*(e3-e1+e4-e2) /(e3-e2) /(e4-e2) * (e-e2)**2 )
				dost   = dost+k
				weight = [ weight[i]+k*m2[i] for i in range(Noofatoms)]
			elif(e<e2 and e>=e1):
				k      = 1./ntetra * 3.*(e-e1)**2 /(e2-e1) /(e3-e1) /(e4-e1)
				dost   = dost+k
				weight = [ weight[i]+k*m2[i] for i in range(Noofatoms)]

	return [dost,weight] 

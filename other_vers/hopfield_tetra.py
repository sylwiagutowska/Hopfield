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
	print ntetra,nbnd
	etetra=[] #etetra[nt][ibnd][0-3]
	wtetra=[]
	for nt in range(ntetra):
		etetra.append([])
		wtetra.append([])
		for ibnd in range(nbnd):
			m  =np.array([ et[ibnd][tetra[i][nt]]  for i in range(4) ])
			ind=m.argsort()
			m  =[m[i] for i in ind]
			etetra[nt].append(m)

			n  =np.array([[weights[na][ibnd][tetra[i][nt]]  for i in ind] for na in range(Noofatoms) ])
			wtetra[nt].append(n)

	return [np.array(etetra),np.array(wtetra)]


#####this modified as it is in matdyn.f90, function dos_gam (the one, which calculates the a2f=integral of gamma/w)
def dos_t(etetra,wtetra,nbnd,ntetra,e,Noofatoms):
  print ntetra,nbnd
  o13=1./3.
  [P1,P2,P3,P4]=[0.,0.,0.,0.]
  G=0.
  dost=0
  weight=[0 for i in range(Noofatoms)]
  k=0
  for nt in range(ntetra):
   m1    = etetra[nt]
   for ibnd in range(nbnd):
    m2 = wtetra[nt][ibnd]
    [e1,e2,e3,e4]=m1[ibnd] #energie - wierzcholki czworoscianow
 #   if e1==e2 or e1==e3 or e1==e4 or e2==e3 or e2==e4 or e3==e4:
 #    print('Tetrahedron '+str(nt)+'ignored')
 #    continue
 
    if  (e<e4 and e>=e3):
        f14 = (e-e4)/(e1-e4)
        f24 = (e-e4)/(e2-e4)
        f34 = (e-e4)/(e3-e4)

        G  =  3.0 * f14 * f24 * f34 / (e4-e)
        P1 =  f14 * o13
        P2 =  f24 * o13
        P3 =  f34 * o13
        P4 =  (3.0 - f14 - f24 - f34 ) * o13
    elif(e<e3 and e>=e2):	
        f13 = (e-e3)/(e1-e3)
        f31 = 1.0 - f13
        f14 = (e-e4)/(e1-e4)
        f41 = 1.0-f14
        f23 = (e-e3)/(e2-e3)
        f32 = 1.0 - f23
        f24 = (e-e4)/(e2-e4)
        f42 = 1.0 - f24
        G   =  3.0 * (f23*f31 + f32*f24)
        P1  =  f14 * o13 + f13*f31*f23 / G
        P2  =  f23 * o13 + f24*f24*f32 / G
        P3  =  f32 * o13 + f31*f31*f23 / G
        P4  =  f41 * o13 + f42*f24*f32 / G
        G   =  G / (e4-e1)
    elif(e<e2 and e>=e1):
        f12 = (e-e2)/(e1-e2)
        f21 = 1.0 - f12
        f13 = (e-e3)/(e1-e3)
        f31 = 1.0 - f13
        f14 = (e-e4)/(e1-e4)
        f41 = 1.0 - f14
        G  =  3.0 * f21 * f31 * f41 / (e-e1)
        P1 =  o13 * (f12 + f13 + f14)
        P2 =  o13 * f21
        P3 =  o13 * f31
        P4 =  o13 * f41
    else:
        G=0.
    if G>1e5 or G<-1e-5: continue
    dost = dost + G /ntetra
    #m2[9][4]
    for i in range(len(weight)):
     weight[i]=weight[i]+ G*(P1*m2[i][0]+P2*m2[i][1]+P3*m2[i][2]+P4*m2[i][3])/ntetra
  return [dost,weight] 

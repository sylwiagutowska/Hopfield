from inputs import *
from math_stuff import *
from multiprocessing import Process,Pool
#############FOR Nb eta=7.627eV/ang^2 [papaconsta.. 2002] or 4.7 eV/ang^2 (allen-dynes)
##mathematica SphericalHarmonicy(l,m,theta,phi). theta in [0,pi), phi in [0,2pi]
##scipy scipy.special.sph_harm(m, l, theta, phi) , theta in [0,2pi], phi in [0,pi]
##sympy.functions.special.spherical_harmonics.Ynm(l,m,theta,phi . theta in [0,pi), phi in [0,2pi]
CIN=1e-22
class hopfield_spherical(inputs):
 def __init__(self,real_structure,band_structure):
  self.so=inputs.so
  self.RADWF=real_structure.RADWF
  self.RADWFsmall=real_structure.RADWFsmall
  self.Vtot=real_structure.Vtot
  self.V=[vi[0] for vi in self.Vtot] #only spehrical part
  self.DR=real_structure.DR
  self.RMESH=real_structure.RMESH
  self.volume=real_structure.volume
  self.Z_of_atoms=real_structure.Z_of_atoms
  self.n_at=band_structure.n_at
  self.dos=band_structure.dos
  self.ldos=band_structure.LDOS
  self.result=[]
 def calc_hopfield_spherical(self):
  return calc_hopfield_spherical(self.ldos,self.dos,self.V,self.RADWF,self.n_at)
 def print_result(self):
  print_result(self.result,self.dos,self.volume)
  if self.so:   
   print('SOC result: the contribution from small component of radwf" ')
   print_result(self.result_so,self.dos,self.volume)


def calc_hopfield_spherical(LDOS,TOTDOS,V,RADWF,na):
	atomic_eta=[]
	print (LDOS,na)
	for i in range(na):
	 eta_i=0
	 for l in range(len(LDOS[i][2:])-1):
	   eta_l=(2.*l+2.)*LDOS[i][l+2]*LDOS[i][l+3]/(2.*l+1.)/(2.*l+3.)/TOTDOS
	   INTi=[]
	   integral=0 #[0,0]
	   for ri in range(1,len(V[i])):
	#    dr=RMESH[i][ri]-RMESH[i][ri-ri]
	    dVdr=(V[i][ri]-V[i][ri-1]) #/dr
	    '''
	    inti=  [RADWF[i][l][ri  ][0]*dVdr*RADWF[i][l+1][ri  ][0],CIN*(RADWF[i][l][ri  ][1]*dVdr*RADWF[i][l+1][ri  ][1])]
	    inti_1=[RADWF[i][l][ri-1][0]*dVdr*RADWF[i][l+1][ri-1][0],CIN*(RADWF[i][l][ri-1][1]*dVdr*RADWF[i][l+1][ri-1][1])]
	    integral=[integral[m]+(inti[m]+inti_1[m])        *dr    /2. for m in range(2)]
	   eta_l=eta_l*(integral[0]*integral[0]+integral[1]*integral[1])
	   eta_i=eta_i+eta_l
	    '''
	    inti=  RADWF[i][l][0][ri  ]*dVdr*RADWF[i][l+1][0][ri  ]+CIN*(RADWF[i][l][1][ri  ]*dVdr*RADWF[i][l+1][1][ri  ])
	    inti_1=RADWF[i][l][0][ri-1]*dVdr*RADWF[i][l+1][0][ri-1]+CIN*(RADWF[i][l][1][ri-1]*dVdr*RADWF[i][l+1][1][ri-1])
	    integral=integral+(inti+inti_1)/2. #*dr
	   eta_l=eta_l*integral*integral
	   eta_i=eta_i+eta_l
	 atomic_eta.append(eta_i)
	print (atomic_eta)
	return atomic_eta


def parallel_job(no_of_pool,function,arguments):
 print ('No of pools:',no_of_pool)
 with Pool(no_of_pool) as pol:
   result=pol.map(function, arguments)
 return result

def print_result(result,dos,volume):
  h=open('hopfield_results.dat','w')
  tab=['big component','small component']
  print(result)
  '''
  print('Hopfield, Hopfield SOC:')
  for i in range(2): 
   suma=(sum([m[i] for m in result]))
   h.write(tab[i]+str(suma))
   print (suma),
  print('')
  print('the same times dos*rytoev/abtoang^2')
  for i in range(2): 
   suma=(sum([m[i]*dos*Ry_to_eV/ab_to_ang/ab_to_ang for m in result]))
   h.write(suma),
   print(suma)
  '''
  h.write('N(EF)='+str(dos/Ry_to_eV)+'(1/eV)\n')
  h.write('Volume of cell='+str(volume)+'Ang^3\n')
  print('Total <I^2> [default units]'),
  I2=sum(result)
  print (I2)
  h.write('<I^2>='+str(I2)+'\n')
  print('Total hopfield eta= <I^2>*N(EF)  [ev/Ang^2]'),
  eta=I2*dos *Ry_to_eV/ab_to_ang/ab_to_ang
  print(eta)
  h.write('eta='+str(eta))
  print('\n')

def blm_integral(l3,m3,l5,m5,ALMBLM_at, weight,n_k):
 results=np.array([0j,0j,0j,0j,0j,0j,0j,0j,0j,0j,0j,0j,0j,0j,0j,0j])
 for k in range(n_k):
  for iband in  range(len(ALMBLM_at[l3][l3+m3][k])): 
   results+=weight[k][iband]*np.array([\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][0])*ALMBLM_at[l5][l5+m5][k][iband][0],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][0])*ALMBLM_at[l5][l5+m5][k][iband][1],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][0])*ALMBLM_at[l5][l5+m5][k][iband][2],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][0])*ALMBLM_at[l5][l5+m5][k][iband][3],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][1])*ALMBLM_at[l5][l5+m5][k][iband][0],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][1])*ALMBLM_at[l5][l5+m5][k][iband][1],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][1])*ALMBLM_at[l5][l5+m5][k][iband][2],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][1])*ALMBLM_at[l5][l5+m5][k][iband][3],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][2])*ALMBLM_at[l5][l5+m5][k][iband][0],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][2])*ALMBLM_at[l5][l5+m5][k][iband][1],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][2])*ALMBLM_at[l5][l5+m5][k][iband][2],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][2])*ALMBLM_at[l5][l5+m5][k][iband][3],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][3])*ALMBLM_at[l5][l5+m5][k][iband][0],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][3])*ALMBLM_at[l5][l5+m5][k][iband][1],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][3])*ALMBLM_at[l5][l5+m5][k][iband][2],\
np.conjugate(ALMBLM_at[l3][l3+m3][k][iband][3])*ALMBLM_at[l5][l5+m5][k][iband][3]])
 return results

def sum_weights(weight,n_k_total):
 suma=sum([ sum(i) for i in weight]) 
 print (suma/n_k_total, 'should equal to DOS in 1/Ry')
 return suma #suma devided by n_k_total should equal to N(EF)

def blm_all_integrals(ALMBLM_at, weight,n_l,n_k,n_k_total,so):
 if so: step=1
 else: step=2
 suma=sum_weights(weight,n_k_total)
 blm_integrals=[[[[ [] for m5 in range(-l5,l5+1) ] for l5 in range(n_l) ] for m3 in range(-l3,l3+1) ] for l3 in range(n_l) ]
 for l3 in range(n_l):
  for m3 in range(-l3,l3+1):
   for l5 in range(n_l):
    for m5 in range(-l5,l5+1):
     blm_integrals[l3][l3+m3][l5][l5+m5]=(blm_integral(l3,m3,l5,m5,ALMBLM_at, weight, n_k))\
                                          /suma/step #I am not sure, but in my experience it is needed to get agreement soc and no soc values. Probably because Almblm are doubled in no/soc case to get DOS per both spins
 return blm_integrals

def r_integral(RADWF_at_l3,RADWF_at_l4,Vtot_at,DVDR1,RMESH_at,Z_of_atom):
 #     norm=np.trapz(np.multiply(RADWF_at_l3,RADWF_at_l4)*(RMESH_at*RMESH_at),x=RMESH_at)
 #     print(norm)
 #     if norm==0: return [0,0]
 #     else:
#       print(len(DVDR1),len(RADWF_at_l3),len(RADWF_at_l4),len(RMESH_at))
       norma=np.trapz(np.multiply(RADWF_at_l3,RADWF_at_l4),x=RMESH_at)
       Ba_and_C_r=np.trapz(np.multiply(np.multiply(Vtot_at,RADWF_at_l3),(RADWF_at_l4)/RMESH_at),x=RMESH_at)
       Aa_r=np.trapz(np.multiply(np.multiply(DVDR1,RADWF_at_l3),RADWF_at_l4),x=RMESH_at)
       '''
      Ba_and_C_r,Aa_r=np.array([0.j,0.j]),np.array([0.j,0.j])
      for nr1,r1 in enumerate(RMESH_at): 
       if nr1==0: continue  
       #if Vr^2 devided by r (so we have V*r)
       dv1dr=(Vtot_at[nr1]-Vtot_at[nr1-1]-Vtot_at[nr1]*dr1)*r1
       # every *r1 -> /r1 and *r1^2 -> 1 because v is multiplied by r^2
       Ba_and_C_r+=Vtot_at[nr1]* (RADWF_at_l3[nr1])*RADWF_at_l4[nr1]*dr1   #*r
       Aa_r+= dv1dr * (RADWF_at_l3[nr1])*RADWF_at_l4[nr1] #*r*
       '''
#       if norma==0: return [0,0]
#       return [Ba_and_C_r/norma,Aa_r/norma]
       return [Ba_and_C_r,Aa_r]

def r_all_integrals(RADWF_at,Vtot_at,DVDR1,RMESH_at,n_l, Z_of_atom):
 r_integrals=[[[ []  for l5 in range(n_l) ] for l3 in range(n_l) ] for l1 in range(len(DVDR1))]
 for l1 in range(len(DVDR1)):
  for l3 in range(n_l):
   for l5 in range(n_l):
      r_integrals[l1][l3][l5]=[\
r_integral(np.conjugate(RADWF_at[l3][0]),RADWF_at[l5][0],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][0]),RADWF_at[l5][1],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][0]),RADWF_at[l5][2],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][0]),RADWF_at[l5][3],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][1]),RADWF_at[l5][0],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][1]),RADWF_at[l5][1],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][1]),RADWF_at[l5][2],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][1]),RADWF_at[l5][3],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][2]),RADWF_at[l5][0],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][2]),RADWF_at[l5][1],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][2]),RADWF_at[l5][2],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][2]),RADWF_at[l5][3],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][3]),RADWF_at[l5][0],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][3]),RADWF_at[l5][1],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][3]),RADWF_at[l5][2],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom),\
r_integral(np.conjugate(RADWF_at[l3][3]),RADWF_at[l5][3],Vtot_at[l1],DVDR1[l1],RMESH_at,Z_of_atom)]
 return r_integrals

def chosen_r_integrals(r_int34):
 B_and_C_and_A_r= np.array([\
r_int34[0],r_int34[0],r_int34[0],\
r_int34[0],r_int34[0],r_int34[0],\
r_int34[0],r_int34[0],r_int34[0],\
r_int34[0],r_int34[0],r_int34[0],\
r_int34[0],r_int34[0],r_int34[0],\
r_int34[0],r_int34[1],r_int34[1],\
r_int34[1],r_int34[1],r_int34[1],\
r_int34[1],r_int34[1],r_int34[1],\
r_int34[1],r_int34[1],r_int34[1],\
r_int34[1],r_int34[1],r_int34[1],\
r_int34[1],r_int34[1],r_int34[2],\
r_int34[2],r_int34[2],r_int34[2],\
r_int34[2],r_int34[2],r_int34[2],\
r_int34[2],r_int34[2],r_int34[2],\
r_int34[2],r_int34[2],r_int34[2],\
r_int34[2],r_int34[2],r_int34[2],\
r_int34[3],r_int34[3],r_int34[3],\
r_int34[3],r_int34[3],r_int34[3],\
r_int34[3],r_int34[3],r_int34[3],\
r_int34[3],r_int34[3],r_int34[3],\
r_int34[3],r_int34[3],r_int34[3],\
r_int34[3],r_int34[4],r_int34[4],\
r_int34[4],r_int34[4],r_int34[4],\
r_int34[4],r_int34[4],r_int34[4],\
r_int34[4],r_int34[4],r_int34[4],\
r_int34[4],r_int34[4],r_int34[4],\
r_int34[4],r_int34[4],r_int34[5],\
r_int34[5],r_int34[5],r_int34[5],\
r_int34[5],r_int34[5],r_int34[5],\
r_int34[5],r_int34[5],r_int34[5],\
r_int34[5],r_int34[5],r_int34[5],\
r_int34[5],r_int34[5],r_int34[5],\
r_int34[6],r_int34[6],r_int34[6],\
r_int34[6],r_int34[6],r_int34[6],\
r_int34[6],r_int34[6],r_int34[6],\
r_int34[6],r_int34[6],r_int34[6],\
r_int34[6],r_int34[6],r_int34[6],\
r_int34[6],r_int34[7],r_int34[7],\
r_int34[7],r_int34[7],r_int34[7],\
r_int34[7],r_int34[7],r_int34[7],\
r_int34[7],r_int34[7],r_int34[7],\
r_int34[7],r_int34[7],r_int34[7],\
r_int34[7],r_int34[7],r_int34[8],\
r_int34[8],r_int34[8],r_int34[8],\
r_int34[8],r_int34[8],r_int34[8],\
r_int34[8],r_int34[8],r_int34[8],\
r_int34[8],r_int34[8],r_int34[8],\
r_int34[8],r_int34[8],r_int34[8],\
r_int34[9],r_int34[9],r_int34[9],\
r_int34[9],r_int34[9],r_int34[9],\
r_int34[9],r_int34[9],r_int34[9],\
r_int34[9],r_int34[9],r_int34[9],\
r_int34[9],r_int34[9],r_int34[9],\
r_int34[9],r_int34[10],r_int34[10],\
r_int34[10],r_int34[10],r_int34[10],\
r_int34[10],r_int34[10],r_int34[10],\
r_int34[10],r_int34[10],r_int34[10],\
r_int34[10],r_int34[10],r_int34[10],\
r_int34[10],r_int34[10],r_int34[11],\
r_int34[11],r_int34[11],r_int34[11],\
r_int34[11],r_int34[11],r_int34[11],\
r_int34[11],r_int34[11],r_int34[11],\
r_int34[11],r_int34[11],r_int34[11],\
r_int34[11],r_int34[11],r_int34[11],\
r_int34[12],r_int34[12],r_int34[12],\
r_int34[12],r_int34[12],r_int34[12],\
r_int34[12],r_int34[12],r_int34[12],\
r_int34[12],r_int34[12],r_int34[12],\
r_int34[12],r_int34[12],r_int34[12],\
r_int34[12],r_int34[13],r_int34[13],\
r_int34[13],r_int34[13],r_int34[13],\
r_int34[13],r_int34[13],r_int34[13],\
r_int34[13],r_int34[13],r_int34[13],\
r_int34[13],r_int34[13],r_int34[13],\
r_int34[13],r_int34[13],r_int34[14],\
r_int34[14],r_int34[14],r_int34[14],\
r_int34[14],r_int34[14],r_int34[14],\
r_int34[14],r_int34[14],r_int34[14],\
r_int34[14],r_int34[14],r_int34[14],\
r_int34[14],r_int34[14],r_int34[14],\
r_int34[15],r_int34[15],r_int34[15],\
r_int34[15],r_int34[15],r_int34[15],\
r_int34[15],r_int34[15],r_int34[15],\
r_int34[15],r_int34[15],r_int34[15],\
r_int34[15],r_int34[15],r_int34[15],\
r_int34[15]  ])
 return np.array([i[0] for i in B_and_C_and_A_r]),np.array([i[1] for i in B_and_C_and_A_r])

def chosen_rp_integrals(r_int56):
 B_and_C_and_A_rp= np.array([\
r_int56[0],r_int56[1],\
r_int56[2],r_int56[3],r_int56[4],\
r_int56[5],r_int56[6],r_int56[7],\
r_int56[8],r_int56[9],r_int56[10],\
r_int56[11],r_int56[12],r_int56[13],\
r_int56[14],r_int56[15],r_int56[0],\
r_int56[1],r_int56[2],r_int56[3],\
r_int56[4],r_int56[5],r_int56[6],\
r_int56[7],r_int56[8],r_int56[9],\
r_int56[10],r_int56[11],r_int56[12],\
r_int56[13],r_int56[14],r_int56[15],\
r_int56[0],r_int56[1],r_int56[2],\
r_int56[3],r_int56[4],r_int56[5],\
r_int56[6],r_int56[7],r_int56[8],\
r_int56[9],r_int56[10],r_int56[11],\
r_int56[12],r_int56[13],r_int56[14],\
r_int56[15],r_int56[0],r_int56[1],\
r_int56[2],r_int56[3],r_int56[4],\
r_int56[5],r_int56[6],r_int56[7],\
r_int56[8],r_int56[9],r_int56[10],\
r_int56[11],r_int56[12],r_int56[13],\
r_int56[14],r_int56[15],r_int56[0],\
r_int56[1],r_int56[2],r_int56[3],\
r_int56[4],r_int56[5],r_int56[6],\
r_int56[7],r_int56[8],r_int56[9],\
r_int56[10],r_int56[11],r_int56[12],\
r_int56[13],r_int56[14],r_int56[15],\
r_int56[0],r_int56[1],r_int56[2],\
r_int56[3],r_int56[4],r_int56[5],\
r_int56[6],r_int56[7],r_int56[8],\
r_int56[9],r_int56[10],r_int56[11],\
r_int56[12],r_int56[13],r_int56[14],\
r_int56[15],r_int56[0],r_int56[1],\
r_int56[2],r_int56[3],r_int56[4],\
r_int56[5],r_int56[6],r_int56[7],\
r_int56[8],r_int56[9],r_int56[10],\
r_int56[11],r_int56[12],r_int56[13],\
r_int56[14],r_int56[15],r_int56[0],\
r_int56[1],r_int56[2],r_int56[3],\
r_int56[4],r_int56[5],r_int56[6],\
r_int56[7],r_int56[8],r_int56[9],\
r_int56[10],r_int56[11],r_int56[12],\
r_int56[13],r_int56[14],r_int56[15],\
r_int56[0],r_int56[1],r_int56[2],\
r_int56[3],r_int56[4],r_int56[5],\
r_int56[6],r_int56[7],r_int56[8],\
r_int56[9],r_int56[10],r_int56[11],\
r_int56[12],r_int56[13],r_int56[14],\
r_int56[15],r_int56[0],r_int56[1],\
r_int56[2],r_int56[3],r_int56[4],\
r_int56[5],r_int56[6],r_int56[7],\
r_int56[8],r_int56[9],r_int56[10],\
r_int56[11],r_int56[12],r_int56[13],\
r_int56[14],r_int56[15],r_int56[0],\
r_int56[1],r_int56[2],r_int56[3],\
r_int56[4],r_int56[5],r_int56[6],\
r_int56[7],r_int56[8],r_int56[9],\
r_int56[10],r_int56[11],r_int56[12],\
r_int56[13],r_int56[14],r_int56[15],\
r_int56[0],r_int56[1],r_int56[2],\
r_int56[3],r_int56[4],r_int56[5],\
r_int56[6],r_int56[7],r_int56[8],\
r_int56[9],r_int56[10],r_int56[11],\
r_int56[12],r_int56[13],r_int56[14],\
r_int56[15],r_int56[0],r_int56[1],\
r_int56[2],r_int56[3],r_int56[4],\
r_int56[5],r_int56[6],r_int56[7],\
r_int56[8],r_int56[9],r_int56[10],\
r_int56[11],r_int56[12],r_int56[13],\
r_int56[14],r_int56[15],r_int56[0],\
r_int56[1],r_int56[2],r_int56[3],\
r_int56[4],r_int56[5],r_int56[6],\
r_int56[7],r_int56[8],r_int56[9],\
r_int56[10],r_int56[11],r_int56[12],\
r_int56[13],r_int56[14],r_int56[15],\
r_int56[0],r_int56[1],r_int56[2],\
r_int56[3],r_int56[4],r_int56[5],\
r_int56[6],r_int56[7],r_int56[8],\
r_int56[9],r_int56[10],r_int56[11],\
r_int56[12],r_int56[13],r_int56[14],\
r_int56[15],r_int56[0],r_int56[1],\
r_int56[2],r_int56[3],r_int56[4],\
r_int56[5],r_int56[6],r_int56[7],\
r_int56[8],r_int56[9],r_int56[10],\
r_int56[11],r_int56[12],r_int56[13],\
r_int56[14],r_int56[15]  ])
 return np.array([i[0] for i in B_and_C_and_A_rp]),np.array([i[1] for i in B_and_C_and_A_rp])

def chosen_k_integrals(k_int36):
 k_integrals= np.array([\
k_int36[0],\
k_int36[1],k_int36[2],k_int36[3],\
k_int36[0],k_int36[1],k_int36[2],\
k_int36[3],k_int36[0],k_int36[1],\
k_int36[2],k_int36[3],k_int36[0],\
k_int36[1],k_int36[2],k_int36[3],\
k_int36[0],k_int36[1],k_int36[2],\
k_int36[3],k_int36[0],k_int36[1],\
k_int36[2],k_int36[3],k_int36[0],\
k_int36[1],k_int36[2],k_int36[3],\
k_int36[0],k_int36[1],k_int36[2],\
k_int36[3],k_int36[0],k_int36[1],\
k_int36[2],k_int36[3],k_int36[0],\
k_int36[1],k_int36[2],k_int36[3],\
k_int36[0],k_int36[1],k_int36[2],\
k_int36[3],k_int36[0],k_int36[1],\
k_int36[2],k_int36[3],k_int36[0],\
k_int36[1],k_int36[2],k_int36[3],\
k_int36[0],k_int36[1],k_int36[2],\
k_int36[3],k_int36[0],k_int36[1],\
k_int36[2],k_int36[3],k_int36[0],\
k_int36[1],k_int36[2],k_int36[3],\
k_int36[4],k_int36[5],k_int36[6],\
k_int36[7],k_int36[4],k_int36[5],\
k_int36[6],k_int36[7],k_int36[4],\
k_int36[5],k_int36[6],k_int36[7],\
k_int36[4],k_int36[5],k_int36[6],\
k_int36[7],k_int36[4],k_int36[5],\
k_int36[6],k_int36[7],k_int36[4],\
k_int36[5],k_int36[6],k_int36[7],\
k_int36[4],k_int36[5],k_int36[6],\
k_int36[7],k_int36[4],k_int36[5],\
k_int36[6],k_int36[7],k_int36[4],\
k_int36[5],k_int36[6],k_int36[7],\
k_int36[4],k_int36[5],k_int36[6],\
k_int36[7],k_int36[4],k_int36[5],\
k_int36[6],k_int36[7],k_int36[4],\
k_int36[5],k_int36[6],k_int36[7],\
k_int36[4],k_int36[5],k_int36[6],\
k_int36[7],k_int36[4],k_int36[5],\
k_int36[6],k_int36[7],k_int36[4],\
k_int36[5],k_int36[6],k_int36[7],\
k_int36[4],k_int36[5],k_int36[6],\
k_int36[7],k_int36[8],k_int36[9],\
k_int36[10],k_int36[11],k_int36[8],\
k_int36[9],k_int36[10],k_int36[11],\
k_int36[8],k_int36[9],k_int36[10],\
k_int36[11],k_int36[8],k_int36[9],\
k_int36[10],k_int36[11],k_int36[8],\
k_int36[9],k_int36[10],k_int36[11],\
k_int36[8],k_int36[9],k_int36[10],\
k_int36[11],k_int36[8],k_int36[9],\
k_int36[10],k_int36[11],k_int36[8],\
k_int36[9],k_int36[10],k_int36[11],\
k_int36[8],k_int36[9],k_int36[10],\
k_int36[11],k_int36[8],k_int36[9],\
k_int36[10],k_int36[11],k_int36[8],\
k_int36[9],k_int36[10],k_int36[11],\
k_int36[8],k_int36[9],k_int36[10],\
k_int36[11],k_int36[8],k_int36[9],\
k_int36[10],k_int36[11],k_int36[8],\
k_int36[9],k_int36[10],k_int36[11],\
k_int36[8],k_int36[9],k_int36[10],\
k_int36[11],k_int36[8],k_int36[9],\
k_int36[10],k_int36[11],k_int36[12],\
k_int36[13],k_int36[14],k_int36[15],\
k_int36[12],k_int36[13],k_int36[14],\
k_int36[15],k_int36[12],k_int36[13],\
k_int36[14],k_int36[15],k_int36[12],\
k_int36[13],k_int36[14],k_int36[15],\
k_int36[12],k_int36[13],k_int36[14],\
k_int36[15],k_int36[12],k_int36[13],\
k_int36[14],k_int36[15],k_int36[12],\
k_int36[13],k_int36[14],k_int36[15],\
k_int36[12],k_int36[13],k_int36[14],\
k_int36[15],k_int36[12],k_int36[13],\
k_int36[14],k_int36[15],k_int36[12],\
k_int36[13],k_int36[14],k_int36[15],\
k_int36[12],k_int36[13],k_int36[14],\
k_int36[15],k_int36[12],k_int36[13],\
k_int36[14],k_int36[15],k_int36[12],\
k_int36[13],k_int36[14],k_int36[15],\
k_int36[12],k_int36[13],k_int36[14],\
k_int36[15],k_int36[12],k_int36[13],\
k_int36[14],k_int36[15],k_int36[12],\
k_int36[13],k_int36[14],k_int36[15]])
 return k_integrals

def chosen_kp_integrals(k_int45):
 kp_integrals=np.array([\
k_int45[0],k_int45[0],k_int45[0],\
k_int45[0],k_int45[1],k_int45[1],\
k_int45[1],k_int45[1],k_int45[2],\
k_int45[2],k_int45[2],k_int45[2],\
k_int45[3],k_int45[3],k_int45[3],\
k_int45[3],k_int45[4],k_int45[4],\
k_int45[4],k_int45[4],k_int45[5],\
k_int45[5],k_int45[5],k_int45[5],\
k_int45[6],k_int45[6],k_int45[6],\
k_int45[6],k_int45[7],k_int45[7],\
k_int45[7],k_int45[7],k_int45[8],\
k_int45[8],k_int45[8],k_int45[8],\
k_int45[9],k_int45[9],k_int45[9],\
k_int45[9],k_int45[10],k_int45[10],\
k_int45[10],k_int45[10],k_int45[11],\
k_int45[11],k_int45[11],k_int45[11],\
k_int45[12],k_int45[12],k_int45[12],\
k_int45[12],k_int45[13],k_int45[13],\
k_int45[13],k_int45[13],k_int45[14],\
k_int45[14],k_int45[14],k_int45[14],\
k_int45[15],k_int45[15],k_int45[15],\
k_int45[15],k_int45[0],k_int45[0],\
k_int45[0],k_int45[0],k_int45[1],\
k_int45[1],k_int45[1],k_int45[1],\
k_int45[2],k_int45[2],k_int45[2],\
k_int45[2],k_int45[3],k_int45[3],\
k_int45[3],k_int45[3],k_int45[4],\
k_int45[4],k_int45[4],k_int45[4],\
k_int45[5],k_int45[5],k_int45[5],\
k_int45[5],k_int45[6],k_int45[6],\
k_int45[6],k_int45[6],k_int45[7],\
k_int45[7],k_int45[7],k_int45[7],\
k_int45[8],k_int45[8],k_int45[8],\
k_int45[8],k_int45[9],k_int45[9],\
k_int45[9],k_int45[9],k_int45[10],\
k_int45[10],k_int45[10],k_int45[10],\
k_int45[11],k_int45[11],k_int45[11],\
k_int45[11],k_int45[12],k_int45[12],\
k_int45[12],k_int45[12],k_int45[13],\
k_int45[13],k_int45[13],k_int45[13],\
k_int45[14],k_int45[14],k_int45[14],\
k_int45[14],k_int45[15],k_int45[15],\
k_int45[15],k_int45[15],k_int45[0],\
k_int45[0],k_int45[0],k_int45[0],\
k_int45[1],k_int45[1],k_int45[1],\
k_int45[1],k_int45[2],k_int45[2],\
k_int45[2],k_int45[2],k_int45[3],\
k_int45[3],k_int45[3],k_int45[3],\
k_int45[4],k_int45[4],k_int45[4],\
k_int45[4],k_int45[5],k_int45[5],\
k_int45[5],k_int45[5],k_int45[6],\
k_int45[6],k_int45[6],k_int45[6],\
k_int45[7],k_int45[7],k_int45[7],\
k_int45[7],k_int45[8],k_int45[8],\
k_int45[8],k_int45[8],k_int45[9],\
k_int45[9],k_int45[9],k_int45[9],\
k_int45[10],k_int45[10],k_int45[10],\
k_int45[10],k_int45[11],k_int45[11],\
k_int45[11],k_int45[11],k_int45[12],\
k_int45[12],k_int45[12],k_int45[12],\
k_int45[13],k_int45[13],k_int45[13],\
k_int45[13],k_int45[14],k_int45[14],\
k_int45[14],k_int45[14],k_int45[15],\
k_int45[15],k_int45[15],k_int45[15],\
k_int45[0],k_int45[0],k_int45[0],\
k_int45[0],k_int45[1],k_int45[1],\
k_int45[1],k_int45[1],k_int45[2],\
k_int45[2],k_int45[2],k_int45[2],\
k_int45[3],k_int45[3],k_int45[3],\
k_int45[3],k_int45[4],k_int45[4],\
k_int45[4],k_int45[4],k_int45[5],\
k_int45[5],k_int45[5],k_int45[5],\
k_int45[6],k_int45[6],k_int45[6],\
k_int45[6],k_int45[7],k_int45[7],\
k_int45[7],k_int45[7],k_int45[8],\
k_int45[8],k_int45[8],k_int45[8],\
k_int45[9],k_int45[9],k_int45[9],\
k_int45[9],k_int45[10],k_int45[10],\
k_int45[10],k_int45[10],k_int45[11],\
k_int45[11],k_int45[11],k_int45[11],\
k_int45[12],k_int45[12],k_int45[12],\
k_int45[12],k_int45[13],k_int45[13],\
k_int45[13],k_int45[13],k_int45[14],\
k_int45[14],k_int45[14],k_int45[14],\
k_int45[15],k_int45[15],k_int45[15],\
k_int45[15] ])
 return kp_integrals








'''
Ciekawostka: w wien2k (lapw5) harmoniki sferyczne sa produkowane iteracyjnie
        Start:
!                        +------+
!                        |   1
!           Y(0,0) =  -+ | -----
!                       \| 4(Pi)
!
!                                   +------+
!                                   |   3
!           Y(1,0) =  cos(Theta) -+ | -----
!                                  \| 4(Pi)
!
!                                     +------+
!                                     |   3    i(Phi)
!           Y(1,1) =  - sin(Theta) -+ | ----- e
!                                    \| 8(Pi)
!
!        Formula 1:
!
!           Y(l,l) =
!                           +--------+
!                           | (2l+1)   i(Phi)
!            -sin(Theta) -+ | ------  e       Y(l-1,l-1)
!                          \|   2l
!
!        Formula 2:
!                                  +---------------+
!                                  |  (2l-1)(2l+1)
!           Y(l,m) = cos(Theta) -+ | -------------- Y(l-1,m)  -
!                                 \|   (l-m)(l+m)
!
!                                    +--------------------+
!                                    |(l-1+m)(l-1-m)(2l+1)
!                              -  -+ |-------------------- Y(l-2,m)
!                                   \|  (2l-3)(l-m)(l+m)
!

oraz z wykorzystaiem wlasnosci
!                                          m	  *
!                            Y(l,-m) = (-1) Y(l,m)
'''

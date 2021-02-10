import numpy as np
import os
import glob
#import hopfield_tetra
#from math import factorial as sil
from scipy.special import sph_harm
from scipy.integrate import quad
from scipy.misc import derivative
from multiprocessing import Process,Pool
sqrtpm1=(np.pi)**(-0.5)
degauss=0.05
prefix=os.getcwd().split('/')[-1]
print(prefix)

SOC=0
CIN=4e-22
if os.stat("DOS/DOS.inso").st_size != 0:  
 CIN=1./(137.0359895*2) #for treating big and small component of wave function
 SOC=1
#############FOR Nb eta=7.627eV/ang^2 [papaconsta.. 2002] or 4.7 eV/ang^2 (allen-dynes)
##mathematica SphericalHarmonicy(l,m,theta,phi). theta in [0,pi), phi in [0,2pi]
##scipy scipy.special.sph_harm(m, l, theta, phi) , theta in [0,2pi], phi in [0,pi]
##sympy.functions.special.spherical_harmonics.Ynm(l,m,theta,phi . theta in [0,pi), phi in [0,2pi]
def round_complex(z,precis):
 return round(z.real,precis)+round(z.imag,precis)*1j

def sil(n):
 if n<0: raise ValueError('w silni n='+str(n)+'<0')
 elif n==0: return 1.
 else:
  a=1
  for i in range(1,int(n)+1):
   a=a*i
  return float(a)

def clebsch_gordon(l1,m1,l2,m2,l3,m3):
# if m3!=(m2+m1): return 0
 if (l3+l1-l2)<0 or (l3-l1+l2)<0 or (l1+l2-l3)<0: return 0 
 fac1=( ( (2*l3+1.)*sil(l3+l1-l2)*sil(l3-l1+l2)*sil(l1+l2-l3) )/ sil(l1+l2+l3+1) )**0.5
 fac2=( sil(l3+m3)*sil(l3-m3)*sil(l1-m1)*sil(l1+m1)*sil(l2-m2)*sil(l2+m2) )**0.5
 fac3=0
 for k in range(0,1000):
  if ((l1+l2-l3-k)<0) or ((l1-m1-k)<0) or ((l2+m2-k)<0) or ((l3-l2+m1+k)<0) or ((l3-l1-m2+k)<0): continue
  fac3+= ((-1)**k)/(sil(k)*sil(l1+l2-l3-k)*sil(l1-m1-k)*sil(l2+m2-k)*sil(l3-l2+m1+k)*sil(l3-l1-m2+k))
 return fac1*fac2*fac3

def three_y(l1,m1,l2,m2,l3,m3):
 return round(( (2*l1+1.)*(2*l2+1.)/(4.*np.pi*(2*l3+1.)) )**0.5\
        *clebsch_gordon(l1,0,l2,0,l3,0)*clebsch_gordon(l1,m1,l2,m2,l3,m3),9)

def ynre(m,l,tet,fi):
 return np.real(sph_harm(m1,l1,fi,tet))
def ynim(m,l,tet,fi):
 return np.imag(sph_harm(m1,l1,fi,tet))

def spec_int(l1,m1,l2,m2,l3,m3):
# if m1==0: return 0.
 out1=quad(lambda tet: np.real(sph_harm(m1,l1,0,tet)*np.conj(sph_harm(m2,l2,0,tet))*sph_harm(m3,l3,0,tet)), 0.0, np.pi)
 out2=quad(lambda tet: np.imag(sph_harm(m1,l1,0,tet)*np.conj(sph_harm(m2,l2,0,tet))*sph_harm(m3,l3,0,tet)), 0.0, np.pi)
 return 2*np.pi*m1*complex(-round(out2[0],6),round(out1[0],6)) #calka*2*pi*m*i


################RADWF is already multiplied by r and Vtotal - by r^2
'''
for l1 in range(2):
 for m1 in range(-l1,l1+1):
  for l2 in range(2):
   for m2 in range(-l2,l2+1):
    for l3 in range(3):
     for m3 in range(-l3,l3+1):
      print l1,m1,l2,m2,l3,m3,three_y(l1,m1,l3,m3,l2,m2),three_y(l1,m1,l2,m2,l3,m3)
exit()
'''


#prefix='Nb' #raw_input('Prefix: ')

Ry_to_eV=13.605693122994
ab_to_ang=0.529177210903

#j_to_ev=1/(1.602176634)*1e19
#hbar=6.582119569e-16 #ev/s
#hbar=hbar/13.6056980659 #ry/s
pi=3.141592653589793238462643
#Na=6.02214076e23



'''
print('Reading EF and Volume...'),
os.system('grep ":FER " DOS/DOS.scf |head -n 1 >>hopfield_tmp')
os.system('grep ":VOL " '+prefix+'.scf0 |head -n 1 >>hopfield_tmp')
os.system('grep ":VOL " '+prefix+'.scf |head -n 1 >>hopfield_tmp')
h=open('hopfield_tmp')
tmp=h.readlines()
h.close()
os.system('rm hopfield_tmp')
EF=float(tmp[0].split()[-1])
volume=float(tmp[-1].split()[-1])
print(volume,EF)

print('Reading Z of atoms...'),
os.system('grep "Z: " '+prefix+'.struct  >hopfield_tmp')
h=open('hopfield_tmp')
tmp=h.readlines()
h.close()
os.system('rm hopfield_tmp')
Z_of_atoms=[float(i.split()[-1]) for i in tmp]
print(Z_of_atoms)
'''

print('Read intrasitial total potential...')
h=open('vtotal/vtotal.vtotal')
tmp=h.readlines()
h.close()
Vtot=[] #[atom][lm][r]
LM=[] #[atom][list of lm numbers]
for i in tmp:
 if 'IN INTERSTITIAL' in i: break
 if 'ATOMNUMBER' in i: 
  Vtot.append([])
  LM.append([])
 if 'VLM(R)' in i:
  Vtot[-1].append([])
  LM[-1].append([int(i.split()[3]),int(i.split()[5])])
 elif len(i.split())==0: continue
 elif len(Vtot)!=0:
   m=i[3:-1]
   for j in range(int(len(m)/19)):
    Vtot[-1][-1].append(float(m[j*19:(j+1)*19]))

n_at=len(Vtot) #number of atoms
nr=len(Vtot[0][0]) #number of r-points
print (n_at,nr, LM)
'''
V_inter=[]
K_inter=[]
print('Read intersitial total potential...')
for ni,i in enumerate(tmp):
 if 'TOTAL POTENTIAL IN INTERSITIAL' in i:
  n_K_intersitial=tmp[ni+2].split()[0]
  tmp=tmp[ni+3:n_K_intersitial+3]
  break
for i in tmp:
 ii=i.replace('-',' -').split()
 K_inter.append([int(m) for m in ii[:3]])
 V_inter.append(complex(float(ii[3]),float(ii[4])))

'''


print('Reading radial wave funtions...')
h=open('DOS/DOS.radwf','r')
tmp=h.readlines()
h.close()

RADWF=[]   #RADWF[i][j][k][0-1] i-atoms, j- l (orbital No), k - r-mesh, [0-1]-large  component of radwf and udot
RADWFsmall=[] #RADWFsmall[i][j][k][0-1] i-atoms, j- l (orbital No), k - r-mesh, [0-1]-  small component of radwf and udot
mesh_info=[]
for i in tmp:
 if str(nr) in i.split(): 
  RADWF.append([])
  RADWFsmall.append([])
  mesh_info.append([float(m) for m in i.split()])
 elif len(i.split())==1: 
  RADWF[-1].append([])
  RADWFsmall[-1].append([])
 else:
  try: 
   RADWF[-1][-1].append(np.array([float(i.split()[0]),float(i.split()[2])]))
   RADWFsmall[-1][-1].append(np.array([CIN*float(i.split()[1]),CIN*float(i.split()[3])]))
  except:continue
n_l=min([ len(i) for i in RADWF])


for i in range(len(RADWF)):
 print('For atom no. '+str(i)+' I found '+str(len(RADWF[i]))+' wave functions u_l(r) at the r-meshes as follows:')
 for j in range(len(RADWF[i])):
  print(len(RADWF[i][j])),
 print(' ')

print('Making r-mesh...')
RMESH=[] #RMESH[i][j] i -atom, j -r-point
#(r_i=r0*exp((i-1)*dx)
#   1 781   0.0000100000   0.0159348926   2.5000000000
# write(23,4645)jatom,jri(jatom),r0(jatom),dx(jatom),rmt(jatom)
for i in mesh_info:
 [r0,dx,rmt]=i[2:]
 RMESH.append([])
 for j in range(nr):
  RMESH[-1].append(r0*np.exp(j*dx))
 print('RMT='+str(rmt)+' should ='+str(RMESH[-1][-1]))

###devide V by Q
###and multiply by 2 to obtain Hartree -> Ry
#oryginalnie w pliku *vtotal jest zapisany V*R^2, wiec dzielimy przez R^2 by dostac czysty V
#Vtot=[ [[Vtot[i][j][k]/(RMESH[i][k]**2)/Z_of_atoms[i]*2 for k in range(len(Vtot[i][j]))] for j in range(len(Vtot[i]))] for i in range(len(Vtot))]   
#Vtot=[ [[Vtot[i][j][k]/(RMESH[i][k]**2) for k in range(len(Vtot[i][j]))] for j in range(len(Vtot[i]))] for i in range(len(Vtot))]  
Vtot=[ [[Vtot[i][j][k]/(RMESH[i][k]) for k in range(len(Vtot[i][j]))] for j in range(len(Vtot[i]))] for i in range(len(Vtot))]  

#RADWF[i][j][k][0-1] i-atoms, j- l (orbital No), k - r-mesh, [0-1]-large  component of radwf and udot
#oryginalnie w pliku *radwf jest zapisany radwf*R, wiec dzielimy przez R by dostac czysty radwf
RADWF=[ [[[RADWF[i][j][k][0]/(RMESH[i][k]),RADWF[i][j][k][1] ] for k in range(len(RADWF[i][j]))] for j in range(len(RADWF[i]))] for i in range(len(RADWF))]   
RADWFsmall=[ [[[RADWFsmall[i][j][k][0]/(RMESH[i][k]),RADWFsmall[i][j][k][1] ] for k in range(len(RADWFsmall[i][j]))] for j in range(len(RADWFsmall[i]))] for i in range(len(RADWFsmall))]   

print('Write spherical V and radwf to files V_i.dat and RADWF_i.dat...')
for i in range(len(Vtot)):
 h=open('V_'+str(i)+'.dat','w')
 for j in range(len(Vtot[i][0])):
  h.write(str(RMESH[i][j])+' ')
  for k in range(len(Vtot[i])):
   h.write(str(Vtot[i][k][j])+' ')
  h.write('\n')
 h.close()


for i in range(len(RADWF)):
 print('aa')
 h=open('RADWF_'+str(i)+'.dat','w')
 h.write('# r , l=0- big component,small component, l=1- big component,small component....\n')
 for k in range(len(RADWF[i][1])):
  h.write(str(RMESH[i][k])+' ')
  for j in range(len(RADWF[i])):
   h.write(str(RADWF[i][j][k][0])+' '+str(RADWF[i][j][k][1])+' ')
  h.write('\n')
 h.close() 

for i in range(len(RADWFsmall)):
 print('aa')
 h=open('RADWFsmall_'+str(i)+'.dat','w')
 h.write('# r , l=0- big component,small component, l=1- big component,small component....\n')
 for k in range(len(RADWFsmall[i][1])):
  h.write(str(RMESH[i][k])+' ')
  for j in range(len(RADWFsmall[i])):
   h.write(str(RADWFsmall[i][j][k][0])+' '+str(RADWFsmall[i][j][k][1])+' ')
  h.write('\n')
 h.close() 



print('Read EF and total DOS...'),

h=open('DOS/tot/DOS.outputt','r')
tmp=h.readlines()
h.close()
sign=0
for j in tmp[-10:]:
 if sign==1:
   [E_f,dos]=[ round(float(m),5) for m in j.split()]
   break
 if ' ******** EF and DOS at fermi level *******' in j:
   sign=1





print('Reading alm coefficients and kpoints...'),
tmp=[]
for i in range(1,64):
 try: 
  h=open('DOS/DOS.almblm_'+str(i),'r')
  tmp.extend([m.split()  for m in h.readlines() if len(m.split())!=0])
  h.close()
 except:
  break
if len(tmp)==0:
  h=open('DOS/DOS.almblm','r')
  tmp.extend([m.split()  for m in h.readlines() if len(m.split())!=0])
  h.close()


#ALM are in file only for occupied (or partially occupied) bands, so the number of bands is smaller here than in *energy files
ALMBLM=[] 
#NONEQ=[] #list of nonequiv kpoints
for i in tmp:
 if 'K-POINT' in i[0]: 
  ALMBLM.append([])
#  NONEQ.append(np.array([round(float(m),3) for m in i[1:4]]))
 elif 'ATOM' in i[1]: ALMBLM[-1].append([])
 elif 'weight' in i[-1]:  
  ALMBLM[-1][-1].append( [ [] for k in range(n_l) ] )
 elif len(i)>10:
  ALMBLM[-1][-1][-1][int(i[0])].append([complex(float(i[3]),float(i[4])),complex(float(i[6]),float(i[7]))]) 


#transpose to get: from ALMBLM[kp][at] to ALMBLM[at][kp]
ALMBLM=[[ALMBLM[i][j] for i in range(len(ALMBLM))] for j in range(n_at)]
#ALMBLM[i][k][j][l][m][0-3]  i-atoms,k-kpoint, j-band, l (orbital No), m, [0: Re[Alm]+j*Im[Alm], 1:  Re[Blm]+j*Im[Blm]]


print('Reading band energies...'),
ENE=[]
tmp=[]
print('Files:'),
kweights=[]
for i in range(1,64):
 try: 
  h=open('DOS/DOS.energy_'+str(i),'r')
  print( ' DOS/DOS.energy_'+str(i)),
  tmp.extend([m.split() for m in h.readlines()[2:]])
  h.close()
 except:
  break
if i==1:
  h=open('DOS/DOS.energy','r')
  print( ' DOS/DOS.energy'),
  tmp.extend([m.split() for m in h.readlines()[2:]])
  h.close()

#ENE[i][j] i-kpoint, j- band
for i in tmp:
  if len(i)>4 and int(i[-4])==len(ENE)+1: 
   kweights.append(float(i[-1]))
   ENE.append([])
  elif len(i)==2: 
   ENE[-1].append(i[1])
#rearrange
#et[ibnd][tetra[i][nt]]
print("; No of noneq kpoints="+str(len(ENE))), 
n_k,n_band=len(kweights),min([ len(i) for i in ENE])
n_k_total=sum(kweights)
print("; Total number of kpoints="+str(n_k_total))
#choose only energies of  which are considered in *almblm (only (fully or partially) occupied)
ENE=[ [float(ENE[k][i]) for i in range(len(ALMBLM[0][k]))] for k in range(len(ENE))]
#ENE=[ [float(ENE[j][i]) for j in range(len(ENE))] for i in range(n_band)]
print("; No of bands="+str(len(ENE[0])))
ENE_weights=[ [] for j in range(n_k)]
# weights for integrals over k (cold smearing)
#  ! cold smearing  (Marzari-Vanderbilt-DeVita-Payne)
#  if (n.eq. - 1) then
#     arg = min (200.d0, (x - 1.0d0 / sqrt (2.0d0) ) **2)
#     w0gauss = sqrtpm1 * exp ( - arg) * (2.0d0 - sqrt ( 2.0d0) * x)
#     return
#           w0g1 = w0gauss ( (ef1 - et (ibnd, ikk) ) / degauss1, ngauss1) &
#                / degauss1

for k in range(n_k):
 for iband in ENE[k]:
          eef=(E_f-iband)/degauss
          arg=min(200., (eef - 2.**(-0.5) ) **2)
          ENE_weights[k].append((kweights[k]*sqrtpm1*np.exp(-arg)*(2.-(2**0.5*eef)))/degauss)



print(" No of noneq kpoints="+str(n_k))
'''
print('Calculate intersitial part of Hopfield...'):
for K1 in range(n_K_intersitial):
 for K2 in range(n_K_intersitial):
   for

'''


print ('Calculate intrasitial part of Hopfield...')
#Vtot=[] #[atom][lm][r]
#LM=[] #[atom][list of lm 

#n_l nieograniczone daje dokladnie te same wyniki co n_l=5; n_l=4 zmienia wynik o 0.01%
n_l=3

dr=[[0]+[ RMESH[at][i]-RMESH[at][i-1] for i in range(1,nr)] for at in range(n_at)]  

def r_integral(RMESH_at,RADWF_at,Vtot_at,dr_at,l3,l4,nlm1):
      Ba_and_C_r,Aa_r=np.array([0.j,0.j]),np.array([0.j,0.j])
      Ba_and_C_r2,Aa_r2=np.array([0.j,0.j]),np.array([0.j,0.j])
      for nr1,r1 in enumerate(RMESH_at): 
       if nr1==0: continue  
       dr1=dr_at[nr1]
       '''  
       #if Vr^2 devided by r^2
       dv1dr=(Vtot_at[nlm1][nr1]-Vtot_at[nlm1][nr1-1])
       B_and_C_r+=r1*Vtot_at[nlm1][nr1]*np.conjugate(RADWF_at[l3][nr1])*RADWF_at[l4][nr1]*dr1 
       A_r+= r1**2*dv1dr *np.conjugate(RADWF_at[l3][nr1])*RADWF_at[l4][nr1] 
       B_and_C_r2+=r1*Vtot_at[nlm1][nr1]*np.conjugate(RADWF_at[l3][nr1])*RADWF_at[l4][nr1][::-1]*dr1
       A_r2+= r1**2*dv1dr *np.conjugate(RADWF_at[l3][nr1])*RADWF_at[l4][nr1][::-1]
       '''
       '''
       #if Vr^2 not devided 
       dv1dr=(Vtot_at[nlm1][nr1]-Vtot_at[nlm1][nr1-1])-2*Vtot_at[nlm1][nr1]/r1*dr1
       # every *r1 -> /r1 and *r1^2 -> 1 because v is multiplied by r^2
       B_and_C_r+=Vtot_at[nlm1][nr1]*np.conjugate(RADWF_at[l3][nr1])*RADWF_at[l4][nr1]*dr1  /r1 #*r**2
       A_r+= dv1dr *np.conjugate(RADWF_at[l3][nr1])*RADWF_at[l4][nr1] #*r**2
       B_and_C_r2+=Vtot_at[nlm1][nr1]*np.conjugate(RADWF_at[l3][nr1])*RADWF_at[l4][nr1][::-1]*dr1 /r1 #*r**2
       A_r2+= dv1dr *np.conjugate(RADWF_at[l3][nr1])*RADWF_at[l4][nr1][::-1]#*r**2
       '''
       #if Vr^2 devided by r (so we have V*r)
       dv1dr=(Vtot_at[nlm1][nr1]-Vtot_at[nlm1][nr1-1]-Vtot_at[nlm1][nr1]*dr1)*r1
       # every *r1 -> /r1 and *r1^2 -> 1 because v is multiplied by r^2
       Ba_and_C_r+=Vtot_at[nlm1][nr1]*np.conjugate(RADWF_at[l3][nr1])*RADWF_at[l4][nr1]*dr1   #*r
       Aa_r+= dv1dr *np.conjugate(RADWF_at[l3][nr1])*RADWF_at[l4][nr1] #*r*
       Ba_and_C_r2+=Vtot_at[nlm1][nr1]*np.conjugate(RADWF_at[l3][nr1])*RADWF_at[l4][nr1][::-1]*dr1 #*r
       Aa_r2+= dv1dr *np.conjugate(RADWF_at[l3][nr1])*RADWF_at[l4][nr1][::-1]#*r
      return Ba_and_C_r,Aa_r,Ba_and_C_r2,Aa_r2

print(max([max(i) for i in ENE_weights]))
def k_integral(ALMBLM_at,ENE_weights,l5,m5,l4,m4,n_k_total):
         Ia_kp,I2a_kp=np.array([0.j,0.j]),np.array([0.j,0.j])
         for kp in range(len(ALMBLM_at)):
            for iband in range(len(ALMBLM_at[kp])):
#             for jband in range(n_band):
              Ia_kp+=ENE_weights[kp][iband]*np.conjugate(ALMBLM_at[kp][iband][l5][m5+l5])*ALMBLM_at[kp][iband][l4][m4+l4]
              I2a_kp+=ENE_weights[kp][iband]*np.conjugate(ALMBLM_at[kp][iband][l5][m5+l5])*ALMBLM_at[kp][iband][l4][m4+l4][::-1]  
         return Ia_kp/n_k_total,I2a_kp/n_k_total


def Hopfield(args):
  [at,nlm1,dr,LM,RMESH,Vtot,RADWF,ALMBLM,n_k_total]=args
#def Hopfield(at,nlm1,dr,LM,RMESH,Vtot,RADWF,ALMBLM):
#  brrr=[0,0,0,0,0,0]
  print (at,nlm1,LM[at])
  Hopfield,Hopfield2=0.j,0.j
  HopfieldS,Hopfield2S=0.j,0.j #with SOC
  B_and_C_rpS,A_rpS,B_and_C_rp2S,A_rp2S,B_and_C_rS,A_rS,B_and_C_r2S,A_r2S=0.j,0.j,0.j,0.j,0.j,0.j,0.j,0.j
  do_b=0
  [l1,m1]=LM[at][nlm1]
  if m1==0: 
   Bx=0.
   do_b=0
  else: do_b=1
  for nlm2 in range(len(LM[at])):
   [l2,m2]=LM[at][nlm2]
   print (l1,': ',l2,Hopfield,Hopfield2)
   if m2==0: 
    Bx=0.
    do_b=0
   else: do_b=1
   for l3 in range(n_l):#len(LM[at])):
    for m3 in range(-l3,l3+1):#len(LM[at])):
     for l4 in range(n_l):#len(LM[at])):
      B_and_C_r,A_r,B_and_C_r2,A_r2=r_integral(RMESH[at],RADWF[at],Vtot[at],dr[at],l3,l4,nlm1)
      if SOC==1:
       B_and_C_rS,A_rS,B_and_C_r2S,A_r2S=r_integral(RMESH[at],RADWFsmall[at],Vtot[at],dr[at],l3,l4,nlm1)
      for m4 in range(-l4,l4+1):#len(LM[at])):
       if m4!=m3-m1: 
        continue 
       Ax1=three_y(l1,m1,l4,m4,l3,m3)
       if do_b: Bx1=spec_int(l1,m1,l3,m3,l4,m4)
       for l5 in range(n_l):#len(LM[at])):
        for m5 in range(-l5,l5+1):#len(LM[at])):
         I_kp,I2_kp=k_integral(ALMBLM[at],ENE_weights,l5,m5,l4,m4,n_k_total)
         for l6 in range(n_l):#len(LM[at])):
          B_and_C_rp,A_rp,B_and_C_rp2,A_rp2=r_integral(RMESH[at],RADWF[at],Vtot[at],dr[at],l5,l6,nlm2)
          if SOC==1:
           B_and_C_rpS,A_rpS,B_and_C_rp2S,A_rp2S=r_integral(RMESH[at],RADWFsmall[at],Vtot[at],dr[at],l5,l6,nlm2)
          for m6 in range(-l6,l6+1):#len(LM[at])):
           if m6!=m5-m2:  
            continue
           Ax=Ax1*three_y(l2,m2,l6,m6,l5,m5)
           if do_b: Bx=Bx1*spec_int(l2,m2,l5,m5,l6,m6)
           Cx=-m1*m2*Ax 
           I_k,I2_k=k_integral(ALMBLM[at],ENE_weights,l3,m3,l6,m6,n_k_total)
 #          print( round_complex(I_k[0],0), round_complex(B_and_C_r[0],0), round_complex(A_r[0],0), round(Bx,0), round(Ax,0))
           Hopfield+=np.dot(B_and_C_r*B_and_C_rp*(Bx+Cx)+A_r*A_rp*Ax, I_k*I_kp)
           Hopfield2+=np.dot(B_and_C_r2*B_and_C_rp2*(Bx+Cx)+A_r2*A_rp2*Ax, I2_k*I2_kp)
           HopfieldS+=np.dot(B_and_C_rS*B_and_C_rpS*(Bx+Cx)+A_rS*A_rpS*Ax, I_k*I_kp)
           Hopfield2S+=np.dot(B_and_C_r2S*B_and_C_rp2S*(Bx+Cx)+A_r2S*A_rp2S*Ax, I2_k*I2_kp)
 #          brrr1=[br.real for br in [B_and_C_r[0],Bx,A_r[0],Ax,I_k[0],I_kp[0]]]
 #          brrr=[max(brrr1[br],brrr[br]) for br in range(6)]
 #     print(l1,l2,':', Hopfield,brrr)
  print('finally: ',at,nlm1,(Hopfield),(Hopfield2))
  print('finally: ',at,nlm1,round_complex(Hopfield,5),round_complex(Hopfield2,5))
  return Hopfield,Hopfield2,HopfieldS,Hopfield2S
#ENE[i][j] i-kpoint, j- band
#RADWF[i][j][k][0-1] i-atoms, j- l (orbital No), k - r-mesh, [0-1]-large  component of radwf and udot
#ALMBLM[i][k][j][l][m][0-3]  i-atoms,k-kpoint, j-band, l (orbital No), m, [0: Re[Alm]+j*Im[Alm], 1:  Re[Blm]+j*Im[Blm]]


print('Calculate eta...')
no_of_pool=sum([len(LM[at]) for at in range(n_at)])
print ('No of pools:',no_of_pool)
if __name__ == '__main__':
 with Pool(no_of_pool) as pol:
  result=pol.map(Hopfield,[[at,nlm1,dr,LM,RMESH,Vtot,RADWF,ALMBLM,n_k_total] for at in range(n_at) for nlm1 in range(len(LM[at])) ])
  print('Hopfield, Hopfield udot, Hopfield SOC, Hopfield udot SOC:')
  for i in range(4): print(sum([m[i] for m in result])),
  print('')
  print('the same times dos*rytoev/abtoang^2')
  for i in range(4): print(sum([m[i]/dos*Ry_to_eV/ab_to_ang/ab_to_ang for m in result])),
  print('Total hopfield'),
  print(sum([sum(i) for i in result]))
  print('Total hopfield times dos*rytoev/abtoang^2'),
  print(sum([sum(i) for i in result])/dos*Ry_to_eV/ab_to_ang/ab_to_ang)

'''
if __name__ == '__main__':
 pool = Pool()
 results = [pool.apply(Hopfield, args=(at,nlm1,dr,LM,RMESH,Vtot,RADWF,ALMBLM))  for at in range(n_at) for nlm1 in range(len(LM[at]))]
 pool.close()
'''
'''
   processes = []
   for at in range(n_at):
     for nlm1 in range(len(LM[at])):
        p = Process(target=Hopfield, args=(at,nlm1,dr,LM,RMESH,Vtot,RADWF,ALMBLM,))
        processes.append(p)
        p.start()        
   for process in processes:
        process.join()
'''

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

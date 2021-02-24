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
 out1=quad(lambda tet: np.real(sph_harm(m1,l1,0,tet)*np.conj(sph_harm(m2,l2,0,tet))*sph_harm(m3,l3,0,tet)), 0.0, np.pi)[0]
 out2=quad(lambda tet: np.imag(sph_harm(m1,l1,0,tet)*np.conj(sph_harm(m2,l2,0,tet))*sph_harm(m3,l3,0,tet)), 0.0, np.pi)[0]
 if np.isnan(out1) or np.isnan(out2): 
  print('spec_int',l1,m1,l2,m2,l3,m3,':infty')
  return 0.j
 return complex(round(out1,6),round(out2,6))


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
n_r=len(Vtot[0][0]) #number of r-points
print (n_at,n_r, LM)
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

RADWF=[]   #RADWF[i][j][k][0-1] i-atoms, j- l (orbital No), k - r-mesh, [0-2]-large  component of radwf and udot, and lo if exists
RADWFsmall=[] #RADWFsmall[i][j][k][0-1] i-atoms, j- l (orbital No), k - r-mesh, [0-1]-  small component of radwf and udot, and lo if exists
mesh_info=[]

if_lo=[]
for i in tmp:
 if str(n_r) in i.split(): 
  RADWF.append([])
  RADWFsmall.append([])
#  if_lo.append([])
  mesh_info.append([float(m) for m in i.split()])
 elif len(i.split())==1: 
  RADWF[-1].append([])
  RADWFsmall[-1].append([])
#  if_lo[-1].append(0)
 else:
  try: 
   RADWF[-1][-1].append(np.array([float(i.split()[0]),float(i.split()[2])]))
   RADWFsmall[-1][-1].append(np.array([CIN*float(i.split()[1]),CIN*float(i.split()[3])]))
  except:continue
 # try: 
#   if_lo[-1]=1
#   RADWF[-1][-1][2]=float(i.split()[4])+float(i.split()[6])+float(i.split()[8])
 #  RADWFsmall[-1][-1][2]=float(i.split()[5])+float(i.split()[7])+float(i.split()[9])
 # except:continue

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
 for j in range(n_r):
  RMESH[-1].append(r0*np.exp(j*dx))
 print('RMT='+str(rmt)+' should ='+str(RMESH[-1][-1]))

RMESH=[ np.array([m for m in i]) for i in RMESH]
DR=[np.diff(RMESH[at]) for at in range(n_at)] #len=n_r-1 (element 0 is excluded in diff) 


#oryginalnie w pliku *vtotal jest zapisany V*R^2, wiec dzielimy przez R by dostac czysty V*R  
Vtot=[ [np.array([Vtot[i][j][k]/(RMESH[i][k]) for k in range(len(Vtot[i][j]))]) for j in range(len(Vtot[i]))] for i in range(len(Vtot))] 
DVDR=[ [  (np.gradient(Vtot[at][nlm1])[1:]-Vtot[at][nlm1][1:]*DR[at])*RMESH[at][1:] for nlm1 in range(len(LM[at]))] for at in range(n_at)]


#RADWF[i][j][k][0-1] i-atoms, j- l (orbital No), k - r-mesh, [0-1]-large  component of radwf and udot
#oryginalnie w pliku *radwf jest zapisany radwf*R, wiec dzielimy przez R by dostac czysty radwf
#dodatkowo z RADWF[at][l][r][u,udot] -> RADWF[at][l][u,udot][r]
RADWFsmall=    [[np.transpose(np.array([ [RADWFsmall[i][j][r][0]/RMESH[i][r],RADWFsmall[i][j][r][1]] for r in range(1,len(RADWF[i][j]))])) for j in range(len(RADWF[i]))] for i in range(len(RADWF))]
RADWF=         [[np.transpose(np.array([ [RADWF[i][j][r][0]/RMESH[i][r],     RADWF[i][j][r][1]]      for r in range(1,len(RADWF[i][j]))])) for j in range(len(RADWF[i]))] for i in range(len(RADWF))]


RMESH=[ i[1:] for i in RMESH]
Vtot=[ [ Vtot[i][j][1:] for j in range(len(Vtot[i]))] for i in range(len(Vtot))] 
Vtot_times_dr=[ [Vtot[i][j]*DR[i]  for j in range(len(Vtot[i]))] for i in range(len(Vtot))]

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

n_k=len(ALMBLM)
#transpose to get: from ALMBLM[kp][at] to ALMBLM[at][kp]
ALMBLM0=[[ALMBLM[i][j] for i in range(len(ALMBLM))] for j in range(n_at)]
#ALMBLM[i][k][j][l][m][0-3]  i-atoms,k-kpoint, j-band, l (orbital No), m, [0: Re[Alm]+j*Im[Alm], 1:  Re[Blm]+j*Im[Blm]]


#n_l nieograniczone daje dokladnie te same wyniki co n_l=5; n_l=4 zmienia wynik o 0.01%
n_l=4

#alm and blm as numpy arrays
#ALMBLM[i][k][j][l][m][0-3] -> ALM[i][l][m][k][j]
ALM=[[[[np.array([ALMBLM0[i][k][j][l][m][0] for j in range(len(ALMBLM0[i][k]))]) for k in range(n_k) ] for m in range(2*l+1) ] for l in range(n_l) ] for i in range(n_at)]
BLM=[[[[np.array([ALMBLM0[i][k][j][l][m][1] for j in range(len(ALMBLM0[i][k]))]) for k in range(n_k) ] for m in range(2*l+1) ] for l in range(n_l) ] for i in range(n_at)]
ALMBLM=[[[[np.array([ALMBLM0[i][k][j][l][m] for j in range(len(ALMBLM0[i][k]))]) for k in range(n_k) ] for m in range(2*l+1) ] for l in range(n_l) ] for i in range(n_at)]


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
n_k2,n_band=len(kweights),min([ len(i) for i in ENE])
n_k_total=sum(kweights)
print("; No of noneq kpoints="+str(n_k))
if n_k2!=n_k: raise ValueError('no of k_points of energies:',n_k,' and of almblm:',str(n_k2),' is not the same')
print("; Total number of kpoints="+str(n_k_total))

print(' From all'+str(max([len(k) for k in ENE]))),
#choose only energies of  which are considered in *almblm (only (fully or partially) occupied)
ENE=[ [float(ENE[k][i]) for i in range(len(ALMBLM0[0][k]))] for k in range(len(ENE))]
print(' only'+str(max([len(k) for k in ENE]))+' are fully or partially occupied')

print('From all '+str(max([len(k) for k in ENE]))),
##which bands nearly cross the FS
who_cross=[]
for k in range(len(ENE)):
 for i in range(len(ENE[k])):
  if i not in who_cross and abs(ENE[k][i]-E_f)<1/13.606: who_cross.append(i)
minn,maxx=min(who_cross),max(who_cross)   
for k in range(len(ENE)):
 ENE[k]=ENE[k][minn:(maxx+1)]
print(' only '+str(max([len(k) for k in ENE]))+' cross EF or are at least 1 ev close')

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






def r_integral(RADWF_at_l3,RADWF_at_l4,Vtot_at,DVDR1,dr_at):
      Ba_and_C_r=np.dot((Vtot_at*RADWF_at_l3),(RADWF_at_l4*dr_at))
      Aa_r=np.dot((DVDR1*RADWF_at_l3)*RADWF_at_l4)
      print ('r integral')
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
      return Ba_and_C_r,Aa_r


def Hopfield(args):
  [at,nlm1,DR,LM,RMESH,Vtot_times_dr,DVDR,RADWF,ALMBLM,n_k,n_k_total]=args
#def Hopfield(at,nlm1,dr,LM,RMESH,Vtot,RADWF,ALMBLM):
  print (at,nlm1,LM[at])
  Hopfield,Hopfield2=0.j,0.j
  HopfieldS,Hopfield2S=0.j,0.j #with SOC
  Ax,Bx,Cx,B_and_C_rp,A_rp,B_and_C_r,A_r=0.j,0.j,0.j,0.j,0.j,0.j,0.j
  [l1,m1]=LM[at][nlm1]
  for nlm2 in range(len(LM[at])):
   [l2,m2]=LM[at][nlm2]
   print (l1,': ',l2,Hopfield,Hopfield2)
   Bx,Cx=0.,0.
   if m1==0 or m2==0:    do_b=0
   else:     do_b=1
   if l1==0 or l2==0:   do_c=0
   else: do_c=1
   for l3 in range(n_l):#len(LM[at])):
    for m3 in range(-l3,l3+1):#len(LM[at])):
     for l4 in range(n_l):#len(LM[at])):
#      for m4 in range(-l4,l4+1):#len(LM[at])):
#       if m4!=m3-m1: 
#        continue 
       m4=m3-m1
       if abs(m4)>abs(l4): continue
       Ax1=three_y(l1,m1,l4,m4,l3,m3)
       R34=np.array([ np.conjugate(RADWF[at][l3][k])*RADWF[at][l4][m] for k in range(len(RADWF[at][l3])) for m in range(len(RADWF[at][l4]))]) #array[4][nr]
       A_r=np.sum(DVDR[at][nlm1]*R34,1)
       print ('sh',R34.shape,A_r,A_r.shape)
       if do_c: 
        Cx=4*np.pi**2*l1*(l1+1)/(((2*l1+1)*(2*l1+3))**0.5)*spec_int(l1+1,m1,l3,m3,l4,m4) 
        if l1!=1 and l1!=m1: 
         Cx+= -4*np.pi**2*l1*(l1-1)/(((2*l1-1)*(2*l1+1))**0.5)*spec_int(l1-1,m1,l3,m3,l4,m4)
       if do_b: Bx1=(-2*np.pi*m1)*spec_int(l1,m1,l3,m3,l4,m4)
       if do_b or do_c:
        B_and_C_r= np.sum(Vtot_times_dr[at][nlm1]*R34,1)
       for l5 in range(n_l):#len(LM[at])):
        for m5 in range(-l5,l5+1):#len(LM[at])):
         print('l6 ended')
         for l6 in range(n_l):#len(LM[at])):
#          for m6 in range(-l6,l6+1):#len(LM[at])):
#           if m6!=m5-m2:  
#            continue
           m6=m5-m2
           if abs(m6)>abs(l6): continue
           Ax=Ax1*three_y(l2,m2,l6,m6,l5,m5)
           R56=np.array([ np.conjugate(RADWF[at][l5][k])*RADWF[at][l6][m] for k in range(len(RADWF[at][l5])) for m in range(len(RADWF[at][l6]))])
           A_rp=      np.sum(DVDR[at][nlm2]*R56,1)
           if do_b: Bx=Bx1*(-2*np.pi*m2)*spec_int(l2,m2,l5,m5,l6,m6)
           if do_c: Cx=Cx*(l2*(l2+1)/(((2*l2+1)*(2*l2+3))**0.5)*spec_int(l2+1,m2,l5,m5,l6,m6) \
              -l2*(l2-1)/(((2*l2-1)*(2*l2+1))**0.5)*spec_int(l2-1,m2,l5,m5,l6,m6))
           if do_b or do_c:
            B_and_C_rp=np.sum(Vtot_times_dr[at][nlm2]*R56,1)     
           ww=0  
           B_and_C_x=Bx+Cx
           for k in range(n_k):
            for iband in  range(len(ENE[k])):
             for kp in range(n_k):
              for ibandp in range(len(ENE[kp])):
#               almblm_r=np.kron(np.conjugate(ALMBLM[at][l3][m3][k][iband]),ALMBLM[at][l4][m4][kp][ibandp])
#               almblm_rp=np.kron(np.conjugate(ALMBLM[at][l5][m5][kp][ibandp]),ALMBLM[at][l6][m6][k][iband])
               almblm_r=np.array([w1*w2 for w1 in np.conjugate(ALMBLM[at][l3][m3][k][iband]) for w2 in ALMBLM[at][l4][m4][kp][ibandp]])
               almblm_rp=np.array([w1*w2 for w1 in np.conjugate(ALMBLM[at][l5][m5][kp][ibandp]) for w2 in ALMBLM[at][l6][m6][k][iband]])
               A=np.dot(almblm_r,A_r) *np.dot(almblm_rp,A_rp)
#               A=sum([almblm_r[i]*A_r[i] for i in range(len(A_r))])*sum([almblm_rp[i]*A_rp[i] for i in range(len(A_rp))])
               B_and_C=0.j
               if do_b or do_c: 
                B_and_C=np.dot(almblm_r,B_and_C_r)*np.dot(almblm_r,B_and_C_rp)
#                B_and_C=sum([almblm_r[i]*B_and_C_r[i] for i in range(len(B_and_C_r))])*sum([almblm_rp[i]*B_and_C_rp[i] for i in range(len(B_and_C_rp))])
                #r_integral(R3,R4,Vtot[at][nlm1][1:],DVDR1,dr[at])
                #r_integral(R5,R6,Vtot[at][nlm2],DVDR2,dr[at])
               Hopfield+=ENE_weights[k][iband]*ENE_weights[kp][ibandp]*(Ax*A+B_and_C*B_and_C_x)
           print('integrals ended') 
#             print('integral ended') 
  print('finally: ',at,nlm1,(Hopfield),(Hopfield2))
  print('finally: ',at,nlm1,round_complex(Hopfield,5),round_complex(Hopfield2,5))
  return Hopfield,Hopfield2,HopfieldS,Hopfield2S
#ENE[i][j] i-kpoint, j- band
#RADWF[i][j][k][0-1] i-atoms, j- l (orbital No), k - r-mesh, [0-1]-large  component of radwf and udot
#ALMBLM[i][k][j][l][m][0-3]  i-atoms,k-kpoint, j-band, l (orbital No), m, [0: Re[Alm]+j*Im[Alm], 1:  Re[Blm]+j*Im[Blm]]
# ALM[at][l][m][k][j]

print('Calculate eta...')
no_of_pool=sum([len(LM[at]) for at in range(n_at)])
print ('No of pools:',no_of_pool)
if __name__ == '__main__':
 with Pool(no_of_pool) as pol:
  result=pol.map(Hopfield,[[at,nlm1,DR,LM,RMESH,Vtot_times_dr,DVDR,RADWF,ALMBLM,n_k,n_k_total] for at in range(n_at) for nlm1 in range(len(LM[at])) ])
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

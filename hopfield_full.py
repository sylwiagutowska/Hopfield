import numpy as np
import os
import glob
import hopfield_tetra
#from math import factorial as sil
from scipy.special import sph_harm
from scipy.integrate import quad
from scipy.misc import derivative
sqrtpm1=(np.pi)**(-0.5)
degauss=0.02
#############FOR Nb eta=7.627eV/ang^2 [papaconsta.. 2002] or 4.7 eV/ang^2 (allen-dynes)

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


################RADWF AND V are already multiplied by r
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

prefix='Al' #raw_input('Prefix: ')

Ry_to_eV=13.605693122994
ab_to_ang=0.529177210903

#j_to_ev=1/(1.602176634)*1e19
#hbar=6.582119569e-16 #ev/s
#hbar=hbar/13.6056980659 #ry/s
pi=3.141592653589793238462643
#Na=6.02214076e23

SOC=0
if SOC==1: CIN=1./137.0359895**2 #for treating big and small component of wave function
else: CIN=4.0*1e-22


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

print('Read total potential...')
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
   for j in range(len(m)/19):
    Vtot[-1][-1].append(float(m[j*19:(j+1)*19]))

n_at=len(Vtot) #number of atoms
nr=len(Vtot[0][0]) #number of r-points
print n_at,nr, LM


print('Reading radial wave funtions...')
h=open('RADWF/RADWF.radwf','r')
tmp=h.readlines()
h.close()

RADWF=[]   #RADWF[i][j][k][0-3] i-atoms, j- l (orbital No), k - r-mesh, [0-1]-large or small component of radwf, [2-3]- large or small component of udot
RADWFsmall=[]
mesh_info=[]
for i in tmp:
 if str(nr) in i.split(): 
  RADWF.append([])
  RADWFsmall.append([])
  mesh_info.append([float(m) for m in i.split()])
 elif len(i.split())==1: RADWF[-1].append([])
 else:
  try: 
   RADWF[-1][-1].append(np.array([float(i.split()[0]),float(i.split()[2])]))
   RADWFsmall[-1][-1].append(np.array([float(i.split()[1]),float(i.split()[3])]))
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
#oryginalnie w pliku *vtotal jest zapisany V*R^2, wiec dzielimy przez R by dostac czysty V
### for L=0 M=0 (czesc sferyczna) Vtot jest rowny V.vsp [zobacz plik V_0.dat]
Vtot=[ [[Vtot[i][j][k]/(RMESH[i][k]**2)/Z_of_atoms[i]*2 for k in range(len(Vtot[i][j]))] for j in range(len(Vtot[i]))] for i in range(len(Vtot))]   

print('Write spherical V and radwf to files V_i.dat and RADWF_i.dat...')
for i in range(len(Vtot)):
 h=open('V_'+str(i)+'.dat','w')
 for j in range(len(Vtot[i][0])):
  h.write(str(RMESH[i][j])+' '+str(Vtot[i][0][j])+' ')
  h.write('\n')
 h.close()

for i in range(len(RADWF)):
 h=open('RADWF_'+str(i)+'.dat','w')
 h.write('# r , l=0- big component,small component, l=1- big component,small component....\n')
 for j in range(len(RADWF[i])):
  h.write(str(RMESH[i][j])+' ')
  for k in range(len(RADWF[i][j])):
   h.write(str(RADWF[i][j][k][0])+' '+str(RADWF[i][j][k][1])+' ')
  h.write('\n')
 h.close() 

print('Read EF and total DOS...'),

h=open('DOS/tot/DOS.outputt','r')
tmp=h.readlines()
h.close()
sign=0
for j in tmp[-10:]:
 if sign==1:
   TOTDOS=[ round(float(m),5) for m in j.split()]
   break
 if ' ******** EF and DOS at fermi level *******' in j:
   sign=1
for j in TOTDOS:
  print(str(j)+' '),
print(' ')

E_f=TOTDOS[0]

print('Reading band energies...'),
ENE=[]
tmp=[]
print('Files:'),
kweights=[]
for i in range(1,64):
 try: 
  h=open('DOS/DOS.energy_'+str(i),'r')
  print ' DOS/DOS.energy_'+str(i),
  tmp.extend([m.split() for m in h.readlines()[2:]])
  h.close()
 except:
  break
if i==1:
  h=open('DOS/DOS.energy','r')
  print ' DOS/DOS.energy',
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
ENE=[ [float(ENE[j][i]) for i in range(n_band)] for j in range(len(ENE))]
#ENE=[ [float(ENE[j][i]) for j in range(len(ENE))] for i in range(n_band)]
print("; No of bands="+str(len(ENE[0])))
ENE_weights=[ [0 for i in range(n_band)] for j in range(n_k)]
# weights for integrals over k (cold smearing)
for k in range(n_k):
 for iband in range(n_band):
          eef=E_f-ENE[k][iband]
          arg=min(200., (eef - 2.**(-0.5) ) **2)
          ENE_weights[k][iband]=(kweights[k]*sqrtpm1*np.exp(-arg)*(2.-(2**0.5*eef)))
'''
ENE2=[ [float(ENE[j][i]) for j in range(len(ENE)) ] for i in range(n_band)]

print('Which bands cross Fermi level: '),
band_f=[]
for ni,i in enumerate(ENE2):
 old=i[0]
 is_fermi=0
 for j in i:
  if (j-E_f)*(old-E_f)<0:
   is_fermi=1
   band_f.append(ni)
   break

n_band=len(band_f)
print n_band,'bands:',band_f
ENE=[ [float(ENE[j][i]) for i in band_f] for j in range(len(ENE))]
'''


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

#ALMBLM[k][i][j][l][m][0-3] k-kpoint, i-atoms, j-band, l (orbital No), m, [0: Re[Alm]+j*Im[Alm], 1:  Re[Blm]+j*Im[Blm]]
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



'''
h=open('kpoints.dat','w')
for i in NONEQ:
 for j in i:
  h.write(str(j)+' ')
 h.write('\n')
h.close()
'''

#rearrange
#from: #ALM[k][i][j][l] k-kpoint, i-atoms, j-band, l (orbital No) (summed over m)
#to: weights[na][ibnd][kpoint]

#ALM=[ [[[ALM[l][i][k][j] for l in range(n_k)] for k in range(n_band)] for j in range(n_l)] for i in range(n_at)]
print(" No of noneq kpoints="+str(n_k))
 





print('Call spherical harmonics')
#YLM[at][lm][nfi][ntet]
n_fi,n_tet=20,10
YLM=[[] for i in range(n_l)]
mesh=np.meshgrid(np.linspace(0.,pi,n_tet),np.linspace(0.,2.*pi,n_fi))
for l in range(n_l):
 for m in range(-l,l+1):
#  YLM[l].append([])
#  for i in range(n_fi):
#   YLM[l][-1].append([])
#   for j in range(n_tet):
    YLM[l].append(sph_harm(m,l,mesh[1],mesh[0]))
print len(YLM[3]), len(YLM[3][0]), (YLM[3][0][0])


#all_spec_int=[[[[[[spec_int(l1,m1,l2,m2,l3,m3,n_fi,n_tet,YLM) for m3 in range(-l3,l3+1)] for l3 in range(n_l)] for m2 in range(-l2,l2+1)] for l2 in range(n_l)] for m1 in range(-l1,l1+1)] for l1 in range(n_l)]


#Vtot=[] #[atom][lm][r]
#LM=[] #[atom][list of lm 
n_l=5
Hopfield,Hopfield2=0.,0.

dtet=pi/n_tet
dfi=2*pi/n_fi
sintet=np.sin(np.linspace(0,pi,n_tet))
dr=[[0]+[ RMESH[at][i]-RMESH[at][i-1] for i in range(1,nr)] for at in range(n_at)]
print len(dr[0])
print('Calculate eta...')
k,kp,iband,jband=0,0,0,0
for at in range(1):#n_at):
 do_b=0
 for nlm1 in range(1,len(LM[at])):
  print nlm1
  [l1,m1]=LM[at][nlm1]
  if m1==0: 
   Bx=0.
   do_b=0
  else: do_b=1
  for nlm2 in range(1,len(LM[at])):
   [l2,m2]=LM[at][nlm2]
   if m2==0: 
    Bx=0.
    do_b=0
   else: do_b=1
   for l3 in range(6):#len(LM[at])):
    for m3 in range(-l3,l3+1):#len(LM[at])):
     for l4 in range(6):#len(LM[at])):
      print l1,l2,l3,l4,Hopfield,Hopfield2
      for m4 in range(-l4,l4+1):#len(LM[at])):

       if m4==m3-m1:  
        Ax1=three_y(l1,m1,l4,m4,l3,m3)
        if do_b: Bx1=spec_int(l1,m1,l3,m3,l4,m4)
       else: 
        continue

       B_and_C_r,A_r=np.array([0.j,0.j]),np.array([0.j,0.j])
       B_and_C_r2,A_r2=np.array([0.j,0.j]),np.array([0.j,0.j])
       for nr1,r1 in enumerate(RMESH[at]):   
        dv1dr=Vtot[at][nlm1][nr1]-Vtot[at][nlm1][nr1-1]
        dr1=dr[at][nr1]
        B_and_C_r+=r1*Vtot[at][nlm1][nr1]*np.conjugate(RADWF[at][l3][nr1])*RADWF[at][l4][nr1]*dr1
        A_r+=r1**2. * dv1dr *np.conjugate(RADWF[at][l3][nr1])*RADWF[at][l4][nr1]
        B_and_C_r2+=r1*Vtot[at][nlm1][nr1]*np.conjugate(RADWF[at][l3][nr1])*RADWF[at][l4][nr1][::-1]*dr1
        A_r2+=r1**2. * dv1dr *np.conjugate(RADWF[at][l3][nr1])*RADWF[at][l4][nr1][::-1]      
       #wsp=np.conjugate(ALMBLM[k][at][iband][l3][m3])*ALMBLM[kp][at][jband][l4][m4]
       #B_and_C_r *=wsp
       #A_r *=wsp

       for l5 in range(6):#len(LM[at])):
        for m5 in range(-l5,l5+1):#len(LM[at])):

         I_kp,I2_kp=np.array([0.j,0.j]),np.array([0.j,0.j])
         for kp in range(len(ALMBLM)):
            for iband in range(len(ALMBLM[kp])):
#             for jband in range(n_band):
              I_kp+=ENE_weights[kp][iband]*np.conjugate(ALMBLM[kp][at][iband][l5][m5+l5])*ALMBLM[kp][at][iband][l4][m4+l4]
              I2_kp+=ENE_weights[kp][iband]*np.conjugate(ALMBLM[kp][at][iband][l5][m5+l5])*ALMBLM[kp][at][iband][l4][m4+l4][::-1]  


         for l6 in range(6):#len(LM[at])):
          for m6 in range(-l6,l6+1):#len(LM[at])):

           if m6==m5-m2:  
            Ax=Ax1*three_y(l2,m2,l6,m6,l5,m5)
	    if do_b: Bx=Bx1*spec_int(l2,m2,l5,m5,l6,m6)
            Cx=-m1*m2*Ax 
           else: 
             continue

	   B_and_C_rp,A_rp=np.array([0.j,0.j]),np.array([0.j,0.j])
	   B_and_C_rp2,A_rp2=np.array([0.j,0.j]),np.array([0.j,0.j])
           for nr1,r1 in enumerate(RMESH[at]):
            dv1dr=Vtot[at][nlm2][nr1]-Vtot[at][nlm2][nr1-1]
            dr1=dr[at][nr1]
            B_and_C_rp+=r1*Vtot[at][nlm2][nr1]*np.conjugate(RADWF[at][l5][nr1])*RADWF[at][l6][nr1]*dr1
            A_rp+=r1**2. * dv1dr *np.conjugate(RADWF[at][l5][nr1])*RADWF[at][l6][nr1]
            B_and_C_rp2+=r1*Vtot[at][nlm2][nr1]*np.conjugate(RADWF[at][l5][nr1])*RADWF[at][l6][nr1][::-1]*dr1
            A_rp2+=r1**2. * dv1dr *np.conjugate(RADWF[at][l5][nr1])*RADWF[at][l6][nr1][::-1]
	    
	   I_k,I2_k=np.array([0.j,0.j]),np.array([0.j,0.j])
           for k in range(len(ALMBLM)):
            for iband in range(len(ALMBLM[k])):
#             for jband in range(n_band):
              I_k +=ENE_weights[k][iband]*np.conjugate(ALMBLM[k][at][iband][l3][m3+l3])*ALMBLM[k][at][iband][l6][m6+l6]
              I2_k+=ENE_weights[k][iband]*np.conjugate(ALMBLM[k][at][iband][l3][m3+l3])*ALMBLM[k][at][iband][l6][m6+l6][::-1]

           Hopfield+=np.dot(B_and_C_r*B_and_C_rp*(Bx+Cx)+A_r*A_rp*Ax, I_k*I_kp)
           Hopfield2+=np.dot(B_and_C_r2*B_and_C_rp2*(Bx+Cx)+A_r2*A_rp2*Ax, I2_k*I2_kp)


#ENE[i][j] i-kpoint, j- band
#RADWF[h][i][j][k][0-1] h-kpoint, i-atoms, j- l (orbital No), k - r-mesh, [0]-radwf, [1]-udot
#ALMBLM[k][i][j][l][m][0-3] k-kpoint, i-atoms, j-band, l (orbital No), m, [0: Re[Alm]+j*Im[Alm], 1:  Re[Blm]+j*Im[Blm]]



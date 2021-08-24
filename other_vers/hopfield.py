import numpy as np
import os
import glob
from scipy.integrate import simps
from scipy.special import hankel1,jv
from scipy.interpolate import interp1d
#import hopfield_tetra
#############FOR Nb eta=7.627eV/ang^2 [papaconsta.. 2002] or 4.7 eV/ang^2 (allen-dynes)

################RADWF AND V are already multiplied by r


def interpolate_V(RMESH0,V):
 nib=80
 nbl=10
 r0=(RMESH0[0])
 x2=(RMESH0[-1]-r0)/(nib*(2.**nbl-1.))
 while x2<r0:
  nbl=nbl-1
  x2=(RMESH0[-1])/(nib*(2.**nbl)-1.) 
 print nbl
 dx=x2
 RMESH=[]
 xx=0
 for i in range(nbl):
  for j in range(nib):
   xx+=dx
   RMESH.append(xx)
  dx=dx+dx
 print RMESH0[-1],RMESH[-1]
 Vf=interp1d(RMESH0,V)
 V2=Vf(RMESH)
 h=open('lipmanV.dat','w')
 for i in range(len(RMESH)):
  h.write(str(RMESH[i])+' '+str(V2[i])+'\n')
 h.close()
 return np.array(RMESH),V2


def lipman(l,RMESH,Rstart,Ef,V2,Z_of_atom): #schrea.f from kkr
 nib=80
 VIN=V2-2*Z_of_atom/RMESH+l*(l+1)/RMESH/RMESH #[V[r]+(-2*Z_of_atom+l*(l+1)/RMESH[r])/(RMESH[r]) for r in range(len(RMESH))]
 A=[1,-2*Z_of_atom/(2*l+2),0]
 A[2]=(-2*Z_of_atom*A[1]+(V2[0]-Ef)*A[0])/(4*l+6)
 i1=6
 SSR=np.zeros(len(RMESH))
 print(len(SSR))
 for i in range(i1+1):
  SSR[i]=RMESH[i]**(l+1)*(A[0]+RMESH[i]*(A[1]+RMESH[i]*A[2]))
 dx=RMESH[i1]-RMESH[i1-1]
 dx=dx*dx/12
 Y=np.array([1+dx*(Ef-VIN[i1-1]),1+dx*(Ef-VIN[i1]),0])
 P=np.array([SSR[i1-1],SSR[i1],0])
 ia=i1+1
 ib=len(RMESH)
 idnint=i1+1
 for i in range(ia,ib):
  Y[2]=1+dx*(Ef-VIN[i])
  P[2]=((12-10*Y[1])*P[1]-Y[0]*P[0])/Y[2]
  Y[0]=Y[1]
  Y[1]=Y[2]
  P[0]=P[1]
  P[1]=P[2]
  SSR[i]=P[2]
  idnint=idnint+1
  if idnint==nib :
   idnint=0
#   dx=(RMESH[i]-RMESH[i-1])*(RMESH[i]-RMESH[i-1])/12
   dx=4*dx
   Y[1]=1+dx*(Ef-VIN[i])
   Y[0]=1+dx*(Ef-VIN[i-2])
   P[1]=SSR[i]
   P[0]=SSR[i-2]   
 SSR=SSR/(simps(SSR*SSR,RMESH)**0.5)


 return SSR

prefix='hopfield_calc' #raw_input('Prefix: ')
rmesh_multiply=1


Ry_to_eV=13.605693122994
ab_to_ang=0.529177210903

#j_to_ev=1/(1.602176634)*1e19
#hbar=6.582119569e-16 #ev/s
#hbar=hbar/13.6056980659 #ry/s
pi=3.141592653589793238462643
#Na=6.02214076e23

SOC=1
if SOC==1: CIN=1./137.0359895**2 #for treating big and small component of wave function
else: CIN=4.0*1e-22


print('Reading EF and Volume...'),
os.system('grep ":FER " '+'hopfield_calc/hopfield_calc.scf2 |head -n 1 >hopfield_tmp')
os.system('grep ":VZERO" '+'hopfield_calc/hopfield_calc.scf0 |head -n 1 >>hopfield_tmp')
os.system('grep ":NOE " '+'hopfield_calc/hopfield_calc.scf2 |head -n 1 >>hopfield_tmp')
os.system('grep ":VOL " '+'hopfield_calc/hopfield_calc.scf0 |head -n 1 >>hopfield_tmp')
os.system('grep ":VOL " '+'hopfield_calc/hopfield_calc.scf |head -n 1 >>hopfield_tmp')
h=open('hopfield_tmp')
tmp=h.readlines()
h.close()
#os.system('rm hopfield_tmp')
EF=float(tmp[0].split()[-1])
vzero=float(tmp[1].split()[1])
NOE=float(tmp[2].split()[-1])
volume=float(tmp[-1].split()[-1])

print(EF,volume)

print('Reading Z of atoms...'),
os.system('grep "Z: " '+'hopfield_calc/'+prefix+'.struct  >hopfield_tmp')
h=open('hopfield_tmp')
tmp=h.readlines()
h.close()
os.system('rm hopfield_tmp')
Z_of_atoms=[float(i.split()[-1]) for i in tmp]
print(Z_of_atoms)


print('Reading potential...')
h=open('hopfield_calc/hopfield_calc.vsp','r') #stored as V*r
tmp=h.readlines()
h.close()

V=[] # V[i][j] i-th atm, j-th r-point
for i in tmp:
  if 'ATOMNUMBER' in i:
   V.append([])
  elif 'VLM(R) FOR L  0   M= 0' in i: continue
  elif 'VLM(R)' in i: break
  else:
   try: V[-1].append(float(i[0:22]))
   except: continue
   try: V[-1].append(float(i[22:41]))
   except: continue
   try: V[-1].append(float(i[41:60]))
   except: continue
   try: V[-1].append(float(i[60:]))
   except: continue

print('Reading potential...')
h=open('hopfield_calc/'+prefix+'.vcoul','r') #stored as V*r
tmp=h.readlines()
h.close()

Vcoul=[] # V[i][j] i-th atm, j-th r-point
for i in tmp:
  if 'ATOMNUMBER' in i:
   Vcoul.append([])
  elif 'VLM(R) FOR L  0   M= 0' in i: continue
  elif 'VLM(R)' in i: break
  else:
   try: Vcoul[-1].append(float(i[0:22]))
   except: continue
   try: Vcoul[-1].append(float(i[22:41]))
   except: continue
   try: Vcoul[-1].append(float(i[41:60]))
   except: continue
   try: Vcoul[-1].append(float(i[60:]))
   except: continue

print('Reading potential...')
h=open('RADWF/RADWF.r2v','r') #stored as V*r
tmp=h.readlines()
h.close()

Vr2v=[] # V[i][j] i-th atm, j-th r-point
for i in tmp:
  if 'ATOMNUMBER' in i:
   Vr2v.append([])
  elif 'VLM(R) FOR L  0   M= 0' in i: continue
  elif 'VLM(R)' in i: break
  else:
   try: Vr2v[-1].append(float(i[0:22]))
   except: continue
   try: Vr2v[-1].append(float(i[22:41]))
   except: continue
   try: Vr2v[-1].append(float(i[41:60]))
   except: continue
   try: Vr2v[-1].append(float(i[60:]))
   except: continue

print('number of atoms:'+str(len(V)))
print('number of r points:'+str(len(V[0])))
na=len(V) #number of atoms
nr=len(V[0]) #number of r-points
###devide V by Q
###and multiply by 2 to obtain Hartree -> Ry
#V=[ [V[i][j]/Z_of_atoms[i]*2 for j in range(len(V[i]))] for i in  range(len(V))]   
print('Reading radial wave funtions...')
h=open('RADWF/RADWF.radwf','r')
tmp=h.readlines()
h.close()

RADWF=[]   #RADWF[i][j][k][0/1] h-kpoint, i-atoms, j- l (orbital No), k - r-mesh, [0/1]-large or small component
mesh_info=[]
for i in tmp:
 if str(nr) in i.split(): 
  RADWF.append([])
  mesh_info.append([float(m) for m in i.split()])
 elif len(i.split())==1: RADWF[-1].append([])
 else:
  try: RADWF[-1][-1].append([float(i.split()[0]),float(i.split()[1])])
  except:continue
lmin=min([ len(i) for i in RADWF])

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
 dx=np.log(rmt/r0)/(nr-1)
 r0=r0*rmesh_multiply
 RMESH.append([])
 for j in range(nr):
  RMESH[-1].append(r0*np.exp(j*dx))
 print('RMT='+str(rmt)+' should ='+str(RMESH[-1][-1]))
# for j in range(nr):
#  RMESH[-1][j]=RMESH[-1][j]/RMESH[-1][-1]*rmt


print ('vzero=',vzero)
print ('Z=',Z_of_atoms)
print ('NOE=',NOE)
nat=len(RMESH)
nstart=[]
for i in range(nat): 
 RMESH[i]=np.array(RMESH[i]) #np.array([ RMESH[i][j] for j in range(len(RMESH[i])) if RMESH[i][j]>0.00177])
 V[i]= np.array(V[i])


RADWF_11= [ [ np.array([m[0] for m in n]) for n in p] for p in RADWF] #large component
RADWF_22= [ [ [m[1] for m in n] for n in p] for p in RADWF] #small component


RADWF_1=[ np.transpose(np.array(m)) for m in RADWF_11]  #RADWF_1[i][j][k] i-atom,j-rmesh, k=l (orbital No)
RADWF_2=[ np.transpose(np.array(m)) for m in RADWF_22]
for i in range(len(RADWF_1)):
 h=open('RADWF_'+str(i)+'.dat','w')
 h.write('# r , l=0- big component,small component, l=1- big component,small component....\n')
 for j in range(len(RADWF_1[i])):
  h.write(str(RMESH[i][j])+' ')
  for k in range(len(RADWF_1[i][j])):
   h.write(str(RADWF_1[i][j][k])+' ')
  h.write('\n')
 h.close() 



V= [ V[i]+2*Z_of_atoms[i] for i in range(nat)]
V= [ V[i]/RMESH[i] for i in range(nat)]

#V=[ V[i]+(Z_of_atoms[i])*2 for i in range(nat)]
print('Write to files V_i.dat and RADWF_i.dat...')
for i in range(len(V)):
 h=open('V_'+str(i)+'.dat','w')
 for j in range(len(V[i])):
  h.write(str(RMESH[i][j])+' '+str(V[i][j])+' ')
  h.write('\n')
 h.close()

#RMESH[0][0]=0
#V[0][0]=0
RMESH[0],V[0]=interpolate_V(RMESH[0],V[0])

for l in range(len(RADWF[0])):
 RADWF[0][l]=lipman(l,RMESH[0],RADWF_11[0][l],EF,V[0],Z_of_atoms[0])
nr=len(V[0])

h=open('lipman.dat','w')
for i in range(len(RMESH[0])):
 h.write(str(RMESH[0][i])+' ')
 for l in range(len(RADWF[0])):
  h.write(str((RADWF[0][l][i]))+' ')
 h.write('\n')
h.close()


DVDR=[   (np.gradient(V[at],RMESH[at],edge_order=2)+2*Z_of_atoms[at]/(RMESH[at]*RMESH[at]))  for at in range(nat)]




'''
V2= [ V[i]/RMESH[i]+2*Z_of_atoms[i]/RMESH[i] for i in range(nat)]
DVDR=[   (np.gradient(V2[at],RMESH[at],edge_order=2)+ 2*Z_of_atoms[at]/RMESH[at]/RMESH[at])  for at in range(nat)]
'''


for i in range(len(V)):
 h=open('DVDR_'+str(i)+'.dat','w')
 for j in range(len(V[i])):
  h.write(str(RMESH[i][j])+' '+str(DVDR[i][j])+' ')
  h.write('\n')
 




print('Read lDOS...')
LDOS=[]
for i in range(na):
 h=open('DOS/ATOM'+str(i+1)+'/DOS.outputt','r')
 tmp=h.readlines()
 h.close()
 sign=0
 for j in tmp[-10:]:
  if sign==1:
   LDOS.append([ round(float(m),5) for m in j.split()])
   break
  if ' ******** EF and DOS at fermi level *******' in j:
   sign=1
for i in range(len(LDOS)):
 print('EF, DOS and LDOS of atom '+str(i)+': '),
 for j in LDOS[i]:
  print(str(j)+' '),
sum_ldos=sum([i[1] for i in LDOS])
print('sum of ldos=',sum_ldos)

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


print('Calculate eta...')
atomic_eta=[]
for i in range(nat):
 eta_i=0
 for l in range(len(LDOS[i][2:])-1):
   print (V[i][-1]-2*Z_of_atoms[i]/RMESH[i][-1])*RADWF[i][l][-1]*RADWF[i][l+1][-1]
   eta_l0=(2.*l+2.)/(2.*l+1)/(2.*l+3)*LDOS[i][l+2]*LDOS[i][l+3]/TOTDOS[1] /2/2*2 #doses per spin
# for l in range(len(DOS_l)-1):
#   eta_l0=(2.*l+2.)/(2.*l+1)/(2.*l+3)*DOS_l[l]*DOS_l[l+1]/sum(DOS_l)
   ###radial integral in def of eta
   integral=simps(np.multiply(np.multiply(DVDR[i],RADWF[i][l]),RADWF[i][l+1]),  RMESH[i])*rmesh_multiply*rmesh_multiply
   #2Z/R * R^2= 2Z*R
   integral+=-(V[i][-1]-2*Z_of_atoms[i]/RMESH[i][-1])*RADWF[i][l][-1]*RADWF[i][l+1][-1]
#   integral+=CIN*CIN*np.trapz(np.multiply(np.multiply(DVDR[i],RADWF_22[i][l]),RADWF_22[i][l+1]),x=RMESH[i])
   '''
   for ri in range(1,len(RMESH[i])):
    dr=RMESH[i][ri]-RMESH[i][ri-1]
    dVdr=(V[i][ri]-V[i][ri-1])-V[i][ri]*dr
    inti=      (RADWF[i][l][ri  ][0]*dVdr*RADWF[i][l+1][ri  ][0]+\
           CIN*(RADWF[i][l][ri  ][1]*dVdr*RADWF[i][l+1][ri  ][1]))  
    integral=integral+(inti) #*dr
   '''
   eta_l=eta_l0*integral*integral
   print (integral,integral*integral,eta_l0, eta_l)
   eta_i=eta_i+eta_l           
 atomic_eta.append(eta_i)


print('atomic eta:'),
print(atomic_eta)
eta=sum(atomic_eta) #j/m^2
eta=eta*Ry_to_eV/(ab_to_ang**2) #/((4./3*3.14*rmt**3)**2)/4*pi*4*pi #ev/A^2
print ('eta= '+str(eta) +' eV/A^2')
 
   
'''
#ten sam wynik
   for ri in range(1,len(RMESH[i])):
    dr=RMESH[i][ri]-RMESH[i][ri-1]
    dVdr=(V[i][ri]-V[i][ri-1])-V[i][ri]*dr
    for rip in range(1,len(RMESH[i])):  
     drp=RMESH[i][rip]-RMESH[i][rip-1]
     dVdrp=(V[i][rip]-V[i][rip-1])-V[i][rip]*drp
     inti=      (RADWF[i][l][ri  ][0]*dVdr*RADWF[i][l+1][ri  ][0]+\
           CIN*(RADWF[i][l][ri  ][1]*dVdr*RADWF[i][l+1][ri  ][1])) *\
                (RADWF[i][l][rip ][0]*dVdrp*RADWF[i][l+1][rip ][0]+\
           CIN*(RADWF[i][l][rip ][1]*dVdrp*RADWF[i][l+1][rip ][1])) 
     integral=integral+(inti) #*dr
   eta_l=eta_l0*integral #*integral
   print integral,eta_l0, eta_l
   eta_i=eta_i+eta_l           
 atomic_eta.append(eta_i)
'''
'''
#calka przeksztalcona by nie robic pochodnej V
print('Calculate eta...')
atomic_eta=[]
for i in range(na):
 eta_i=0
 for l in range(len(LDOS[i][2:])-1):
   eta_l0=(2.*l+2.)/(2.*l+1)/(2.*l+3)*LDOS[i][l+2]*LDOS[i][l+3]/sum_ldos /2/2*2 #doses per spin
   ###radial integral in def of eta
 #  integral=(RADWF[i][l][-1 ][0]*RADWF[i][l+1][-1 ][0]*V[i][-1][0]) #[0,0]
   #int x * y' = x*y - int x'*y
   #x=r^2R_l*R_l1, x'=r^2*(R_l'R_l1 + R_l*R'_l1)+2r*R_l*R_l1
   for ri in range(1,len(RMESH[i])):
    dr=RMESH[i][ri]-RMESH[i][ri-1] 
    dRldr =RADWF[i][l  ][ri  ][0]-RADWF[i][l  ][ri-1][0]
    dRl1dr=RADWF[i][l+1][ri  ][0]-RADWF[i][l+1][ri-1][0]
    if ri!=1: 
     dRldr1 =RADWF[i][l  ][ri-1 ][0]-RADWF[i][l  ][ri-2][0]
     dRl1dr1=RADWF[i][l+1][ri-1 ][0]-RADWF[i][l+1][ri-2][0]
    else:
     dRldr1=dRldr
     dRl1dr1=dRl1dr
    integral=-(dRldr *RADWF[i][l+1][ri][0]+\
               dRl1dr*RADWF[i][l  ][ri][0]+\
               2*RADWF[i][l+1][ri][0]*RADWF[i][l][ri][0]/RMESH[i][ri]*dr)\
              *V[i][ri][0]/2.
    integral=-(dRldr1 *RADWF[i][l+1][ri-1][0]+\
               dRl1dr*RADWF[i][l  ][ri-1][0]+\
               2*RADWF[i][l+1][ri-1][0]*RADWF[i][l][ri][0]/RMESH[i][ri-1]*dr)\
              *V[i][ri-1][0]/2.
   print integral
   eta_l=eta_l0*integral*integral
   print integral,eta_l0, eta_l
   eta_i=eta_i+eta_l           
 atomic_eta.append(eta_i)
'''



'''
print('Read k-mesh...')

h=open('DOS/DOS.outputkgen')
tmp=h.readlines()
h.close()
for numi,i in enumerate(tmp):
 if 'DIVISION' in i: no_of_kpoints=int(i.split()[-1])
 if 'relation' in i:
  tmp=tmp[numi+1:numi+1+(no_of_kpoints+1)**3]
  break
EQUIV_orig=[ [int(m) for m in i.split()[:5]] for i in tmp if int(i.split()[1])!=no_of_kpoints and int(i.split()[2])!=no_of_kpoints and int(i.split()[3])!=no_of_kpoints]
NONEQ_map=[i[0] for i in (EQUIV_orig) if i[0]==i[4]]


#NONEQ_map=[[i,numi] for numi,i in enumerate(NONEQ_map)]
EQUIV=[]
for i in EQUIV_orig:
 for numj,j in enumerate(NONEQ_map):
  if i[4]==j:
   EQUIV.append(numj)
   break
#EQUIV=[ int(i.split()[4])-1 for i in tmp if int(i.split()[1])!=no_of_kpoints and int(i.split()[2])!=no_of_kpoints and int(i.split()[3])!=no_of_kpoints]
print("No_of noneq="+str(len(NONEQ_map))+". No of all kpoints="+str(len(EQUIV))+' should be = '+str(no_of_kpoints**3))
  

print('Reading band energies...'),
ENE=[]
tmp=[]
print('Files:'),
weights=[]
for i in range(1,64):
 try: 
  h=open('DOS/DOS.energy_'+str(i),'r')
  print ' DOS/DOS.energy_'+str(i),
  tmp.extend([m.split() for m in h.readlines()[2:]])
  h.close()
 except:
  break

#ENE[i][j] i-kpoint, j- band
for i in tmp:
  if len(i)>4 and int(i[-4])==len(ENE)+1: 
   weights.append(i[-1])
   ENE.append([])
  elif len(i)==2: 
   ENE[-1].append(i[1])
#rearrange
#et[ibnd][tetra[i][nt]]
print("; No of noneq kpoints="+str(len(ENE))), 
n_band=min([ len(i) for i in ENE])
ENE=[ [float(ENE[j][i]) for j in range(len(ENE))] for i in range(n_band)]
print("; No of bands="+str(len(ENE)))


print('Reading alm coefficients and kpoints...'),
tmp=[]
for i in range(1,64):
 try: 
  h=open('DOS/DOS.almblm_'+str(i),'r')
  tmp.extend([m.split()  for m in h.readlines() if len(m.split())!=0])
  h.close()
 except:
  break
#ALM[k][i][j][l] k-kpoint, i-atoms, j-band, l (orbital No) (summed over m)
ALM=[] #[ [[] for l in range(len(RADWF[i]))] for i in range(na)]
NONEQ=[] #list of nonequiv kpoints
for i in tmp:
 if 'K-POINT' in i[0]: 
  ALM.append([])
  NONEQ.append(np.array([round(float(m),3) for m in i[1:4]]))
 elif 'ATOM' in i[1]: ALM[-1].append([])
 elif 'weight' in i[-1]:  ALM[-1][-1].append( [ 0 for k in range(lmin) ] )
 elif len(i)>10:
  ALM[-1][-1][-1][int(i[0])] += (float(i[3])**2+float(i[4])**2) #alm* * alm = |alm|^2=re^2+im^2

h=open('kpoints.dat','w')
for i in NONEQ:
 for j in i:
  h.write(str(j)+' ')
 h.write('\n')
h.close()

#rearrange
#from: #ALM[k][i][j][l] k-kpoint, i-atoms, j-band, l (orbital No) (summed over m)
#to: weights[na][ibnd][kpoint]
[n_k,n_at,n_l]=[len(ALM),1,lmin]
ALM=[ [[[ALM[l][i][k][j] for l in range(n_k)] for k in range(n_band)] for j in range(n_l)] for i in range(n_at)]
print(" No of noneq kpoints="+str(n_k))
 
tetra=hopfield_tetra.tetrahedra(no_of_kpoints,no_of_kpoints,no_of_kpoints,EQUIV) #the matrix with 6*len(noneq) tetrahedrons. tetra[i][j] i-no of wierzcholek (i=1:4), j-no of tetrahedron  
print('No of tetrahedrons='+str(len(tetra[0])))

[etetra,wtetra]=\
  hopfield_tetra.e_tetra(ENE, len(ENE),len(tetra[0]),tetra,ALM[0],len(ALM[0]))
print np.transpose(tetra)[100:110]

#etetra - energies  arranged into tetrahedron 
#(1 energy for 1 vertex of tetrahedron = 4energies per tetrahedron)
#wtetra - weights arranged into tetrahedron and avaraged over tetrahedron
# so we have 1 weight for 1 tetrahedron
#def dos_t(etetra,wtetra,nbnd,ntetra,e,Noofatoms):

[DOS,DOS_l] = hopfield_tetra.dos_t(etetra,wtetra,len(ENE),len(tetra[0]),EF,len(ALM[0]))
print 'DOS',(DOS),DOS_l,sum(DOS_l)
'''






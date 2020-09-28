import numpy as np
import os
import glob
import hopfield_tetra
#############FOR Nb eta=7.627eV/ang^2 [papaconsta.. 2002] or 4.7 eV/ang^2 (allen-dynes)

################RADWF AND V are already multiplied by r


prefix=raw_input('Prefix: ')
no_of_kpoints=12


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


print('Reading Volume...'),
os.system('grep ":VOL " '+prefix+'.scf0 |head -n 1 >>hopfield_tmp')
os.system('grep ":VOL " '+prefix+'.scf |head -n 1 >>hopfield_tmp')
h=open('hopfield_tmp')
tmp=h.readline()
h.close()
os.system('rm hopfield_tmp')
volume=float(tmp.split()[-1])
print(volume)

print('Reading potential...')
h=open(prefix+'.vsp','r') #stored as V*r
tmp=h.readlines()
h.close()

V=[]
for i in tmp:
  if 'ATOMNUMBER' in i:
   V.append([])
  elif 'VLM(R)' in i:
   V[-1].append([])
  else:
   try: V[-1][-1].append(float(i[0:22]))
   except: continue
   try: V[-1][-1].append(float(i[22:41]))
   except: continue
   try: V[-1][-1].append(float(i[41:60]))
   except: continue
   try: V[-1][-1].append(float(i[60:]))
   except: continue


V=[ np.transpose(np.array(i)) for i in V] # V[i][j][k] i-th atm, j-th r-point,
print('number of atoms:'+str(len(V)))
print('number of r points:'+str(len(V[0])))
na=len(V) #number of atoms
nr=len(V[0]) #number of r-points
 


print('Reading radial wave funtions...')
h=open('RADWF/RADWF.radwf','r')
tmp=h.readlines()
h.close()

RADWF=[]   #RADWF[h][i][j][k][0/1] h-kpoint, i-atoms, j- l (orbital No), k - r-mesh, [0/1]-large or small component
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
 RMESH.append([])
 for j in range(nr):
  RMESH[-1].append(r0*np.exp(j*dx))
 print('RMT='+str(rmt)+' should ='+str(RMESH[-1][-1]))

#renormalizacja f-cji radialnych do 1
#niepotrzebne - fcje sa unormowane, a zreszta to sie skraca w ecie
'''
for i in range(len(RADWF)):
 for j in range(len(RADWF[i])):
  norma=np.sqrt(sum([(RADWF[i][j][k][0]*RADWF[i][j][k][0]+CIN*RADWF[i][j][k][1]*RADWF[i][j][k][1])*(RMESH[i][k]-RMESH[i][k-1]) for k in range(1,len(RADWF[i][j]))]))
  RADWF[i][j]=[k/norma for k in RADWF[i][j]]
'''

#oryginalnie w pliku *vsp jest zapisany V*R, wiec dzielimy przez R by dostac czysty V
V=[ [V[i][j]/RMESH[i][j]*(4*pi)**0.5 for j in range(len(V[i]))] for i in range(len(V))]   



print('Write to files V_i.dat and RADWF_i.dat...')
for i in range(len(V)):
 h=open('V_'+str(i)+'.dat','w')
 for j in range(len(V[i])):
  h.write(str(RMESH[i][j])+' ')
  for k in range(len(V[i][j])):
   h.write(str(V[i][j][k])+' ')
  h.write('\n')
 h.close()

RADWF_1= [ [ [m[0] for m in n] for n in p] for p in RADWF] #large component
RADWF_2= [ [ [m[1] for m in n] for n in p] for p in RADWF] #small component
RADWF_1=[ np.transpose(np.array(m)) for m in RADWF_1]  #RADWF_1[i][j][k] i-atom,j-rmesh, k=l (orbital No)
RADWF_2=[ np.transpose(np.array(m)) for m in RADWF_2]
for i in range(len(RADWF_1)):
 h=open('RADWF_'+str(i)+'.dat','w')
 h.write('# r , l=0- big component,small component, l=1- big component,small component....\n')
 for j in range(len(RADWF_1[i])):
  h.write(str(RMESH[i][j])+' ')
  for k in range(len(RADWF_1[i][j])):
   h.write(str(RADWF_1[i][j][k])+' '+str(RADWF_2[i][j][k])+' ')
  h.write('\n')
 h.close() 




print('Reading alm coefficients and kpoints...')
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
#weights[na][ibnd][tetra[i][nt]]
[n_k,n_at,n_band,n_l]=[len(ALM),len(ALM[0]),min([len(i[0]) for i in ALM]),lmin]
ALM=[ [[[ALM[l][i][k][j] for l in range(n_k)] for k in range(n_band)] for j in range(n_l)] for i in range(n_at)]
print(" No of noneq kpoints="+str(n_k))

print('Read k-mesh...')
h=open('DOS/DOS.outputkgen')
tmp=h.readlines()
h.close()
for numi,i in enumerate(tmp):
 if 'relation' in i:
  tmp=tmp[numi+1:numi+1+(no_of_kpoints+1)**3]
  break
EQUIV=[ int(i.split()[4])-1 for i in tmp]

  

print('Reading band energies...'),
ENE=[]
tmp=[]
weights=[]
for i in range(1,64):
 try: 
  h=open('DOS/DOS.energy_'+str(i),'r')
  print 'Reading file DOS/DOS.energy_'+str(i)
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
print("No of noneq kpoints="+str(len(ENE))+'='+str(n_k))
ENE=[ [ENE[j][i] for j in range(len(ENE))] for i in range(n_band)]



tetra=hopfield_tetra.tetrahedra(no_of_kpoints,no_of_kpoints,no_of_kpoints,EQUIV) #the matrix with 6*len(noneq) tetrahedrons. tetra[i][j] i-no of wierzcholek (i=1:4), j-no of tetrahedron  
[etetra,wtetra]=hopfield_tetra.e_tetra(ENE, len(ENE),len(tetra),tetra,ALM[0],len(ALM[0]))
#etetra - energies  arranged into tetrahedron 
#(1 energy for 1 vertex of tetrahedron = 4energies per tetrahedron)
#wtetra - weights arranged into tetrahedron and avaraged over tetrahedron
# so we have 1 weight for 1 tetrahedron
EF=0.4809652770
[DOSofE,Atomweight] = hopfield_tetra.dos_t(etetra,wtetra,len(ENE),len(tetra),EF,len(ALM[0]))
print 'DOS',(DOSofE)




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

for i in range(na):
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
for i in range(na):
 eta_i=0
 for l in range(len(LDOS[i][2:])-1):
   eta_l=(2.*l+2.)/(2.*l+1)/(2.*l+3)*LDOS[i][l+2]*LDOS[i][l+3]/TOTDOS[1] /2/2*2 #doses per spin
   ###radial integral in def of eta
   INTi=[]
   integral=0 #[0,0]
   for ri in range(1,len(RMESH[i])):
#    dr=RMESH[i][ri]-RMESH[i][ri-ri]
    dVdr=(V[i][ri][0]-V[i][ri-1][0]) #/dr #do not devide here by dr, because then the integral is multiplied by dr, so dr is cancelled out           
    inti=      (RADWF[i][l][ri  ][0]*dVdr*RADWF[i][l+1][ri  ][0]+CIN*(RADWF[i][l][ri  ][1]*dVdr*RADWF[i][l+1][ri  ][1])) #*RMESH[i][ri]*RMESH[i][ri]
    inti=inti+ (RADWF[i][l][ri-1][0]*dVdr*RADWF[i][l+1][ri-1][0]+CIN*(RADWF[i][l][ri-1][1]*dVdr*RADWF[i][l+1][ri-1][1])) #*RMESH[i][ri-1]*RMESH[i][ri-1]
    integral=integral+(inti/2.) #*dr
   print integral,eta_l
   eta_l=eta_l*integral*integral
   eta_i=eta_i+eta_l           
 atomic_eta.append(eta_i)

print('atomic eta:'),
print(atomic_eta)
eta=sum(atomic_eta) #j/m^2
eta=eta*Ry_to_eV/(ab_to_ang**2) #ev/A^2
print ('eta=',eta,'eV/A^2')
 
   









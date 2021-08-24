import numpy as np
import os
import scipy.special as sp
#############FOR Nb eta=7.627eV/ang^2 [papaconsta.. 2002] or 4.7 eV/ang^2 (allen-dynes)

################RADWF AND V are already multiplied by r


prefix=raw_input('Prefix: ')

Ry_to_eV=13.605693122994
ab_to_ang=0.529177210903

#j_to_ev=1/(1.602176634)*1e19
#hbar=6.582119569e-16 #ev/s
#hbar=hbar/13.6056980659 #ry/s
#pi=3.141592653589793238462643
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


V=[ np.transpose(np.array(i)) for i in V]
print('number of atoms:'+str(len(V)))
print('number of r points:'+str(len(V[0])))
na=len(V) #number of atoms
nr=len(V[0]) #number of r-points



print('Reading radial wave funtions...')
h=open('RADWF/RADWF.radwf','r')
tmp=h.readlines()
h.close()

RADWF=[]   #RADWF[i][j][k][0/1] i-atoms, j- l (orbital No), k - r-mesh, [0/1]-large or small component
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
RMESH=[]
#(r_i=r0*exp((i-1)*dx)
#   1 781   0.0000100000   0.0159348926   2.5000000000
# write(23,4645)jatom,jri(jatom),r0(jatom),dx(jatom),rmt(jatom)
for i in mesh_info:
 [r0,dx,rmt]=i[2:]
 RMESH.append([])
 for i in range(nr):
  RMESH[-1].append(r0*np.exp(i*dx))
 print('RMT='+str(rmt)+' should ='+str(RMESH[-1][-1]))




#oryginalnie w pliku *vsp jest zapisany V*R, wiec dzielimy przez R by dostac czysty V
V=[ [V[i][j]/RMESH[i][j] for j in range(len(V[i]))] for i in range(len(V))]   



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
 print(' ')
sum_ldos=sum([i[1] for i in LDOS])

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
TOTDOS[1]=TOTDOS[1]/2 #we need totdos per spin

print('Calculate phase shifts...'),
krs=TOTDOS[0]**0.5*RMESH[i][-1] #sqrt(ef)*RMT
DELTA=[[] for i in range(na)]
for i in range(na):
 for l in range(len(LDOS[i][2:])):
  L=(RADWF[i][l][-1][0]/RMESH[i][-1]-RADWF[i][l][-2][0]/RMESH[i][-2])/(RMESH[i][-1]-RMESH[i][-2])/(RADWF[i][l][-1][0]/RMESH[i][-1])
  tgd= (sp.spherical_jn(l,krs,derivative=True)-sp.spherical_jn(l,krs)*L)/\
       (sp.spherical_yn(l,krs,derivative=True)-sp.spherical_yn(l,krs)*L)
  DELTA[-1].append(np.arctan(tgd))
print DELTA


#renormalizacja f-cji radialnych do 1
#niepotrzebne - fcje sa unormowane, a zreszta to sie skraca w ecie
#renormalizacja zgodnie z pickettem

print('Renormalize wave function...'),
for i in range(na):
 for j in range(len(LDOS[i][2:])):
  coeff=(sp.spherical_jn(j,krs)*np.cos(DELTA[i][j])-sp.spherical_yn(j,krs)*np.sin(DELTA[i][j]))/RADWF[i][j][-1][0]
  print(coeff),
  RADWF[i][j]=[[k[0]*coeff,k[1]] for k in RADWF[i][j]]


print('\nCalculate free scattering dos...'),
N1=[[0 for j in range(len(LDOS[i][2:]))] for i in range(na)]
for i in range(na):
 for l in range(len(LDOS[i][2:])):
   n1l=0 #integral f(x)dx= (f(x)+f(x-1))*dx/2
   for ri in range(1,len(RMESH[i])):
    dr=RMESH[i][ri]-RMESH[i][ri-1]
    n1l= n1l+ RADWF[i][l][ri][0]*RADWF[i][l][ri][0]*dr #RMESH[i][ri]*RMESH[i][ri] #small comp.
    n1l= n1l+ CIN*RADWF[i][l][ri][1]*RADWF[i][l][ri][1]*dr #*RMESH[i][ri]*RMESH[i][ri] #big comp.
    n1l= n1l+ RADWF[i][l][ri-1][0]*RADWF[i][l][ri-1][0]*dr #*RMESH[i][ri-1]*RMESH[i][ri-1] #small comp.
    n1l= n1l+ CIN*RADWF[i][l][ri-1][1]*RADWF[i][l][ri-1][1]*dr #*RMESH[i][ri-1]*RMESH[i][ri-1] #big comp.
   N1[i][l]=(n1l/2.)*(2*l+1) #/volume
print(N1)

print('Calculate eta...')
atomic_eta=[]
for i in range(na):
 eta_i=0
 for l in range(len(LDOS[i][2:])-1):
   eta_l=(2.*l+2.)*LDOS[i][l+2]*LDOS[i][l+3]/sum_ldos/N1[i][l]/N1[i][l+1]*np.sin(DELTA[i][l]-DELTA[i][l+1])*np.sin(DELTA[i][l]-DELTA[i][l+1])
   eta_i=eta_i+eta_l            #papalous- nie mnozyc #*2*LDOS[0][0] *mass[i]/hbar/hbar/pi/pi #LDOS[0][0]=Ef
 atomic_eta.append(eta_i)

print('atomic eta:')
print(atomic_eta)
eta=sum(atomic_eta) #j/m^2
eta=eta*Ry_to_eV/(ab_to_ang**2) #ev/A^2
print ('eta=',eta,'eV/A^2')
 
   









from inputs import *

class real_structure(inputs):
 def __init__(self,EF):
  self.so=inputs.so
  if self.so: self.CIN=(1/137.) #**2
  else: self.CIN= 1e-22
  self.Vtot=[]
  self.Vtot_interstitial=[]
  self.kpoints_interstitial=[]
  self.LM=[]
  self.RADWF=[]
  self.RADWFsmall=[]
  self.RMESH=[]
  self.DR=[]
  self.DVDR=[]
  self.n_r=0
  self.n_at=0
  self.volume=0
  self.Z_of_atoms=[]
  self.EF=EF
 def read_potential(self):
  self.Vtot,self.LM,self.n_r,self.n_at,self.Vtot_interstitial,self.kpoints_interstitial=read_potential(inputs.pot_file)
 def read_radwf(self):
  self.RADWF,self.RADWFsmall,self.RMESH,self.DR=read_radwf(inputs.radwf_file,self.n_r,self.n_at,self.CIN)
 def div_pot_and_radwf_by_r(self):
  self.Vtot,self.DVDR,self.RADWF,self.RADWFsmall,self.RMESH=div_pot_and_radwf_by_r(self.Vtot,self.RADWF,self.RADWFsmall,self.RMESH, self.DR,self.n_at, self.LM,self.Z_of_atoms,self.EF)
  self.n_r=[len(self.Vtot[at][0]) for at in range(self.n_at)]
#  self.LM=[[self.LM[at][0]] for at in range(self.n_at)]
#  self.Vtot=[[self.Vtot[at]] for at in range(self.n_at)]
#  self.DVDR=[[self.DVDR[at]] for at in range(self.n_at)]
 def read_volume(self):
  self.volume,self.Z_of_atoms=read_volume(inputs.scf_file)

def read_volume(scf_file):
 print('Read volume...')
 h=open(scf_file)
 tmp=h.readlines()
 h.close()
 for i in tmp: 
  if ':VOL' in i: 
   volume=float(i.split()[-1])
   break
 print('Reading Z of atoms...'),
 os.system('grep "Z: " hopfield_calc/hopfield_calc.struct  >hopfield_tmp')
 h=open('hopfield_tmp')
 tmp=h.readlines()
 h.close()
 os.system('rm hopfield_tmp')
 Z_of_atoms=[float(i.split()[-1]) for i in tmp]
 print('Z=',Z_of_atoms)
 return volume,Z_of_atoms

def read_potential(pot_file):
 print('Read intrasitial total potential...')
 h=open(pot_file)
 tmp=h.readlines()
 h.close()
 Vtot=[] #[atom][lm][r]
 LM=[] #[atom][list of lm numbers]
 Vtot_interstitial=[]
 kpoints_interstitial=[]
 for ni,i in enumerate(tmp):
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
 print('Read intersitial total potential...')
 no_of_vk=int(tmp[ni+2].split()[0])
 print(no_of_vk)
 for i in tmp[ni+3:ni+3+no_of_vk]:
  kpoints_interstitial.append([int(k) for k in i[:18].split()])
  Vtot_interstitial.append(float(i[18:37])+1j*float(i[37:]))

 print (n_at,n_r, LM)
 return Vtot,LM,n_r,n_at,Vtot_interstitial,kpoints_interstitial

def read_radwf(radwf_file,n_r,n_at,CIN):
 print('Reading radial wave funtions...')
#     IF(MODUS.EQ.'ALM  ') then
#        write(23,4645)jatom,jri(jatom),r0(jatom),dx(jatom),rmt(jatom)
#        do l=0,lmax2
#           write(23,*) l
#           if (l.le.lomax) then
#              write(23,4646) (RRAD1(jrj,l),RRAD2(jrj,l),RADE1(jrj,l),RADE2(jrj,l), &
#                   a1lo(jrj,1,l),b1lo(jrj,1,l),a1lo(jrj,2,l),b1lo(jrj,2,l), &
#                   a1lo(jrj,3,l),b1lo(jrj,3,l),jrj=1,jri(jatom))
#           else
#              write(23,4647) (RRAD1(jrj,l),RRAD2(jrj,l),RADE1(jrj,l),RADE2(jrj,l), &
#                    jrj=1,jri(jatom))
#           endif
#        enddo
#     endif

 h=open(radwf_file,'r')
 tmp=h.readlines()
 h.close()

 RADWF=[]   #RADWF[i][j][k][0-1] i-atoms, j- l (orbital No), k - r-mesh, [0-2]-large  component of radwf and udot, and lo if exists
 RADWFsmall=[] #RADWFsmall[i][j][k][0-1] i-atoms, j- l (orbital No), k - r-mesh, [0-1]-  small component of radwf and udot, and lo if exists
 mesh_info=[]


 #read radwf
 #factor i^l/sqrt(4pi) comes from the formula from main.F in SRC_lapw7. by RMT we will multiply later
 #! psi(r) = e^ikR Sum(lm) w_lm,a(|r-R-R_a|) Y*(*)_lm(T_a^-1(r-R-R_a))   with
 #!
 #!  w_lm,a(r) = 4pi*i^l [ A_l,m,a *      u_l,a(r,E _l,a) +
 #!                       B_l,m,a * d/dE u_l,a(r,E _l,a) +
 #!                       C_l,m,a *      u_l,a(r,E'_l,a) ] * Rmt_a^2
 #where R_a=pos of atom, R=pos of considered unit cell (here R=0, because we consider only 1 unit cell 
 #--not needed, lapw2 already does it
 for i in tmp:
  if str(n_r) in i.split(): 
   RADWF.append([])
   RADWFsmall.append([])
   mesh_info.append([float(m) for m in i.split()])
  elif len(i.split())==1: 
   l=int(i.split()[0])
   RADWF[-1].append([])
   RADWFsmall[-1].append([])
  else:
   try: 
    RADWF[-1][-1].append(np.array([float(i.split()[0]),float(i.split()[2]),float(i.split()[4]),float(i.split()[6])   ]))  #with RLO
    RADWFsmall[-1][-1].append(np.array([CIN*float(i.split()[1]),CIN*float(i.split()[3]),CIN*float(i.split()[5]),CIN*float(i.split()[7])    ]))
   except:
    try: 
     RADWF[-1][-1].append(np.array([float(i.split()[0]),float(i.split()[2]),0.,0.]))  #w/o RLO
     RADWFsmall[-1][-1].append(np.array([CIN*float(i.split()[1]),CIN*float(i.split()[3]),0.,0.]))
    except:continue

 for i in range(len(RADWF)):
  print('For atom no. '+str(i)+' I found '+str(len(RADWF[i]))+' wave functions u_l(r) at the r-meshes as follows:')
#  if not np.any([sum(j) for j in
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
 return RADWF, RADWFsmall, RMESH, DR

def div_pot_and_radwf_by_r(Vtot,RADWF,RADWFsmall,RMESH,DR, n_at,LM,Z_of_atoms,EF):
 #oryginalnie w pliku *vtotal jest zapisany V*R^2, wiec dzielimy przez R by dostac czysty V*R
 '''
 #renormalization of potential
 for i in range(len(Vtot)):
  for j in range(len(Vtot[i])):
   suma=sum([Vtot[i][j][k]*DR[i][k-1] for k in range(1,len(RMESH[i]))])
   Vtot[i][j]=np.array(Vtot[i][j])/suma
 '''
# for i in range(n_at):
#  Vtot[i][0]=np.array(Vtot[i][0])/((4*np.pi)**0.5)

 '''
 Vtot=[ [np.array([(Vtot[i][j][k]/(RMESH[i][k])+2*Z_of_atoms[i])/RMESH[i][k] for k in range(len(Vtot[i][j]))]) for j in range(len(Vtot[i]))] for i in range(len(Vtot))]  
 DVDR=[ [  (np.gradient(Vtot[at][nlm1],RMESH[at],edge_order=2)+2*Z_of_atoms[at]/RMESH[at]/RMESH[at]) for nlm1 in range(len(LM[at]))] for at in range(n_at)] #in fact it is r^2 * dV/dr  =(d(Vr^2)/(dr) - 2Vr)
 Vtot=[ [np.array([Vtot[i][j][k]-2*Z_of_atoms[i]/RMESH[i][k] for k in range(len(Vtot[i][j]))]) for j in range(len(Vtot[i]))] for i in range(len(Vtot))]  
 '''

 Vtot=[ [np.array([(Vtot[i][j][k]/(RMESH[i][k]))/RMESH[i][k]//(4*np.pi)**0.5 for k in range(len(Vtot[i][j]))]) for j in range(len(Vtot[i]))] for i in range(len(Vtot))]  
 DVDR=[ [  np.gradient(Vtot[at][nlm1],RMESH[at],edge_order=2) for nlm1 in range(len(LM[at]))] for at in range(n_at)] 

#RADWF[i][j][k][0-1] i-atoms, j- l (orbital No), k - r-mesh, [0-1]-large  component of radwf and udot

#oryginalnie w pliku *radwf jest zapisany radwf*R, wiec dzielimy przez R by dostac czysty radwf
#dodatkowo z RADWF[at][l][r][u,udot] -> RADWF[at][l][u,udot][r]
 RADWF=    np.transpose(RADWF,axes=(0,1,3,2))
 RADWFsmall=    np.transpose(RADWFsmall,axes=(0,1,3,2))


# Vtot=[Vtot[at][0] for at in range(n_at)]
# DVDR=[DVDR[at][0] for at in range(n_at)]

 h=open('hopfield_calc/hopfield_calc.output1_1')
 tmpp=h.readlines()
 h.close()
 for ni,i in enumerate(tmpp):
  if 'POTENTIAL PARAMETERS' in i:
   EN=[float(m.split()[1]) for m in tmpp[ni+2:ni+12]]
   break


 RMESH0=RMESH.copy()
 RADWF0=RADWF.copy()
 RADWF2=[[[[] for k in range(len(RADWF0[at][l]))] for l in range(len(RADWF0[at]))] for at in range(n_at)]
 RADWFsmall2=[[[[] for k in range(len(RADWF0[at][l]))] for l in range(len(RADWF0[at]))] for at in range(n_at)]
 RMESH2=[[] for at in range(n_at)]
 Vtot2=[[[] for nlm1 in range(len(LM[at]))] for at in range(n_at)]
 for at in range(n_at):
  for nlm1 in range(len(LM[at])):
   RMESH2[at],Vtot2[at][nlm1]=interpolate_V(RMESH0[at],Vtot[at][nlm1])
#  RADWF2[at][0][0]=lipman(0,RMESH2[at],/RADWF0[at][0][0],EF-0.2,Vtot2[at],Z_of_atoms[at])
  for l in range(len(RADWF[at])):
   for k in range(len(RADWF[at][l])):
    RMESH2[at],RADWFsmall2[at][l][k]=interpolate_V(RMESH0[at],RADWFsmall[at][l][k])
#    if l==0 and k==0: RADWF2[at][l][0]=lipman(0,RMESH2[at],RADWF0[at][l][0],EN[l],Vtot2[at][0],Z_of_atoms[at])
    if k==0: RADWF2[at][l][0]=lipman(l,RMESH2[at],RADWF0[at][l][0],EF,Vtot2[at][0],Z_of_atoms[at])
#    if l==0 and k==0: continue
    else:
#     RADWF2[at][l][k]=np.zeros(len(RMESH2[at]))
     RMESH2[at],RADWF2[at][l][k]=interpolate_V(RMESH0[at],RADWF0[at][l][k])
#    RADWF[at][l][k]=RADWF[at][l][k]/(simps(RADWF[at][l][k]*RADWF[at][l][k],RMESH[at])**0.5)
 DVDR2=[ [  np.gradient(Vtot2[at][nlm1],RMESH2[at],edge_order=2) for nlm1 in range(len(LM[at]))] for at in range(n_at)] 



 '''
 print('Write spherical V and radwf to files V_i.dat and RADWF_i.dat...')
 for i in range(len(Vtot)):
  h=open('V_'+str(i)+'.dat','w')
  for j in range(len(Vtot[i][0])):
   h.write(str(RMESH[i][j])+' ')
   for k in range(len(Vtot[i])):
#    h.write(str((Vtot[i][k][j]/((4*np.pi)**0.5)+2*82)/RMESH[i][j])+' ')
    h.write(str(Vtot[i][k][j])+' ')
   h.write('\n')
  h.close()

 for i in range(len(Vtot)):
  h=open('DVDR_'+str(i)+'.dat','w')
  for j in range(len(Vtot[i][0])):
   h.write(str(RMESH[i][j])+' ')
   for k in range(len(Vtot[i])):
    h.write(str(DVDR[i][k][j])+' ')
   h.write('\n')
  h.close()
 '''

 print('Write spherical V and radwf to files V_i.dat and RADWF_i.dat...')
 for i in range(len(Vtot2[0])):
  h=open('V_'+str(i)+'.dat','w')
  for j in range(len(Vtot2[0][i])):
   h.write(str(RMESH2[0][j])+' '+str(Vtot2[0][i][j])+'\n')
  h.close()


#RADWF[at][l][u,udot][r]
 for i in range(len(RADWF2)):
  print('aa')
  h=open('RADWF_'+str(i)+'.dat','w')
  h.write('# r , l=0- big component,small component, l=1- big component,small component....\n')
  for k in range(len(RADWF2[i][1][1])):
   h.write(str(RMESH2[i][k])+' ')
   for j in range(len(RADWF2[i])):
    h.write(str(RADWF2[i][j][0][k])+' '+str(RADWF2[i][j][1][k])+' '+str(RADWF2[i][j][2][k])+' ')
   h.write('\n')
  h.close() 

 return Vtot2,DVDR2,RADWF2,RADWFsmall2,RMESH2

# return Vtot,DVDR,RADWF,RADWFsmall,RMESH

def interpolate_V(RMESH0,V):
 nib=80
 nbl=20
 r0= (RMESH0[0])
 x2= (RMESH0[-1]-r0)/(nib*(2.**nbl-1.))+1e-3
 while x2<r0:
  nbl=nbl-1
  x2=(RMESH0[-1]-r0)/(nib*(2.**nbl)-1.) +1e-3
 dx=x2
 RMESH=[]
 xx=0
 for i in range(nbl):
  for j in range(nib):
   xx+=dx
   if  xx<RMESH0[-1]:
    RMESH.append(xx)
  dx=dx+dx
 Vf=interp1d(RMESH0,V)
 V2=Vf(RMESH)
 h=open('lipmanV.dat','w')
 for i in range(len(RMESH)):
  h.write(str(RMESH[i])+' '+str(V2[i])+'\n')
 h.close()
 return np.array(RMESH),V2


def lipman(l,RMESH,Rstart,Ef,V2,Z_of_atom): #schrea.f from kkr
 #Ef=0.3
 print('EF,Z,nr=',Ef,Z_of_atom,len(RMESH))
 nib=80
 V3=V2 +2*Z_of_atom/RMESH
 VIN=V2 +l*(l+1)/RMESH/RMESH #-2*Z_of_atom/RMESH#[V[r]+(-2*Z_of_atom+l*(l+1)/RMESH[r])/(RMESH[r]) for r in range(len(RMESH))]
 A=[1,-2*Z_of_atom/(2*l+2),0]
 A[2]=(-2*Z_of_atom*A[1]+(V3[0]-Ef)*A[0])/(4*l+6)
 i1=6
 SSR=np.zeros(len(RMESH))
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

import numpy as np
import os
import glob
#import hopfield_tetra
#from math import factorial as sil
from scipy.special import sph_harm
from scipy.integrate import quad,simps
from scipy.misc import derivative
from multiprocessing import Process,Pool
#from scipy.special import hankel1,jv
from scipy.interpolate import interp1d
Ry_to_eV=13.605693122994
ab_to_ang=0.529177210903

#j_to_ev=1/(1.602176634)*1e19
#hbar=6.582119569e-16 #ev/s
#hbar=hbar/13.6056980659 #ry/s
pi=3.141592653589793238462643
#Na=6.02214076e23
sqrtpm1=(np.pi)**(-0.5)


#if os.stat("DOS/DOS.inso").st_size != 0:  
# CIN=1./(137.0359895*2) #for treating big and small component of wave function
# SOC=1
def if_so():
 so=''
 sofile=glob.glob('run_lapw/*inso')
 for i in sofile:
  if os.stat(i).st_size: 
   so=' -so '
   break
 return so

class inputs:
  so=if_so()
  dos_file='hopfield_calc/hopfield_calc.outputt'
  almblm_file='hopfield_calc/hopfield_calc.almblm'
  ene_file='hopfield_calc/hopfield_calc.energy'
  if so: ene_file+='so'
  pot_file='hopfield_calc/hopfield_calc.vtotal'
  radwf_file='hopfield_calc/hopfield_calc.radwf'
  scf_file='hopfield_calc/hopfield_calc.scf0'
  prefix=os.getcwd().split('/')[-1]
  n_l=5 #n_l nieograniczone daje dokladnie te same wyniki co n_l=5; n_l=4 zmienia wynik o 0.01%
  nk=48
  def run_calculations(self):
   yn='n'
   yn=input('should I run calc? y/[n]: ')
   if yn=='y': 
    run_calc(inputs)







def run_calc(inputs):
 print(inputs.so)
# if os.stat(inputs.prefix+'.inso').st_size!=0: so=' -so '
 os.system('rm -r hopfield_calc_old; mv hopfield_calc hopfield_calc_old; cp -r run_lapw  hopfield_calc')
 os.system('cp .machines hopfield_calc')
 os.chdir('hopfield_calc')
# os.system('pwd')
# exit()
 os.system('rename_files '+inputs.prefix+' hopfield_calc')
 os.system("x kgen <<'EOF'\n0\n+"+str(inputs.nk)+" "+str(inputs.nk)+" "+str(inputs.nk)+" \n0\nEOF")
 os.system("sed -i 's/NR2V/R2V/' hopfield_calc.in0")
# os.system("sed -i 's/1      (GLOBAL/0      (GLOBAL/' hopfield_calc.in1")
# os.system("sed -i 's/CONT 1/CONT 0/' hopfield_calc.in1")
# os.system("sed -i 's/STOP 1/STOP 0/' hopfield_calc.in1")
 '''
 h=open("hopfield_calc.in1","r")
 tmp=h.readlines()
 h.close()
 change=0
 tmp2=[]
 for ni in range(len(tmp)):
  if ni<len(tmp)-1 and tmp[ni].split()[0]==tmp[ni+1].split()[0] and tmp[ni].split()[1]=='0.30':
   change+=1
   continue
  elif ni>0 and tmp[ni].split()[0]==tmp[ni-1].split()[0] and tmp[ni].split()[1]=='0.30':
   change+=1
   continue
  else: tmp2.append(tmp[ni])
 aa=str(int(tmp2[2].split()[1])-change)+' '
 tmp2[2]=tmp2[2][:9]+aa+tmp2[2][12:]
 h=open("hopfield_calc.in1","w")
 for i in tmp2: h.write(i)
 h.close()
 '''
# os.system("run_lapw  -i 100  -p")
# os.system("rm -r run_nrel_lapw; save_lapw -d run_nrel_lapw")
 os.system("run_lapw -ec 0.00000001 -cc 0.0001 -i 100 "+inputs.so+" -p")
 os.system("rm -r run_lapw; save_lapw -d run_lapw")
 os.system("x lapw1 -p")
 if len(inputs.so): os.system("x lapwso -p")
# os.system("sed -i 's/TOT /FERMI /' hopfield_calc.in2")
 os.system("x lapw2 "+inputs.so+" -alm  -p")
 os.system("x lapw2 "+inputs.so+" -qtl  -p")
 os.system("x tetra ; x tetra "+inputs.so+"-p")
 os.chdir('..')

#so=if_so()
#run_calc(inputs)

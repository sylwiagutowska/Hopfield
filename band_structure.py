from inputs import *

class band_structure(inputs):
 def __init__(self):
  self.E_f=0
  self.dos=0
  self.n_k=0
  self.so=inputs.so
  self.n_k_total=0
  self.n_at=0
  self.ALMBLM=[]
  self.ENE=[] 
  self.kweights=[]
  self.ENE_kweights=[] 
  self.degauss=0
  self.LDOS=[]
 def read_ef_and_dos(self):
  self.E_f,self.dos=read_ef_and_dos(inputs.dos_file,self.so)
 def read_almblm(self):
  self.ALMBLM,self.No_of_bands,self.n_k,self.n_at=read_almblm(inputs.almblm_file,inputs.n_l)
 def read_ene(self):
  self.ENE,self.kweights,self.n_k_total=read_ene(inputs.ene_file,self.No_of_bands,self.n_k,self.n_at,self.E_f)
 def which_bands_cross_ef(self):
  self.ENE,self.ALMBLM=which_bands_cross_ef_and_treat_so(self.ENE,self.ALMBLM,self.E_f,self.so)
 def calc_kweights(self):
  self.degauss=adjust_degauss(self.kweights,self.ENE,self.E_f,self.n_k,self.n_k_total,self.dos)
  print('degauss equal ',self.degauss,' gives the best estimation of dos')
  self.ENE_kweights=calc_kweights(self.kweights,self.ENE,self.E_f,self.n_k,self.degauss)
 def read_ldos(self):
  self.LDOS=read_ldos(self.n_at)

def read_ef_and_dos(dos_file,so):
 print('Read EF and total DOS from '+dos_file+'...')
 E_f,dos=0,0
 h=open(dos_file,'r')
 tmp=h.readlines()[-10:]
 h.close()
 for nj,j in enumerate(tmp):
  if 'EF and DOS at fermi level' in j:
    [E_f,dos]=[ round(float(m),5) for m in tmp[nj+1].split()]
    break
 dos=dos/2.
 print (E_f,dos,' 1/Ry per spin')
 return E_f,dos


#ENE[i][j] i-kpoint, j- band
#ALM[i][l][m][k][j]
#ENE_weights[k][iband]
def which_bands_cross_ef_and_treat_so(ENE,ALMBLM,E_f,so):
 chosen_bands=[]
 nbnd=max([len(k) for k in ENE])
 band_ranges=[[1,-1] for i in range(nbnd)]
 for k in range(1,len(ENE)):
  for i in range(len(ENE[k])):
   if ENE[k][i]<band_ranges[i][0]: band_ranges[i][0]=ENE[k][i]
   elif ENE[k][i]>band_ranges[i][1]: band_ranges[i][1]=ENE[k][i]
 print (band_ranges)

 if so: step=2
 else: step=1

 for ni,i in enumerate(band_ranges):
   if i[0]<E_f+0.01 and i[1]>E_f-0.01: chosen_bands.append(ni)
 mmin=min(chosen_bands)
 mmax=max(chosen_bands)+1
 print('Only ',len(chosen_bands)/step,'bands are in [-0.01,0.01] range around E_F')
 #print(np.array(ENE[:][mmin:mmax]).shape)

 ENE= ([ k[mmin:mmax:step] for k in ENE])
# ENE_kweights= ([ k[mmin:mmax] for k in ENE_kweights])
 ALMBLM=([[[[ k[mmin:mmax:step] for k in m] for m in l] for l in i] for i in ALMBLM] )
# ENE_kweights=[np.array(ENE_kweights)[:][mmin:mmax]
# ALMBLM=np.array(ALMBLM)[:][:][:][:][mmin:mmax]

 return ENE,ALMBLM

def read_almblm(almblm_file,n_l):
 #write(24,4893)l,m,index,alm(INDEX,NUM),blm(INDEX,NUM),(clm(INDEX,NUM,jlo),jlo=1,3)
 print('Reading alm coefficients and kpoints from '+almblm_file+'...')
 tmp=[]
 for i in range(1,64):
  try: 
   h=open(almblm_file+'_'+str(i),'r')
   tmp.extend([m.split()  for m in h.readlines() if len(m.split())!=0])
   h.close()
  except:
   break
 if len(tmp)==0:
   h=open(almblm_file,'r')
   tmp.extend([m.split()  for m in h.readlines() if len(m.split())!=0])
   h.close()
 #ALM are in file only for occupied (or partially occupied) bands, so the number of bands is smaller here than in *energy files
 #in fact, ALM are equal to ! X_l,m,a = Sum(K) c_K/sqrt(V) exp(i(K+k)R_a) Y(*)_lm(T_a^-1(K+k)) X_l,a(|K+k|), (c_k = gaunt coeff.), to liczy augpw.frc, ale w almblm mamy to juz gotowe
 ALMBLM=[] 
 for i in tmp:
  if 'K-POINT' in i[0]: 
   ALMBLM.append([])
  elif 'ATOM' in i[1]: ALMBLM[-1].append([])
  elif 'weight' in i[-1]:  
   ALMBLM[-1][-1].append( [ [] for k in range(n_l) ] )
  elif len(i)>10:
   if int(i[0])<n_l: 
    ALMBLM[-1][-1][-1][int(i[0])].append([ complex(float(i[3]),float(i[4])),complex(float(i[5]),float(i[6])),complex(float(i[7]),float(i[8])),complex(float(i[9]),float(i[10])) ]) 


 n_k=len(ALMBLM)
 n_at=len(ALMBLM[0])
 No_of_bands=[ len(i[0]) for i in ALMBLM]
 #transpose to get: from ALMBLM[kp][at] to ALMBLM[at][kp]
 ALMBLM0=[[ALMBLM[i][j] for i in range(len(ALMBLM))] for j in range(n_at)]
 #ALMBLM[i][k][j][l][m][0-3]  i-atoms,k-kpoint, j-band, l (orbital No), m, [0: Re[Alm]+j*Im[Alm], 1:  Re[Blm]+j*Im[Blm]]
 

#alm and blm as numpy arrays
#ALMBLM[i][k][j][l][m][0-3] -> ALM[i][l][m][k][j]
 ALMBLM=[[[[np.array([ALMBLM0[i][k][j][l][m] for j in range(len(ALMBLM0[i][k]))]) for k in range(n_k) ] for m in range(2*l+1) ] for l in range(n_l) ] for i in range(n_at)]
 return ALMBLM,No_of_bands,n_k,n_at

def read_ene(ene_file,No_of_bands,n_k,n_at,E_f):
 print('Reading band energies...'),
 ENE=[]
 tmp=[]
 print('Files:'),
 kweights=[]
 for i in range(1,64):
  try: 
   h=open(ene_file+'_'+str(i),'r')
   print(ene_file+'_'+str(i)),
   tmp.extend([m.split() for m in h.readlines()[2:]])
   h.close()
  except:
   break
 if i==1:
   h=open(ene_file,'r')
   print(ene_file),
   tmp.extend([m.split() for m in h.readlines()[2:]])
   h.close()

#ENE[i][j] i-kpoint, j- band
 for i in tmp:
   if len(i)>4 and int(i[-4])==len(ENE)+1: 
    kweights.append(float(i[-1]))
    ENE.append([])
   elif len(i)==2: 
    ENE[-1].append(i[1])

 n_k2,n_band=len(kweights),min([ len(i) for i in ENE])
 n_k_total=sum(kweights)
 print("; No of noneq kpoints="+str(n_k))
 if n_k2!=n_k: raise ValueError('no of k_points of energies:',n_k,' and of almblm:',str(n_k2),' is  not the same')
 print("; Total number of kpoints="+str(n_k_total))

 print(' From all '+str(max([len(k) for k in ENE]))+' bands')
#choose only energies of  which are considered in *almblm (only (fully or partially) occupied)
 ENE=[ [float(ENE[k][i]) for i in range(No_of_bands[k])] for k in range(len(ENE))]
 print(' only '+str(max([len(k) for k in ENE]))+' are fully or partially occupied')

 return ENE, kweights,n_k_total

def calc_kweights(kweights,ENE,E_f,n_k,degauss):
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
          eef=(iband-E_f)/degauss
          arg=min(200., (eef - 2.**(-0.5) ) **2)
          ENE_weights[k].append((kweights[k]*sqrtpm1*np.exp(-arg)*(2.-(2**0.5*eef)))/degauss)
 return ENE_weights

def adjust_degauss(kweights,ENE,E_f,n_k,n_k_total,dos):
 min_diff_dos=[1000,1000]
 for i in range(1,200):
  degauss=0.0001*i
  weight=calc_kweights(kweights,ENE,E_f,n_k,degauss)
  dos2=sum([ sum(m) for m in weight])/n_k_total
  if abs(dos2-dos)<min_diff_dos[1]: min_diff_dos=[degauss,abs(dos2-dos)]
 return min_diff_dos[0]



def read_ldos(na):
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
	return LDOS

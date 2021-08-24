#class math_stuff:
from numpy import real,imag,isnan,pi,conj,sin,trapz,cos,array,transpose,flipud,gradient
from scipy.integrate  import quad,dblquad
from scipy.misc import derivative 
from scipy.special  import sph_harm,lpmn


def generate_Plm(lmax):
 #Plm[l][m][cos(tet)]
# scipy.special.lpmn(m,l,cos(theta)) return P_lm and dP_lm/d(cos(theta))=-sin(theta)*dP_lm/d/theta
 [Plm,dPlm]=transpose(array([lpmn(lmax,lmax,cos(pi*x/500.)) for x in range(500)]),axes=(1,3,2,0))
 Plm=[[ Plm[l][m] for m in range(l+1)] for l in range(lmax)]
 dPlm=[[ dPlm[l][m] for m in range(l+1)] for l in range(lmax)]
 #Plm[l][-m]= (-1)**m*sil(l-m)/sil(l+m)*Plm
 #reverse to get an order from -l to 0, not from 0 to -l
 for l in range(1,lmax):
  pl_m=[ (-1)**m*sil(l-m)/sil(l+m)*Plm[l][m] for m in range(l,0,-1)]
  Plm[l]=pl_m + Plm[l]
  dpl_m=[ (-1)**m*sil(l-m)/sil(l+m)*dPlm[l][m]  for m in range(l,0,-1)]
  dPlm[l]= dpl_m + dPlm[l]
 return Plm,dPlm


def generate_Ylm_0(Plm):
#Ylm_0[l][m][tet]=Y_lm(tet,fi=0)
# scipy.special.lpmn(m,l,cos(theta)) return P_lm and dP_lm/d(cos(theta))=-sin(theta)*dP_lm/d/theta
 Ylm=array([[ ( (2*l+1)/(4*pi)* sil(l-(m-l))/sil(l+(m-l)) )**0.5*Plm[l][m]  for m in range(len(Plm[l]))] for l in range(len(Plm))]) #litle bit faster than sph_harm and gives the same result (checked)
# Ylm=[[ array([ sph_harm(m,l,0,tet/500.*pi) for tet in range(500)]) for m in range(-l,l+1)] for l  in range(10)]
 return Ylm



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
  fac1=( ( (2*l3+1.)*sil(l3+l1-l2)*sil(l3-l1+l2)*sil(l1+l2-l3) )/ sil(l1+l2+l3+1) )**0.5
  fac2=( sil(l3+m3)*sil(l3-m3)*sil(l1-m1)*sil(l1+m1)*sil(l2-m2)*sil(l2+m2) )**0.5
  fac3=0
  for k in range(0,500):
   if ((l1+l2-l3-k)<0) or ((l1-m1-k)<0) or ((l2+m2-k)<0) or ((l3-l2+m1+k)<0) or ((l3-l1-m2+k)<0): continue
   fac3+= ((-1)**k)/(sil(k)*sil(l1+l2-l3-k)*sil(l1-m1-k)*sil(l2+m2-k)*sil(l3-l2+m1+k)*sil(l3-l1-m2+k))
  return fac1*fac2*fac3


def threejsymbol(l1,m1,l2,m2,l3,m3):
  fac1=(-1)**(-l1+l2+m3)*(sil(l1+l2-l3)*sil(l1-l2+l3)*sil(l2+l3-l1)*sil(l1+m1)*sil(l1-m1)*sil(l2+m2)*sil(l2-m2)*sil(l3+m3)*sil(l3-m3))**0.5
  fac2=(sil(l1+l2+l3+1))**0.5
  fac3=0
  for k in range(0,500):
   if ((l1+l2-l3-k)<0) or ((l1-m1-k)<0) or ((l2+m2-k)<0) or ((l3-l2+m1+k)<0) or ((l3-l1-m2+k)<0): continue
   fac3+= ((-1)**k)/(sil(k)*sil(l1+l2-l3-k)*sil(l1-m1-k)*sil(l2+m2-k)*sil(l3-l2+m1+k)*sil(l3-l1-m2+k))
  return fac1*fac3/fac2

def three_y(l1,m1,l2,m2,l3,m3):
  if m1+m2-m3!=0: return 0
  elif abs(l1-l2)>l3 or l1+l2<l3: return 0
#  elif l2==0:
#   if l1==l3 and m1==m3:
#    threey=( (2*l1+1.)*(2*l2+1.)/(4.*pi*(2*l3+1.)) )**0.5 
#   else: threey=0
#  else:
  threey= ( (2*l1+1.)*(2*l2+1.)*(2*l3+1.)/(4.*pi ))**0.5\
        *threejsymbol(l1,0,l2,0,l3,0)*threejsymbol(l1,m1,l2,m2,l3,m3) #*(-1)**m3 #-1^m comes from conjugation of latter harmonic 
  '''
  elif m3<0 and l1<l2:
   threey= (-1)**(l3-l2-l1)*( (2*l1+1.)*(2*l2+1.)/(4.*pi*(2*l3+1.)) )**0.5\
        *clebsch_gordon(l2,0,l1,0,l3,0)*clebsch_gordon(l2,-m2,l1,-m1,l3,-m3)   
  elif m3<0: 
   threey= (-1)**(l3-l1-l2)*( (2*l1+1.)*(2*l2+1.)/(4.*pi*(2*l3+1.)) )**0.5\
        *clebsch_gordon(l1,0,l2,0,l3,0)*clebsch_gordon(l1,-m1,l2,-m2,l3,-m3)
  elif l1<l2: 
   threey= ( (2*l1+1.)*(2*l2+1.)/(4.*pi*(2*l3+1.)) )**0.5\
        *clebsch_gordon(l2,0,l1,0,l3,0)*clebsch_gordon(l2,m2,l1,m1,l3,m3)
  '''
  return round(threey,9)

def ynre(m,l,tet,fi):
  return real(sph_harm(m1,l1,fi,tet))
def ynim(m,l,tet,fi):
  return imag(sph_harm(m1,l1,fi,tet))


def spec_int_1(l1,m1,l2,m2,l3,m3,Ylm0):
 # if m1==0: return 0.
  tet=[t*pi/500 for t in range(500)]
  out1=trapz(cos(tet)*sin(tet)*real(Ylm0[l1][m1+l1]*conj(Ylm0[l2][m2+l2])*Ylm0[l3][m3+l3]),x=tet)
  out2=trapz(cos(tet)*sin(tet)*imag(Ylm0[l1][m1+l1]*conj(Ylm0[l2][m2+l2])*Ylm0[l3][m3+l3]),x=tet)
#  out1=quad(lambda tet: real(Ylm0[l1][m1+l1][int(tet)]*conj(Ylm0[l2][m2+l2][int(tet)])*Ylm0[l3][m3+l3][int(tet)]), 0.0, 500)[0]*pi/500
#  out2=quad(lambda tet: imag(Ylm0[l1][m1+l1][int(tet)]*conj(Ylm0[l2][m2+l2][int(tet)])*Ylm0[l3][m3+l3][int(tet)]), 0.0, 500)[0]*pi/500

  if isnan(out1): # or isnan(out2): 
   print('spec_int',l1,m1,l2,m2,l3,m3,':infty')
   return 0.j
  return complex(round(out1,6),round(out2,6))

def spec_int_2(l1,m1,l2,m2,l3,m3,Ylm0):
 # if m1==0: return 0.
  tet=[t*pi/500 for t in range(500)]
  dY=gradient(Ylm0[l1][m1+l1],tet)
  out1=trapz(sin(tet)**2.*real(dY*conj(Ylm0[l2][m2+l2])*Ylm0[l3][m3+l3]),x=tet)
  out2=trapz(sin(tet)**2.*imag(dY*conj(Ylm0[l2][m2+l2])*Ylm0[l3][m3+l3]),x=tet)
#  out1=quad(lambda tet: real(Ylm0[l1][m1+l1][int(tet)]*conj(Ylm0[l2][m2+l2][int(tet)])*Ylm0[l3][m3+l3][int(tet)]), 0.0, 500)[0]*pi/500
#  out2=quad(lambda tet: imag(Ylm0[l1][m1+l1][int(tet)]*conj(Ylm0[l2][m2+l2][int(tet)])*Ylm0[l3][m3+l3][int(tet)]), 0.0, 500)[0]*pi/500

  if isnan(out1): # or isnan(out2): 
   print('spec_int',l1,m1,l2,m2,l3,m3,':infty')
   return 0.j
  return complex(round(out1,6),round(out2,6))

def spec_int(l1,m1,l2,m2,l3,m3,Ylm0):
 # if m1==0: return 0.
  tet=[t*pi/500 for t in range(500)]
  out1=trapz(real(Ylm0[l1][m1+l1]*conj(Ylm0[l2][m2+l2])*Ylm0[l3][m3+l3]),x=tet)
  out2=trapz(imag(Ylm0[l1][m1+l1]*conj(Ylm0[l2][m2+l2])*Ylm0[l3][m3+l3]),x=tet)
#  out1=quad(lambda tet: real(Ylm0[l1][m1+l1][int(tet)]*conj(Ylm0[l2][m2+l2][int(tet)])*Ylm0[l3][m3+l3][int(tet)]), 0.0, 500)[0]*pi/500
#  out2=quad(lambda tet: imag(Ylm0[l1][m1+l1][int(tet)]*conj(Ylm0[l2][m2+l2][int(tet)])*Ylm0[l3][m3+l3][int(tet)]), 0.0, 500)[0]*pi/500

  if isnan(out1): # or isnan(out2): 
   print('spec_int',l1,m1,l2,m2,l3,m3,':infty')
   return 0.j
  return complex(round(out1,6),round(out2,6))


def spec_int2(l1,m1,l2,m2,l3,m3,Plm,dPlm):
 # if m1==0: return 0.
# we  use the def. of sph harm: Y_lm=sqrt( (2l+1)/(4pi) (l-m)!/(l+m)! )*P_lm(costheta) *e(im*fi) and make integral over P_lm(costheta)
# the integral over fi is equal 2*pi for m3=m2-1 and 0 otherwise, so we do not need to calc. it.
# we multiply by sin^2(theta) (not sin(tet)) because scipy.special.lpmn(m,l,cos(theta)) return P_lm and dP_lm/d(cos(theta)) and dP_lm/dtet=-sin(theta)*dP_lm/(d(cos(tet)))
  fac=-2*pi*( (2*l1+1)/(4*pi)* sil(l1-m1)/sil(l1+m1) )**0.5 *( (2*l2+1)/(4*pi)* sil(l2+m2)/sil(l2+m2) )**0.5*( (2*l3+1)/(4*pi)* sil(l3-m3)/sil(l3+m3) )**0.5
  tet=[t*pi/500 for t in range(500)]
  out1=fac*trapz(real(array([(sin(t))**2 for t in tet])*dPlm[l1][m1+l1]*conj(Plm[l2][m2+l2])*Plm[l3][m3+l3]),x=tet)
  out2=fac*trapz(imag(array([(sin(t))**2 for t in tet])*dPlm[l1][m1+l1]*conj(Plm[l2][m2+l2])*Plm[l3][m3+l3]),x=tet) 
#  out1=fac*quad(lambda tet: real((sin(tet*pi/500))**2*dPlm[l1][m1+l1][int(tet)]*conj(Plm[l2][m2+l2][int(tet)])*Plm[l3][m3+l3][int(tet)]), 0,500)[0]*pi/500
#  out2=fac*quad(lambda tet: imag((sin(tet*pi/500))**2*dPlm[l1][m1+l1][int(tet)]*conj(Plm[l2][m2+l2][int(tet)])*Plm[l3][m3+l3][int(tet)]), 0.0, 500)[0]*pi/500
  if isnan(out1) or isnan(out2): 
   print('spec_int2',l1,m1,l2,m2,l3,m3,':infty')
   return 0.j
  return complex(round(out1,6),round(out2,6))





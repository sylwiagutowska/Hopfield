from math_stuff import *
'''
Plm,dPlm=generate_Plm(10)
Ylm0=generate_Ylm_0(Plm)
for l1 in range(1,5):
  for m1 in range(-l1,l1+1):
   for l2 in range(2):
    for m2 in range(-l2,l2+1):
     for l3 in range(2):
      for m3 in range(-l3,l3+1):
       if m3!=m2-m1 or l1==0 or abs(m1)>l1-1: b,a=0,0 #(for l=0 d Ylm/d theta= 0) 
       else: 
        b=2*pi*(l1*(l1+1)/(((2*l1+1)*(2*l1+3))**0.5)*spec_int(l1+1,m1,l2,m2,l3,m3,Ylm0) \
              -l1*(l1-1)/(((2*l1-1)*(2*l1+1))**0.5)*spec_int(l1-1,m1,l2,m2,l3,m3,Ylm0)) 
        a=spec_int2(l1,m1,l2,m2,l3,m3,Plm,dPlm)
       if m3==m2-m1: print  (l1,m1,l2,m2,l3,m3,a,b)
'''

'''
for l1 in range(1):
  for m1 in range(-l1,l1+1):
   for l2 in range(5):
    for m2 in range(-l2,l2+1):
     for l3 in range(5):
      for m3 in range(-l3,l3+1):
       a=three_y(l1,m1,l2,m2,l3,m3)
       if a!=0: print l1,m1,l2,m2,l3,m3,a
'''

'''
import numpy as np
b_vec=np.array([\
   [0.000000,  0.160357  ,0.160357],\
   [0.160357,  0.000000  ,0.160357],\
   [0.160357,  0.160357  ,0.000000],\
])


kp=np.array([0,-1,-1])
kp=np.array([0   , 0  , -2])
kp=np.array([1   , -1  , -2])
print np.matmul(np.linalg.inv(b_vec),np.transpose(kp))
'''

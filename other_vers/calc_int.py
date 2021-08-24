from sympy import *
from sympy import Function, Symbol
l1,m1,l2,m2,l3,m3 = symbols('l1 m1 l2 m2 l3 m3', integer=True)
l2 = Symbol('l2', integer=True)
theta = Symbol("theta")
phi = Symbol("phi")

#integrate (sin(theta)*exp(-I*phi)*Ynm(l1, m1, theta, phi)*conjugate(Ynm(l2, m2, theta, phi))*Ynm(l3, m3, theta, phi), (theta,0,pi),(phi,0,2*pi))
#print Ynm(0,0,0,0)
#print integrate (sin(theta)*Ynm(l1, m1, theta, phi)*conjugate(Ynm(l2, m2, theta, phi)), (phi,0,2*pi),(theta,0,pi))

print diff((phi**2-1)**l1,phi,l1)
print diff(2*phi*(phi**2 - 1)**l1*(l1*log(phi**2 - 1) + 1)/(phi**2 - 1),phi,l1)

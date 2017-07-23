#
# cp6: Solution to Computing Project 6: Quantum Wire Project
#

from vpython import *
from numpy import *
from vpython import rate

canvas()
g1 = gcurve(color=color.black)

hbar=1.05e-34                    # Js
m=0.067*9.1e-31                  # and m_eff=kg
NA=80                            # how many arrows?
NA2=int(NA/2)                    # half of the arrows
a=5.34e-9                        # width for 5 bound states and 1eV well depth.

x = linspace(-2.0*a, 2.0*a, NA)  # NA locations from -2a to 2a

z0 = 2.1*pi                      # Not a great choice. Why not? Pick a better one
k0 = z0/a                        # get k0

#
# numerical solutions for z when z0 = 2.1*pi, you should find a better z0 and
# find solutions for that choice.
#

z1 = 1.36274
z2 = 2.71711

k1 = k0*z1/z0
k2 = k0*z2/z0
kap1 = sqrt(k0**2 - k1**2)
kap2 = sqrt(k0**2 - k2**2)

E1 = -(hbar*kap1)**2/(2.0*m)
E2 = -(hbar*kap2)**2/(2.0*m)

w1 = E1/hbar
w2 = E2/hbar
wn=[w1,w2]

t = 0.0
dt = (2*pi/(w2-w1))/200.0

psis = zeros((2,NA),double)

def f1(x):
    return cos(k1*x)
    
def f2(x):
    return sin(k2*x)
    
def f3(x):
    return f1(a)*exp(-abs(kap1*x))/exp(-abs(kap1*a))
    
def f4(x):
    return sign(x)*f2(a)*exp(-abs(kap2*x))/exp(-abs(kap2*a))
    
psis[0] = piecewise(x, [x<-a, (x>=-a)&(x<=a), x>a], [f3, f1, f3])
psis[0] = psis[0]/sqrt((abs(psis[0])**2).sum())
psis[1] = piecewise(x, [x<-a, (x>=-a)&(x<=a), x>a], [f4, f2, f4])
psis[1] = psis[1]/sqrt((abs(psis[1])**2).sum())

#
# Decide how much of each bound state to use in the superposition state:
#
# Edit the "if" statement to pick the one you want. The first
# clause with a '1' wins.
#

if 1:
    # Equal parts 1 and 2
    c1 = 1.0/sqrt(2)
    c2 = 1.0/sqrt(2)
elif 0:
    # all 2 and no 1
    c1=0.0
    c2=1.0
elif 0:
    # all 1 and no 2
    c1=1.0
    c2=0.0

cn=[c1, c2]                               # array of amplitudes
t = 0.0                                   # start at t=0

psi = zeros(NA, complex)                  # construct psi at time '0'
for i in range(len(cn)):
    psi = psi + cn[i]*psis[i]

arrowScale = a/psis[0][NA2]              # scale to make the middle of psis[0] about 3a high

def SetArrowFromCN( cn, a):
    """
    SetArrowWithCN takes a complex number  cn  and an arrow object  a .
    It sets the  y  and  z  components of the arrow s axis to the real 
    and imaginary parts of the given complex number. 
    
    Just like Computing Project 1, except y and z for real/imag.
    """
    a.axis.y = cn.real
    a.axis.z = cn.imag
    a.axis.x = 0.0

alist = []
for i in range(NA):
    alist.append(arrow(pos=vec(x[i],0,0), axis=vec(0,a,0), color=color.red))
    SetArrowFromCN(arrowScale*psi[i],alist[i])

rate(1)  # put a "fake" rate command in to render what we've got so far.

#
# Now, all the arrows are made, and the basis functions and coefficients are 
# set. Create a loop that produces the correct time evolutions.
#


# In[ ]:

while True:
    rate(30)

    t = t+dt
    psi = zeros(NA, complex)
    for i in range(2):
        psi = psi + cn[i]*psis[i]*exp(-1j*wn[i]*t)

    for i in range(NA):
        SetArrowFromCN(arrowScale*psi[i],alist[i])

    xexp=(x*abs(psi*psi)).sum()
    g1.plot(pos=(t,xexp))





#!/Library/Frameworks/Python.Framework/Versions/2.6/bin/python
#
# 468/cp1: Solution for Computing Project 1: 2 Particles in a 1-D box
#

from visual import *
from visual.graph import *
from safcnP1 import SetArrowFromCN  # defined in a file safcnP8.py 

from random import random, seed

scene.background=color.white

gd = gdisplay(title="Prob(x1<a/2 AND x2<a/2,t)", xtitle="t", ytitle="P", 
                foreground=color.black, background=color.white,
                ymax=10.1)

gr = gcurve(color=color.black)
gra = gcurve(color=color.red)

hbar=1.0                       # use units where hbar = 1
m=1.0                          # and m=1.0
NA=20                          # how many arrows per dimension NAxNA grid
N1=[1,2,3,4,5,6,7,8,9,10]      # which fourier terms in for particle 1
N2=[1,2,3,4,5,6,7,8,9,10]      # which fourier terms in for particle 2
a=10.0                         # size of box
dt=0.1                         # time step
t=0.0                          # start time at zero...
arrowScale = 20.0              # just make the arrows long enough to see..
eigenstates = {}               # dictionary for precomputed eigenstates
coefs = {}                     # dictionary for precomputed fourier coef.
omegas = {}                    # dictionary for energies

seed(0)

BOSONS=0
FERMIONS=1

x1, x2 = mgrid[0:NA:NA*1j, 0:NA:NA*1j]*(a/NA) # trick to make x1 and x2 arrays

r0 = vector(-a/2, -a/2, 0)   # place origin of arrows

#
# build arrows.... in 2-D space, store them in a set of nested lists
#

alist = []
for i in range(NA):
    sublist = []
    alist.append(sublist)
    for j in range(NA):
        r = r0 + vector(x1[i,j], x2[i,j], 0)
        sublist.append(cylinder(pos=r, axis=(0,0,1), color=color.red))

#
# compute the eigenstates and store them in a dictionary 'eigenstates'
#

for n1 in N1:
    for n2 in N2:
        if FERMIONS:
            psinm = sin(n1*pi*x1/a)*sin(n2*pi*x2/a) - sin(n2*pi*x1/a)*sin(n1*pi*x2/a)   # compute the n1,n2 energy eigenstate
        elif BOSONS:
            psinm = sin(n1*pi*x1/a)*sin(n2*pi*x2/a) + sin(n2*pi*x1/a)*sin(n1*pi*x2/a)   # compute the n1,n2 energy eigenstate
        else:
            psinm = sin(n1*pi*x1/a)*sin(n2*pi*x2/a)
            
        if (FERMIONS and (n1 != n2)) or (not FERMIONS):
            psinm = psinm/sqrt((abs(psinm)**2).sum())  # normalize it.
            eigenstates[(n1,n2)] = psinm


omega0 = hbar*pi**2/(2*m*a**2)

sumSq = 0.0
for nmPair in eigenstates.keys():
    n1, n2 = nmPair
    psinm = eigenstates[nmPair]                          # get nth basis
    cn1n2 = 2.0*(random()-1)/(n1*n2)    
    if BOSONS or FERMIONS:
        if n1 <= n2:
            sumSq += abs(cn1n2)**2
            coefs[nmPair] = cn1n2                            # save it.
        else:
            coefs[nmPair] = 0.0
    else:
        sumSq += abs(cn1n2)**2
        coefs[nmPair] = cn1n2                            # save it.
        
    omega = omega0*(n1**2+n2**2)                         # get omega for nmPair, multiple of omega0
    omegas[nmPair] = omega
    
for nmPair in eigenstates.keys():
    coefs[nmPair] = coefs[nmPair]/sqrt(sumSq)

#
# build up psi via fourier series
#
        
psi = zeros((NA,NA), complex) # initialize psi

for nmPair in eigenstates.keys():
    psi += coefs[nmPair]*eigenstates[nmPair]

for i in range(NA):
    for j in range(NA):
        SetArrowFromCN(arrowScale*psi[i,j],alist[i][j])

updateScreen = False
npts = 0
esum = 0.0

while True:
    rate(20.0/dt)

    if scene.kb.keys: # event waiting to be processed?
        s = scene.kb.getkey() # get keyboard info
        if s == ' ':
            updateScreen ^= 1
    
    if updateScreen:
        psi = zeros((NA,NA), complex) # initialize psi
        
        #
        # Here's where you put in your code to compute the 
        # wavefunction psi, at later times
        #

        for nmPair in eigenstates.keys():
            psi += coefs[nmPair]*eigenstates[nmPair]*exp(-1j*omegas[nmPair]*t)
                
        for i in range(NA):
            for j in range(NA):
                SetArrowFromCN(arrowScale*psi[i,j],alist[i][j])

        psi = psi/(sqrt((abs(psi)**2).sum()))       
        expectation = (abs(psi**2*abs(x1-x2)).sum())
        esum += expectation
        npts += 1
        gr.plot(pos=(t, expectation))
        gra.plot(pos=(t, esum/npts))
        t += dt

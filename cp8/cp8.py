#!/Library/Frameworks/Python.framework/Versions/2.6/bin/python

#
# cp8: Solution for Computing Project 8: Particle in a 2-D box
#

from visual import *
from visual.graph import *
from safcnP8 import SetArrowFromCN  # defined in a file safcnP8.py 

scene.background=color.white

gd = gdisplay(title="Prob(x<a/2 AND y<a/2,t)", xtitle="t", ytitle="P", 
                foreground=color.black, background=color.white,
                ymax=1.1)

gr = gcurve(color=color.black)

hbar=1.0                       # use units where hbar = 1
m=1.0                          # and m=1.0
NA=20                          # how many arrows per dimension NAxNA grid
NX=[1,2,3,5,6,7,9,10]          # which fourier terms in x direction
NY=[1,2,3,5,6,7,9,10]          # which fourier terms in y direction
a=10.0                         # size of box
dt=0.1                         # time step
t=0.0                          # start time at zero...
arrowScale = 20.0              # just make the arrows long enough to see..
eigenstates = {}               # dictionary for precomputed eigenstates
coefs = {}                     # dictionary for precomputed fourier coef.
omegas = {}                    # dictionary for energies

x, y = mgrid[0:NA:NA*1j, 0:NA:NA*1j]*(a/NA) # trick to make x and y arrays

r0 = vector(-a/2, -a/2, 0)   # place origin of arrows

#
# build arrows.... in 2-D space, store them in a set of nested lists
#

alist = []
for i in range(NA):
    sublist = []
    alist.append(sublist)
    for j in range(NA):
        r = r0 + vector(x[i,j], y[i,j], 0)
        sublist.append(cylinder(pos=r, axis=(0,0,1), color=color.red))

#
# compute the eigenstates and store them in a dictionary 'eigenstates'
#

for nx in NX:
    for ny in NY:
        psinxmy = sin(nx*pi*x/a)*sin(ny*pi*y/a)      # compute the n,m energy eigenstate
        psinxmy = psinxmy/sqrt((abs(psinxmy)**2).sum())  # normalize it.
        eigenstates[(nx,ny)] = psinxmy
        
psi0 = zeros((NA,NA),complex)
psi0[:NA/2,:NA/2] = 1.0
psi0 = psi0/sqrt((abs(psi0)**2).sum())             # get psi at t=0, normalized

omega0 = hbar*pi**2/(2*m*a**2)

for nmPair in eigenstates.keys():
    nx, ny = nmPair
    psinxmy = eigenstates[nmPair]                          # get nth basis
    cnm = ((psi0[:NA/2,:NA/2]*psinxmy[:NA/2,:NA/2]).sum()) # compute fourier coef.
    coefs[nmPair] = cnm                                  # save it.
    omega = omega0*(nx**2+ny**2)                         # get omega for nmPair, multiple of omega0
    omegas[nmPair] = omega                               # save it.

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

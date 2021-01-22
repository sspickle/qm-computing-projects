"""

Copyright (c) 2011, Steve Spicklemire

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
from vpython import *
from numpy import *
from vpython import rate

from ThreeD_QM_Models import TwoStateModelDemo

"""
A 3D representation of a superposition of two 1D quantum states using the y and z directions
 to represent the phasor wave function.
"""
omega = 1.0
L = 1.0
N = 100 # number of phasors to model line..
xarray = linspace(0, L, N)
t=0.0
dt=0.01
relPhase = 0
gndPhase = 0.0  # phase of gnd state, useful for "rotating" with the gnd state

def state(n, x, t):
    """
    return the value of the (not normalized) psi(x,t) for the nth quantum state
    """
        
    if relPhase:
        return sin(n*pi*(x+L/2)/L)*exp(-1j*((n*n - 1)*omega*t + gndPhase))
    else:
        return sin(n*pi*(x+L/2)/L)*exp(-1j*(n*n)*omega*t)

#print "Use 1-4 keys to toggle display of 1:psi1, 2:psi2, 3:psi1+psi2, 4:psi*psi."
#print "Use 5 key to rotate with the ground state phasor."
            
TSM = TwoStateModelDemo(xarray, state)
updateScreen = True

while True:
    rate(100)
    if updateScreen:
        TSM.update(t)
        t+=dt
        
    if scene.kb.keys: # event waiting to be processed?
        s = scene.kb.getkey() # get keyboard info
        if s in ('1','2','3','4'):
            TSM.toggleVisible(int(s))
        if s == '5':
            relPhase ^= 1
            if relPhase:
                gndPhase = omega*t
            else:
                t=gndPhase/omega
        if s==' ':
            updateScreen ^= 1

            


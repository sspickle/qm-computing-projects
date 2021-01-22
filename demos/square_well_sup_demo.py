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
from __future__ import print_function

from visual import *
scene = display(title='ISW Superposition Coefs', width=640, height=480)

from ThreeD_QM_Models import SqWellSupDemo

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
gndPhase = 0.0  # phase of gnd state, useful for "rotating" with the gnd state

def state(n, x, t):
    """
    return the value of the (not normalized) psi(x,t) for the nth quantum state
    """
    return (1.0/(pi*n))*(1.0-cos(n*pi/2))*sin(n*pi*(x+L/2)/L)*exp(-1j*(n*n)*omega*t)
            
TSM = SqWellSupDemo(xarray, state, Nterms=7)
updateScreen = False
TSM.update(0.0)

print("'space' to start/stop time evolution.")
print("Type 1-7 for first seven eigenstates.")
print("Type 's' for the superposition wavefunction of all displayed eigenstates.")
print("Type 'p' for a representation of the probability density of the superposition.")
print("right and left arrow to single step forward and backward in time.")

while True:
    rate(100)
    if updateScreen:
        TSM.update(t)
        t+=dt
    if scene.kb.keys: # event waiting to be processed?
        s = scene.kb.getkey() # get keyboard info
        if s in ('1','2','3','4','5','6','7','s','p'):
            if s=='s':
                s=8
            if s=='p':
                s=9
            TSM.toggleVisible(int(s))
        elif s=='r':
            t=0.0
            TSM.update(t)
            updateScreen=False
        elif s==' ':
            updateScreen ^= 1
        elif s=='right':
            t+=dt
            TSM.update(t)
            updateScreen=False
        elif s=='left':
            t-=dt
            TSM.update(t)
            updateScreen=False

            


            


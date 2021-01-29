"""

Copyright (c) 2011, Steve Snp.picklemire

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell conp.pies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all conp.pies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
import vpython as vp
import numpy as np

class PhasorModel:
    """
    Model a superposition of two states in a square well potential
    """

    def __init__(self, xarray, state, k, phaseToColor=False):
        L=xarray[-1]-xarray[0]
        N=len(xarray)
        self.state = state
        self.k = k
        self.phaseToColor = phaseToColor
        self.slist = [(vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.001, color=vp.color.red, axis=vp.vec(0,.1,0),opacity=1, visible=True),
                       vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.001, color=vp.color.green, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.001, color=vp.color.blue, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       ) for x in xarray]
        
    def update(self, t):
        for s1,s2,s3 in self.slist:
            psi = self.state(self.k,s1.pos.x,t)/5
            s1.axis.y = psi.real
            s1.axis.z = psi.imag
            s2.axis.y = psi.real
            s2.axis.z = 0.0
            s3.axis.y = 0.0
            s3.axis.z = psi.imag
            if self.phaseToColor:
                phase1 = np.arctan2(psi.imag, psi.real)
                if phase1 < 0.0:
                    phase1 += 2*np.pi
                s1.color = vp.color.hsv_to_rgb((1.0 - phase1/(2*np.pi),1.0,1.0))
            
    def toggleVisible(self, index):
        for s in self.slist:
            s[index-1].visible ^= True

class TravellingWaveModel:
    """
    Model a superposition of two states in a square well potential
    """

    def __init__(self, xarray, state, k, phaseToColor=False):
        L=xarray[-1]-xarray[0]
        N=len(xarray)
        self.state = state
        self.k = k
        self.phaseToColor = phaseToColor
        self.slist = [(vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.001, color=vp.color.red, axis=vp.vec(0,.1,0),opacity=1, visible=True),
                       vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.001, color=vp.color.green, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.001, color=vp.color.blue, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.cylinder(pos=vp.vec(x,0,0), radius=.1, color=vp.color.green, axis=vp.vec(L/N,0,0),opacity=1, visible=False),
                       
                       ) for x in xarray]
        
    def update(self, t):
        for s1,s2,s3,s4 in self.slist:
            psi1 = self.state(self.k,s1.pos.x,t)/5
            psi2 = self.state(-self.k,s1.pos.x,t)/5
            psi = psi1 + psi2
            s1.axis.y = psi1.real
            s1.axis.z = psi1.imag
            s2.axis.y = psi2.real
            s2.axis.z = psi2.imag
            s3.axis = s1.axis + s2.axis
            if s1.visible and s2.visible:
                s4.radius = abs(psi)**2
            elif s2.visible:
                s4.radius = abs(psi2)**2
            elif s1.visible:
                s4.radius = abs(psi1)**2
            if self.phaseToColor:
                phase1 = arctan2(psi1.imag, psi1.real)
                if phase1 < 0.0:
                    phase1 += 2*np.pi
                phase2 = arctan2(psi2.imag, psi2.real)
                if phase2 < 0.0:
                    phase2 += 2*np.pi
                phase = arctan2(psi.imag, psi.real)
                if phase < 0.0:
                    phase += 2*np.pi
                s1.color = vp.color.hsv_to_rgb((1.0 - phase1/(2*np.pi),1.0,1.0))
                s2.color = vp.color.hsv_to_rgb((1.0 - phase2/(2*np.pi),1.0,1.0))
                s3.color = vp.color.hsv_to_rgb((1.0 - phase/(2*np.pi),1.0,1.0))
            
    def toggleVisible(self, index):
        for s in self.slist:
            s[index-1].visible ^= True


class TwoStateModel:
    """
    Model a superposition of two states in a square well potential
    """

    def __init__(self, xarray, state):
        L=xarray[-1]-xarray[0]
        N=len(xarray)
        self.state = state
        self.slist = [(vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.001, color=vp.color.red, axis=vp.vec(0,.1,0),opacity=1, visible=True),
                       vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.001, color=vp.color.green, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.001, color=vp.color.blue, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.cylinder(pos=vp.vec(x,0,0), radius=.1, color=vp.color.green, axis=vp.vec(L/N,0,0),opacity=1, visible=False),
                       
                       ) for x in xarray]
        
    def update(self, t):
        for s1,s2,s3,s4 in self.slist:
            psi1 = self.state(1,s1.pos.x,t)/5
            psi2 = self.state(2,s1.pos.x,t)/5
            psi = psi1 + psi2
            s1.axis.y = psi1.real
            s1.axis.z = psi1.imag
            s2.axis.y = psi2.real
            s2.axis.z = psi2.imag
            s3.axis = s1.axis + s2.axis
            if s4.visible:
                if (s1.visible and s2.visible) or s3.visible:
                    s4.radius = abs(psi)**2
                elif s1.visible:
                    s4.radius = abs(psi1)**2
                elif s2.visible:
                    s4.radius = abs(psi2)**2
                else:
                    s4.radius = 0.0
            
    def toggleVisible(self, index):
        for s in self.slist:
            s[index-1].visible ^= True

class SHOSupDemo:
    """
    Model a superposition of 7 states in a SHO potential
    
    0-6: quantum states 1-7
    7:   superposition
    8:   prob distr
    """

    def __init__(self, xarray, state, Nterms=7):
        self.L = L=xarray[-1]-xarray[0]
        N=len(xarray)
        self.xarray = xarray
        self.state = state
        self.Nterms=Nterms
        self.slist = [(vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.01, color=vp.color.red, axis=vp.vec(0,.1,0),opacity=1, visible=True),
                       vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.01, color=vp.color.green, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.01, color=vp.color.blue, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.01, color=vp.color.yellow, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.01, color=vp.color.cyan, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.01, color=vp.color.gray(0.7), axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.01, color=vp.color.magenta, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.01, color=vp.color.orange, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.cylinder(pos=vp.vec(x,0,0), radius=.1, color=vp.color.yellow, axis=vp.vec(L/N,0,0),opacity=1, visible=False)
                       ) for x in xarray]
                       
    def SetArrowFromCN(self, cn, a):
        a.axis.y=cn.real
        a.axis.z=cn.imag
        
    def update(self, t):
        
        self.visible = []
        for j in range(len(self.slist)):
            stuple=self.slist[j]
            for i in range(7):
                if j==0:
                    self.visible.append(stuple[i].visible)
                self.SetArrowFromCN(self.state(i,j,t), stuple[i])
                
        psiNet = np.zeros(len(self.xarray),complex)
        for i in range(self.Nterms):
            if i<len(self.visible) and self.visible[i]:
                psiNet += self.state(i, -1, t)
            
        C=sqrt((abs(psiNet)**2).sum())
        if C:
            psiNetNorm = 3*psiNet/C
        else:
            psiNetNorm = psiNet
        for i in range(len(self.slist)):
            stuple =self.slist[i]
            self.SetArrowFromCN(psiNet[i], stuple[7])
            stuple[8].radius = abs(psiNetNorm[i])**2
            
    def toggleVisible(self, index):
        for stuple in self.slist:
            stuple[index-1].visible ^= True
            
    def bumpVisible(self, amount):
        stuple = self.slist[0]
        j=0
        while j<7 and not stuple[j].visible:
            j+=1
        if j<7: # we found one!
            k=j+amount
            if k < 0 or k > 6:
                pass
            else:
                #
                # we have valid values!
                for stuple in self.slist:
                    for i in range(7):
                        if i==k:
                            stuple[i].visible = True
                        else:
                            stuple[i].visible = False

                

class SqWellSupDemo:
    """
    Model a superposition of 7 states in a square well potential
    
    0-6: quantum states 1-7
    7:   superposition
    8:   prob distr
    """

    def __init__(self, xarray, state, Nterms=7):
        self.L = L=xarray[-1]-xarray[0]
        N=len(xarray)
        self.xarray = xarray
        self.state = state
        self.Nterms=Nterms
        self.slist = [(vp.arrow(pos=vp.vec(x-L/2,0,0), shaftwidth=0.001, color=vp.color.red, axis=vp.vec(0,.1,0),opacity=1, visible=True),
                       vp.arrow(pos=vp.vec(x-L/2,0,0), shaftwidth=0.001, color=vp.color.green, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x-L/2,0,0), shaftwidth=0.001, color=vp.color.blue, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x-L/2,0,0), shaftwidth=0.001, color=vp.color.yellow, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x-L/2,0,0), shaftwidth=0.001, color=vp.color.cyan, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x-L/2,0,0), shaftwidth=0.001, color=vp.color.gray(0.7), axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x-L/2,0,0), shaftwidth=0.001, color=vp.color.magenta, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x-L/2,0,0), shaftwidth=0.001, color=vp.color.orange, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.cylinder(pos=vp.vec(x-L/2,0,0), radius=.1, color=vp.color.yellow, axis=vp.vec(L/N,0,0),opacity=1, visible=False)
                       ) for x in xarray]
                       
    def SetArrowFromCN(self, cn, a):
        a.axis.y=cn.real
        a.axis.z=cn.imag
        
    def update(self, t):
        
        for stuple in self.slist:
            for i in range(7):
                self.SetArrowFromCN(self.state(i+1,stuple[0].pos.x,t), stuple[i])
                
        psiNet = np.zeros(len(self.xarray),'d')*(1.0+0j)
        for i in range(self.Nterms):
            if self.slist[0][i].visible:
                psiNet += self.state(i+1, self.xarray-self.L/2.0, t)
            
        for i in range(len(self.slist)):
            stuple =self.slist[i]
            self.SetArrowFromCN(psiNet[i], stuple[7])
            stuple[8].radius = abs(psiNet[i])**2
            
    def toggleVisible(self, index):
        for stuple in self.slist:
            stuple[index-1].visible ^= True

class TwoStateModelDemo:
    """
    Model a superposition of two states in a square well potential
    """

    def __init__(self, xarray, state):
        L=xarray[-1]-xarray[0]
        N=len(xarray)
        self.state = state
        self.slist = [(vp.arrow(pos=vp.vec(x-L/2,0,0), shaftwidth=0.001, color=vp.color.red, axis=vp.vec(0,.1,0),opacity=1, visible=True),
                       vp.arrow(pos=vp.vec(x-L/2,0,0), shaftwidth=0.001, color=vp.color.green, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.arrow(pos=vp.vec(x-L/2,0,0), shaftwidth=0.001, color=vp.color.blue, axis=vp.vec(0,.1,0),opacity=1, visible=False),
                       vp.cylinder(pos=vp.vec(x-L/2,0,0), radius=.1, color=vp.color.green, axis=vp.vec(L/N,0,0),opacity=1, visible=False),
                       
                       ) for x in xarray]
        
    def update(self, t):
        for s1,s2,s3,s4 in self.slist:
            psi1 = self.state(1,s1.pos.x,t)/5
            psi2 = self.state(2,s1.pos.x,t)/5
            psi = psi1 + psi2
            s1.axis.y = psi1.real
            s1.axis.z = psi1.imag
            s2.axis.y = psi2.real
            s2.axis.z = psi2.imag
            s3.axis = s1.axis + s2.axis
            if s4.visible:
                if (s1.visible and s2.visible) or s3.visible:
                    s4.radius = abs(psi)**2
                elif s1.visible:
                    s4.radius = abs(psi1)**2
                elif s2.visible:
                    s4.radius = abs(psi2)**2
                else:
                    s4.radius = 0.0
            
    def toggleVisible(self, index):
        for s in self.slist:
            s[index-1].visible ^= True


class WaveFunctionRep:
    """
    Model a wave-function in 3D as phasors.
    """

    def __init__(self, xarray):
        L=xarray[-1]-xarray[0]
        N=len(xarray)
        self.slist = [(vp.arrow(pos=vp.vec(x,0,0), shaftwidth=0.0001, color=vp.color.red, axis=vp.vec(0,.01,0),opacity=1),
                       vp.cylinder(pos=vp.vec(x,0,0), radius=.01, color=vp.color.green, axis=vp.vec(L/N,0,0),opacity=1),
                       ) for x in xarray]
        
    def update(self, psi):
        for i in range(len(self.slist)):
            s1,s2 = self.slist[i]
            pval = psi[i]
            s1.axis.y = pval.real*20
            s1.axis.z = pval.imag*20
            s2.radius = abs(pval)**2
            

#
# Gabbiani & Cox,  Mathematics for Neuroscientists
#
# Solve the Passive Membrane Equation subject to a
# current pulse, via Backward Euler
#
#  usage   beps(dt,Tfin)
#
#  e.g.    beps(0.01,40)
#

from scipy import *
from matplotlib.pyplot import *

close('all')

def beps(dt,Tfin):

    VCl = -68	    # mV
    A = 4*pi*1e-6	# cm^2 patch area
    Cm = 1         # micro F/cm^2
    gCl = 0.3      # mS/cm^2
    tau = Cm/gCl   # ms
    
    Nt = round(1+Tfin/dt)                  # number of time steps
    v = zeros(Nt)  
    t = zeros(Nt) 
    
    v[0] = VCl                 # initialize v
    
    for j in range(1,Nt):
       t[j] = j*dt
       Istim = (t[j]>2)*(t[j]<22)*1e-5    # 10 pA 20 ms pulse
       v[j] = (v[j-1] + dt*(VCl/tau + Istim/A/Cm))/(1+dt/tau)  # backward euler
    
    plot(t,v, linewidth=2)
    xlabel('t  (ms)',fontsize=16)
    ylabel('V  (mV)',fontsize=16)
    
    

#
# Gabbiani & Cox,  Mathematics for Neuroscientists
# 
# Solve the Passive Membrane Equation subject to a
# synaptic square wave of period T, via the trapezoid rule
#
#  usage   curvssyn(dt,gbar,T,Tfin)
#
#  e.g.    curvssyn(0.01,0.1,2,40)
#          curvssyn(0.01,0.1,4,40)
#

from scipy import *
from matplotlib.pyplot import *

close('all')

def curvssyn(dt,gbar,T,Tfin):

    VCl = -68	    # mV
    A = 4*pi*1e-6	# cm^2 patch area
    Cm = 1         # micro F/cm^2
    gCl = 0.3      # mS/cm^2
    tau = Cm/gCl   # ms
    E = 60
    z = 2*Cm/dt
    
    #Tfin = 10*T
    Nt = round(1+Tfin/dt)                  # number of time steps
    v = zeros(Nt)  
    w = zeros(Nt) 
    t = zeros(Nt)  
    sq = zeros(Nt) # preallocate space
    
    g0 = 0
    
    for j in range(1,Nt):
    
       t[j] = j*dt
    
       sq[j] = (mod(t[j],T)<1) 
       g1 = gbar*sq[j]
       
       v[j] = ( (z-gCl-g0)*v[j-1] + (g0+g1)*E ) / (z+gCl+g1)
       
       w[j] = ( (z-gCl)*w[j-1] + (g0+g1)*E )/ (z+gCl)
       
       g0 = g1
    
    
    plot(t,v,'k')
    plot(t,w,'r')
    legend(['v (syn input)','w (current approx)'],loc='best')
    plot(t,sq)
    xlabel('t  (ms)',fontsize=16)
    ylabel('(mV)',fontsize=16)
    

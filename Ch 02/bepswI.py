#
# Gabbiani & Cox,  Mathematics for Neuroscientists
# 
# Solve the Passive Membrane Equation subject to a
# current pulse, via Backward Euler and plot V and both
# currents
#
#  usage   bepswI(dt,Tfin)
#
#  e.g.    bepswI(0.01,39)
#

from scipy import *
from matplotlib.pyplot import *

def bepswI(dt,Tfin):    # function definition

    close('all')
    
    VCl = -68	    # mV
    A = 4*pi*1e-6   	# cm^2 patch area
    Cm = 1         # micro F/cm^2
    gCl = 0.3      # mS/cm^2
    tau = Cm/gCl   # ms
    
    Nt = round(1+Tfin/dt)                  # number of time steps
    v = zeros(Nt)  
    t = zeros(Nt)
    ICl = zeros(Nt) 
    IC = zeros(Nt)  
    
    v[0] = VCl                 # initialize v
    
    for j in range(1,Nt):
       t[j] = j*dt
       Istim = (t[j]>2)*(t[j]<22)*1e-5    # 10 pA 20 ms pulse
       v[j] = (v[j-1] + dt*(VCl/tau + Istim/A/Cm))/(1+dt/tau)  # backward euler
       ICl[j] = gCl*(v[j]-VCl)  # chloride current
       IC[j] = Istim/A - ICl[j] # capacitive current
       
    subplot(1,2,1)
    plot(t, v)
    xlabel('t  (ms)',fontsize=14)
    ylabel('V  (mV)',fontsize=14)
    
    subplot(1,2,2)
    plot(t, IC, 'r', t, ICl, 'k')
    legend(['$I_C$','$I_{Cl}$'], loc='best')
    xlabel('t  (ms)',fontsize=14)
    ylabel('I  $(\mu A/cm^2)$',fontsize=14)




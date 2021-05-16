#
# Gabbiani & Cox,  Mathematics for Neuroscientists
# 
# Solve the Passive Membrane Equation subject to a
# synaptic input, via the trapezoid rule
#

from scipy import *
from matplotlib.pyplot import *

def trapsyn(dt,Tfin,gsyn):

    VCl = -68	    # mV
    Cm = 1         # micro F/cm^2
    gCl = 0.3      # mS/cm^2
    z = 2/dt
    
    Nt = round(1+Tfin/dt)                    # number of time steps
    V = zeros(Nt)               # preallocate space
    t = zeros(Nt)  	 
    g = zeros((Nt,size(gsyn['gmax'])))
    
    V[0] = VCl
    a0 = gCl/Cm
    b0 = gCl*VCl/Cm
    
    for j in range(1,Nt):
       t[j] = j*dt
       tmp = (t[j]-gsyn['t1'])/gsyn['taua']
       g[j,:] = gsyn['gmax']*tmp*exp(1-tmp)*(t[j]>gsyn['t1'])
       a1 = (gCl + sum(g[j,:]))/Cm
       b1 = ( gCl*VCl + dot(g[j,:], gsyn['Vsyn']) )/Cm
       V[j] = ( (z-a0)*V[j-1] + b0 + b1 ) / (z + a1)
       a0 = a1
       b0 = b1
      
    return(t,V,g) 

close('all')   # now drive this function
    
gsyn = {     #   single excitatory conductance
   'gmax': 0.2, 
   'taua': 2, 
   't1': 5, 
   'Vsyn': 0
}
(t, V1, g1) = trapsyn(0.01, 35, gsyn)

gsyn = {        # an ihibitory and an excitatory conductance
   'gmax': 0.2*array([1, 1]), 
   'taua': 2*array([1, 1]), 
   't1': array([4, 5]), 
   'Vsyn': array([-68, 0])
}
(t, V2, g2) = trapsyn(0.01, 35, gsyn)

figure(1)
plot(t,g1,'k')
plot(t,g2[:,0],'r',t,g2[:,1],'r--')

legend(['excitatory, solo','inhibitory, paired','excitatory, paired'], loc='best')
xlabel('t  (ms)', fontsize=14)
ylabel('$g_{syn}$  $(\mu S/cm^2)$', fontsize=14)

figure(2)
plot(t,V1,'k')
plot(t,V2,'r')
legend(['excitatory only','inhibitory and excitatory'], loc='best')
xlabel('t  (ms)', fontsize=14)
ylabel('V  (mV)', fontsize=14)






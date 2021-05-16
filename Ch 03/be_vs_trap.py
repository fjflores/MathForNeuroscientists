#
# Gabbiani & Cox,  Mathematics for Neuroscientists
#
# Contrast the accuracy of Backward Euler and Trapezoid
# in the solution of the passive isopotential cell
#
#   be_vs_trap
#

from scipy import *
from matplotlib.pyplot import *

close('all')

VCl = -68      # mV        
A = 4*pi*1e-6  # cm^2 patch area
Cm = 1         # micro F/cm^2
gCl = 0.3      # mS/cm^2
tau = Cm/gCl

# Backward Euler on a Passive Sphere with alpha-like input

def bepsa(dt,Tfin):
    
    j = 1
    v = zeros(1+int(Tfin/dt))
    v[0] = VCl
    t = 0
    
    while (t < Tfin):
       t = j*dt
       Istim = (2e-5)*t*exp(-t/tau)/tau
       v[j] = (v[j-1] + dt*(VCl/tau + Istim/A/Cm))/(1+dt/tau)
       j = j + 1
   
    return v

# Trapezoid Rule on a Passive Sphere with alpha-like input

def trapsa(dt,Tfin):

    j = 1
    t = 0
    v = zeros(1+int(Tfin/dt))
    Istim = zeros(1+int(Tfin/dt))
    v[0] = VCl
    Istim[0] = 0
    
    while (t < Tfin):
       t = j*dt
       Istim[j] = (2e-5)*t*exp(-t/tau)/tau
       v[j] = ( (1-dt/tau/2)*v[j-1] + dt*( VCl/tau + (Istim[j]+Istim[j-1])/2/A/Cm ) )/(1+dt/tau/2)
       j = j + 1
      
    return v

dt = array([.5, .1, .05, .01, .005, .001])   # compare methods a many dt
be_err = zeros(6)
cn_err = zeros(6)
Tfin = 20

for k in range(6):

    vb = bepsa(dt[k],Tfin)

    cnv = trapsa(dt[k],Tfin)

    t = linspace(0,20,1+int(Tfin/dt[k]))

    vex = VCl + (2e-5)*exp(-t/tau)*(t**2)/2/A/Cm/tau

    be_err[k] = amax( abs(vb-vex) )

    cn_err[k] = amax( abs(cnv-vex) )

loglog(dt, be_err, 'kx-')
loglog(dt, cn_err, 'ro-')
xlabel('dt  (ms)',fontsize=14)
ylabel('maximum absolute error  (mV)', fontsize=14)
legend(['Backward Euler','Trapezoid'], loc='best')
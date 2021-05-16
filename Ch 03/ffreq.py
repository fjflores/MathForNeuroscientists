#
# Gabbiani & Cox, Mathematics for Neuroscientists
#
# frequency response of the passive isopotential cell
#

from scipy import *
from matplotlib.pyplot import *

close('all')

A = 4*pi*(1e-6)    # (cm)^2
Cm = 1             # muF/(cm)^2
GCl = 0.3          # mS/(cm)^2

tau = Cm/GCl       #in msec
tau_s = tau*1e-3   #converted to sec
f = linspace(0,100,1000)      #Hz

#because Cm is in muF/cm2 and A is in cm2, rin is in MOhms

rin = 1/(A*Cm*sqrt( (2*pi*f)**2 + 1/tau_s**2 ))

subplot(1,2,1)
plot(f,rin,'k')
xlabel('$\omega$  (Hz)',fontsize=14)
ylabel('$R_{in}$  $(M\Omega)$',fontsize=14)

ph = -arctan(2*pi*f*tau_s)
ph_deg = ph*180/pi
subplot(1,2,2)
plot(f,ph_deg,'k')
xlabel('$\omega$ (Hz)',fontsize=14)
ylabel('phase (deg)',fontsize=14)





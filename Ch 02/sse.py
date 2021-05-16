#
# sse.py
#
# Gabbiani & Cox, Mathematics for Neuroscientists
#
# steady-state depolarization as a function of normalized conductance for 
# one excitatory synapse
#

from scipy import *
from matplotlib.pyplot import *

close('all')

ce = linspace(0,15,150)     # normalized conductance

Vss = -68/(1+ce)      # steady-state response

plot(ce, Vss, 'k')
xlabel('$c_{syn}$',fontsize=14)
ylabel('$V_{ss}$ (mV)',fontsize=14)



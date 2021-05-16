#
# Gabbiani & Cox,  Mathematics for Neuroscientists
# 
# ss2e.py
#

from scipy import *
from matplotlib.pyplot import *

close('all')

VCl = -68
ce = linspace(0,15,150)

Vsse = VCl - (ce/(1+ce))*VCl

Vss2e = VCl - (2*ce/(1+2*ce))*VCl

Vsse2 = VCl - (2*ce/(1+ce))*VCl

plot(ce,Vsse,'k')
plot(ce,Vss2e,'r')
plot(ce,Vsse2,'k--')

legend(['$V_{ss,e}$','$V_{ss,2e}$','$V_{ss,e2}$'], loc='best')
xlabel('$c_e$', fontsize=14)
ylabel('mV', fontsize=14)

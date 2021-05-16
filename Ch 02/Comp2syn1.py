#
# Gabbiani & Cox,  Mathematics for Neuroscientists
# 
# Comp2syn1.py
#

from scipy import *
from matplotlib.pyplot import *

close('all')

ce = logspace(-2,3,100)

ci = array([0, 1, 2])
sty = ['k', 'r', 'k--']

for i in range(3):
    vp = (100/9)*ce/((ci[i]+10/9)*(ce+10/9)-1/81)
    semilogx(ce, vp, sty[i])

legend([ '$c_i$ = 0', '$c_i = 1$','$c_i = 2$' ], loc='best')
xlabel('$c_e$', fontsize=14)
ylabel('$v_p$', fontsize=14)




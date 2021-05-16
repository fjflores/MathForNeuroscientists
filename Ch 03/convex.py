#
# Gabbiani & Cox,  Mathematics for Neuroscientists
# 
# explicit example of convolution
#

from scipy import *
from matplotlib.pyplot import *

close('all')

t = linspace(0,10,1000)
t1 = 1
t2 = 2
f = exp(-t/t1)
g = exp(-t/t2)
fcg = (t1*t2)*(exp(-t/t1)-exp(-t/t2))/(t1-t2)

plot(t,f,'k')
plot(t,g,'k--')
plot(t,fcg,'r')
legend(['exp(-t)','exp(-t/2)','their convolution'])
xlabel('t',fontsize=14)
%
% ghk.m
%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% plot the driving force Phi
%

V = -100:.1:50;
tmp = exp(2*V/25.8);
u = 5e-5;
Phi = V.*(1-u*tmp)./(1-tmp);
plot(V,Phi,'color','k')
grid
xlabel('V   (mV)','fontsize',14)
ylabel('\Phi  (mV)','fontsize',14)
box off
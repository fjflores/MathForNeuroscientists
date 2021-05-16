%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  Rinxs.m
%
%  Graph the input resistance as a function of stimulus location
%

g_L = 1/15; %0.3;     		% mS / cm^2
R_2 = 0.3; % 0.034;		% k Ohm cm
a = 1e-4;
ell = 1e-1;
dx = 1e-4;
x = dx/2:dx:ell-dx/2;       % vector of patch midpoints

lambda = sqrt(a/(2*R_2*g_L));

Rin = 1./(2*pi*a*lambda*g_L*(tanh(x/lambda)+tanh((ell-x)/lambda)));

plot(x,Rin*1e-3,'k')
box off
xlabel('x_s  (cm)','fontsize',14)
ylabel('R_{in}  (M\Omega)','fontsize',14)

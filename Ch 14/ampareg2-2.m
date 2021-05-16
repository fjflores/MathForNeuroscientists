%
% Gabbiani and Cox, Mathematics for Neuroscientists
%
% ampareg2.m    with calcium now dependent on spine potential, V, and
%       presynaptic frequency, f, compute
%       dependence of AMPA conductance on V and f
%
%  From Castellini et al 2001
%
%  usage:    ampareg2
%

close all

[V,f] = meshgrid(-80:0,1:50);

Vr = 130;

c = -(V-Vr).*f./(1+exp(-0.062*V)/3.57)/200;   % 200 puts c in phys range

km1 = 1+30*c.^2./(1+c.^2);
km2 = km1;
k1 = 1+100*c.^2./(64+c.^2);
k2 = k1;

den = km1.*km2 + k1.*km2 + k2.*km1 + k1.*k2;

a0 = km1.*km2./den;
a1 = k1.*km2./den;
a2 = k2.*km1./den;
a3 = k1.*k2./den;

g = a0 + 2*(a1+a2) + 4*a3;

figure(1)
mesh(V,f,c)
xlabel('V  (mV)','fontsize',14)
ylabel('f  (Hz)','fontsize',14)
zlabel('c','fontsize',14)
text(-75,50,25,'(A)','fontsize',20)
colormap([0 0 0])

figure(2)
mesh(V,f,g)
xlabel('V  (mV)','fontsize',14)
ylabel('f  (Hz)','fontsize',14)
zlabel('g_{AMPA}','fontsize',14)
text(-75,50,3,'(B)','fontsize',20)
colormap([0 0 0])

f = 0:0.1:50;
V = 4*f/5 - 60;
c = -(V-Vr).*f./(1+exp(-0.062*V)/3.57)/200;   % 200 puts c in phys range

km1 = 1+30*c.^2./(1+c.^2);
km2 = km1;
k1 = 1+100*c.^2./(64+c.^2);
k2 = k1;

den = km1.*km2 + k1.*km2 + k2.*km1 + k1.*k2;

a0 = km1.*km2./den;
a1 = k1.*km2./den;
a2 = k2.*km1./den;
a3 = k1.*k2./den;

g = a0 + 2*(a1+a2) + 4*a3;

figure(3)
plot(f,g,'k','linewidth',1.5)
hold on

c = -3*(V-Vr).*f./(1+exp(-0.062*V)/3.57)/200;   % 200 puts c in phys range

km1 = 1+30*c.^2./(1+c.^2);
km2 = km1;
k1 = 1+100*c.^2./(64+c.^2);
k2 = k1;

den = km1.*km2 + k1.*km2 + k2.*km1 + k1.*k2;

a0 = km1.*km2./den;
a1 = k1.*km2./den;
a2 = k2.*km1./den;
a3 = k1.*k2./den;

g = a0 + 2*(a1+a2) + 4*a3;

plot(f,g,'r','linewidth',1.5)
xlabel('f  (Hz)','fontsize',14)
ylabel('g_{AMPA}','fontsize',14)
legend('g_{NMDA} = 0.005','g_{NMDA} = 0.015')
text(5,3,'(C)','fontsize',20)
hold off
box off


%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% htrack.m
%
%  Track V_r and QA eigenvalues of the HH+I_h system over a
%  a range of g_h values
%

close all

figure(1)
V = -140:.1:0;
qinf = 1./(1+exp((V+69)/7.1));
tauq = 1000./(exp((V+66.4)/9.3)+exp(-(V+81.6)/13));
[ax,h1,h2]=plotyy(V,qinf,V,tauq);
set(h1,'color','k')
set(ax(1),'ycolor','k')
set(ax(1),'ytick',[0 0.2 0.4 0.6 0.8 1.0])
set(h2,'color','r')
set(ax(2),'ycolor','r')
set(ax(2),'ylim',[0 1000])
set(ax(2),'ytick',[0 200 400 600 800 1000])
xlabel('V  (mV)','fontsize',14)
ylabel(ax(1),'q_\infty','fontsize',14)
ylabel(ax(2),'\tau_q  (ms)','fontsize',14)
box off

G = struct('K', 36, 'Na', 120, 'Cl', 0.3, 'h', 1);

for g = 0:.04:2,
    G.h = g;
    [Vr,B] = hhsymh(G);
    figure(2)
    plot(g,Vr,'.','color',[g/2 0 0])
    hold on
    figure(3)
    e = eig(B);
    plot(real(e),imag(e),'.','color',[g/2 0 0])
    hold on
end
figure(2)
xlabel('g_h  (mS/cm^2)','fontsize',14)
ylabel('V_r  (mV)','fontsize',14)
box off
figure(3)
xlim([-.2 0.01])
grid
xlabel('real(z)','fontsize',14)
ylabel('imag(z)','fontsize',14)

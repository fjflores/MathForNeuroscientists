%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
%  stEcabQandA.m
%
%  Contrast the active and quasi-active responses to small random input
%
function stEcabQandA

cab = struct('rad',1e-4,'ell',1e-1,'dx',1e-3,'dt',0.03);
g = struct('K', 36, 'Na', 120, 'Cl', 1/15);

stim = struct('t1',1,'t2',2,'amp',20e-6,'loc',0.03,'Tfin',50);

stim.t1 = rand(10,1)*30;
stim.t2 = stim.t1 + 1;
stim.loc = rand(10,1)*cab.ell;

[t,VA] = stEcab(cab,g,stim,0);
plot(t,VA,'k')

[t,VQ] = stEQcab(cab,g,stim,0);

close all
plot(t,VA,'k')
hold on
plot(t,VQ,'r')
hold off
box off
legend('Active','Quasi-Active')
xlabel('t  (ms)','fontsize',14)
ylabel('V_{soma}  (mV)','fontsize',14)



%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% hyEcabCa3traindrive
%
% turn gKCa on and off in a cell with periodic input and
% show that spikes are missed
%

function hyEcabCa3traindrive

cab = struct('rad',.0001,'ell',.1,'N',200,'dt',.02);
stim = struct('pw',1,'per',20,'Tfin',100,'Iapp',3e-4,'BT',5e2,...
              'knaca',100,'Jpmca',2e-6,'gKCa',0);

[v0, t] = hyEcabCa3train(cab,stim,0);

stim.gKCa = 10;

[v10, t] = hyEcabCa3train(cab,stim,0);

plot(t,v0,'k')
hold on
plot(t,v10,'r--')
xlim([0 90])
ylim([-80 50])
box off
hold off
xlabel('t  (ms)','fontsize',14,'color','k')
ylabel('V  (mV)','fontsize',14,'color','k')


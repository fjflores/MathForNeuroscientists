%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
%  hyEcabCa2drive
%
function hyEcabCa2drive

cab = struct('rad',.0001,'ell',.1,'N',200,'dt',.02);
stim = struct('t1',1,'t2',3,'Tfin',10,'Iapp',3e-4,'BT',5e2,'knaca',100,'Jpmca',2e-6);
pinc = 10;

[cafull,t] = hyEcabCa2(cab,stim,pinc);
figure(5)
plot(t,cafull,'k')
hold on
pinc = 0;

stim = struct('t1',1,'t2',3,'Tfin',10,'Iapp',3e-4,'BT',5e2,'knaca',0,'Jpmca',0);
[caB,t] = hyEcabCa2(cab,stim,pinc);
plot(t,caB,'k--')

stim = struct('t1',1,'t2',3,'Tfin',10,'Iapp',3e-4,'BT',5e2,'knaca',100,'Jpmca',0);
[caBx,t] = hyEcabCa2(cab,stim,pinc);
plot(t,caBx,'r')

stim = struct('t1',1,'t2',3,'Tfin',10,'Iapp',3e-4,'BT',0,'knaca',0,'Jpmca',0);
[canull,t] = hyEcabCa2(cab,stim,pinc);
plot(t,canull,'r--')
legend('full','B only','B and Ex','null','location','best')
box off
hold off


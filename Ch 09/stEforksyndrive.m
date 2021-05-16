%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% stEforksyndrive
%

g = struct('K', 36, 'Na', 120, 'Cl', 1/15);
stim = struct('t1',[1 1],'tau',[1 1],'Gsyn',5e-7*[1 1],'Vsyn',[0 0],'loc',[100 350],'Tfin',30);
[t,vtogether]=stEforksyn(stim,g,0);

stim = struct('t1',1,'tau',1,'Gsyn',5e-7,'Vsyn',0,'loc',100,'Tfin',30);
[t,v100]=stEforksyn(stim,g,0);

stim = struct('t1',1,'tau',1,'Gsyn',5e-7,'Vsyn',0,'loc',350,'Tfin',30);
[t,v350]=stEforksyn(stim,g,0);

figure(1)
plot(t,v100(2,:)+v350(2,:)-v100(2,1)-v350(2,1),'r')
hold on
plot(t,vtogether(end,:)-vtogether(end,1),'k')
ax = axis;
box off
legend('sum of single responses','response to simultaneous input')
legend boxoff
xlabel('t  (ms)','fontsize',14)
ylabel('Relative Soma Potential  (mV)','fontsize',14)
hold off
axes('position',[1/5 1/5 1/5 1/5])
plot([.9 .9 ],[.1 .9],'k','linewidth',2)
hold on
plot([.5 .9],[.5 .5],'k','linewidth',2)
plot(.9,.05,'ko','markersize',10)
plot(.7,.5,'rd','markersize',8)
plot(.9,.7,'rd','markersize',8)
axis equal
axis off
hold off

%print -depsc stEforksyndriveA

stim = struct('t1',[1 1],'tau',[1 1],'Gsyn',5e-7*[1 2],'Vsyn',[0 -70],'loc',[625 700],'Tfin',20);
[t,vidist]=stEforksyn(stim,g,0);

stim = struct('t1',[1 1],'tau',[1 1],'Gsyn',5e-7*[2 1],'Vsyn',[-70 0],'loc',[550 625],'Tfin',20);
[t,viprox]=stEforksyn(stim,g,0);

figure(2)
plot(t,viprox(2,:),'r')
hold on
plot(t,vidist(2,:),'k')
box off
legend('distal inhibition','proximal inhibition')
legend boxoff
xlabel('t (ms)','fontsize',14)
ylabel('Soma Potential (mV)','fontsize',14)

hold off
axes('position',[1/5 1/5 1/5 1/5])
plot([.9 .9 ],[.1 .9],'k','linewidth',2)
hold on
plot([.5 .9],[.5 .5],'k','linewidth',2)
plot(.9,0.05,'ko','markersize',10)
plot(.9,.4,'kx','markersize',8)
plot(.9,.2,'kx','markersize',8)
plot(.9,.3,'rd','markersize',8)
axis equal
axis off
hold off

%print -depsc stEforksyndriveB

stim = struct('t1',[1 3],'tau',[1 1],'Gsyn',5e-7*[1 1],'Vsyn',[0 0],'loc',[550 700],'Tfin',20);
[t,vtosoma]=stEforksyn(stim,g,0);

stim = struct('t1',[3 1],'tau',[1 1],'Gsyn',5e-7*[1 1],'Vsyn',[0 0],'loc',[550 700],'Tfin',20);
[t,vfromsoma]=stEforksyn(stim,g,0);

figure(3)
plot(t,vtosoma(2,:),'r')
hold on
plot(t,vfromsoma(2,:),'k')
box off
legend('distal then proximal','proximal then distal')
legend boxoff
xlabel('t (ms)','fontsize',14)
ylabel('Soma Potential (mV)','fontsize',14)
hold off
axes('position',[1/5 1/5 1/5 1/5])
plot([.9 .9 ],[.1 .9],'k','linewidth',2)
hold on
plot([.5 .9],[.5 .5],'k','linewidth',2)
plot(.9,0.05,'ko','markersize',10)
plot(.9,.4,'rd','markersize',8)
plot(.9,.2,'rd','markersize',8)
axis equal
axis off
hold off

%print -depsc stEforksyndriveC


%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% inject a current at a fixed location and of fixed duration but
% of variable amplitude, and track the max potential (and associated time)
% as a function of amplitude
%
%  stEcabthresh
%

function stEcabthresh

close all
stim = struct('t1',1,'t2',2,'amp',4e-4,'loc',0.05,'Tfin',8);
cab = struct('rad',1e-4,'ell',1e-1,'dx',1e-4,'dt',0.02);
g = struct('K', 36, 'Na', 120, 'Cl', 1/15);

Vr = cabrest(cab,g);

for j = 1:250

    I0(j) = 50e-6 + (j-1)*1e-6;
    stim.amp = I0(j);
    [t,Vhot] = stEcabwVr(cab,g,stim,0,Vr);
    [vmax(j),tj] = max(Vhot);
    tmax(j) = t(tj);
    j

end

[ax,h1,h2] = plotyy(I0*10^6,vmax,I0*10^6,tmax)
set(ax,'xlim',[50 300])
set(ax(1),'ycolor','k','ylim',[-80 60],'ytick',[-80 -60 -40 -20 0 20 40 60])
set(ax(2),'ycolor','r','ylim',[0 10],'ytick',[0 2 4 6 8 10])
set(h1,'color','k')
set(h2,'color','r')
xlabel('I_{0} (pA)','fontsize',14)
ylabel(ax(1),'V_{max} (mV)','fontsize',14)
ylabel(ax(2),'t_{max} (ms)','fontsize',14)
box off

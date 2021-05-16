%
% Gabbiani and Cox, Mathematics for Neuroscientists
%
%  gcnstdp.m    illustrate stdp curves
%
%  usage:   gcnstdp
%

Dt = 0:.1:100;   % ms
tau = 6;    % ms

figure(1)
plot(Dt,1.2*exp(-Dt/(4*tau)),'k','linewidth',1.5)
hold on
plot(-Dt,-exp(-Dt/(3*tau)),'k','linewidth',1.5)
plot([-100 100],[0 0],'k')
plot([0 0],[-1.2 1.4],'k')
text(8,1,'1.2exp(-\Deltat/24)')
text(-40,-0.8,'-exp(-\Deltat/18)')
hold off
box off
xlabel('\Deltat  (ms)','fontsize',14)

figure(2)
plot(Dt,-2*exp(-Dt/(2*tau)),'k','linewidth',1.5)
hold on
plot(Dt,-7*exp(-Dt/(2*tau)),'r','linewidth',1.5)
plot(-Dt,exp(-Dt/(4*tau)),'k','linewidth',1.5)
plot(-Dt,3.5*exp(-Dt/(4*tau)),'r','linewidth',1.5)
plot([-100 100],[0 0],'k')
plot([0 0],[-7.5 4],'k')
text(1,1,'exp(-\Deltat/24)')
text(-60,2,'3.5exp(-\Deltat/24)','color','r')
text(-35,-2,'-2exp(-\Deltat/12)')
text(10,-4,'-7exp(-\Deltat/12)','color','r')
hold off
box off
xlabel('\Delta t  (ms)','fontsize',14)

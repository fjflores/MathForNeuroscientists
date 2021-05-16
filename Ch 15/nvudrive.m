%
%  Gabbiani and Cox, Mathematics for Neuroscientists
%
%  nvudrive
%
%  drive nvu and display results 
%
dt = 0.005;
tfin = 300;
Jip3lum = 0.18;
[t,kstim,kplo,Vmlo,Jcamlo,camlo,radlo] = nvu(dt,tfin,Jip3lum);
Jip3lum = 0.4;
[t,kstim,kphi,Vmhi,Jcamhi,camhi,radhi] = nvu(dt,tfin,Jip3lum);

figure(1)
plot(t,kstim*1e5/(4e-2),'k','linewidth',1.5)   % same scaling as astrocyte2
xlabel('t  (s)','fontsize',14)
ylabel('J_{stim}   (10^{-5} cm mM/s)','fontsize',14)
text(10,6,'(A)','fontsize',20)
box off

figure(2)
plot(t,kplo,'k','linewidth',1.5)
hold on
plot(t,kphi,'r','linewidth',1.5)
legend('J_{lum}=0.18','J_{lum}=0.4','location','best')
xlabel('t  (s)','fontsize',14)
ylabel('[K^+]_p   (mM)','fontsize',14)
text(10,11,'(B)','fontsize',20)
hold off
box off

figure(3)
plot(t,Vmlo,'k','linewidth',1.5)
hold on
plot(t,Vmhi,'r','linewidth',1.5)
legend('J_{lum}=0.18','J_{lum}=0.4','location','best')
xlabel('t  (s)','fontsize',14)
ylabel('V_m   (mV)','fontsize',14)
text(50,-20,'(C)','fontsize',20)
hold off
box off

figure(4)
plot(t,Jcamlo,'k','linewidth',1.5)
hold on
plot(t,Jcamhi,'r','linewidth',1.5)
legend('J_{lum}=0.18','J_{lum}=0.4','location','best')
xlabel('t  (s)','fontsize',14)
ylabel('J_{Ca,m}   (\mu M/s)','fontsize',14)
text(50,-0.05,'(D)','fontsize',20)
hold off
box off

figure(5)
plot(t,camlo,'k','linewidth',1.5)
hold on
plot(t,camhi,'r','linewidth',1.5)
legend('J_{lum}=0.18','J_{lum}=0.4','location','best')
xlabel('t  (s)','fontsize',14)
ylabel('[Ca^{2+}]_m   (\mu M)','fontsize',14)
text(50,0.8,'(E)','fontsize',20)
hold off
box off

figure(6)
plot(t,radlo,'k','linewidth',1.5)
hold on
plot(t,radhi,'r','linewidth',1.5)
legend('J_{lum}=0.18','J_{lum}=0.4','location','best')
xlabel('t  (s)','fontsize',14)
ylabel('\rho   (\mu m)','fontsize',14)
text(250,26,'(F)','fontsize',20)
hold off
box off

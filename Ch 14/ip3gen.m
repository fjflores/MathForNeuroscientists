%
% ip3gen.m
%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% generate IP3 via mgluR stimulation
%
% rates are from Destexhe, Mainen and Sejnowski (except last two)
%
% usage: ip3gen
%

function ip3gen

Mrate = [6.6e5*1e-6*1e-3 20e-3 5.3e-3 17e-3 8.3e-8 7.9e-3 10e-3 5e-3];
t1 = 10; t2 = 12; glu = 1e3; Tfin = 600;
dt = 0.1;
M0 = [-Mrate(2) Mrate(3) 0 0;
       Mrate(4) -Mrate(3) 0 0;
       Mrate(4) 0 -Mrate(6) 0;
       0 0 Mrate(7) -Mrate(8)];
M = M0;
nM = M0;
I = eye(4);
gam = [0 0 0 0]';
ngam = gam;

t = 0;

Nt = ceil(Tfin/dt);         % number of time steps
m = zeros(4,Nt);            % preallocate space

for j=1:Nt-1,

     t = j*dt;
     
     stim = Mrate(1)*glu*(mod(t,50)>t1)*(mod(t,50)<t2)*(t<500); 
     
     ngam(1) = stim;
     
     nM(1,1) = M0(1,1) - stim;
     nM(1,2) = M0(1,2) - stim;
     
     m(:,j+1) = (2*I/dt - nM)\ ( (2*I/dt + M)*m(:,j) + gam + ngam );
     
     gam = ngam;
     
     M = nM;
     
end

tim = linspace(0,Tfin,Nt);
plot(tim,m(1,:),'r')
hold on
plot(tim,m(2,:),'k')
plot(tim,m(3,:),'r--')
plot(tim,m(4,:),'k--')
hold off
box off

legend('mA','mI','G \muM','IP3 \muM','location','NW')
legend('boxoff')
xlabel('t  (ms)','fontsize',14)

     

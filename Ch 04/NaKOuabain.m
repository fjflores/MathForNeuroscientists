%
% Gabbiani and Cox Mathematics for Neuroscientists
%
% NaKOuabain, adapted from NaKdyn.m
%
% usage:   NaKOuabain
%
function NaKOubain
Ke = 4; % mM
Ne = 116; % mM
Cle = Ke+Ne; % mM
S = 308;  % mM
P = struct('K',5e-5,'Na',5e-7,'Cl',5e-3); % cm/s

q = 1.55e-6;   % micro mol, amount of impermeants
A = 0.8e-5;    % cm^2 about 3 times that of the equivalent sphere
q = q/A;       % w to u converter

options = optimset('disp','off');

lam = 1e-5; % cm/s
vbar = fsolve(@(v) veq(v,lam,P,Ke,Ne,S),-1,options)
ev = exp(vbar);   
Nibar = Ne*vbar./((lam/P.Na).*(ev-1) + vbar.*ev)
Kibar = ( Ke*vbar*P.K + Nibar.*(2*lam/3).*(ev-1) )/vbar/ev/P.K
Clibar = ev*Cle;   % mM
ubar = q./(Nibar + Kibar - Clibar);   % cm

jfin = 5000;
dt = 1;  % s
u = zeros(jfin,1);
t = u; N = u; v = u; K = N; Cl = K;
u(1) = ubar;
N(1) = Nibar;
K(1) = Kibar;
v(1) = vbar;
Cl(1) = Nibar + Kibar - q/ubar;
t(1) = 0;
uN = u(1)*N(1);
uK = u(1)*K(1);
Lp = 5e-8;   % specific water permeability through passive water channels
             % Ostby supp 
             
lam = 0;   % Oubain stops pump

for j=1:jfin
    
    t(j+1) = j*dt;    
    ev = exp(v(j));
    uN = uN + dt*( P.Na*v(j)*(Ne-N(j)*ev)/(ev-1) - lam*N(j));
    uK = uK + dt*( P.K*v(j)*(Ke-K(j)*ev)/(ev-1) + 2*lam*N(j)/3);
    u(j+1) = u(j) + Lp*dt*( 2*(N(j)+K(j)) - S );
    N(j+1) = uN/u(j+1);
    K(j+1) = uK/u(j+1);
    Cl(j+1) = N(j+1) + K(j+1) - q/u(j+1);
    v(j+1) = fsolve(@(v) veqdyn(v,lam,P,Cle,Ne,N(j+1),Ke,K(j+1),u(j+1),q),v(j),options);
    
end

figure(1)
[hax,hL1,hL2] = plotyy(t,v*25.8,t,A*u);
xlabel('t  (s)','fontsize',14)
ylabel(hax(1),'V  (mV)','fontsize',14)
ylabel(hax(2),'w  (cm^3)','fontsize',14)
text(400,1.9e-8,'(A)','fontsize',20,'parent',hax(2))

figure(2)
plot(t,N,'k')
hold on
plot(t,K,'r')
plot(t,Cl,'g')
hold off
xlabel('t  (s)','fontsize',14)
ylabel('(mM)','fontsize',14)
legend('[Na]_i','[K]_i','[Cl]_i','location','best')
text(400,135,'(B)','fontsize',20)

return

function val = veq(v,lam,P,Ke,Ne,S)
a = exp(v);
b = 1-1/a;
val = a*(S/2) - Ne*(v + (2*lam/P.K/3)*b)./(v + (lam/P.Na)*b) - Ke;

function val = veqdyn(v,lam,P,Cle,Ne,Ni,Ke,Ki,u,q)
a = lam*Ni/3;
val = exp(v)*((P.Na*Ni+P.K*Ki+P.Cl*Cle)*v + a) - ...
    ((P.Na*Ne+P.K*Ke+P.Cl*(Ni+Ki-q/u))*v + a);




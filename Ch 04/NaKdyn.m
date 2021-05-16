%
% Gabbiani and Cox Mathematics for Neuroscientists
%
% NaKdyn.m
%
% Dynamic solution of the sodium/potassium pump model
%
%  usage: NaKdyn 
%
function NaKdyn
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

% now propagate

jfin = 2400;
dt = 0.25;  % s
u = zeros(jfin,1);
t = u; N = u; v = u; K = N;
u(1) = 1.05*ubar;
N(1) = 1.05*Nibar;
K(1) = 0.95*Kibar;

v(1) = fsolve(@(v) veqdyn(v,lam,P,Cle,Ne,N(1),Ke,K(1),u(1),q),vbar,options);
t(1) = 0;
uN = u(1)*N(1);
uK = u(1)*K(1);
Lp = 5e-8;   % specific water permeability through passive water channels
             % Ostby supp

for j=1:jfin
    
    t(j+1) = j*dt;    
    ev = exp(v(j));
    uN = uN + dt*( P.Na*v(j)*(Ne-N(j)*ev)/(ev-1) - lam*N(j));
    uK = uK + dt*( P.K*v(j)*(Ke-K(j)*ev)/(ev-1) + 2*lam*N(j)/3);
    u(j+1) = u(j) + Lp*dt*( 2*(N(j)+K(j)) - S );
    N(j+1) = uN/u(j+1);
    K(j+1) = uK/u(j+1);
    v(j+1) = fsolve(@(v) veqdyn(v,lam,P,Cle,Ne,N(j+1),Ke,K(j+1),u(j+1),q),v(j),options);
    
end

Naload = A*(u(1)*N(1)-ubar*Nibar)
J = lam*A*(N-Nibar);
Naremoved = trapz(t,J)

figure(1)
[hax,hL1,hL2] = plotyy(t,v*25.8,t,A*u);
xlabel('t  (s)','fontsize',14)
ylabel(hax(1),'V  (mV)','fontsize',14)
ylabel(hax(2),'w  (cm^3)','fontsize',14)
text(500,1.08e-8,'(A)','fontsize',20,'parent',hax(2))

figure(2)
hax = plotyy(t,N,t,K);
xlabel('t  (s)','fontsize',14)
ylabel(hax(1),'[Na]_i  (mM)','fontsize',14)
ylabel(hax(2),'[K]_i  (mM)','fontsize',14)
text(500,127,'(B)','fontsize',20,'parent',hax(2))

return

function val = veq(v,lam,P,Ke,Ne,S)
a = exp(v);
b = 1-1/a;
val = a*(S/2) - Ne*(v + (2*lam/P.K/3)*b)./(v + (lam/P.Na)*b) - Ke;

function val = veqdyn(v,lam,P,Cle,Ne,Ni,Ke,Ki,u,q)
a = lam*Ni/3;
val = exp(v)*((P.Na*Ni+P.K*Ki+P.Cl*Cle)*v + a) - ...
    ((P.Na*Ne+P.K*Ke+P.Cl*(Ni+Ki-q/u))*v + a);




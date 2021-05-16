%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
%  stEfork.m 
%
% solve the active fork problem via staggered Euler
%
% usage:	stEfork(fork,g,stim,pinc)
%
%		cab.rad = fiber radius (cm)
%		cab.ell = fiber length (cm)
%       cab.dx = space step (cm)
%       cab.dt = timestep (ms)
%		stim.t1 = start of current pulse (ms)
%		stim.t2 = end of current pulse (ms)
%
%       g = struct('K', 36, 'Na', 120, 'Cl', 1/15);
%
%       stim.amp = amplitude of current pulse (micro amps)
%       stim.loc = location of current pulse (cm)
%       stim.Tfin = stopping time (ms)
%
%       pinc = number of time steps between plots
%
% example:	
% fork = struct('rad',1e-4*[1 1 1],'ell',250e-4*[1 1 1],'dx',5e-4,'dt',0.05);
% stim = struct('t1',.5,'t2',1.5,'amp',4e-4,'loc',10,'Tfin',7);
% pinc = 2 
%

function [t,Vhot] = stEfork(fork,g,stim,pinc)

a = fork.rad;
ell = fork.ell; 
dx = fork.dx;
dt = fork.dt;
N = ell/dx;
A3 = 2*pi*a(3)*dx;
As = 4*pi*1e-6;
rho = A3/As;
%As = A3;
R2 = 0.3; 
Cm = 1;
lam = a/(2*R2)/dx^2;   % lambda^2/dx^2
r = a/a(3);
Hd = [2*lam(1)*ones(1,N(1)) 2*lam(2)*ones(1,N(2)) 2*lam(3)*ones(1,N(3)+1)];
Hd(1) = lam(1);
Hd(N(1)+1) = lam(2);
Hd(N(1)+N(2)+1) =  lam*r';
Hd(end) = rho*lam(3);
Hlen = length(Hd);

Hu = [-lam(1)*ones(1,N(1)-1) 0 -lam(2)*ones(1,N(2)) -lam(3)*ones(1,N(3))];
Hl = [-lam(1)*ones(1,N(1)-1) 0 -lam(2)*ones(1,N(2)-1) -r(2)*lam(2) -lam(3)*ones(1,N(3))];
Hl(end) = rho*Hl(end);

B = spdiags( [[Hl 0]' Hd' [0 Hu]'], -1:1, Hlen, Hlen);

B(N(1)+N(2)+1,N(1)) = -r(1)*lam(1);
B(N(1),N(1)+N(2)+1) = -lam(1);

dB = diag(B);
Nx = Hlen;
A = A3;
Nt = ceil(stim.Tfin/dt);
e = ones(Nx,1);

E = struct('K', -77, 'Na', 56, 'Cl', -68);
options = optimset('Jacobian','on');
V = fsolve(@(V) Iss(V,E,g,B),-70*e,options);    % initial conditions

eloc = stim.loc;
e1 = zeros(Nx,1);
e1(eloc) = 1;

if eloc == Hlen,
   Iapp = stim.amp*e1/As;
else
   Iapp = stim.amp*e1/A3;
end

Vhot = zeros(Nt,1);

Vhot(1) = V(eloc);

n = an(V)./(an(V)+bn(V)); 
m = am(V)./(am(V)+bm(V)); 
h = ah(V)./(ah(V)+bh(V)); 

x3 = 0:dx:ell(3);
x1 = ell(3):dx:ell(3)+ell(1);
x2 = ell(3):dx:ell(3)+ell(2);

t = 0;
if pinc
  figure(1)
  v1 = fliplr(V(1:N(1))');
  v2 = fliplr(V(N(1)+1:N(1)+N(2))');
  v3 = fliplr(V(N(1)+N(2)+1:end)');
  plot3(x3,t*ones(size(x3)),v3,'r')
  hold on
  plot3(x2,t*ones(size(x2)),[v3(end) v2],'r')
  plot3(x1,t*ones(size(x1)),[v3(end) v1],'k')
end

for j=2:Nt,

      t = (j-1)*dt;

      I = Iapp*(t-dt/2>stim.t1)*(t-dt/2<stim.t2);

      a = an(V);  c = (a+bn(V))/2;
      n = ( (1/dt-c).*n + a) ./ (1/dt + c); n4 = n.^4;

      a = am(V);  c = (a+bm(V))/2;
      m = ( (1/dt-c).*m + a) ./ (1/dt + c);

      a = ah(V);  c = (a+bh(V))/2;
      h = ( (1/dt-c).*h + a) ./ (1/dt + c); m3h = m.^3.*h;

      d = g.Na.*m3h + g.K.*n4 + g.Cl;

      f = g.Na.*m3h*E.Na + g.K.*n4*E.K + g.Cl.*E.Cl + I;

      B(1:Nx+1:end) = dB + d + 2*Cm/dt;         % update the diagonal

      Vmid = B\(2*Cm*V/dt + f);
      
      V = 2*Vmid - V;

      if mod(j,pinc) == 0
         v1 = fliplr(V(1:N(1))');
         v2 = fliplr(V(N(1)+1:N(1)+N(2))');
         v3 = fliplr(V(N(1)+N(2)+1:end)');
         plot3(x3,t*ones(size(x3)),v3,'r')
         plot3(x2,t*ones(size(x2)),[v3(end) v2],'r')
         plot3(x1,t*ones(size(x1)),[v3(end) v1],'k')
      end
  
      Vhot(j) = V(eloc);

end

if pinc
    xlabel('x  (cm)','fontsize',16)
    ylabel('t  (ms)','fontsize',16)
    zlabel('V  (mV)','fontsize',16)
    hold off
end

t = linspace(0,stim.Tfin,Nt)';

return

function [val, jac] = Iss(V,E,g,B)
Nx = length(V);
a = am(V); b = bm(V);
m = a./(a+b);
dm = (dam(V).*b - a.*dbm(V))./(a+b).^2;

a = ah(V); b = bh(V);
h = a./(a+b);
dh = (dah(V).*b - a.*dbh(V))./(a+b).^2;

a = an(V); b = bn(V);
n = a./(a+b);
dn = (dan(V).*b - a.*dbn(V))./(a+b).^2;

m3h = m.^3.*h;
n4 = n.^4;

val = B*V + g.Na.*m3h.*(V-E.Na) + g.K.*n4.*(V-E.K) + g.Cl.*(V-E.Cl);

dj = g.Na.*((3*dm.*m.^2.*h + m.^3.*dh).*(V-E.Na) + m3h) + ...
     g.K.*(4*dn.*n.^3.*(V-E.K) + n4) + g.Cl;
jac = B + spdiags(dj,0,Nx,Nx);

function val = an(v)
val = .01*(10-(v+71))./(exp(1-(v+71)/10)-1);

function val = dan(v)
tmp = exp(-(61+v)/10);
val = -( tmp.*(71+v) - 10 )./(tmp-1).^2/1000;

function val = bn(v)
val = .125*exp(-(v+71)/80);

function val = dbn(v)
val = -exp(-(v+71)/80)/640;

function val = am(v)
val = .1*(25-(v+71))./(exp(2.5-(v+71)/10)-1);

function val = dam(v)
tmp = exp(-(46+v)/10);
val = -( tmp.*(56+v) - 10 )./(tmp-1).^2/100;

function val = bm(v)
val = 4*exp(-(v+71)/18);

function val = dbm(v)
val = -(2/9)*exp(-(v+71)/18);

function val = ah(v)
val = 0.07*exp(-(v+71)/20);

function val = dah(v)
val = -(7/2000)*exp(-(v+71)/20);

function val = bh(v)
val = 1./(exp(3-(v+71)/10)+1);

function val = dbh(v)
tmp = exp(-(v+41)/10);
val = tmp./(tmp+1).^2/10;



%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
%  forkrest.m  
%
% find rest pot of active fork via fsolve with jacobian
%
% usage:  forkrest(fork,g)
%
%		cab.rad = fiber radius (cm)
%		cab.ell = fiber length (cm)
%       cab.dx = space step (cm)
%       cab.dt = timestep (ms)
%
%       g = struct('K', 36, 'Na', 120, 'Cl', 1/15);
%
% example:	
% fork = struct('rad',1e-4*[1 1 1],'ell',250e-4*[1 1 1],'dx',5e-4,'dt',0.05);
%    g = struct('K', 36, 'Na', 120, 'Cl', 1/15);
%

function Vr = forkrest(fork,g)

a = fork.rad;
ell = fork.ell; 
dx = fork.dx;
N = ell/dx;
A3 = 2*pi*a(3)*dx;
As = 4*pi*1e-6;
%As = A3;
rho = A3/As;
R2 = 0.3; 
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

E = struct('K', -77, 'Na', 56, 'Cl', -68);
e = ones(Hlen,1);

options = optimset('Jacobian','on');
Vr = fsolve(@(V) Iss(V,E,g,B),-70*e,options);    % initial conditions

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

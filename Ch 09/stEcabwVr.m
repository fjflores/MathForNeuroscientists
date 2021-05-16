
%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
%  stEcabwVr.m   
%
% solve the active cable problem via staggered Euler, with Vr (rest)
% precomputed. This saves time when conducting parameter sweeps across
% various inputs, e.g., by changing input amplitude, duration or location
%
%  Return time and voltage at site of stimulation
%
% usage:	[t, Vhot] = stEcabwVr(cab,g,stim,pinc,Vr)
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
%       pinc = number of time steps between plots
%
% example:	cab = struct('rad',1e-4,'ell',1e-1,'dx',1e-4,'dt',0.02)
%           stim = struct('t1',1,'t2',2,'amp',4e-4,'loc',0.03,'Tfin',8)
%
%  or   x = cab.dx/2:cab.dx:cab.ell-cab.dx/2; hot = ((x>0.005).*(x<0.01))';
%       g = struct('K', 40-20*hot, 'Na', 44+560*hot, 'Cl', 1/15);
%           pinc = 1 
%

function [t,Vhot] = stEcabwVr(cab,g,stim,pinc,Vr)

Cm = 1;		% micro F / cm^2
R2 = 0.3; %0.034;		% k Ohm cm
dx = cab.dx;
dt = cab.dt;
Nx = cab.ell/dx;		% patch length
A = 2*pi*cab.rad*dx;		% patch surface area
x = dx/2:dx:cab.ell-dx/2;	% vector of patch midpoints

Nt = ceil(stim.Tfin/dt);

E = struct('K', -77, 'Na', 56, 'Cl', -68);

e1 = zeros(Nx,1);
eloc = round(Nx*stim.loc/cab.ell);
e1(eloc) = 1;
Iapp = stim.amp*e1/A;

Vhot = zeros(Nt,1);

e = ones(Nx,1);
B = spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx/dx;
B(1,1) = 1/dx/dx;
B(Nx,Nx) = 1/dx/dx;
B = (cab.rad/2/R2)*B;
dB = diag(B);

V = Vr;    % initial conditions
Vhot(1) = V(eloc);

n = an(V)./(an(V)+bn(V)); 
m = am(V)./(am(V)+bm(V)); 
h = ah(V)./(ah(V)+bh(V)); 

t = 0;
if pinc > 0
   plot3(x,t*ones(Nx,1),V,'k')
   hold on
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
         plot3(x,t*ones(Nx,1),V,'k')
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

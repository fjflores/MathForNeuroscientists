%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
%  stEforksyn.m   
%
% solve the active fork problem via staggered Euler
%
% usage:	stEforksyn(cab,stim,pinc)
%
%		stim.t1 = start of current pulse (ms)
%       stim.amp = amplitude of current pulse (micro amps)
%       stim.loc = location of current pulse (cm)
%       stim.Tfin = stopping time (ms)
%       pinc = number of time steps between plots
%
% e.g.:    stim = struct('t1',1,'tau',1,'Gsyn',5e-7,'Vsyn',0,'loc',100,'Tfin',20)
%           pinc = 5 
%

function [t,Vrec] = stEforksyn(stim,g,pinc)

a = 1e-4*[1 1 1];
ell = [2.5 2.5 2.5]/100;
dx = .0001;
dt = 0.05;
N = ell/dx;
A3 = 2*pi*a(3)*dx;
As = 4*pi*1e-6;
rho = A3/As; 
R2 = 0.3; %0.034;
gL = 1/15; %0.3;
Cm = 1;
tau = Cm/gL;
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

H = spdiags( [[Hl 0]' Hd' [0 Hu]'], -1:1, Hlen, Hlen);

H(N(1)+N(2)+1,N(1)) = -r(1)*lam(1);
H(N(1),N(1)+N(2)+1) = -lam(1);

%full(H(Hlen-3:Hlen,Hlen-3:Hlen))
%pause

dH = diag(H);

eloc = round(Hlen*1e-4*stim.loc/sum(ell));
bloc = (eloc-1)*(Hlen+1) + 1;
dHe = H(bloc);

x3 = 0:dx:ell(3);
x1 = ell(3):dx:ell(3)+ell(1);
x2 = ell(3):dx:ell(3)+ell(2);

Nt = ceil(stim.Tfin/dt)+1;

E = struct('K', -77, 'Na', 56, 'Cl', -68);
e = ones(Hlen,1);

options = optimset('Jacobian','on');
V = fsolve(@(V) Iss(V,E,g,H),-77*e,options);     % initial conditions
%figure(1)
%plot(V,'k')
%box off
%xlabel('x  (cm)','fontsize',14)
%ylabel('V_{rest}  (mV)','fontsize',14)

n = an(V)./(an(V)+bn(V)); 
m = am(V)./(am(V)+bm(V)); 
h = ah(V)./(ah(V)+bh(V)); 

t = 0;
x = linspace(0,ell(3)+max(ell(1:2)),Hlen);
if pinc > 0
   V1 = fliplr(V(1:N(1))');
   V2 = fliplr(V(N(1)+1:N(1)+N(2))');
   V3 = fliplr(V(N(1)+N(2)+1:end)');
   plot3(x3,t*ones(size(x3)),V3,'r')
   hold on
   plot3(x2,t*ones(size(x2)),[V3(end) V2],'r')
   plot3(x1,t*ones(size(x1)),[V3(end) V1],'k')
end

gsyn = 0*e;
Iapp = stim.Gsyn/(2*pi*dx*a(1));

if eloc == Hlen
   Iapp = Iapp*2*pi*dx*a(1)/As;
end
%pause

nloc = length(stim.loc);

Vrec = zeros(nloc+1,Nt); 
Vrec(1:nloc,1) = V(stim.loc); Vrec(nloc+1,1) = V(end);

for j=2:Nt,

      t = (j-1)*dt;

      t2 = t-dt/2;
      gsyn(eloc) = Iapp.*((t2-stim.t1)./stim.tau).*exp(1-(t2-stim.t1)./stim.tau).*(t2>stim.t1);

      a = an(V);  c = (a+bn(V))/2;
      n = ( (1/dt-c).*n + a) ./ (1/dt + c); n4 = n.^4;

      a = am(V);  c = (a+bm(V))/2;
      m = ( (1/dt-c).*m + a) ./ (1/dt + c);

      a = ah(V);  c = (a+bh(V))/2;
      h = ( (1/dt-c).*h + a) ./ (1/dt + c); m3h = m.^3.*h;

      d = g.Na.*m3h + g.K.*n4 + g.Cl + gsyn;

      f = g.Na.*m3h*E.Na + g.K.*n4*E.K + g.Cl*E.Cl;
      f(eloc) = f(eloc) + gsyn(eloc).*stim.Vsyn';

      H(1:Hlen+1:end) = dH + d + 2*Cm/dt;         % update the diagonal

      Vmid = H\(2*Cm*V/dt + f);
      
      V = 2*Vmid - V;
    
      Vrec(1:nloc,j) = V(stim.loc); Vrec(nloc+1,j) = V(end);

      if mod(j,pinc) == 0
         V1 = fliplr(V(1:N(1))');
         V2 = fliplr(V(N(1)+1:N(1)+N(2))');
         V3 = fliplr(V(N(1)+N(2)+1:end)');
         plot3(x3,t*ones(size(x3)),V3,'r')
         plot3(x2,t*ones(size(x2)),[V3(end) V2],'r')
         plot3(x1,t*ones(size(x1)),[V3(end) V1],'k')
      end

end

t = linspace(0,stim.Tfin,Nt)';

if pinc
xlabel('x  (cm)','fontsize',14)
ylabel('t  (ms)','fontsize',14)
zlabel('V  (mV)','fontsize',14)
hold off
end

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


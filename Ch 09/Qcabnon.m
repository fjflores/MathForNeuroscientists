%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
%  Qcabnon.m   
%
% solve the nonuniform quasi-active cable problem 
%
% usage:	Qcabnon(cab)
%
%		cab.rad = fiber radius (cm)
%		cab.ell = fiber length (cm)
%               cab.dx = space step (cm)
%
% example:	cab = struct('rad',1e-4,'ell',0.1,'dx',1e-3)
%

function Qcabnon(cab)

Cm = 1;		% micro F / cm^2
R2 = 0.3; %0.034;		% k Ohm cm
dx = cab.dx;
Nx = cab.ell/dx;		% patch length
A = 2*pi*cab.rad*dx;		% patch surface area
x = dx/2:dx:cab.ell-dx/2;	% vector of patch midpoints

E = struct('K', -77, 'Na', 56, 'Cl', -68);
g = struct('K', 36, 'Na', 120, 'Cl', 1/15);
hot = ((x>0.005).*(x<0.01))';
g = struct('K', 40-20*hot, 'Na', 44+560*hot, 'Cl', 1/15);

e = ones(Nx,1);
S = spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx/dx;
S(1,1) = 1/dx/dx;
S(Nx,Nx) = 1/dx/dx;
S = (cab.rad/2/R2)*S;
Z = 0*speye(Nx);

options = optimset('Jacobian','on');
V = fsolve(@(V) Iss(V,E,g,S),-77*e,options);    % initial conditions

a = am(V); b = bm(V);
taum = 1./(a+b);
m = a.*taum;
dm = (dam(V).*b - a.*dbm(V)).*taum.^2;

a = ah(V); b = bh(V);
tauh = 1./(a+b);
h = a.*tauh;
dh = (dah(V).*b - a.*dbh(V)).*tauh.^2;

a = an(V); b = bn(V);
taun = 1./(a+b);
n = a.*taun;
dn = (dan(V).*b - a.*dbn(V)).*taun.^2;

Q11 = -(S + spdiags(g.Cl+g.Na.*m.^3.*h+g.K.*n.^4,0,Nx,Nx))/Cm;
Q12 = -spdiags(3*g.Na.*m.^2.*h.*(V-E.Na)/Cm,0,Nx,Nx);
Q13 = -spdiags(g.Na.*m.^3.*(V-E.Na)/Cm,0,Nx,Nx);
Q14 = -spdiags(4*g.K.*n.^3.*(V-E.K)/Cm,0,Nx,Nx);

Q21 = spdiags(dm./taum,0,Nx,Nx);
Q22 = spdiags(-1./taum,0,Nx,Nx);

Q31 = spdiags(dh./tauh,0,Nx,Nx);
Q33 = spdiags(-1./tauh,0,Nx,Nx);

Q41 = spdiags(dn./taun,0,Nx,Nx);
Q44 = spdiags(-1./taun,0,Nx,Nx);

B = [Q11 Q12 Q13 Q14;
     Q21 Q22  Z   Z;
     Q31  Z  Q33  Z;
     Q41  Z   Z  Q44];

[V,D] = eig(full(B));
d = diag(D);
figure(1)
plot(real(d),imag(d),'ko','markersize',6)
xlim([-0.35 -.05])
xlabel('real','fontsize',14)
ylabel('imaginary','fontsize',14')
text(-.4,.165,'(A)','fontsize',18)
%print -depsc QcabnonA

figure(2)
ind = find(imag(d)>1e-2);

plot(x,real(V(1:Nx,ind(1))),'k')
hold on
plot(x,imag(V(1:Nx,ind(1))),'r')
plot(x,real(V(1:Nx,ind(2))),'k--')
plot(x,imag(V(1:Nx,ind(2))),'r--')

box off
legend('q_0 (real)','q_0 (imag)','q_1 (real)','q_1 (imag)','location','best')
legend boxoff
axis tight
xlabel('x  (cm)','fontsize',14)
text(-0.01,.22,'(B)','fontsize',18)
%print -depsc QcabnonB

hold off

q0non = real(V(1:Nx,ind(1)));

save vec0 q0non

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

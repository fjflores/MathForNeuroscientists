%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
%  stEcabspine.m 
%
% solve the active nonuniform cable problem via staggered Euler
%
% usage:	stEcabspine(cab,stim,pinc)
%
%		cab.rad = fiber radius (cm)
%		cab.ell = fiber length (cm)
%       cab.dx = space step (cm)
%       cab.dt = timestep (ms)
%		stim.t1 = start of current pulse (ms)
%		stim.t2 = end of current pulse (ms)
%       stim.amp = amplitude of current pulse (micro amps)
%       stim.loc = location of current pulse (cm)
%       stim.Tfin = stopping time (ms)
%       pinc = number of time steps between plots
%
% example:	cab = struct('rad',1e-4,'ell',0.1,'dx',1e-3,'dt',0.01)
%           stim = struct('t1',1,'t2',2,'dose',1,'loc',0.04,'Tfin',12)
%           pinc = 10 
%

function stEcabspine(cab,stim,pinc)

Cm = 1;		% micro F / cm^2
R2 = 0.3; %0.034;		% k Ohm cm
dx = cab.dx;
dt = cab.dt;
Nx = cab.ell/dx;		% patch length
A = 2*pi*cab.rad*dx;		% patch surface area
x = dx/2:dx:cab.ell-dx/2;	% vector of patch midpoints

Nt = ceil(stim.Tfin/dt)+1;

E = struct('K', -77, 'Na', 56, 'Cl', -68);
hot = ((x>0.005).*(x<0.01))';
g = struct('K', 40-20*hot, 'Na', 30+570*hot, 'Cl', 1/15);

ampa = struct('kp',1.1,'km',0.19,'gbar',2e2);
nmda = struct('kp',0.072,'km',0.0066,'gbar',1e2);

e1 = zeros(Nx,1);
eloc = round(Nx*stim.loc/cab.ell);
e1(eloc) = 1;

e = ones(Nx,1);
B = spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx/dx;
B(1,1) = 1/dx/dx;
B(Nx,Nx) = 1/dx/dx;
B = (cab.rad/2/R2)*B;
dB = diag(B);

ellss = 1e-4; ass = 1e-5; Ash = 1e-8;
Rss = ellss*R2/pi/ass^2;
gam1 = 1/Rss/Ash;
gam2 = 1/Rss/2/pi/cab.rad;
Vsyn = 20;

co =  gam2*g.Cl/(g.Cl+gam1)/dx;
options = optimset('Jacobian','on');
V = fsolve(@(V) Issp(V,E,g,B,co,eloc),-70*e,options);    % initial conditions
figure(1)
plot(x,V,'k')
box off
xlabel('x  (cm)','fontsize',14)
ylabel('V_{rest}  (mV)','fontsize',14)

W = zeros(Nt,1); Vsp = W;
W(1) = (g.Cl*E.Cl + gam1*V(eloc))/(g.Cl+gam1);
Vsp(1) = V(eloc);

figure(2)

n = an(V)./(an(V)+bn(V)); 
m = am(V)./(am(V)+bm(V)); 
h = ah(V)./(ah(V)+bh(V)); 

rampa = 0;
rnmda = 0;

t = 0;
if pinc > 0
   plot3(x,t*ones(Nx,1),V,'k')
   hold on
end

Iampa = zeros(Nt,1);
Inmda = zeros(Nt,1);

for j=2:Nt,

      t = (j-1)*dt;

      T = stim.dose*(t>stim.t1)*(t<stim.t2);
      rampa = ( rampa + dt*ampa.kp*T ) / (1 + dt*(ampa.kp*T + ampa.km) );
      rnmda = ( rnmda + dt*nmda.kp*T ) / (1 + dt*(nmda.kp*T + nmda.km) );
      gampa = ampa.gbar*rampa;
      gnmda = nmda.gbar*rnmda;

      tmp = gampa + gnmda*W(j-1)/(1+2*exp(-0.062*W(j-1))/3.57);

      Wnum = (Cm/dt)*W(j-1) + g.Cl*E.Cl + tmp*Vsyn + gam1*V(eloc);
      Wden = (Cm/dt) + g.Cl + tmp + gam1;
      W(j) = Wnum/Wden;

      a = an(V);  
      n = (n + dt*a) ./ (1 + dt*(a+bn(V)) ); n4 = n.^4;

      a = am(V);
      m = (m + dt*a) ./ (1 + dt*(a+bm(V)) );

      a = ah(V);  
      h = (h + dt*a) ./ (1 + dt*(a+bh(V)) ); m3h = m.^3.*h;

      d = g.Na.*m3h + g.K.*n4 + g.Cl;
      d(eloc) = d(eloc)  + gam2/dx;

      f = g.Na.*m3h*E.Na + g.K.*n4*E.K + g.Cl*E.Cl;
      f(eloc) = f(eloc)  + gam2*W(j)/dx;

      B(1:Nx+1:end) = dB + d + Cm/dt;         % update the diagonal

      V = B\(Cm*V/dt + f);

      Vsp(j) = V(eloc);
      
      Iampa(j) = gampa*(W(j)-Vsyn);
      Inmda(j) = gnmda*(W(j)-Vsyn)/(1+2*exp(-0.062*W(j))/3.57);

      if mod(j,pinc) == 0
         plot3(x,t*ones(Nx,1),V,'k')
      end

end

xlabel('x  (cm)','fontsize',16)
ylabel('t  (ms)','fontsize',16)
zlabel('V  (mV)','fontsize',16)
ylim([0 16])
zlim([-80 60])

hold off

tim = linspace(0,stim.Tfin,Nt)';

figure(3)
plot(tim,Vsp,'k')
hold on
plot(tim,W,'r')
box off
legend('V(x_s,t)','W(t)','location','best')
xlabel('t  (ms)','fontsize',14)
ylabel('mV','fontsize',14)
hold off

figure(4)
plot(tim,1e6*Iampa*Ash,'k')
hold on
plot(tim,100*1e6*Inmda*Ash,'r')
box off
legend('I_{ampa}','I_{nmda}*100','location','best')
xlabel('t  (ms)','fontsize',14)
ylabel('pA','fontsize',14)
hold off

return

function [val, jac] = Issp(V,E,g,B,co,eloc)
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
val(eloc) = val(eloc) - co*(E.Cl-V(eloc));

dj = g.Na.*((3*dm.*m.^2.*h + m.^3.*dh).*(V-E.Na) + m3h) + ...
     g.K.*(4*dn.*n.^3.*(V-E.K) + n4) + g.Cl;
jac = B + spdiags(dj,0,Nx,Nx);
jac(eloc,eloc) = jac(eloc,eloc) + co;

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


%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  stEcabextrapop.m  
%
% solve for the extracellular potential in the neighborhood of 
% three active cables
%
% usage:	stEcabextrapop
%

function stEcabextrapop

cab = struct('rad',1e-4,'ell',1e-1,'dx',1e-4,'dt',0.01);
g = struct('K', 36, 'Na', 120, 'Cl', 0.3);
stim = struct('amp',2e-4,'loc',cab.ell/2,'Tfin',8);

close all

Cm = 1;		% micro F / cm^2
R2 = 0.3; %0.034;		% k Ohm cm
dx = cab.dx;
dt = cab.dt;
Nx = cab.ell/dx;		% patch length
A = 2*pi*cab.rad*dx;		% patch surface area
x = dx/2:dx:cab.ell-dx/2;	% vector of patch midpoints

Nt = ceil(stim.Tfin/dt);

E = struct('K', -77, 'Na', 56, 'Cl', -68);

eloc = round(Nx*stim.loc/cab.ell);
Iapp = stim.amp/A;

e = ones(Nx,1);
B = spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx/dx;
B(1,1) = 1/dx/dx;
B(Nx,Nx) = 1/dx/dx;
B = (cab.rad/2/R2)*B;
dB = diag(B);

options = optimset('Jacobian','on');
V = fsolve(@(V) Iss(V,E,g,B),-70*e,options);    % initial conditions

n = an(V)./(an(V)+bn(V)); 
m = am(V)./(am(V)+bm(V)); 
h = ah(V)./(ah(V)+bh(V)); 

t = 0;
I = zeros(Nx,1);
Imem = zeros(Nx,Nt);

for j=2:Nt,

      t = (j-1)*dt;

      I(eloc) = Iapp.*t.*exp(1-t);

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
  
      % collect the membrane currents

    INa = g.Na*m3h.*(V-E.Na);
    IK = g.K*n4.*(V-E.K);
    ICl = g.Cl*(V-E.Cl);
    ICap = 2*Cm*(V-Vmid)/dt;
    Imem(:,j) = A*(INa + IK + ICl - I + ICap);
    
end

% compute contributions at each electrode

ncab = 3;
t1 = [1 6 12];
tfin = 20;
xloc = [-.005 0.0025 0.005];
yloc = [-.005 0.005 0];
ycenters = ((1:Nx)'-1/2)*cab.ell/Nx;
phi = zeros(9,ceil(tfin/dt));
sig = 0.003; % S/cm 

for c = 1:ncab
    xcenter = xloc(c); %0.005;
    jinc = ceil(t1(c)/dt);
    for j=1:Nt,
        for iy = 1:9
            yrec = 0.01*iy;
            dscale = 4*pi*sig*sqrt(xcenter^2+(yrec-(ycenters+yloc(c))).^2);
            phi(iy,j+jinc) = phi(iy,j+jinc) + sum(Imem(:,j)./dscale);
        end
    end
end

% draw thee electrode/cables configuration

figure(1)

plot(3*[-1 1 1 -1 -1]-50,[0 0 1e3 1e3 0]-50,'k')   % cell 1
hold on
text(-52,-100,'1','fontsize',14)
plot(3*[-1 1 1 -1 -1]+25,[0 0 1e3 1e3 0]+50,'k')   % cell 2
text(23,0,'2','fontsize',14)
plot(3*[-1 1 1 -1 -1]+50,[0 0 1e3 1e3 0],'k')   % cell 3
text(48,-50,'3','fontsize',14)
plot([-10 -10 0 10 10],[1.1e3 0 -100 0 1.1e3],'r')   % shank
for j=1:9,
    fill(10*[-1 1 1 -1],j*100 + 5*[-1 -1 1 1],'r')
end
text(-35,1100,'(A)','fontsize',20)
hold off
xlabel('x  (\mum)','fontsize',16)
ylabel('y  (\mum)','fontsize',16)
box off

% graph extracellular potentials

figure(2)
tim = dt:dt:tfin;
for iy=1:9,
    plot(tim,phi(iy,:)+12*(iy-1),'k','linewidth',1.5)
    hold on
end
plot([1 3],92*[1 1],'k')
text(1.25,88.25,'2 ms','fontsize',10)
plot([6 6],[87 93],'k')
text(6.2,90,'6 \muV','fontsize',10)
text(1,105,'(B)','fontsize',20)
hold off
axis tight
axis off
    
    
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

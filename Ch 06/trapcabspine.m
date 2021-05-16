%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  trapcabspine.m
%
% solve the passive fiber with alpha func synaptic input via the trapezoid rule
%
% usage:        trapcabspine(cab,stim,pinc)
%
%               cab.rad = cable radius (cm)
%               cab.ell = cable length (cm)
%    	        cab.dx = space step (cm)
%               cab.dt = timestep (ms)
%               stim.t1 = stim start time (ms)
%               stim.tau = time constant of alpha function (ms)
%               stim.Gsyn = amplitude of synaptic conductance (mS)
%               stim.loc = location of current pulse (cm)
%               stim.Tfin = stopping time (ms)
%               pinc = number of time steps between plots
%
% example:      cab = struct('rad',1e-4,'ell',1e-1,'dx',1e-3,'dt',0.05)
%               stim = struct('t1',1,'tau',1,'Gsyn',100,'loc',0.05,'Tfin',10)
%               pinc = 1 
%
% or, for multiple stimuli
%
%               stim = struct('t1',[1 3],'tau',[1 1]/2,'Gsyn',100*[1 1],'loc',[0.06 0.04],'Tfin',10)
%

function trapcabsynspine(cab,stim,pinc)

Cm = 1;		% micro F / cm^2
g_L = 1/15; %0.3;     		% mS / cm^2
R_2 = .3; %0.034;		% k Ohm cm
dx = cab.dx;
dt = cab.dt;
Nx = cab.ell/dx;                % patch length
A = 2*pi*cab.rad*dx;            % patch surface area
x = dx/2:dx:cab.ell-dx/2;       % vector of patch midpoints

Nt = ceil(stim.Tfin/dt);

v0 = zeros(Nx,1);		% initial conditions

e = ones(Nx,1);
S = spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx/dx;
S(1,1) = 1/dx/dx;
S(Nx,Nx) = 1/dx/dx;

tau = Cm/g_L;
lambda = sqrt(cab.rad/(2*R_2*g_L))
A = 2*pi*cab.rad*dx;

e1 = zeros(Nx,1);
eloc = round(Nx*stim.loc/cab.ell);
vsyn = 70;
Iapp = stim.Gsyn/g_L;

B = (speye(Nx)+lambda^2*S)/tau;
B = (2/dt)*speye(Nx) + B;
bloc = (eloc-1)*(Nx+1)+1;
dBe = B(bloc);

ellss = 1e-4; ass = 1e-5; Ash = 1e-8;
Rss = ellss*R_2/pi/ass^2; Gsh = g_L*Ash;
g1 = 1/Rss/Gsh;
g2 = 1/Rss/A/Cm;

t = 0;
subplot(1,2,1)
figure(1)
plot3(x,t*e,v0)
hold on

w = zeros(Nt,length(eloc));
vhot = zeros(length(eloc),Nt);

c0 = Iapp.*(t./stim.tau).*exp(1-t./stim.tau).*(t>stim.t1);   % j=1
xi0 = g2 - g1*g2./(2*tau/dt + 1 + g1 + c0);
w(1,:) = 0;

t = dt;
c1 = Iapp.*(t./stim.tau).*exp(1-t./stim.tau).*(t>stim.t1);  % j=2
f = g2*( (4*tau/dt + c1 - c0)*w(1) + vsyn*(c0+c1) )./( 2*tau/dt + 1 + c1 + g1); 
xi1 = g2 - g1*g2./(2*tau/dt + 1 + g1 + c1)
dBe
size(dBe)

r = zeros(Nx,1);
r(eloc) = f;
B(bloc) = dBe + xi1;
v1 = B\r;
w(2,:) = ( (2*tau/dt - 1 - c0 - g1).*w(1,:) + vsyn*(c0+c1) + (v0(eloc) + v1(eloc))'*g1 )./ ( 2*tau/dt + 1 + c1 + g1);

for j=3:Nt,

    t = (j-1)*dt;

    c0 = c1; xi0 = xi1; v0 = v1;

    c1 = Iapp.*((t-stim.t1)./stim.tau).*exp(1-(t-stim.t1)./stim.tau).*(t>stim.t1);  % j=2
    f = g2*( (4*tau/dt + c1 - c0).*w(j-1,:) + vsyn*(c0+c1) )./( 2*tau/dt + 1 + c1 + g1);
    xi1 = g2 - g1*g2./(2*tau/dt + 1 + g1 + c1);

    r = (4/dt)*v0 - r;
    r(eloc) = r(eloc) + (xi0-xi1)'.*v0(eloc) + f';
    B(bloc) = dBe + xi1;
    v1 = B\r;
    w(j,:) = ( (2*tau/dt - 1 - c0 - g1).*w(j-1,:) + vsyn*(c0+c1) + (v0(eloc) + v1(eloc))'*g1 ) ./ ( 2*tau/dt + 1 + c1 + g1);

    vhot(:,j) = v1(eloc);

    if mod(j,pinc) == 0
        plot3(x,t*e,v1)
    end

end

xlabel('x (cm)','fontsize',14)
ylabel('t (ms)','fontsize',14)
zlabel('v (mV)','fontsize',14)

hold off

subplot(1,2,2)
t = linspace(0,stim.Tfin,Nt);
plot(t',w,'--')
hold on
plot(t,vhot)
hold off
box off

xlabel('t (ms)','fontsize',14)
ylabel('(mV)','fontsize',14)


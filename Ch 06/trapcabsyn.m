%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  trapcabsyn.m
%
% solve the passive fiber with alpha func synaptic input via the trapezoid rule
%
% usage:        trapcabsyn(cab,stim,pinc)
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
% example:      cab = struct('rad',1e-4,'ell',.1,'dx',1e-3,'dt',0.05)
%               stim = struct('t1',1,'tau',1,'Gsyn',1e-4,'loc',0.05,'Tfin',10)
%               pinc = 4 
%
% or, for multiple stimuli
%
%   stim = struct('t1',[1 3],'tau',[1 1]/2,'Gsyn',1e-6*[1 1],'loc',[0.06 0.04],'Tfin',10)
%

function trapcabsyn(cab,stim,pinc)

Cm = 1;		% micro F / cm^2
g_L = 1/15; %0.3;     		% mS / cm^2
R_2 = 0.3; % 0.034;		% k Ohm cm
dx = cab.dx;
dt = cab.dt;
Nx = cab.ell/dx;                % patch length
A = 2*pi*cab.rad*dx;            % patch surface area
x = dx/2:dx:cab.ell-dx/2;       % vector of patch midpoints

Nt = ceil(stim.Tfin/dt);

v = zeros(Nx,1);		% initial conditions

e = ones(Nx,1);
S = spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx/dx;
S(1,1) = 1/dx/dx;
S(Nx,Nx) = 1/dx/dx;

tau = Cm/g_L;
lambda = sqrt(cab.rad/(2*R_2*g_L));
A = 2*pi*cab.rad*dx;

e1 = zeros(Nx,1);
eloc = round(Nx*stim.loc/cab.ell);
vsyn = 70;
Iapp = stim.Gsyn*(dt/2)/A/Cm;

B = (dt/2)*(speye(Nx)+lambda^2*S)/tau;
B = speye(Nx) + B;
bloc = (eloc-1)*(Nx+1)+1;
dBe = B(bloc);

t = 0;
figure(1)
plot3(x,t*e,v,'k')
hold on

c0 = Iapp.*(t./stim.tau).*exp(1-t./stim.tau).*(t>stim.t1);
t = dt;
c1 = Iapp.*(t./stim.tau).*exp(1-t./stim.tau).*(t>stim.t1);

r = zeros(Nx,1);
r(eloc) = vsyn*(c0 + c1)';

vhot = zeros(length(eloc),Nt);

for j=2:Nt,
    
    B(bloc) = dBe + c1;

    v = B\r; 

    vhot(:,j) = v(eloc);

    if mod(j,pinc) == 0
        plot3(x,t*e,v,'k')
    end

    t = j*dt;
    c0 = c1;
    c1 = Iapp.*((t-stim.t1)./stim.tau).*exp(1-(t-stim.t1)./stim.tau).*(t>stim.t1);

    r = 2*v - r;
    r(eloc) = r(eloc) + vsyn*(c0 + c1)';

end

xlabel('x (cm)','fontsize',14)
ylabel('t (ms)','fontsize',14)
zlabel('v (mV)','fontsize',14)

hold off

figure(2)
t = linspace(0,stim.Tfin,Nt);
plot(t,vhot(1,:),'k') 
hold on
plot(t,vhot(2,:),'r')
hold off
box off
xlabel('t (ms)','fontsize',14)
ylabel('v (mV)','fontsize',14)


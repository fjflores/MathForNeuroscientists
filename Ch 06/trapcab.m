%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  trapcab.m
%
% solve the passive cable with current injection via the trapezoid rule
%
% usage:        trapcab(cab,stim,pinc)
%
%               cab.rad = cable radius (cm)
%               cab.ell = cable length (cm)
%    	        cab.dx = space step (cm)
%               cab.dt = timestep (ms)
%               stim.t1 = start of current pulse (ms)
%               stim.t2 = end of current pulse (ms)
%               stim.amp = amplitude of current pulse (micro amps)
%               stim.loc = location of current pulse (cm)
%               stim.Tfin = stopping time (ms)
%               pinc = number of time steps between plots
%
% example:      cab = struct('rad',1e-4,'ell',.1,'dx',1e-3,'dt',0.05)
%               stim = struct('t1',1,'t2',2,'amp',1e-4,'loc',0.06,'Tfin',10)
%               pinc = 4 
%
% or, for multiple stimuli
%
%     stim = struct('t1',[1 2],'t2',[2 3],'amp',1e-4*[1 1],'loc',[0.06 0.04],'Tfin',10)
%

function trapcab(cab,stim,pinc)

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

B = (speye(Nx)+lambda^2*S)/tau;

[L,U] = lu(speye(Nx)+B*dt/2);

e1 = zeros(Nx,1);
eloc = round(Nx*stim.loc/cab.ell);
Iapp = stim.amp*(dt/2)/A/Cm;

vhot = zeros(length(eloc),Nt);

t = 0;
subplot(1,2,1)
plot3(x,t*e,v)
hold on

f0 = Iapp.*(t>stim.t1).*(t<stim.t2);
t = dt;
f1 = Iapp.*(t>stim.t1).*(t<stim.t2);

r = zeros(Nx,1);
r(eloc) = (f0 + f1)';

for j=2:Nt,
    
    v = U \ ( L \ r );

    vhot(:,j) = v(eloc);

    if mod(j,pinc) == 0
        plot3(x,t*e,v)
    end

    t = j*dt;
    f0 = f1;
    f1 = Iapp.*(t>stim.t1).*(t<stim.t2);

    r = 2*v - r;
    r(eloc) = r(eloc) + (f0 + f1)';

end

xlabel('x (cm)','fontsize',14)
ylabel('t (ms)','fontsize',14)
zlabel('v (mV)','fontsize',14)

hold off

subplot(1,2,2)
t = linspace(0,stim.Tfin,Nt);
plot(t,vhot(1,:),'k')
hold on
plot(t,vhot(2,:),'r')
hold off
box off
xlabel('t (ms)','fontsize',14)
ylabel('v (mV)','fontsize',14)

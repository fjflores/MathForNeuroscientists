%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% eigcab.m
%
% compute and display the semi-analytical solution of the passive
% cable equation. The "semi" stem from use of eig
%
% usage:        eigcab(cab,stim,pinc)
%
%               cab.rad = cable radius (cm)
%               cab.ell = cable length (cm)
%               cab.dx = space step (cm)
%               cab.dt = timestep (ms)
%               stim.t1 = start of current pulse (ms)
%               stim.t2 = end of current pulse (ms)
%               stim.amp = amplitude of current pulse (micro amps)
%               stim.loc = location of current pulse (cm)
%               stim.Tfin = stopping time (ms)
%               pinc = number of time steps between plots
%
% example:      cab = struct('rad',1e-4,'ell',0.1,'dx',1e-3,'dt',0.05)
%               stim = struct('t1',1,'t2',2,'amp',1e-4,'loc',0.06,'Tfin',10)
%               pinc = 4 
%
function eigcab(cab,stim,pinc)

C_m = 1;		% micro F / cm^2
G_L = 1/15; % 0.3;     		% mS / cm^2
R_2 = 0.3; %0.034;		% k Ohm cm

dx = cab.dx;
Nx = cab.ell/dx;
x = dx/2:dx:cab.ell-dx/2;	% vector of patch midpoints

S = (2*eye(Nx)-diag(ones(Nx-1,1),1)-diag(ones(Nx-1,1),-1))/dx/dx;
S(1,1) = 1/dx/dx;
S(Nx,Nx) = 1/dx/dx;

tau = C_m/G_L;
lambda = sqrt(cab.rad/(2*R_2*G_L));

[Q,th] = eig(S);

th = diag(th)';

z = (1+lambda^2*th)/tau;

Nt = ceil(stim.Tfin/cab.dt);

t = 0;
e = ones(Nx,1);
plot3(x,t*e,0*e)
hold on

sc = stim.amp*Nx/(2*pi*cab.rad*cab.ell);
eloc = round(Nx*stim.loc/cab.ell);

for j=1:Nt,

    t = j*cab.dt;

    et = exp(min(stim.t2-t,0)*z) - exp(min(stim.t1-t,0)*z);
    P = Q*diag(Q(eloc,:).*et./z);
    v = sc*sum(P,2);

    if mod(j,pinc) == 0
        plot3(x,t*e,v)
    end

end

hold off
box off
xlabel('x  (cm)','fontsize',14)
ylabel('t  (ms)','fontsize',14)
zlabel('v  (mV)','fontsize',14)

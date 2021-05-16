%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  cnpfib.m
%
% solve the passive fiber problem with distributed current
%  injection via trap and be
%
% usage:        cnpfib(a,ell,N,dt,T,pinc)
%
%               a = fiber radius (cm)
%               ell = fiber length (cm)
%               N = number of compartments
%               dt = timestep (ms) 
%               T = end time
%               pinc = number of time steps between plots
%
% example:     cnpfib(.001,1,100,.03,10,10)
%

function [verrbe,verr] = cnpfib(a,ell,N,dt,T,pinc)

C_m = 1;		% micro F / cm^2
g_L = 0.3;     		% mS / cm^2
R_2 = 0.034;		% k Ohm cm

dx = ell/N;		% patch length
x = dx/2:dx:ell-dx/2;	% vector of patch midpoints
q1 = cos(pi*x/ell)'*sqrt(2/ell);
q1 = q1/(10*2*pi*a);

v = zeros(N,1);		% initial conditions
vbe = v;

e = ones(N,1);
S = spdiags([-e 2*e -e], -1:1, N, N)/dx/dx;
S(1,1) = 1/dx/dx;
S(N,N) = 1/dx/dx;

tau = C_m/g_L;
lambda = sqrt(a/(2*R_2*g_L));
zeta = (1+lambda^2*4*(N/ell)^2*sin(pi/2/N)^2)/tau;

B = (speye(N)+lambda^2*S)/tau;

[Lbe,Ube] = lu(speye(N) + dt*B);

[L,U] = lu(speye(N)+B*dt/2);

t = 0;
tcnt = 0;
if pinc
   tcnt = 0;
   plot3(x,t*e,v)
   hold on
end

fpre = q1*I(0);
t = dt;
fcur = q1*I(t);

rhs = (speye(N)-B*dt/2)*v + (fpre+fcur)*dt/2;

verr = 0;
verrbe = verr;

while (t <= T)
    
      rhsbe = vbe + fcur*dt;
      vbe = Ube \ ( Lbe \ rhsbe );
      
      v = U \ ( L \ rhs );
      
      vex = q1.*(exp(-t)*(zeta-2)+exp(-2*t)*(1-zeta)+exp(-zeta*t));
      vex = vex./(zeta-1)./(zeta-2);
      
      verr = max(verr,max(abs(vex-v)));
      verrbe = max(verrbe,max(abs(vex-vbe)));
      
      tcnt = tcnt + 1;
      if mod(tcnt,pinc) == 0
         plot3(x,t*e,abs(v-vex))
         %plot3(x,t*e,vex,'r')
      end
      
      t = t+dt;
      fpre = fcur;
      fcur = q1*I(t);
      
      rhs = 2*v - rhs + (fpre+fcur)*dt/2;

end

if pinc
   xlabel('x (cm)','fontsize',16)
   ylabel('t (ms)','fontsize',16)
   zlabel('v (mV)','fontsize',16)
end

hold off

function val = I(t)
val = (exp(-t)-exp(-2*t));

%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% trapezoid on the fork
%
%  usage:   trapfork(dt,stim,pinc)
%
%               stim.t1 = start of current pulse (ms)
%               stim.t2 = end of current pulse (ms)
%               stim.amp = amplitude of current pulse (micro amps)
%               stim.loc = location of current pulse (compartment num)
%               stim.Tfin = stopping time (ms)
% example:      
%       stim = struct('t1',1,'t2',2,'amp',1e-4,'loc',50,'Tfin',10)
%        dt = 0.1,    pinc = 4 
%
% or, for a pair of stimuli
%
% stim = struct('t1',[1 2],'t2',[2 3],'amp',1e-4*[1 1],...
%               'loc',[100 350],'Tfin',10)
%
%  or for say s=10 inputs over time Tfin=60
%  t1 = round(Tfin*rand(1,s));   t2 = t1+1;
%  loc = round(750*rand(1,10));    % our fork has 751 compartments
%  stim = struct('t1',t1,'t2',t2,'amp',1e-4*ones(1,s),'loc',loc,'Tfin',Tfin)
%

function trapfork(dt,stim,pinc)

a = 1e-4*[1 1 1];
ell = [2.5 2.5 2.5]/100;
dx = .0001;
N = ell/dx;
A3 = 2*pi*a(3)*dx;
As = 4*pi*1e-6;
rho = A3/As;
R2 = 0.3;
gL = 1/15;
Cm = 1;
tau = Cm/gL;
lam = a/(2*R2*gL)/dx^2;   % lambda^2
r = a/a(3);
Hd = [2*lam(1)*ones(1,N(1)) 2*lam(2)*ones(1,N(2)) 2*lam(3)*ones(1,N(3)+1)];
Hd(1) = lam(1);
Hd(N(1)+1) = lam(2);
Hd(N(1)+N(2)+1) =  lam*r';
Hd(end) = rho*lam(3);
Hlen = length(Hd)

Hu = [-lam(1)*ones(1,N(1)-1) 0 -lam(2)*ones(1,N(2)) -lam(3)*ones(1,N(3))];
Hl = [-lam(1)*ones(1,N(1)-1) 0 -lam(2)*ones(1,N(2)-1) -r(2)*lam(2) -lam(3)*ones(1,N(3))];
Hl(end) = rho*Hl(end);

H = spdiags( [[Hl 0]' Hd' [0 Hu]'], -1:1, Hlen, Hlen);

H(N(1)+N(2)+1,N(1)) = -r(1)*lam(1);
H(N(1),N(1)+N(2)+1) = -lam(1);

I = speye(Hlen);

Bb = I+(I+H)*(dt/tau/2); 
[L,U] = lu(Bb); 

x3 = 0:dx:ell(3);
x1 = ell(3):dx:ell(3)+ell(1);
x2 = ell(3):dx:ell(3)+ell(2);

v = zeros(Hlen,1);         % initial conditions
rhs = zeros(Hlen,1);

t = 0;
tcnt = 0;
x = linspace(0,ell(3)+max(ell(1:2)),Hlen);
figure(1)
plot3(x,t*ones(Hlen,1),v)
hold on

eloc = round(Hlen*1e-4*stim.loc/sum(ell));

Nt = ceil(stim.Tfin/dt);

vhot = zeros(1,Nt);

stim.amp = stim.amp*(dt/2)/A3/Cm;

f0 = stim.amp.*(t>stim.t1).*(t<stim.t2);
t = dt;
f1 = stim.amp.*(t>stim.t1).*(t<stim.t2);

r = zeros(Hlen,1);
r(eloc) = (f0 + f1)';

for j=2:Nt,

    v = U \ ( L \ r );

    vhot(j) = v(Hlen);

    if mod(j,pinc) == 0
         v1 = fliplr(v(1:N(1))');
         v2 = fliplr(v(N(1)+1:N(1)+N(2))');
         v3 = fliplr(v(N(1)+N(2)+1:end)');
         plot3(x3,t*ones(size(x3)),v3,'r')
         plot3(x2,t*ones(size(x2)),[v3(end) v2],'r')
         plot3(x1,t*ones(size(x1)),[v3(end) v1],'k')
    end

    t = j*dt;
    f0 = f1;
    f1 = stim.amp.*(t>stim.t1).*(t<stim.t2);

    r = 2*v - r;
    r(eloc) = r(eloc) + (f0 + f1)';

end

xlabel('x (cm)','fontsize',14)
ylabel('t (ms)','fontsize',14)
zlabel('v (mV)','fontsize',14)
axis tight
hold off

figure(2)
t = linspace(0,stim.Tfin,Nt);
plot(t,vhot,'k')
xlabel('t (ms)','fontsize',14)
ylabel('v (mV)','fontsize',14)
box off

figure(3)
plot([250 250],[0 500],'k','linewidth',2)
hold on
plot([1 250],[250 250],'k','linewidth',2)
[seloc,sind] = sort(eloc);
t1s = stim.t1(sind);
for j=1:length(eloc)
    p = seloc(j);
    if p<250,
       text(240+15*(-1)^j,500-p,num2str(t1s(j)))
    elseif p<500
       text(p-250,245+15*(-1)^j,num2str(t1s(j)))
    else
       text(240+15*(-1)^j,750-p,num2str(t1s(j)))
    end
end
axis equal
axis off
hold off

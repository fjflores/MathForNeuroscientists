%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% contrast the resonance curves for the uniform and nonuniform cables
%
% stEcabResdrive.m
%
function stEcabResdrive

cab = struct('rad',1e-4,'ell',.1,'dx',1e-3,'dt',0.05);
g = struct('K', 36, 'Na', 120, 'Cl', 1/15);
stim = struct('I0',10e-6,'Tfin',160);

[omHz,RinM] = stEcabRes(cab,g,stim,0);

x = cab.dx/2:cab.dx:cab.ell-cab.dx/2;
hot = ((x>0.005).*(x<0.01))';
g = struct('K', 40-20*hot, 'Na', 44+560*hot, 'Cl', 1/15);

[omHznon,RinMnon] = stEcabRes(cab,g,stim,0);

plot(omHz,RinM,'k')
hold on
plot(omHznon,RinMnon,'r')
box off
legend('Uniform','Nonuniform','location','best')
xlabel('\omega  (Hz)','fontsize',14)
ylabel('R_{in}  (M\Omega)','fontsize',14)

%
%  stEcabRes.m  
%
% solve the active cable problem with distributed input via staggered Euler
%
% usage:	stEcabD(cab,stim,pinc)
%
%		cab.rad = fiber radius (cm)
%		cab.ell = fiber length (cm)
%       cab.dx = space step (cm)
%       cab.dt = timestep (ms)
%
%	x = cab.dx/2:cab.dx:cab.ell-cab.dx/2;       % vector of patch midpoints
%	hot = ((x>0.005).*(x<0.01))';
%	g = struct('K', 40-20*hot, 'Na', 44+560*hot, 'Cl', 1/15);
%
%		stim.I0 = mag of current stim
%       stim.Tfin = stopping time (ms)
%
%       pinc = number of time steps between plots
%
% example:	cab = struct('rad',1e-4,'ell',.1,'dx',1e-3,'dt',0.05)
%           stim = struct('I0',10e-6,'Tfin',160)
%           pinc = 0 
%

function [omHz,RinM] = stEcabRes(cab,g,stim,pinc)

Cm = 1;		% micro F / cm^2
R2 = 0.3;		% k Ohm cm
dx = cab.dx;
dt = cab.dt;
Nx = cab.ell/dx;		% patch length
A = 2*pi*cab.rad;		% patch surface area
x = dx/2:dx:cab.ell-dx/2;	% vector of patch midpoints

Nt = ceil(stim.Tfin/dt);

E = struct('K', -77, 'Na', 56, 'Cl', -68);
%g = struct('K', 36, 'Na', 120, 'Cl', 1/15);

e = ones(Nx,1);
B = spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx/dx;
B(1,1) = 1/dx/dx;
B(Nx,Nx) = 1/dx/dx;
B = (cab.rad/2/R2)*B;
dB = diag(B);

options = optimset('Jacobian','on');
Vr = fsolve(@(V) Iss(V,E,g,B),-70*e,options);    % initial conditions

a = am(Vr); b = bm(Vr);
taum = 1./(a+b);
m = a.*taum;
dm = (dam(Vr).*b - a.*dbm(Vr)).*taum.^2;

a = ah(Vr); b = bh(Vr);
tauh = 1./(a+b);
h = a.*tauh;
dh = (dah(Vr).*b - a.*dbh(Vr)).*tauh.^2;

a = an(Vr); b = bn(Vr);
taun = 1./(a+b);
n = a.*taun;
dn = (dan(Vr).*b - a.*dbn(Vr)).*taun.^2;

Q11 = -(B + spdiags(g.Cl+g.Na.*m.^3.*h+g.K.*n.^4,0,Nx,Nx))/Cm;
Q12 = -spdiags(3*g.Na.*m.^2.*h.*(Vr-E.Na)/Cm,0,Nx,Nx);
Q13 = -spdiags(g.Na.*m.^3.*(Vr-E.Na)/Cm,0,Nx,Nx);
Q14 = -spdiags(4*g.K.*n.^3.*(Vr-E.K)/Cm,0,Nx,Nx);

Q21 = spdiags(dm./taum,0,Nx,Nx);
Q22 = spdiags(-1./taum,0,Nx,Nx);

Q31 = spdiags(dh./tauh,0,Nx,Nx);
Q33 = spdiags(-1./tauh,0,Nx,Nx);

Q41 = spdiags(dn./taun,0,Nx,Nx);
Q44 = spdiags(-1./taun,0,Nx,Nx);

Z = 0*speye(Nx);

Q = [Q11 Q12 Q13 Q14;
     Q21 Q22  Z   Z;
     Q31  Z  Q33  Z;
     Q41  Z   Z  Q44];

[evec,eval] = eig(full(Q));

[val,ind] = max(real(diag(eval)));

q0 = real(evec(1:Nx,ind));
Ix = stim.I0*q0';

%figure(1)
%plot(q0)
%pause

for k=1:45
    
    omega(k) = (k+10)/1000;
    
    k
    
    V = Vr;
    Vmax(k) = 0;
    
    n = an(V)./(an(V)+bn(V));
    m = am(V)./(am(V)+bm(V));
    h = ah(V)./(ah(V)+bh(V));
    
    t = 0;
    if pinc
        plot3(x,t*ones(Nx,1),V)
        hold on
    end
    
    for j=2:Nt,
        
        t = (j-1)*dt;
        
        I = Ix*sin(2*pi*(t-dt/2)*omega(k))/A;
        
        a = an(V);  c = (a+bn(V))/2;
        n = ( (1/dt-c).*n + a) ./ (1/dt + c); n4 = n.^4;
        
        a = am(V);  c = (a+bm(V))/2;
        m = ( (1/dt-c).*m + a) ./ (1/dt + c);
        
        a = ah(V);  c = (a+bh(V))/2;
        h = ( (1/dt-c).*h + a) ./ (1/dt + c); m3h = m.^3.*h;
        
        d = g.Na.*m3h + g.K.*n4 + g.Cl;
        
        f = g.Na.*m3h*E.Na + g.K.*n4*E.K + g.Cl.*E.Cl + I';
        
        B(1:Nx+1:end) = dB + d + 2*Cm/dt;         % update the diagonal
        
        Vmid = B\(2*Cm*V/dt + f);
        
        V = 2*Vmid - V;
        
        if mod(j,pinc) == 0
            plot3(x,t*ones(Nx,1),V)
        end
        
        if j>round(Nt/2)
            Vmax(k) = max(Vmax(k),max(V-Vr));
        end
        
    end
    
    if pinc
        xlabel('x  (cm)','fontsize',16)
        ylabel('t  (ms)','fontsize',16)
        zlabel('V  (mV)','fontsize',16)
        drawnow
        hold off
    end
    
end % k

omHz = omega*1000;
RinM = Vmax/stim.I0/1000;

%figure(2)
%I0 = max(Ix);
%plot(omega*1000,Vmax/I0/1000)
%grid
%xlabel('\omega  (Hz)','fontsize',14)
%ylabel('R_{in}  (M\Omega)','fontsize',14)

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



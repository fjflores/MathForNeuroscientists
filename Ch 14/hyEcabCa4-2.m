%
% hyEcabCa4.m  
%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% solve the active cable via sparse hybrid trapezoid
%
% with 3 calcium channels AND NaCa exchange
%
%  and solution of the u=[c;b] reaction diffusion system
%
%  and gKCa  and ER with Ryanodine Receptors
%
% usage: [ca, t] = hyEcabCa4(cab,stim,pinc)
%
%		cab.rad = cable radius (cm)
%		cab.ell = cable length (cm)
%       cab.N = number of compartments
%       cab.dt = timestep (ms)
%		stim.t1 = start of current pulse (ms)
%		stim.t2 = end of current pulse (ms)
%       stim.Tfin = stopping time (ms)
%		stim.Iapp = size of current pulse (micro amps)
%       stim.cstim = size of cstim
%       pinc = number of time steps between plots
%
% I_stim example: (fig 13.11)
%
% cab = struct('rad',.0001,'ell',.1,'N',200,'dt',.02)
% stim = struct('t1',3,'t2',4,'Tfin',20,'Iapp',3e-4,'cstim',0)
% pinc = 20  
%
% Ca_stim example: (fig 13.12)
%
% cab = struct('rad',.0001,'ell',.1,'N',400,'dt',.02)
% stim = struct('t1',3,'t2',4,'Tfin',100,'Iapp',0,'cstim',1)
% pinc = 100
%
%

function [ca, t] = hyEcabCa4(cab,stim,pinc)

Cm = 1;         % micro F / cm^2
R2 = 0.3;       % 0.034;  % k Ohm cm
N = cab.N;
ell = cab.ell;
dx = ell/N;		% patch length
dt = cab.dt;
A = 2*pi*cab.rad*dx;		% patch surface area
x = dx/2:dx:cab.ell-dx/2;	% vector of patch midpoints
Nmid = round(N/2);

Nt = ceil(stim.Tfin/cab.dt)+1;
vmid = zeros(Nt+1,1); ca = vmid; bca = ca; INaCa = ca;
wmid = ca; IKCa = ca; smid = ca;

E = struct('K', -77, 'Na', 56, 'Cl', -68);
g = struct('K', 36, 'Na', 120, 'Cl', 1/15, ...
           'CaT', .25, 'CaN', 2.5, 'CaL', 2.5, 'KCa',10);
        
e1 = zeros(N,1); e1(Nmid) = 1;
e2 = e1;
Iapp = stim.Iapp*e1/A;
 
e = ones(N,1);

cr = .05*e;     % micro molar
co = 1e3;       % micro molar

buff = struct('k1',1.5e-3,'k2',0.3e-3,'BT',500,...
              'k1e',1.5e-3,'k2e',0.3e-3,'QT',500);

Diff = struct('c',220e-11,'b',110e-11,'s',220e-11,'q',110e-11);

pmca = struct('max', 2e-6, 'K', 1);

naca = struct('max', 100, 'Na3', (50/440)^3);

serca = struct('max', 2e-4, 'K', 2);

ryr = struct('nu',1e-6,'Ka',0.372,'Kb',0.636, 'Kc',0.057, 'tm',1e4);

sr = 0.5*1e3*e;
qr = buff.QT*sr./(sr+buff.k2e/buff.k1e);

rade = cab.rad/2;

br = buff.BT*cr./(cr + buff.k2/buff.k1);
F = 96485.3399;     % coulombs / mole   Faraday's constant
VT = 25.8;  % mV

S = spdiags([-e 2*e -e], -1:1, N, N)/dx/dx;
S(1,1) = 1/dx/dx;
S(N,N) = 1/dx/dx;
lam2 = cab.rad/2/R2;
B = lam2*S;

eN = speye(N);
Ri = [-Diff.c*S-buff.k1*buff.BT*eN  buff.k2*eN; buff.k1*buff.BT*eN -Diff.b*S-buff.k2*eN];
Re = [-Diff.s*S-buff.k1e*buff.QT*eN  buff.k2e*eN; buff.k1e*buff.QT*eN -Diff.q*S-buff.k2e*eN];
R = spalloc(4*N,4*N,16*N-8);
R(1:2*N,1:2*N) = Ri;
R(2*N+1:4*N,2*N+1:4*N) = Re;
Rp = 2*speye(4*N) + dt*R;
Rm = 2*speye(4*N) - dt*R;
[LRm,URm] = lu(Rm);


ci = cr;
Vr = fsolve(@(V) Iss(V,E,g,B,ci,co,VT),-70*e);  % initial conditions
%max(Vr)

n1 = an(Vr)./(an(Vr)+bn(Vr)); n14 = n1.^4;
m1 = am(Vr)./(am(Vr)+bm(Vr)); 
h1 = ah(Vr)./(ah(Vr)+bh(Vr)); m13h1 = m1.^3.*h1;

m1T = amT(Vr)./(amT(Vr)+bmT(Vr)); 
h1T = ahT(Vr)./(ahT(Vr)+bhT(Vr)); m12h1T = m1T.^2.*h1T;

m1N = amN(Vr)./(amN(Vr)+bmN(Vr)); 
h1N = ahN(Vr)./(ahN(Vr)+bhN(Vr)); m12h1N = m1N.^2.*h1N;

m1L = amL(Vr)./(amL(Vr)+bmL(Vr)); m12L = m1L.^2;

inaca = naca.max*(naca.Na3*exp(Vr/2/VT)-(ci/co).*exp(-Vr/2/VT));

m1KCa = amKCa(Vr,cr,VT)./(amKCa(Vr,cr,VT)+bmKCa(Vr,cr,VT));

w1 = winf(cr,ryr);

ca(1) = cr(Nmid);
bca(1) = br(Nmid);
vmid(1) = Vr(Nmid);
wmid(1) = w1(Nmid);
smid(1) = sr(Nmid);

u = [cr; br; sr; qr];           % initially at rest

ICa1 = (g.CaT*m12h1T + g.CaN*m12h1N + g.CaL*m12L).*Phi(u(1:N),co,Vr,VT);

dt = cab.dt;

v = Vr;
a = an(v);  b = bn(v);
n2 = ( (2/dt-a-b).*n1 + 2*a) ./ (2/dt + a + b); n24 = n2.^4;

a = am(v);  b = bm(v);
m2 = ( (2/dt-a-b).*m1 + 2*a) ./ (2/dt + a + b);

a = ah(v);  b = bh(v);
h2 = ( (2/dt-a-b).*h1 + 2*a) ./ (2/dt + a + b); m23h2 = m2.^3.*h2;

a = amT(v);  b = bmT(v);
m2T = ( (2/dt-a-b).*m1T + 2*a) ./ (2/dt + a + b);
a = ahT(v);  b = bhT(v);
h2T = ( (2/dt-a-b).*h1T + 2*a) ./ (2/dt + a + b); m22h2T = m2T.^2.*h2T;

a = amN(v);  b = bmN(v);
m2N = ( (2/dt-a-b).*m1N + 2*a) ./ (2/dt + a + b);
a = ahN(v);  b = bhN(v);
h2N = ( (2/dt-a-b).*h1N + 2*a) ./ (2/dt + a + b); m22h2N = m2N.^2.*h2N;

a = amL(v);  b = bmL(v);
m2L = ( (2/dt-a-b).*m1L + 2*a) ./ (2/dt + a + b); m22L = m2L.^2;

a = amKCa(v,cr,VT);  b = bmKCa(v,cr,VT);
m2KCa = ( (2/dt-a-b).*m1KCa + 2*a) ./ (2/dt + a + b);

wi = winf(cr,ryr); tw = ryr.tm*wi;
w2 = ( (2*tw-dt).*w1 + 2*wi*dt) ./ (2*tw + dt);

inaca = naca.max*(naca.Na3*exp(v/2/VT)-(u(1:N)/co).*exp(-v/2/VT));

ca(2) = u(Nmid);
bca(2) = u(Nmid+N);
vmid(2) = Vr(Nmid);
wmid(2) = w2(Nmid);
smid(2) = u(2*N+Nmid);

t = 0;
if pinc 
   figure(2)
   plot3(x,t*ones(N,1),u(1:N),'color','k')
   hold on
   figure(3)
   plot3(x,t*ones(N,1),u(2*N+1:3*N),'color','k')
   hold on
end

i1 = 0;
t = cab.dt;
i2 = Iapp*(t>stim.t1)*(t<stim.t2);

dB = diag(B);

tf = 2*Cm/cab.dt;

d2 = g.Na*m23h2 + g.K*n24 + g.Cl + g.KCa*m2KCa;

ICa2 = (g.CaT*m22h2T + g.CaN*m22h2N + g.CaL*m22L).*Phi(u(1:N),co,v,VT);
   
f = g.Na*(m23h2 + m13h1)*E.Na + g.K*(n24+n14)*E.K + 2*g.Cl*E.Cl ...
      - 2*inaca + i2 + i1 - (ICa1 + ICa2) + g.KCa*(m2KCa+m1KCa)*E.K;

r = (tf - d2).*v - B*v + f;

B(1:N+1:end) = dB + d2 + tf;         % update the diagonal

v = B\r;   % v at dt

if pinc
   figure(2)
   plot3(x,t*ones(N,1),u(1:N),'color','k')
   figure(3)
   plot3(x,t*ones(N,1),u(2*N+1:3*N),'color','k')
end

i1 = i2; n1 = n2; m1 = m2; h1 = h2; d1 = d2; m13h1 = m23h2; n14 = n24; 
h1N = h2N; m1N = m2N; h1T = h2T; m1T = m2T; m1L = m2L; ICa1 = ICa2;
m1KCa = m2KCa; w1 = w2;

for j=2:Nt,
    
    er = Jryr(u(1:N),u(2*N+1:3*N),w1,ryr) - Jserca(u(1:N),serca);

    pin = 2*cab.rad*((-ICa2+inaca)/(2*F)-Jpmca(u(1:N),pmca))/(cab.rad^2-rade^2) +...
          2*rade*er/(cab.rad^2-rade^2) + (t>stim.t1)*(t<stim.t2)*stim.cstim*e2;

    tmp = buff.k1*u(1:N).*u(N+1:2*N);

    pe = -2*er/rade;

    tmpe = buff.k1e*u(2*N+1:3*N).*u(3*N+1:4*N);

    q = 2*[tmp + pin; -tmp; tmpe + pe; -tmpe];

    u = URm\(LRm\(Rp*u + dt*q));

    t = j*dt;

    i2 = Iapp*(t>stim.t1)*(t<stim.t2);

    a = an(v);  b = bn(v);
    n2 = ( (2/dt-a-b).*n1 + 2*a) ./ (2/dt + a + b); n24 = n2.^4;

    a = am(v);  b = bm(v);
    m2 = ( (2/dt-a-b).*m1 + 2*a) ./ (2/dt + a + b);    

    a = ah(v);  b = bh(v);
    h2 = ( (2/dt-a-b).*h1 + 2*a) ./ (2/dt + a + b); m23h2 = m2.^3.*h2;

    a = amT(v);  b = bmT(v);
    m2T = ( (2/dt-a-b).*m1T + 2*a) ./ (2/dt + a + b);
    a = ahT(v);  b = bhT(v);
    h2T = ( (2/dt-a-b).*h1T + 2*a) ./ (2/dt + a + b); m22h2T = m2T.^2.*h2T;

    a = amN(v);  b = bmN(v);
    m2N = ( (2/dt-a-b).*m1N + 2*a) ./ (2/dt + a + b);
    a = ahN(v);  b = bhN(v);
    h2N = ( (2/dt-a-b).*h1N + 2*a) ./ (2/dt + a + b); m22h2N = m2N.^2.*h2N;

    a = amL(v);  b = bmL(v);
    m2L = ( (2/dt-a-b).*m1L + 2*a) ./ (2/dt + a + b); m22L = m2L.^2;

    a = amKCa(v,u(1:N),VT);  b = bmKCa(v,u(1:N),VT);
    m2KCa = ( (2/dt-a-b).*m1KCa + 2*a) ./ (2/dt + a + b);

    wi = winf(u(1:N),ryr); tw = ryr.tm*wi;
    w2 = ( (2*tw-dt).*w1 + 2*wi*dt) ./ (2*tw + dt);

    d2 = g.Na*m23h2 + g.K*n24 + g.Cl + g.KCa*m2KCa;

    ICa2 = (g.CaT*m22h2T + g.CaN*m22h2N + g.CaL*m22L).*Phi(u(1:N),co,v,VT);
    
    inaca = naca.max*(naca.Na3*exp(v/2/VT)-(u(1:N)/co).*exp(-v/2/VT));

    f = g.Na*(m23h2 + m13h1)*E.Na + g.K*(n24+n14)*E.K + 2*g.Cl*E.Cl ...
        - 2*inaca + i2 + i1 - (ICa1 + ICa2) + g.KCa*(m2KCa+m1KCa)*E.K;

    r = (2*tf - d2 + d1).*v - r + f;

    B(1:N+1:end) = dB + d2 + tf;         % update the diagonal

    v = B\r;

    if mod(j,pinc) == 0
        figure(2)
        plot3(x,t*ones(N,1),u(1:N),'color','k')
        figure(3)
        plot3(x,t*ones(N,1),u(2*N+1:3*N),'color','k')
    end

    vmid(j+1) = v(Nmid);
    ca(j+1) = u(Nmid);
    bca(j+1) = u(Nmid+N);
    wmid(j+1) = w2(Nmid);
    smid(j+1) = u(2*N+Nmid);

    i1 = i2; n1 = n2; m1 = m2; h1 = h2; d1 = d2; m13h1 = m23h2; n14 = n24;
    h1N = h2N; m1N = m2N; h1T = h2T; m1T = m2T; m1L = m2L; ICa1 = ICa2;
    m1KCa = m2KCa;

end

if pinc > 0
   figure(2)
   xlabel('x  (cm)','fontsize',14)
   ylabel('t  (ms)','fontsize',14)
   zlabel('c  (\muM)','fontsize',14)
   hold off
   figure(3)
   xlabel('x  (cm)','fontsize',14)
   ylabel('t  (ms)','fontsize',14)
   zlabel('s   (\muM)','fontsize',14)
end

t = linspace(0,stim.Tfin,length(vmid))';

figure(1)

subplot(4,1,1)
plot(t,vmid,'linewidth',1,'color','k')
axis tight
ylabel('V (mV)','fontsize',14)
set(gca,'xticklabel',[])
box off

subplot(4,1,4)
plot(t,smid,'linewidth',1,'color','k')
axis tight
xlabel('t  (ms)','fontsize',14)
ylabel('s (\muM)','fontsize',14)
box off

subplot(4,1,2)
plot(t,ca,'linewidth',1,'color','k')
axis tight
ylabel('c (\muM)','fontsize',14)
set(gca,'xticklabel',[])
box off

subplot(4,1,3)
plot(t,bca,'linewidth',1,'color','k')
axis tight
ylabel('b (\muM)','fontsize',14)
set(gca,'xticklabel',[])
box off

return

function val = Iss(V,E,g,B,ci,co,VT)
val = B*V + g.Na*(am(V)./(am(V)+bm(V))).^3.*...
            (ah(V)./(ah(V)+bh(V))).*(V-E.Na) + ...
            g.K*(an(V)./(an(V)+bn(V))).^4.*(V-E.K) + g.Cl.*(V-E.Cl) + ...
    g.KCa*(amKCa(V,ci,VT)./(amKCa(V,ci,VT)+bmKCa(V,ci,VT))).*(V-E.K) + ...
   (g.CaT*(amT(V)./(amT(V)+bmT(V))).^2.*(ahT(V)./(ahT(V)+bhT(V))) + ...
    g.CaN*(amN(V)./(amN(V)+bmN(V))).^2.*(ahN(V)./(ahN(V)+bhN(V))) + ...
    g.CaL*(amL(V)./(amL(V)+bmL(V))).^2) .*Phi(ci,co,V,VT);

function val = an(v)
val = .01*(10-(v+71))./(exp(1-(v+71)/10)-1);

function val = bn(v)
val = .125*exp(-(v+71)/80);

function val = am(v)
val = .1*(25-(v+71))./(exp(2.5-(v+71)/10)-1);

function val = bm(v)
val = 4*exp(-(v+71)/18);

function val = ah(v)
val = 0.07*exp(-(v+71)/20);

function val = bh(v)
val = 1./(exp(3-(v+71)/10)+1);

function val = Phi(ci,co,v,VT)
e = exp(2*v/VT);
val = v.*(1-(ci/co).*e)./(1-e);

function val = amL(V)
a_L = 15.69; b_L = 81.5; c_L = 0.29; d_L = 10.86;
val = a_L*(b_L-V)./(exp((b_L-V)/10)-1);

function val = bmL(V)
a_L = 15.69; b_L = 81.5; c_L = 0.29; d_L = 10.86;
val = c_L*exp(-V/d_L);

function val = amN(V)
a_N = 0.19; b_N = 19.88; c_N = 0.046; d_N = 20.73;
e_N = 1.6e-4; f_N=48.46; g_N=39;
val = a_N*(b_N-V)./(exp((b_N-V)/10)-1);

function val = bmN(V)
a_N = 0.19; b_N = 19.88; c_N = 0.046; d_N = 20.73;
e_N = 1.6e-4; f_N=48.46; g_N=39;
val = c_N*exp(-V/d_N);

function val = ahN(V)
a_N = 0.19; b_N = 19.88; c_N = 0.046; d_N = 20.73;
e_N = 1.6e-4; f_N=48.46; g_N=39;
val = e_N*exp(-V/f_N);

function val = bhN(V)
a_N = 0.19; b_N = 19.88; c_N = 0.046; d_N = 20.73;
e_N = 1.6e-4; f_N=48.46; g_N=39;
val  = 1./(1+exp((g_N-V)/10));

function val = amT(V)
a_T = 0.2; b_T = 19.26; c_T = 0.009; d_T = 22.03;
e_T = 1e-6; f_T=16.26; g_T=29.79;
val = a_T*(b_T-V)./(exp((b_T-V)/10)-1);

function val = bmT(V)
a_T = 0.2; b_T = 19.26; c_T = 0.009; d_T = 22.03;
e_T = 1e-6; f_T=16.26; g_T=29.79;
val = c_T*exp(-V/d_T);

function val = ahT(V)
a_T = 0.2; b_T = 19.26; c_T = 0.009; d_T = 22.03;
e_T = 1e-6; f_T=16.26; g_T=29.79;
val = e_T*exp(-V/f_T);

function val = bhT(V)
a_T = 0.2; b_T = 19.26; c_T = 0.009; d_T = 22.03;
e_T = 1e-6; f_T=16.26; g_T=29.79;
val = 1./(1+exp((g_T-V)/10));

function val = amKCa(v,c,VT)
val = 0.28*c./(c+.48*exp(-1.7*v/VT));

function val = bmKCa(v,c,VT)
val = 0.48./(1+(c/.13e-3).*exp(2*v/VT));

function val = winf(c,ryr)
stim.t1 = (ryr.Ka./c).^4;
stim.t2 = (c./ryr.Kb).^3;
top = 1 + stim.t1 + stim.t2;
val = top./(top+1/ryr.Kc);

function val = Jryr(c,s,w,ryr)
tmp1 = (ryr.Ka./c).^4;
tmp2 = (c./ryr.Kb).^3;
top = 1 + tmp2;
val = ryr.nu*w .*top.*(s-c)./(top+tmp1);

function val = Jpmca(c,pmca)
val = pmca.max*( c./(c+pmca.K) );

function val = Jserca(c,serca)
val = serca.max*( c.^2./(c.^2+(serca.K).^2) );


  

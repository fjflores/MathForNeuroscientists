
%
% hyEcabCa1.m   
%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% solve the active fiber via sparse hybrid trapezoid
%
% with 3 calcium channels 
%
% usage:	hyEcabCa1(cab,stim,gCa,pinc,pcur)
%
%		cab.a = cable radius (cm)
%		cab.ell = cable length (cm)
%       cab.N = number of compartments
%       cab.dt = timestep (ms)
%		stim.t1 = start of current pulse (ms)
%		stim.t2 = end of current pulse (ms)
%       stim.Tfin = stopping time (ms)
%		stim.Iapp = size of current pulse (micro amps)
%       pinc = number of time steps between plots
%       pcur = 1 for plot of Ca current, 0 otherwise
%
% example:	
%           cab = struct('rad',.0001,'ell',.1,'N',200,'dt',.01)
%           stim = struct('t1',1,'t2',2,'Tfin',15,'Iapp',.001)
%           gCa = struct('T',.25,'N',2.5,'L',2.5)
%           pinc = 20
%           pcur = 1
%

function [vmid, t] = hyEcabCa1(cab,stim,gCa,pinc,pcur)

Cm = 1;		% micro F / cm^2
R2 = 0.3;		% k Ohm cm
dx = cab.ell/cab.N;		% patch length
A = 2*pi*cab.rad*dx;		% patch surface area
x = dx/2:dx:cab.ell-dx/2;	% vector of patch midpoints
dt = cab.dt;
Nmid = round(cab.N/2);

Nt = ceil(stim.Tfin/cab.dt)+1;
ICaT = zeros(Nt+1,1); ICaN = ICaT; ICaL = ICaT; vmid = ICaT;
INa = ICaT; 

E = struct('K', -77, 'Na', 56, 'Cl', -68);
g = struct('K', 36, 'Na', 120, 'Cl', 1/15, ...
           'CaT', gCa.T, 'CaN', gCa.N, 'CaL', gCa.L);
        
ci = .05;
co = 1e3;

e1 = zeros(cab.N,1);
e1(Nmid/2) = 1;
Iapp = stim.Iapp*e1/A;

tau = Cm;
lam2 = cab.rad/2/R2;

e = ones(cab.N,1);
B = spdiags([-e 2*e -e], -1:1, cab.N, cab.N)/dx/dx;
B(1,1) = 1/dx/dx;
B(cab.N,cab.N) = 1/dx/dx;
B = lam2*B;

Vr = fsolve(@(V) Iss(V,E,g,B,ci,co),-70*e); % initial conditions

n1 = an(Vr)./(an(Vr)+bn(Vr)); n14 = n1.^4;
m1 = am(Vr)./(am(Vr)+bm(Vr)); 
h1 = ah(Vr)./(ah(Vr)+bh(Vr)); m13h1 = m1.^3.*h1;

m1T = amT(Vr)./(amT(Vr)+bmT(Vr)); 
h1T = ahT(Vr)./(ahT(Vr)+bhT(Vr)); m12h1T = m1T.^2.*h1T;

m1N = amN(Vr)./(amN(Vr)+bmN(Vr)); 
h1N = ahN(Vr)./(ahN(Vr)+bhN(Vr)); m12h1N = m1N.^2.*h1N;

m1L = amL(Vr)./(amL(Vr)+bmL(Vr)); m12L = m1L.^2;

ICaT(1) = g.CaT*m12h1T(Nmid)*Phi(ci,co,Vr(Nmid));
ICaN(1) = g.CaN*m12h1N(Nmid)*Phi(ci,co,Vr(Nmid));
ICaL(1) = g.CaL*m12L(Nmid)*Phi(ci,co,Vr(Nmid));
INa(1) = g.Na*m13h1(Nmid)*(Vr(Nmid)-E.Na);
vmid(1) = Vr(Nmid);

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

ICaT(2) = g.CaT*m22h2T(Nmid)*Phi(ci,co,Vr(Nmid));
ICaN(2) = g.CaN*m22h2N(Nmid)*Phi(ci,co,Vr(Nmid));
ICaL(2) = g.CaL*m22L(Nmid)*Phi(ci,co,Vr(Nmid));
INa(2) = g.Na*m23h2(Nmid)*(Vr(Nmid)-E.Na);
vmid(2) = Vr(Nmid);

t = 0;
if pinc > 0
   figure(1)
   plot3(x,t*ones(cab.N,1),v,'k')
   hold on
end

i1 = 0;
t = dt;
i2 = Iapp*(t>stim.t1)*(t<stim.t2);

dB = diag(B);

tf = 2*tau/dt;

d2 = g.Na*m23h2 + g.K*n24 + g.Cl;

f = g.Na*(m23h2 + m13h1)*E.Na + g.K*(n24+n14)*E.K + 2*g.Cl*E.Cl +...
    i2 + i1 - (g.CaT*(m22h2T+m12h1T)+g.CaN*(m22h2N+m12h1N)+...
    g.CaL*(m22L+m12L)).*Phi(ci,co,v);

r = (tf - d2).*v - B*v + f;

B(1:cab.N+1:end) = dB + d2 + tf;         % update the diagonal

v = B\r;

if pinc > 0
   plot3(x,t*ones(cab.N,1),v,'color','k')
   hold on
   box off
end

i1 = i2; n1 = n2; m1 = m2; h1 = h2; d1 = d2; m13h1 = m23h2; n14 = n24; 
m12h1T = m22h2T; m12h1N = m22h2N; m12L = m22L;
h1N = h2N; m1N = m2N; h1T = h2T; m1T = m2T; m1L = m2L;

for j=2:Nt,

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
      h2T = ( (2/dt-a-b).*h1T + 2*a) ./...
          (2/dt + a + b); m22h2T = m2T.^2.*h2T;

      a = amN(v);  b = bmN(v);
      m2N = ( (2/dt-a-b).*m1N + 2*a) ./ (2/dt + a + b);
      a = ahN(v);  b = bhN(v);
      h2N = ( (2/dt-a-b).*h1N + 2*a) ./...
          (2/dt + a + b); m22h2N = m2N.^2.*h2N;

      a = amL(v);  b = bmL(v);
      m2L = ( (2/dt-a-b).*m1L + 2*a) ./ (2/dt + a + b); m22L = m2L.^2;

      d2 = g.Na*m23h2 + g.K*n24 + g.Cl;
      
      f = g.Na*(m23h2 + m13h1)*E.Na + g.K*(n24+n14)*E.K + 2*g.Cl*E.Cl + ...
           i2+i1 - (g.CaT*(m22h2T+m12h1T)+...
        g.CaN*(m22h2N+m12h1N)+g.CaL*(m22L+m12L)).*Phi(ci,co,v);

      r = (2*tf - d2 + d1).*v - r + f;

      B(1:cab.N+1:end) = dB + d2 + tf;         % update the diagonal

      v = B\r;

      if mod(j,pinc) == 0
         plot3(x,t*ones(cab.N,1),v,'color','k')
         box off
      end
      
      ICaT(j+1) = g.CaT*m22h2T(Nmid)*Phi(ci,co,v(Nmid));
      ICaN(j+1) = g.CaN*m22h2N(Nmid)*Phi(ci,co,v(Nmid));
      ICaL(j+1) = g.CaL*m22L(Nmid)*Phi(ci,co,v(Nmid));
      INa(j+1) = g.Na*m23h2(Nmid)*(v(Nmid)-E.Na);
      vmid(j+1) = v(Nmid);

      i1 = i2; n1 = n2; m1 = m2; h1 = h2; d1 = d2; m13h1 = m23h2;n14 = n24;
      m12h1T = m22h2T; m12h1N = m22h2N; m12L = m22L;
      h1N = h2N; m1N = m2N; h1T = h2T; m1T = m2T; m1L = m2L;

end

if pinc > 0
    xlabel('x  (cm)','fontsize',14,'color','k')
    ylabel('t  (ms)','fontsize',14,'color','k')
    zlabel('v  (mV)','fontsize',14,'color','k')
    hold off
end

t = linspace(0,stim.Tfin,length(ICaT))';

if pcur
subplot(2,1,1)
plot(t,vmid,'k')
set(gca,'xticklabel',[]);
ylabel('V  (mV)','fontsize',14,'color','k')
ylim([-80 50])
box off

subplot(2,1,2)
plot(t,ICaL,'-','linewidth',1,'color','r')
hold on
plot(t,ICaN,'--','linewidth',1,'color','k')
plot(t,ICaT,'linewidth',1,'color','k')
legend('I_{Ca,L}','I_{Ca,N}','I_{Ca,T}','location','best')
box off
hold off
ylim([-50 1])
xlabel('t  (ms)','fontsize',14,'color','k')
ylabel('\muA/cm^2','fontsize',14,'color','k')
end

return

function val = Iss(V,E,g,B,ci,co)
val = B*V + g.Na*(am(V)./(am(V)+bm(V))).^3.*(ah(V)./...
    (ah(V)+bh(V))).*(V-E.Na) + g.K*(an(V)./(an(V)+bn(V))).^4.*(V-E.K) +...
    g.Cl.*(V-E.Cl) + (g.CaT*(amT(V)./(amT(V)+bmT(V))).^2.*(ahT(V)./...
    (ahT(V)+bhT(V))) + g.CaN*(amN(V)./...
    (amN(V)+bmN(V))).^2.*(ahN(V)./(ahN(V)+bhN(V))) +...
    g.CaL*(amL(V)./(amL(V)+bmL(V))).^2) .*Phi(ci,co,V);

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

function val = Phi(ci,co,v)
F = 96485;      % C/mol
R = 8.314;      % J/deg K/mole
T = 300;        % deg K
e = exp(2*F*v*1e-3/R/T);
val = v.*(1-(ci/co)*e)./(1-e);

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


  

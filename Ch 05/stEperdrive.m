%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% drive stEper
%
%  stEperdrive
%

function stEperdrive

figure
[t,V]=stEper(0.01,200,35,1000/30); 
plot(t,V,'k')
xlabel('t  (ms)','fontsize',14)
ylabel('V  (mV)','fontsize',14)
box off
ylim([-80 45])

figure
[t,V]=stEper(0.01,200,35,1000/45); 
plot(t,V,'k')
xlabel('t  (ms)','fontsize',14)
ylabel('V  (mV)','fontsize',14)
box off
ylim([-80 45])

figure
[t,V]=stEper(0.01,200,35,1000/60); 
plot(t,V,'k')
xlabel('t  (ms)','fontsize',14)
ylabel('V  (mV)','fontsize',14)
box off
ylim([-80 45])

%
% stEper.m, 
%
% staggered Euler on the active sphere
%
% usage:  [t,v] = stEper(dt,Tfin,I0,per)
%
% example:  [t,v] = stEper(0.01,50,20,20)
%

function [t,V] = stEper(dt,Tfin,I0,per)

A = 4*pi*(1e-6);% (cm)^2

E = struct('K', -77, 'Na', 56, 'Cl', -68); % reversal potentials, mV

G = struct('K', 36, 'Na', 120, 'Cl', 0.3);  % channel conductances, mS/cm^2

Vr = fsolve(@(V) Iss(V,E,G),-71)   % find rest potential 

Nt = ceil(Tfin/dt);
t = zeros(Nt,1); V = t;
t(1) = 0;
V(1) = Vr;
n = an(Vr)/(an(Vr)+bn(Vr));  
m = am(Vr)/(am(Vr)+bm(Vr));  
h = ah(Vr)/(ah(Vr)+bh(Vr));  

td = 1/dt;

for j = 2:Nt;

      t(j) = (j-1)*dt;

      %Istim = I0*(t(j)-dt/2>5)*(t(j)-dt/2<7)*(1e-6);
      Istim = I0*(mod(t(j)-dt/2,per)<2)*(1e-6);
      
      a = an(V(j-1));  b = bn(V(j-1)); c = (a+b)/2;
      n = ( (td-c)*n + a ) / (td + c);
      
      cK = G.K*n^4;
      
      a = am(V(j-1));  b = bm(V(j-1)); c = (a+b)/2;
      m = ( (td-c)*m + a ) / (td + c);
      
      a = ah(V(j-1));  b = bh(V(j-1)); c = (a+b)/2;
      h = ( (td-c)*h + a ) / (td + c);
      
      cNa = G.Na*m^3*h;
      
      top = 2*V(j-1)*td + cK*E.K + cNa*E.Na + G.Cl*E.Cl + Istim/A;
      
      bot = 2*td + cK + cNa + G.Cl;
      
      Vmid = top/bot;

      V(j) = 2*Vmid - V(j-1); 

end

function val = Iss(V,E,G)
val = G.Na*(am(V)/(am(V)+bm(V)))^3*(ah(V)/(ah(V)+bh(V)))*(V-E.Na) + ...
      G.K*(an(V)/(an(V)+bn(V)))^4*(V-E.K) + G.Cl*(V-E.Cl);

function val = an(V)
val = .01*(10-(V+71))./(exp(1-(V+71)/10)-1);

function val = bn(V)
val = .125*exp(-(V+71)/80);

function val = am(V)
val = .1*(25-(V+71))./(exp(2.5-(V+71)/10)-1);

function val = bm(V)
val = 4*exp(-(V+71)/18);

function val = ah(V)
val = 0.07*exp(-(V+71)/20);

function val = bh(V)
val = 1./(exp(3-(V+71)/10)+1);


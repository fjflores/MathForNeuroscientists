%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% staggered Euler on the active and qusi active sphere with h current
%
% usage:  stEqah(dt,Tfin)
%
% example:  stEqah(0.01,40)
%

function stEqah(dt,Tfin)

close all

A = 4*pi*(1e-6);% (cm)^2

%E = struct('K', -77, 'Na', 56, 'Cl', -68, 'h', -55); % reversal potentials, mV
E = struct('K', -77, 'Na', 56, 'Cl', -68, 'h', -40); % reversal potentials, mV

%G = struct('K', 36, 'Na', 120, 'Cl', 0.3, 'h', 10);  % channel conductances, mS/cm^2
G = struct('K', 36, 'Na', 120, 'Cl', 0.3, 'h', 1);  % channel conductances, mS/cm^2

[Vr,B] = hhsymh(G)

Nt = ceil(Tfin/dt);

stc = 0;

%for I0 = 10:10:40,
for I0 = -4:-4:-16,
 
stc = stc + 1;

t = zeros(Nt,1); V = t; Ih = t;
V(1) = Vr;
n = an(Vr)/(an(Vr)+bn(Vr));  
m = am(Vr)/(am(Vr)+bm(Vr));  
h = ah(Vr)/(ah(Vr)+bh(Vr));  
q = qinf(Vr);
Ih(1) = G.h*q^2*(Vr-E.h);

y = zeros(5,Nt);
f = zeros(5,1);

td = 1/dt;
B1 = 2*td*eye(5) - B;
B2 = 2*td*eye(5) + B;

Istim0 = 0;
t1 = 2; t2 = 20;

for j=2:Nt,

      t(j) = (j-1)*dt;

      Istimid = I0*(t(j)-dt/2>t1)*(t(j)-dt/2<t2)*(1e-6);
      
      a = an(V(j-1));  b = bn(V(j-1)); c = (a+b)/2;
      n = ( (td-c)*n + a ) / (td + c);

      cK = G.K*n^4;

      a = am(V(j-1));  b = bm(V(j-1)); c = (a+b)/2;
      m = ( (td-c)*m + a ) / (td + c);

      a = ah(V(j-1));  b = bh(V(j-1)); c = (a+b)/2;
      h = ( (td-c)*h + a ) / (td + c);

      cNa = G.Na*m^3*h;

      temp = tauq(V(j-1));
      q = ( (2*temp - dt)*q + 2*qinf(V(j-1))*dt ) / (2*temp + dt);

      ch = G.h*q^2;

      top = 2*V(j-1)*td + cK*E.K + cNa*E.Na + ch*E.h + G.Cl*E.Cl + Istimid/A;

      bot = 2*td + cK + cNa + ch + G.Cl;

      Vmid = top/bot;

      V(j) = 2*Vmid - V(j-1);

      Ih(j) = G.h*q^2*(V(j)-E.h);

      Istim1 = I0*(t(j)>t1)*(t(j)<t2)*(1e-6);

      f(1) = (Istim0 + Istim1)/A;

      y(:,j) = B1 \ ( B2*y(:,j-1) + f );

      Istim0 = Istim1;

end


subplot(2,2,stc)
plot(t,V,'k')
hold on
plot(t',Vr+y(1,:),'r')
text(5,max(V)-0.2*(max(V)-min(V))/2,['I_0 = ' num2str(I0) ' pA'])
axis tight
box off
hold off


end

subplot(2,2,3)
xlabel('t  (ms)','fontsize',14)
ylabel('V  (mV)','fontsize',14)


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

function val = qinf(V)
val = 1/(1+exp((V+69)/7.1));

function val = tauq(V)
val = 100*10/(exp((V+66.4)/9.3)+exp(-(V+81.6)/13));

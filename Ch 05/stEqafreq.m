
%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% staggered Euler on the active cell, spike frequency
%
% usage:    stEqafreq(dt,Tfin,I0)
%
% where:    dt = timestep        
%           Tfin = final time
%           I0 = current stimulus
%
% example:  stEqafreq(0.1,100,1)
%

function stEqafreq(dt,Tfin,I0)

close all               % closes all open figures

A = 4*pi*(1e-6);        % (cm)^2

E = struct('K', -77, 'Na', 56, 'Cl', -68);  % reversal potentials, mV

G = struct('K', 36, 'Na', 120, 'Cl', 0.3);  % channel conductances, mS/cm^2

[Vr,B] = hhsym;         % find rest potential and quasiactive matrix

Nt = ceil(Tfin/dt);     % number of time steps

for k=1:95,

t = zeros(Nt,1); V = t; qv = t;     % preallocate space
V(1) = Vr;
n = an(Vr)/(an(Vr)+bn(Vr));  
m = am(Vr)/(am(Vr)+bm(Vr));  
h = ah(Vr)/(ah(Vr)+bh(Vr));  

y = zeros(4,Nt);            % preallocate space
f = zeros(4,1);             % preallocate space

td = 1/dt;
B1 = 2*td*eye(4) - B;
B2 = 2*td*eye(4) + B;

Istim0 = 0;

omega(k) = (k+4)/1000;

for j=2:Nt,

      t(j) = (j-1)*dt;

      Istimid = I0*sin(2*pi*(t(j)-dt/2)*omega(k))*(1e-6);
      
      a = an(V(j-1));  b = bn(V(j-1)); c = (a+b)/2;
      n = ( (td-c)*n + a ) / (td + c);

      cK = G.K*n^4;

      a = am(V(j-1));  b = bm(V(j-1)); c = (a+b)/2;
      m = ( (td-c)*m + a ) / (td + c);

      a = ah(V(j-1));  b = bh(V(j-1)); c = (a+b)/2;
      h = ( (td-c)*h + a ) / (td + c);

      cNa = G.Na*m^3*h;

      top = 2*V(j-1)*td + cK*E.K + cNa*E.Na + G.Cl*E.Cl + Istimid/A;

      bot = 2*td + cK + cNa + G.Cl;

      Vmid = top/bot;          % trapezoid scheme

      V(j) = 2*Vmid - V(j-1);

      Istim1 = I0*sin(2*pi*t(j)*omega(k))*(1e-6)/A;

      f(1) = Istim0 + Istim1;

      y(:,j) = B1 \ ( B2*y(:,j-1) + f );

      qv(j) = y(1,j);

      Istim0 = Istim1;

end


V = V-Vr;

vamp(k) = max(abs(V(ceil(end/2):end))); 

qvamp(k) = max(abs(qv(ceil(end/2):end))); 

end

figure(2)
plot(omega*1000,vamp/(I0*1e-6)/1000,'k')
hold on
plot(omega*1000,qvamp/(I0*1e-6)/1000,'r')
xlim([10 100])
xlabel('\omega  (Hz)','fontsize',14)
ylabel('R_{in}  (M\Omega)','fontsize',14)
box off

return


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

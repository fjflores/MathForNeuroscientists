%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% ml1.m
%
% 1 Morris Lecar cells with and without tonic inhibition
%
%  usage   ml1(dt, Tfin, IsA)      IsA  = nA/cm^2
%
%  example    ml1(.001,8,150,0)
%

function ml1(dt, Tfin, IsA, gsyn)

Cm = 1;
g = struct('K',20,'Ca',15,'Cl',5);
E = struct('K',-80,'Ca',100,'Cl',-50,'syn',-80);

Nt = ceil(Tfin/dt);
icnt = 0;

%for IsA = 2300:5:2900
for IsA = 10:1:700

up = 0;
vth = -20;
cnt = 0;
spst = [];

icnt = icnt + 1;

V = -50;
w = minf(V);

for j=2:Nt,

    tmpt = tauw(V);
    tmpm = minf(V);
  
    w = (tmpt.*w + tmpm*dt)./(dt + tmpt);

    top = Cm*V + dt*(g.Ca*tmpm*E.Ca + g.K*w*E.K + g.Cl*E.Cl + gsyn*E.syn + IsA);
    bot = Cm + dt*(g.Ca*tmpm + g.K*w + g.Cl + gsyn);

    V = top./bot;
 
    if V > vth
       if up == 0
           up = 1;           % mark upward crossing
           cnt = cnt + 1;
           spst(cnt) = (j-1)*dt;
       end
   else
       up = 0;
   end

end

isi = inf;
if cnt > 1
   isi = diff(spst);
   isi = isi(end)
end

Iinj(icnt) = IsA/1000;
freq(icnt) = 1000/isi;

end % for IsA

plot(Iinj,freq,'k')
box off
set(gca,'tickdir','out')
xlabel('I_{stim}   (\mu A/cm^2)','fontsize',14)
ylabel('Spike Frequency  (Hz)','fontsize',14)

return

function val = minf(V)
val = (1 + tanh(V/15))/2;
return

function val = tauw(V)
val = 1./cosh(V/30);
return

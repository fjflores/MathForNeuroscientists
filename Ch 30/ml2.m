%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% ml2.m
%
% Voltage traces and phase plane for
% 2 Morris Lecar cells that inhibit one another
%
%  usage   ml2(dt, Tfin, IsA)      IsA  = nA/cm^2
%
%  examples    ml2(.001,8,550)     % fig 27.21
%              ml2(.001,8,1550)    % fig 27.22
%              ml2(.001,8,2550)    % fig 27.23
%              ml2(.001,8,3050)    % fig 27.24
%

function ml2(dt, Tfin, IsA)

Cm = 1;
g = struct('K',20,'Ca',15,'Cl',5,'syn',30);
E = struct('K',-80,'Ca',100,'Cl',-50,'syn',-80);

Nt = ceil(Tfin/dt);

V = zeros(Nt,2); w = V;

V(1,:) = [-10 -20];
w(1,:) = minf(V(1,:));

for j=2:Nt,

    tmpt = tauw(V(j-1,:));
    tmpm = minf(V(j-1,:));
    tmps = sinf(fliplr(V(j-1,:)));
  
    w(j,:) = (tmpt.*w(j-1,:) + tmpm*dt)./(dt + tmpt);

    top = Cm*V(j-1,:) + dt*(g.Ca*tmpm*E.Ca + g.K*w(j,:)*E.K + ...
                            g.Cl*E.Cl + g.syn*tmps*E.syn + IsA*[1 1]);
    bot = Cm + dt*(g.Ca*tmpm + g.K*w(j,:) + g.Cl + g.syn*tmps);

    V(j,:) = top./bot;

end

t = linspace(0,Tfin,Nt);

figure(1)
plot(t',V(:,1),'k')
hold on
plot(t',V(:,2)-140,'r')
plot(t,zeros(1,Nt),'k--')
plot(t,zeros(1,Nt)-140,'r--')
text(Tfin-1,20,'V_1')
text(Tfin-1,20-180,'V_2')
set(gca,'tickdir','out')
hold off
axis tight
axis off

figure(2)
plot(V(:,1),w(:,1),'k')
hold on
plot(V(:,2),w(:,2),'r')
v = -80:.1:80;
plot(v,minf(v),'k--')
w = (IsA - g.Ca*minf(v).*(v-E.Ca) - g.Cl*(v-E.Cl))./(g.K*(v-E.K));
plot(v,w,'k:')
w = (IsA - g.Ca*minf(v).*(v-E.Ca) - g.Cl*(v-E.Cl) - ...
     g.syn*(1+0*sinf(v)).*(v-E.syn))./(g.K*(v-E.K));
plot(v,w,'k-.')
ylim([-.1 1.1])
xlabel('V  (mV)','fontsize',14)
ylabel('n','fontsize',14)
set(gca,'tickdir','out')
box off
hold off

return

function val = minf(V)
val = (1 + tanh(V/15))/2;
return

function val = tauw(V)
val = 1./cosh(V/30);
return

function val = sinf(V)
val = (1 + tanh(V/15))/2;
return

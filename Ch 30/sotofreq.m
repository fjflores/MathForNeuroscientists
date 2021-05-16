%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  sotofreq.m
%
%  Determine the dependence of oscillator frequency on the
%  value of the Calcium target is the Soto model.
%

function sotofreq

C = 8000:250:10000;
freq = zeros(1,length(C));
for cnt=1:length(C)
    freq(cnt) = sotoloc(1,4e5,[1 2],C(cnt));
end
figure(5)
plot(C,freq,'k')
box off
set(gca,'tickdir','out')
xlabel('Calcium target, C  (\muC/cm^2)','fontsize',14)
ylabel('Frequency  (Hz)','fontsize',14)
return     

function freq = sotoloc(dt, Tfin, w0, C)

Cm = 1;
g = struct('K',15,'Ca',20,'Cl',5);
E = struct('K',-80,'Ca',100,'Cl',-50,'syn',-80);

Nt = ceil(Tfin/dt);

V = zeros(Nt,2); c = V; n = V; w = V;

V(1,:) = [-80 -20];
n(1,:) = ninf(V(1,:));
w(1,:) = w0;
c(1,:) = [4000 5000];

tauw = 35e3;
tauCa = 10e3;

for j=2:Nt,

    tmptaun = taun(V(j-1,:));
    tmpm = minf(V(j-1,:));
    tmpn = ninf(V(j-1,:));
    tmps = sinf(fliplr(V(j-1,:)));
  
    n(j,:) = (tmptaun.*n(j-1,:) + tmpn*dt)./(dt + tmptaun);

    c(j,:) = (c(j-1,:) - 1e-3*dt*g.Ca*tmpm.*(V(j-1,:)-E.Ca))/(1+dt/tauCa);
    
    w(j,:) = w(j-1,:)./(1-(dt/tauw)*(c(j,:)-C)/C);

    top = Cm*V(j-1,:) + dt*(g.Ca*tmpm*E.Ca + g.K*n(j,:)*E.K + g.Cl*E.Cl + w(j,:).*tmps*E.syn);
    bot = Cm +          dt*(g.Ca*tmpm +      g.K*n(j,:) +     g.Cl +      w(j,:).*tmps);

    V(j,:) = top./bot;  

end

V1max = max(V(end-2e3:end,1));
hi = find(V(end-2e3:end,1)>floor(V1max)-1);
higap = diff(hi);
freq = 1000/mean(higap(higap>10));

if 1==0
t = linspace(0,Tfin,Nt);
t = t'/1000;

figure(1)
plot(t,V(:,1),'k',t,V(:,2),'r')
legend('Cell 1','Cell 2','orientation','horizontal')
box off
ylim([-80 80])
set(gca,'tickdir','out')
xlabel('t  (s)','fontsize',14)
ylabel('V  (mV)','fontsize',14)

figure(2)
plot(V(:,1),n(:,1),'k',V(:,2),n(:,2),'r')
hold on
v = -80:.1:80;
plot(v,ninf(v),'k--')
n = -(g.Ca*minf(v).*(v-E.Ca) + g.Cl*(v-E.Cl))./(g.K*(v-E.K));
plot(v,n,'k:')
ylim([-.1 1.1])
hold off
box off
set(gca,'tickdir','out')
xlabel('V  (mV)','fontsize',14)
ylabel('n','fontsize',14)

figure(3)
plot(t,w(:,1),'k',t,w(:,2),'r')
legend('Cell 1','Cell 2')
box off
set(gca,'tickdir','out')
xlabel('t  (s)','fontsize',14)
ylabel('w  (mS/cm^2)','fontsize',14)

figure(4)
plot(t,c(:,1),'k',t,c(:,2),'r')
legend('Cell 1','Cell 2')
box off
set(gca,'tickdir','out')
xlabel('t  (s)','fontsize',14)
ylabel('c  (\muC/cm^2)','fontsize',14)
end

return

function val = minf(V)
val = (1 + tanh((V+10)/20))/2;
return

function val = taun(V)
val = 125./cosh(V/30);
return

function val = ninf(V)
val = (1 + tanh((V+10)/5))/2;
return

function val = sinf(V)
val = 1./(1+exp(-(V+58)/10));
return

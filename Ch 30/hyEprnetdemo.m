%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% hyEprnetdemo.m
%
% solve the pinsky-rinzel 2 compartment 2 cell CA3 net
%
%  usage   hyEprnetdemo(dt, Tfin) 
%
%  e.g.,   hyEprnetdemo(.01,100)
%

function A = hyEprnetdemo(dt,Tfin)

Nt = ceil(Tfin/dt);

N = 2;
A = [0 0; 1 0];

%N = 30;		Simple ring
%A = zeros(N);
%A(1,N) = 1;
%for i=2:N,
%    A(i,i-1) = 1;
%end

Vd = zeros(N,Nt);
Vs = zeros(N,Nt);
c = zeros(N,Nt);
n = c;
m = c;
q = c;
h = c;
s = c;
r = c;
e = ones(N,1);

Vs(:,1) = -4.6;
Vd(:,1) = -4.5;
h(:,1) = 0.999*e;
n(:,1) = 0.001*e;
s(:,1) = 0.009*e;
r(:,1) = 0.007*e;
q(:,1) = 0.01*e;
c(:,1) = 0.2;
x = 0*e;
y = 0*e;
Is = -.5*e;
Id = 0*e;
gsyn = zeros(N,Nt);

g = struct('L', 0.1,'Na', 30, 'K', 15, 'Ca', 10, 'AHP', 0.8, 'KC', 15, 'AMPA', 0.0045, 'NMDA', 0.014);
E = struct('Na', 120, 'Ca', 140, 'K', -15, 'L', 0, 'syn', 60);

gc = 2.1;
p = 0.5;
Cm = 3;
oodt = 1/dt;

spk = zeros(1,Nt);
Iampa = spk;
Inmda = spk;

for j=2:Nt,

    a = ah(Vs(:,j-1));  abm = (a+bh(Vs(:,j-1)))/2;
    h(:,j) = ((oodt-abm).*h(:,j-1) + a)./(oodt + abm); 

    a = an(Vs(:,j-1));  abm = (a+bn(Vs(:,j-1)))/2;
    n(:,j) = ((oodt-abm).*n(:,j-1) + a)./(oodt+abm);

    a = am(Vs(:,j-1));  
    m(:,j) = a./(a + bm(Vs(:,j-1)));

    a = as(Vd(:,j-1));  abm = (a + bs(Vd(:,j-1)))/2;
    s(:,j) = ((oodt - abm).*s(:,j-1) + a)./(oodt + abm);

    a = ar(Vd(:,j-1));  abm = (a + br(Vd(:,j-1)))/2;
    r(:,j) = ((oodt - abm).*r(:,j-1) + a)./(oodt + abm);

    a = aq(c(:,j-1));  abm = (a + bq(c(:,j-1)))/2;
    q(:,j) = ((oodt - abm).*q(:,j-1) + a)./(oodt + abm);

    chi = min(c(:,j-1)/250,1);

    e20 = Vs(:,j-1)>20;
    x = (x + dt*A*e20)/(1 + dt/2);
    spk(j) = sum(e20);

    e10 = Vs(:,j-1)>10;
    y = (y + dt*A*e10)/(1 + dt/150);
    y = min(y,125);

    gsyn(:,j) = (g.AMPA*x + g.NMDA*y.*M(Vd(:,j-1)))/(1-p);

    t = (j-1)*dt;

    for i=1:N

        A2 = [Cm + dt*(g.L + g.Na*m(i,j)^2*h(i,j) + g.K*n(i,j) + gc/p)  -dt*gc/p
             -dt*gc/(1-p)  Cm + dt*(g.L + g.Ca*s(i,j)^2 + g.KC*chi(i)*r(i,j) + g.AHP*q(i,j) + gc/(1-p) + gsyn(i,j))];

   
        rhs = [Cm*Vs(i,j-1) + dt*(g.L*E.L + g.Na*m(i,j)^2*h(i,j)*E.Na + g.K*n(i,j)*E.K + (Is(i)+10*(i==1)*(t>10)*(t<13))/p)  
        Cm*Vd(i,j-1) + dt*(g.L*E.L + g.Ca*s(i,j)^2*E.Ca + (g.KC*chi(i)*r(i,j)+g.AHP*q(i,j))*E.K + Id(i)/(1-p) + gsyn(i,j)*E.syn)];

        V = A2 \ rhs;

        Vs(i,j) = V(1);    Vd(i,j) = V(2);

    end

    ICa = g.Ca*(s(:,j).^2).*(Vd(:,j)-E.Ca);

    c(:,j) = (c(:,j-1) - dt*0.13*ICa)/(1 + dt*0.075);

    Iampa(j) = g.AMPA*x(2)*(Vd(2,j)-E.syn);
    Inmda(j) = g.NMDA*y(2)*M(Vd(2,j))*(Vd(2,j) - E.syn);

end

t = linspace(0,Tfin,Nt);
%figure(2)
%plot(t,spk)

figure(1)
plot(t,Vs(1,:),'k')
hold on
plot(t,Vs(2,:),'r')
hold off
box off
set(gca,'tickdir','out')
legend('Cell 1','Cell 2')
xlabel('time  (ms)','fontsize',16)
ylabel('Soma Potential  (mV)','fontsize',16)

figure(2)
plot(t,Iampa,'k')
hold on
plot(t,Inmda,'r')
hold off
box off
set(gca,'tickdir','out')
legend('I_{AMPA}','I_{NMDA}')
xlabel('time  (ms)','fontsize',16)
ylabel('\mu A/cm^2','fontsize',16)

if 1==0
figure(2)
plot(t,c)

figure(3)
plot(t,gsyn)

figure(4)
plot(t,r)

figure(5)
plot(t,s)
end

return

function val = am(v)
val = 0.32*(13.1-v)./(exp((13.1-v)/4)-1);

function val = bm(v)
val = 0.28*(v-40.1)./(exp((v-40.1)/5)-1);

function val = an(v)
val = 0.016*(35.1-v)./(exp((35.1-v)/5)-1);

function val = bn(v)
val = 0.25*exp(0.5-0.025*v);

function val = ah(v)
val = 0.128*exp((17-v)/18);

function val = bh(v)
val = 4./(1+exp((40-v)/5));

function val = as(v)
val = 1.6./(1+exp(-0.072*(v-65)));

function val = bs(v)
val = 0.02*(v-51.1)./(exp((v-51.1)/5)-1);

function val = ar(v)
val = (v<=50).*exp((v-10)/11-(v-6.5)/27)/18.975 + (v>50)*2.*exp((6.5-v)/27);

function val = br(v)
val = (v<50).*(2*exp((6.5-v)/27) - exp((v-10)/11-(v-6.5)/27)/18.975);

function val = aq(Ca)
val = min((0.00002)*Ca,0.01);

function val = bq(v)
val = 0.001;

function val = M(v)
val = 1./(1+0.28*exp(-0.062*(v-60)));

%
% Gabbiani and Cox, Mathematics for Neuroscientists
%
% NaKrest
%
%  resting levels of the sodium/potassium pump
%
% usage:  KaKrest
%
function NaKrest
Ke = 4; % mM
Ne = 116; % mM
S = 308; % mM
PK = 5e-5; % cm/s
PNa = 5e-7; % cm/s
rho = 3/2;
ind = 0;
v = -1;
for j=-7:.1:-3,
    ind = ind + 1;
    lam(ind) = 10^j;
    v = fsolve(@(v) veq(v,lam(ind),PK,PNa,Ke,Ne,S),v);
    vbar(ind) = v;
end

Ni = Ne*vbar./((lam/PNa).*(exp(vbar)-1) + vbar.*exp(vbar));
Ki = Ke./exp(vbar) + Ni.*(lam/PK/rho).*(exp(vbar)-1)./vbar./exp(vbar);
Cli = exp(vbar)*(Ne+Ke);   % mM
w = (1.55e-6)./(Ni + Ki - Cli);   % mL

figure(1)
[hax,hL1,hL2] = plotyy(lam,vbar*25.8,lam,w,'semilogx')
set(hax(1),'ytick',-100:20:0)
%set(hax(2),'ytick',0:5:30)
xlabel('\lambda  (cm/s)','fontsize',14)
ylabel(hax(1),'V  (mV)','fontsize',14)
ylabel(hax(2),'w  (cm^3)','fontsize',14)
text(1e-4,1.8e-8,'(A)','fontsize',20,'parent',hax(2))

figure(2)
semilogx(lam,Ni,'k')
hold on
semilogx(lam,Ki,'r')
semilogx(lam,Cli)
hold off
text(1e-4,120,'(B)','fontsize',20)
legend('[Na]_i','[K]_i','[Cl]_i','location','best')
xlabel('\lambda  (cm/s)','fontsize',14)
ylabel('mM','fontsize',14)

return

function val = veq(v,lam,PK,PNa,Ke,Ne,S)
a = exp(v);
b = 1-1/a;
rho = 3/2;
val = a*(S/2) - Ne*(v + (lam/PK/rho)*b)./(v + (lam/PNa)*b) - Ke;




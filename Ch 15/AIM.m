%
% Gabbiani & Cox
%
%  AIM.m   The Arterial Impulse Model of Kim & Ress, 
%          NeuroImage 124(2016) 394-408
%
%  this code is adapted from their   AIMboldModel.m
%  and requires their data file      Fit.mat
% 

function AIM

close all
load Fit

% Calculate state-state solution for initial condition

P = 80; 
Qb00 = 8.56; 
sigma = 0.1;  
E = 0.428; 
alpha = 0.0178; 
U0 = 7e-4; 
k1 = -1.0430e+3;   
k2 = 197.7358;  
zetac = 0.0135*0.024/0.0135;  % fix a mismatch in original
Hc = 0.2;
Nzc = 80*8;  % number of points in capillary

nuCo = Hc/(1-Hc);
options = optimset('disp','off');
Qp0 = fsolve(@(q) Hill(q,Qb00),0.5,options);
QpL = fsolve(@(q) Hill(q,Qb00*(1-E)),0.5,options);
Gamma0 = sigma*(Qp0*(1-Hc)*P*zetac); % *0.024/0.0135
Lc = -U0*zetac*(1-Hc)*( nuCo*k1*(QpL^2-Qp0^2)+(1+nuCo*k2)*(QpL-Qp0) )/Gamma0;
dzc = Lc /( (Nzc - 1));   % spatial step size in capillary
z = 0:dzc:Lc;              % space points in capillary
% acoef = nuCo*k1;
% bcoef = 1 + nuCo*k2;
% ccoef = -Gamma0/(U0*zetac*Hc*k1);
% mQp0 = -bcoef/(2*acoef)-sqrt(ccoef*z+(Qp0+bcoef/(2*acoef))^2); 
aa = Hc*k1;
bb = 1+Hc*(k2-1);
cc = (Gamma0/U0/zetac)*z - aa*Qp0^2 - bb*Qp0;
mQp0 = (-bb+sqrt(bb^2-4*aa*cc))/2/aa;
mQe0 = mQp0 - Gamma0/(zetac*(1-Hc)*P);

% Extract the experimental data

vMeasure= Fitdata.goal;
vTime = Fitdata.time;
vupperCI = Fitdata.upperCI;
vlowerCI = Fitdata.lowerCI;
maxT = ceil(max(vTime));

% step calculation

dh = 0.2; 
dt = dh*dzc/U0; % Set up time step so that dt << dz/U0
time = 0:dt:maxT+0.1; % add 0.1 to maxT to cover full range of time period for interpolation.
Nt = length(time);

% Gamma variate function for CMRO2

gamma2 = 0.0043*0.024/0.0135; %*0.0135/0.024;
Fa = 2.089;
Fb = 1.3236;
stimDuration = 2;

nGamma = round(stimDuration/dt);
vGamma2 = ones(nGamma,1);
bGamma = zeros(Nt,1);
bGamma(1:nGamma) = vGamma2;
kernel = time.^(Fa-1).*exp(-time/Fb);
Gamma1 = gamma2/sum(kernel);
vGamma = Gamma0*ones(Nt,1);
vGamma2F = Gamma1*conv(bGamma, kernel');
vGamma2F = vGamma2F(1:Nt);
vGammaF = vGamma2F + vGamma;
CMRO2percent = (vGammaF/vGammaF(1)-1)*100;

% Underdamped sinusoidal function for CBF

freq = 0.0570;
tau = 3.7816;
amp = 0.2089;
vU = sin(2*pi*freq.*time).*exp(-time/tau); 
nt = round(2/dt);
vU = conv(vU, ones(1,nt)/nt);
vU = vU(1:Nt);
U1 = amp*1e-3/max(vU);
vU = U1*vU; 
CBFpercent= vU/(U0)*100;

Qe = zeros(Nzc,Nt);
Qe(:,1) = mQe0; 
Qe(1,:) = mQe0(1);

Lvein = 0.002; 
Nzv = Nzc*4;
dzv = Lvein /(Nzv - 1);
Nz = Nzc + Nzv - 1;
z = [z Lc+dzv:dzv:Lc+Lvein]*10^3;  % length in mm

Qp = zeros(Nz,Nt);
Qp(:,1) = [mQp0'; mQp0(end)*ones(Nzv-1,1)];
Qp(1,:) = mQp0(1);

vU = U0 + vU;
Uvein = 0.0025; 
Uvein = vU*(Uvein/U0);

fac1 = alpha*dt*P*(1-Hc);
fac2 = dt/dzv;
fac3 = P*(1-Hc)*dzc*dt;

figure(1)
SO2 = 100*Hill(Qp(:,1),0)/9.2;
plot3(z,zeros(Nz,1),SO2,'k')
hold on
pinc = 1000;

for j=2:Nt,   % solve the transport equations
    
    Qe(2:end,j) = (Qe(2:end,j-1) + fac1*Qp(2:Nzc,j-1) - alpha*dt*vGammaF(j)/zetac)/(1+fac1);
    
    for k=2:Nzc,   %  Capillary
        
        beta = fac3/( 1+ Hc*(2*k1*Qp(k,j-1) + k2 - 1) );   
        num  = Qp(k,j-1)*dzc + vU(j)*dt*Qp(k-1,j) + beta*Qe(k,j);
        den = dzc + vU(j)*dt + beta;
        Qp(k,j) = num/den;
        
    end
    
    for k=Nzc+1:Nz   %  Vein
        
        Qp(k,j) = (Qp(k,j-1) + fac2*Uvein(j)*Qp(k-1,j))/(1 + fac2*Uvein(j));
        
    end
    
    if mod(j,pinc) == 0
        SO2 = 100*Hill(Qp(:,j),0)/9.2;
        plot3(z,dt*(j-1)*ones(Nz,1),SO2,'k')
    end
           
end

hold off
xlabel('z  (mm)','fontsize',14)
ylabel('t  (s)','fontsize',14)
zlabel('SO_2  (%)','fontsize',14)

mQp = Qp(1:Nzc,:);
mQe = Qe;
mQpVein = Qp(Nzc:end,:);
vModelTime = time;
vModelT = time;

% Sampling and Interpolation 
mQpNew = mQp(:,1:100:end);
mQeNew = mQe(:,1:100:end);
mQpNewVein = mQpVein(:,1:100:end);
vNewMtime = vModelTime(1:100:end);
CBF = interp1(vModelT,CBFpercent, vNewMtime,'linear');
CMRO2 = interp1(vModelT, CMRO2percent, vNewMtime, 'linear');

% Initialize BOLD parameters

CstarCap = 102.7;
AstarCap = 17.66;
vc = 0.4*0.06;   % zetac above!
CstarVein = 174.7;
AstarVein = 21.23;
vv = 0.4*0.06;
dchi = 2.64e-7;
HtVein = 0.44;
dg = 2.68e8;
B0 = 3;
SO2ref = 0.95;
TE = 0.025;   % s
lambda = 1.15;
R2sExt = 25.1;   % 1/s

[zl,tl] = size(mQpNew);
TotalBOLD = zeros(1,tl);

% Resting-state parameters

vSO2rest = Hill(mQpNew(:,1),0)./9.2;
SigExt0 = exp(-TE*R2sExt);
ec = lambda*exp(-TE*(AstarCap+CstarCap*(1-vSO2rest).^2))/SigExt0;
ev = lambda*exp(-TE*(AstarVein+CstarVein*(1-vSO2rest(end)).^2))/SigExt0;
mSO2actVein1 = Hill(mQpNewVein(:,1),0)./9.2;

for it = 1:tl
    
  mSO2act = Hill(mQpNew(:,it),0)./9.2;
  mSO2actVein = Hill(mQpNewVein(:,it),0)./9.2;
  
  % Intravascular capillaries
  deltaR2Ic = CstarCap*((1-mSO2act).^2-(1-vSO2rest).^2);
  vIntraC = (ec*vc)./(1-vc+ec*vc).*exp(-TE*deltaR2Ic);
  IntraCap = nanmean(vIntraC);
  
  % Intravascular veins
  deltaR2Iv = CstarVein*((1-mSO2actVein).^2-(1-mSO2actVein1).^2);
  IntraV = (ev*vv)./(1-vv+ev*vv).*exp(-TE*deltaR2Iv);
  IntraVein = nanmean(IntraV);
  
  % Extravascular capillaries
  deltaR2Ec = 0.04*vc*(dchi*Hc*dg*B0)^2*((abs(SO2ref-mSO2act)).^2-(abs(SO2ref-vSO2rest)).^2);
  vExtraC = 1./(1+(ec*vc)./(1-vc)).*exp(-deltaR2Ec.*TE);
  ExtraCap = nanmean(vExtraC);
  
  % Extravascular vein
  deltaR2Ev = 4/3*pi*dchi*HtVein*dg*B0*vc*(abs(SO2ref-mSO2actVein)-abs(SO2ref-mSO2actVein1));
  ExtraVc = 1./(1+((ev*vv)./(1-vv))).*exp(-deltaR2Ev.*TE);
  ExtraVein = nanmean(ExtraVc);
  
  % Total BOLD signal
  TotalBOLD(it) = (IntraCap + IntraVein + ExtraCap + ExtraVein) - 2;
  
end

Bolddata = TotalBOLD*100;

% Interpolation for measured BOLD HRF

vMeasureInterp = interp1(vTime,vMeasure, vNewMtime,'linear');
vupperCIinterp = interp1(vTime,vupperCI, vNewMtime,'linear');
vlowerCIinterp = interp1(vTime,vlowerCI, vNewMtime,'linear');

% Plots

figure(2)
plot(vNewMtime, vMeasureInterp, 'k', vNewMtime, Bolddata, 'r', 'LineWidth', 1.5);
legend('Experimental Data','Arterial Impulse Model','Location','NorthEast');
hold on;
plot(vNewMtime, vupperCIinterp, 'k--', 'LineWidth', 1.5);
plot(vNewMtime, vlowerCIinterp, 'k--', 'LineWidth', 1.5);
xlabel('time (s)','fontsize',14)
ylabel('BOLD Signal (%)','fontsize',14)
box off
axis([0 27 -0.4 1.2])

figure(3)
plot(vNewMtime, CBF, 'k', vNewMtime, CMRO2, 'r', 'LineWidth', 1.5);
legend('CBF','CMRO_2','Location','NorthEast');
xlabel('time (s)','fontsize',14)
ylabel('Response (%)','fontsize',14)
box off
axis([0 27 -5 35])

end

function val = Hill(q,rhs)
h = 2.73;
c = (26*1.39e-3)^h;
q = q.^h;
val = 9.2*q./(q+c) - rhs;
end

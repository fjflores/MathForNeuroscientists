%% Define parameters
%gyromagnetic ratio of H+
gamma = 42.58; %MHz/T
gamma1 = gamma*1e6; %convert to Hz/T
B0 = 3; %assume 3T
Hred = 0.4; %hematocrit level

%difference between oxy and deoxy susceptibility
chi_t = 3.4306e-6;

Y = 0.6; %blood oxygenation level

%assume 4 percent volume fraction for blood 
xi = 0.04;

%measurement time assume 80 ms
te = 80e-3; %s

%Compute omega_c for current parameters
omega_c = gamma1*B0*Hred*(1-Y)*chi_t/2; %1/s

%maximal percent change in signal
a = omega_c *te*xi;

%%

%value of Ances et al. 2008
%a = 0.06; 

%cbf from baseline (=1) to 100 percent increase or 2
f = 1:0.01:2;

%slope of relation between relative CBF and CMRO2 changes
sl = 2;
%assume (f-1)/(m-1) = sl
m = 1 + (f-1)/sl;

%parameters of the Davies model. Note: usually beta = 1.5. Here beta = 1
%match equation for BOLD
alpha = 2/5;
beta = 1.0;
b = a*(1- (m.^beta).*(f.^(alpha-beta)));

hf = figure; 
ha = axes;
%in plot convert to percentage change in CBF and percentage BOLD
hp1 = plot((f-1)*100,b*100,'k');
hold on;

%repeat for a slope of 3
sl = 3;
m = 1 + (f-1)/sl;
b = a*(1- (m.^beta).*(f.^(alpha-beta)));
hp2 = plot((f-1)*100, b*100,'r');
set(ha,'TickDir','out');
legend([hp1 hp2],{'n = 2','n = 3'});
xlabel('CBF (% change)');
ylabel('BOLD (% change)');
%print(hf,'figs_fg/bold_curve.eps','-depsc');
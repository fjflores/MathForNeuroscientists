%% Speed tuning of integrated response low-pass filter case

%low-pass filter time constants in s
tau_s = 20e-3;
tau_f = 2e-3;

%angular separation of the facets in deg
dphi = 2; 

%scaling factor for pulse; in our case no scaling
alpha = 1;

%speed range in deg/s
v = 0.1:0.1:10e4;

%theoretically derived formula
R = (alpha^2/(tau_f + tau_s))*(exp(-dphi./(v*tau_s))-exp(-dphi./(v*tau_f)));

Rpeak = max(R);

%show effect of doubling slow time constant
tau_s2 = 2* tau_s;
Rs2 = (alpha^2/(tau_f + tau_s2))*(exp(-dphi./(v*tau_s2))-exp(-dphi./(v*tau_f)));

%rescale so as to have the same peak 
Rs2rs = (Rpeak*Rs2)/max(Rs2);

tau_f2 = tau_f/2;
Rf2 = (alpha^2/(tau_f2 + tau_s))*(exp(-dphi./(v*tau_s))-exp(-dphi./(v*tau_f2)));

%rescale so as to have the same peak 
Rf2rs = (Rpeak*Rf2)/max(Rf2);

h_f = figure; 
h_a = axes;
h_p1 = semilogx(v,R,'Color','k');
hold on;
h_p2 = semilogx(v,Rs2rs,'Color','r','LineStyle','--');
h_p3 = semilogx(v,Rf2rs,'Color','r','LineStyle',':');
title('Standard correlator response to translating pulse stimuli');
xlabel('translation speed (deg/s)');
ylabel('integrated correlator output (s)');
set(h_a,'TickDir','out');
set(h_a,'XLim', [10 1e4]);
legend([h_p1 h_p2 h_p3],{'standard', 't_s 2', 't_f/2'});

%print(h_f,'figures/stdhrc_lp_d_tuning.eps','-depsc');

%% Speed tuning of integrated response high-pass filter case

%high pass filter time constant in s 
tau_h = 50e-3;

%theoretically derived formula
R2 = ( (alpha^2/(tau_s + tau_h))*exp(-dphi./(v*tau_h)) ) + ...
    ( (alpha^2*tau_h/(tau_s*(tau_s+tau_h)))*exp(-dphi./(v*tau_s)) );

h_f = figure; 
h_a = axes;
semilogx(v,R2,'Color','k');
xlabel('translation speed (deg/s)');
ylabel('integrated correlator output (s)');
set(h_a,'TickDir','out');
set(h_a,'XLim', [10 1e4]);

%print(h_f,'figures/stdhrc_hp_d_tuning.eps','-depsc');

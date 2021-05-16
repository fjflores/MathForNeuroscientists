function  std_corr_sw
%Computes the correlation model output for a translating sinewave stimulus. 
%Compares numerical stimulation with computed steady-state theoretical solution 

%phase of the cosine (see below) 
phi = 0;

alpha = 0.8;
omega_x = 0.05; %c/deg, corresponding to a 20 deg period

t_start = 0;
t_end = 2.0; %s
dt = 0.05e-3; %s
t_vect = t_start:dt:t_end;

v = 10; %deg/s

%position of the first sampling station
x0 = 0;

lum_vect0 = alpha*cos(2*pi*(omega_x*x0-omega_x*v*t_vect+phi));

%distance to the second sampling station
dx = 2; %deg
x1 = x0+dx;
lum_vect1 = alpha*cos(2*pi*(omega_x*x1-omega_x*v*t_vect+phi));

%compute responses of the 1st order low-pass filters
tau_s = 20e-3; %20 ms
tau_f = 2e-3; %2 ms
lv0_s_num = lp1filt_fn(lum_vect0,dt,tau_s);
lv0_f_num = lp1filt_fn(lum_vect0,dt,tau_f);
lv1_s_num = lp1filt_fn(lum_vect1,dt,tau_s);
lv1_f_num = lp1filt_fn(lum_vect1,dt,tau_f);

%half correlator responses
r_hc1 = lv0_s_num.*lv1_f_num;
r_hc2 = lv1_s_num.*lv0_f_num;

%full opponent response
r_hr = r_hc1 - r_hc2;

%theoretical steady-state response
g_s = (1/tau_s)/sqrt(1/tau_s^2 + (2*pi*omega_x*v)^2);
g_f = (1/tau_f)/sqrt(1/tau_f^2 + (2*pi*omega_x*v)^2);
psi_s = atan(2*pi*omega_x*v*tau_s);
psi_f = atan(2*pi*omega_x*v*tau_f);
r_ss_th = alpha^2*g_s*g_f*sin(2*pi*omega_x*dx)*sin(psi_s-psi_f);

%temporal frequency tuning curve
omega_t = 0:0.1:1e3;
R_ss = alpha^2*2*pi*omega_t*(tau_s-tau_f)*sin(2*pi*omega_x*dx)./...
    ((1+(2*pi*omega_t*tau_s).^2).*(1+(2*pi*omega_t*tau_f).^2));

%plot half-correlator responses
h_f = figure;
h_a = axes;
hp1 = plot(t_vect,r_hc1,'k');
hold on;
hp2 = plot(t_vect,r_hc2,'r');
legend([hp1 hp2],{'preferred direction','null'});
xlabel('time (s)');
ylabel('half-correlator responses');
set(h_a,'TickDir','out');
%print(h_f,'figures/std_corr_sw1.eps','-depsc');

%plot full correlator response
h_f = figure;
h_a = axes;
hp3 = plot(t_vect,r_hr,'k');
hold on;
hp4 = plot([t_start t_end],r_ss_th*[1 1],'r--');
legend([hp3 hp4],{'numerical solution','steady-state theor.'});
xlabel('time (s)');
ylabel('correlation model response');
set(h_a,'TickDir','out');
%print(h_f,'figures/std_corr_sw2.eps','-depsc');

%plot temporal frequency tuning curve
h_f = figure;
h_a = axes;
semilogx(omega_t,R_ss,'Color','k');
xlabel('temporal frequency \omega_t (c/s)');
ylabel('steady-state correlation model response');
set(h_a,'TickDir','out');
%print(h_f,'figures/std_corr_sw3.eps','-depsc');

end

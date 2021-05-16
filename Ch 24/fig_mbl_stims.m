%% Stimulus generation code for 1-d microbalanced stimuli

%time interval, 0.1 ms
dt = 0.1e-3; %s

%time vector; could start at zero, but then tick labels in figures below
%would not be round numbers  
t_start = dt; 
t_end = 5; %s
t_vect = t_start:dt:t_end;
n_t = length(t_vect);

% partial randomization in time; in dt_rstep
% multiples of dt 
dt_rstep = 50; %this corresponds to 5 ms

%number of times we can generate dt_rstep identical copies of a single
%spatially randomized vect
rrep = floor(n_t/dt_rstep);
%remaining time steps to be filled
rrem = n_t-rrep*dt_rstep;

%stimulus parameters - time 
omega_t = 10; %temporal frequency, c/s

%stimulus parameters - space
lambda_s = 10; %spatial wavelength (deg)
omega_s = 1/lambda_s; %spatial frequency (c/deg)

%spatial sampling discretization
dx = 0.1;

%spacing of correlator inputs
dx_corr = 2; %deg
%spacing in terms of spatial sampling index step
di_corr = round(dx_corr/dx);

%spatial vector
x_max = 90;
x_vect = dx:dx:x_max; %same remark as above, starting at dx, not zero
n_x = length(x_vect);

%stimulus matrix: each line is a different time point
%and each column a different spatial position
stim = zeros(n_t,n_x);

%partial randomization in time, complete randomization in space
phi_m2 = zeros(n_t,n_x);
for rep_cnt = 1:rrep
    phi_int = randi([0 1],1,n_x);
    phi_m2(1+(rep_cnt-1)*dt_rstep:dt_rstep+(rep_cnt-1)*dt_rstep,1:n_x) = repmat(phi_int,dt_rstep,1);
end

%fill in the remaining time points
phi_int = randi([0 1],1,n_x);
phi_m2(1+rrep*dt_rstep:rrem+rrep*dt_rstep,1:n_x) = repmat(phi_int,rrem,1);
phi_m2(phi_m2==0) = -1;

%Generate Demonstration 2 stimulus of Chubb and Sperling, JOSAA 5:1986-2007, 1988 [CS88]. 
%Use addition rule to compute g_m(t) = cos(2 pi (alpha m/M - beta n/N)),
%see p. 1994 or [CS88]. cos(alpha-beta) = cos(alpha)cos(beta)+sin(alpha)sin(beta)

%space vectors 1 x n_x
c_x = cos(2*pi*omega_s*x_vect);
s_x = sin(2*pi*omega_s*x_vect);

%time vectors n_t x 1
c_t = cos(2*pi*omega_t*t_vect)';
s_t = sin(2*pi*omega_t*t_vect)';

%replicate n_t lines
c_x2 = repmat(c_x,n_t,1);
s_x2 = repmat(s_x,n_t,1);

%replicate over n_x columns
c_t2 = repmat(c_t,1,n_x);
s_t2 = repmat(s_t,1,n_x);

%g_m(t) according to def above
c_mat = c_x2.*c_t2+s_x2.*s_t2;

%final contrast matrix based eq. 3 of [CS88]  
stim = phi_m2.*(c_mat+1)/2;

%plot a fraction of the stimulus
h_f1 = figure;
ha = axes;
%maximal subimage that can be plotted at full resolution
imagesc(stim(1:800,1:800));
colormap(gray)
set(ha,'TickDir','out');
colorbar(ha,'TickDirection','out');

% Decorations
%clear large accessory data structures used in previous code section
clear c_x c_t s_x s_t c_x2 s_x2 c_t2 s_t2 c_mat;

%these are the indices in the stim matrix
tick_x = get(ha,'XTick');
pos_x = x_vect(tick_x); %convert to position
%convert to cell array of strings, one per tick
label_x = cellstr(num2str(pos_x'));
set(ha,'XTickLabel',label_x);
xlabel('position (deg)');

tick_y = get(ha,'YTick');
pos_y = t_vect(tick_y); %convert to position
%convert to cell array of strings, one per tick
label_y = cellstr(num2str(pos_y'));
set(ha,'YTickLabel',label_y);
ylabel('time (s)');

%print(h_f1,'figures/mbl_stims1.eps','-depsc','-r300','-painters');
%print(h_f1,'figures/mbl_stims1.tif','-dtiff','-r300','-painters');

%% Single correlation model output

ind0 = 1; %this corresponds to location dx
lum_vect0 = stim(:,ind0); % at x = dx
ind1 = ind0 + di_corr; %next sampling at index corresponding to 2 deg away 
lum_vect1 = stim(:,ind1); % at x = 2+dx

%compute responses of the 1st order low-pass filters
tau_s = 20e-3; % 20 ms
tau_f = 2e-3; % 2 ms
lv0_s_num = lp1filt_fn(lum_vect0,dt,tau_s);
lv0_f_num = lp1filt_fn(lum_vect0,dt,tau_f);
lv1_s_num = lp1filt_fn(lum_vect1,dt,tau_s);
lv1_f_num = lp1filt_fn(lum_vect1,dt,tau_f);

%half correlator responses
r_hc1 = lv0_s_num.*lv1_f_num;
r_hc2 = lv1_s_num.*lv0_f_num;

%full opponent response
r_hr_single = r_hc1 - r_hc2;

%time averaged response
m_r_hr_single = mean(r_hr_single);

%plot the correlator response
h_f2 = figure; 
h_a2 = axes;
h_p1 = plot(t_vect,r_hr_single,'k');
hold on;
h_p1a = plot([t_vect(1) t_vect(end)],m_r_hr_single*ones(1,2),'r');

xlabel('time (s)');
ylabel('correlation model response');

%% variance estimation

x_1 = x_vect(ind0);
x_2 = x_vect(ind1); %spacing of 2 deg

m1 = (1 + cos(2*pi*omega_s*x_1-2*pi*omega_t*t_vect))/2;
m2 = (1 + cos(2*pi*omega_s*x_2-2*pi*omega_t*t_vect))/2;

%delay estimated from pulse tuning of correlation model
dt_hr = 10e-3;

n_hr = dt_hr/dt;

t1 = m1(1:n_t-n_hr).^2.*m2(1+n_hr:n_t).^2;
t2 = m1(1+n_hr:n_t).^2.*m2(1:n_t-n_hr).^2;
ttot = t1 + t2;

scale_factor = max(r_hr_single)/max(ttot);
ttot_sf = scale_factor*ttot;

figure(h_f2);

%plot the estimated variance
plot(t_vect(1+n_hr:n_t),ttot_sf,'r:');
plot(t_vect(1+n_hr:n_t),-ttot_sf,'r:');

legend([h_p1 h_p1a],{'single','time-average'});
set(h_a2,'XLim',[1 3]);
set(h_a2,'TickDir','out');

%print(h_f2,'figures/mbl_stims2.eps','-depsc','-r300','-painters');

%% Full-wave rectified stimulus

%compute full-wave rectified response
stim_a = abs(stim);

%plot a fraction of the stimulus
h_f2a = figure;
ha2 = axes;
imagesc(stim_a(1:800,1:800));
colormap(gray)
set(ha2,'TickDir','out');
colorbar(ha2,'TickDirection','out');

%these are the indices in the stim matrix
tick_x = get(ha2,'XTick');
pos_x = x_vect(tick_x); %convert to position
%convert to cell array of strings, one per tick
label_x = cellstr(num2str(pos_x'));
set(ha2,'XTickLabel',label_x);
xlabel('position (deg)');

tick_y = get(ha2,'YTick');
pos_y = t_vect(tick_y); %convert to position
%convert to cell array of strings, one per tick
label_y = cellstr(num2str(pos_y'));
set(ha2,'YTickLabel',label_y);
ylabel('time (s)');

%print(h_f2a,'figures/mbl_stims3.eps','-depsc','-r300','-painters');
%print(h_f2a,'figures/mbl_stims3.tif','-dtiff','-r300','-painters');

%% Compute the response of a single correlation model
lum_vect0 = stim_a(:,ind0); % at x = 0
lum_vect1 = stim_a(:,ind1); % at x = 2

%compute responses of the 1st order low-pass filters
lv0_s_num = lp1filt_fn(lum_vect0,dt,tau_s);
lv0_f_num = lp1filt_fn(lum_vect0,dt,tau_f);
lv1_s_num = lp1filt_fn(lum_vect1,dt,tau_s);
lv1_f_num = lp1filt_fn(lum_vect1,dt,tau_f);

%half correlator responses
r_hc1 = lv0_s_num.*lv1_f_num;
r_hc2 = lv1_s_num.*lv0_f_num;

%full opponent response
r_hr_a = r_hc1 - r_hc2;

%time averaged response
m_r_hr_a = mean(r_hr_a);

%plot the rectified single correlator response
h_f3 = figure;
h_a3 = axes;
h_p3 = plot(t_vect,r_hr_a,'k');
hold on;
h_p3a = plot([t_vect(1) t_vect(end)],m_r_hr_a*ones(1,2),'r');

xlabel('time (s)');
ylabel('correlation model response');

%plot the averaged correlator response
figure(h_f3); 
legend([h_p3 h_p3a],{'single', 'time average'});
set(h_a3,'XLim',[1 3]);
set(h_a3,'TickDir','out');

%print(h_f3,'figures/mbl_stims4.eps','-depsc','-r300','-painters');

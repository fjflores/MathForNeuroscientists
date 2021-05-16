%% Generate figure with output of extended Barlow Levick inhibitory model
%
%Compute extended Barlow-Levick model outputs for a translating sinewave stimulus.
%Plot tuning of model with frequency

%set default value of the phase of the cosine wave, this gives a sine
phi = 1/4;

%peak contrast
alpha = 0.5;

omega_x = 0.05; %c/deg, corresponding to a 20 deg period

t_start = 0;
t_end = 5; %s
dt = 0.05e-3; %s
t_vect = t_start:dt:t_end;

%translation speed
v = 10; %deg/s

%temporal frequency
omega_t = omega_x*v;

%position of the first sampling station
x0 = 0;

lum_vect0 = alpha*cos(2*pi*(omega_x*x0-omega_t*t_vect+phi));

%distance to the second sampling station
dx = 2; %deg
x1 = x0+dx;
lum_vect1 = alpha*cos(2*pi*(omega_x*x1-omega_t*t_vect+phi));

plot_debug = 0;
if (plot_debug == 1)
    figure;
    hp1 = plot(t_vect,lum_vect0);
    hold on;
    hp2 = plot(t_vect,lum_vect1,'r:');
    legend([hp1 hp2],{'position 1','position 2'});
    xlabel('time (s)');
    ylabel('luminance');
    title('input signals at facets 1 and 2');
end

%% Responses of 1st order low-pass filters

%compute responses of the 1st order low-pass filters to standard stims
tau_s = 20e-3; %20 ms
tau_f = 2e-3; %2 ms
lv0_s_num = lp1filt_fn(lum_vect0,dt,tau_s);
lv0_f_num = lp1filt_fn(lum_vect0,dt,tau_f);
lv1_s_num = lp1filt_fn(lum_vect1,dt,tau_s);
lv1_f_num = lp1filt_fn(lum_vect1,dt,tau_f);

%theoretical steady-state response of correlator model
r_ss_th = sin(2*pi*omega_x*dx)*alpha^2*2*pi*omega_t*(tau_s-tau_f)/ ...
    ((1 +(2*pi*omega_t*tau_s)^2)*(1+(2*pi*omega_t*tau_f)^2));

%% And-Not response based on the de Polavieja Neural computation paper
%
%Baseline variant
%
%Model output should be a continuous approximation of 
%  [B - (B .~ A_lp)] - [A - (A .~ B_lp)]
% where B is input 2, A is input 1, _lp denotes the low-pass filtered
% implementation of the delay.

%extended bl-model half responses using and-not (an) gates
% 1 .~ 1 = 0, 1 .~ 0 = 1, 0 .~ 0 = 0, 0 .~ 1 = 0
% This an gate is approaximated by the H(A-B-1/2), where H is the 
% Heaviside function. This will be 1 only when A = 1, B = 0. 
% The continous approximation is with a smoothed Heaviside as below, based
% on an hyperbolic tangent.

% This scale factor implements how steep the smoothed Heaviside is, and
% thus how close it is to the real one. Value in paper is 5. Note: for a 
% value of 5, the following ang1/2 formula  can be different from zero
% for lv1_f_num-lv0_s_num < 0. This is not the case with a factor of 10.
heavi_scale = 10; 

% This implements B .~ A_lp
ang1 = (1+tanh(heavi_scale*((lv1_f_num-lv0_s_num)-0.5)))/2;

% This implements A .~ B_lp
ang2 = (1+tanh(heavi_scale*((lv0_f_num-lv1_s_num)-0.5)))/2;

if (plot_debug == 1)
    figure;
    h_p1 = plot(t_vect,ang1);
    hold on;
    h_p2 = plot(t_vect,ang2,'r');
    xlabel('time (s)');
    ylabel('Smoothed Heaviside output');
    legend([h_p1 h_p2], {'B .~ A_{lp}', 'A .~ B_{lp}'});
    title('Smoothed approximation to and-not');
end

%This implements B - (B .~ A_lp)
r_bl1_an = lv1_f_num - ang1;

%This implements A - (A .~ B_lp)
r_bl2_an = lv0_f_num - ang2;

%full opponent response
r_bl_an = r_bl1_an - r_bl2_an;

%t_min is used to get rid of potential initial transients. t_max is set so
%as to integrate over 1 period. Since omega_t = 0.5 c/s, this is 2 s.
t_min = 1;
t_max = 3;

%compute steady-state signal amplitude
sss_amp_an = mean(r_bl_an( (t_vect>=t_min) & (t_vect<t_max) ) );

%plot AN output
h_f = figure; 
h_ax = axes;
hp3 = plot(t_vect,r_bl_an,'k');
hold on;
hp4 = plot([t_start t_end],r_ss_th*[1 1],'r--');
hp5 = plot(t_vect,ones(size(t_vect))*sss_amp_an,'k--');
legend([hp3 hp4 hp5],{'eBL','HR steady-state','eBL time-averaged'});
xlabel('time (s)');
ylabel('correlator responses');
set(h_ax,'TickDir','out');

%print(h_f,'figures/bl_response1.eps','-depsc');

%% Plot the temporal frequency tuning of HR
omega_tv = 0:0.1:1e3;
R_ss = alpha^2*2*pi*omega_tv*(tau_s-tau_f)*sin(2*pi*omega_x*dx)./...
    ((1+(2*pi*omega_tv*tau_s).^2).*(1+(2*pi*omega_tv*tau_f).^2));

if ( plot_debug == 1)
    h_f1 = figure;
    h_a = axes;
    semilogx(omega_tv,R_ss,'Color','k');
    xlabel('temporal frequency \omega_t (c/s)');
    ylabel('steady-state correlation model response');
    set(h_a,'TickDir','out');
end

%% Performance of HR at select frequencies

t_end = 30; %s
t_vect = t_start:dt:t_end;

%select a temporal frequency to sample
omega_td = [0.2 0.4 1 2 4 8 10 14 20 30 40 60 100]; %c/s

%theoretical steady-state response of correlator model
rssth_v = sin(2*pi*omega_x*dx)*alpha^2*2*pi*omega_td*(tau_s-tau_f)./ ...
    ((1 +(2*pi*omega_td*tau_s).^2).*(1+(2*pi*omega_td*tau_f).^2));

if (plot_debug == 1)
    figure;
    plot(omega_td,rssth_v);
end

%% And-Not response based on the de Polavieja Neural computation paper
%
%Model output should be a continuous approximation of 
%  [B - (B .~ A_lp)] - [A - (A .~ B_lp)]
% where B is input 2, A is input 1, _lp denotes the low-pass filtered
% implementation of the delay.

%extended bl-model half responses using and-not (an) gates
% 1 .~ 1 = 0, 1 .~ 0 = 1, 0 .~ 0 = 0, 0 .~ 1 = 0
% This an gate is approaximated by the H(A-B-1/2), where H is the 
% Heaviside function. This will be 1 only when A = 1, B = 0. 
% The continous approximation is with a smoothed Heaviside as below, based
% on an hyperbolic tangent.

% This scale factor implements how steep the smoothed Heaviside is, and
% thus how close it is to the real one. Value in paper is 5. Note: for a 
% value of 5, the following ang1/2 formula  can be different from zero
% for lv1_f_num-lv0_s_num < 0. This is not the case with a factor of 10.
heavi_scale = 10; 

%t_min is used to get rid of potential initial transients in temporal 
%transients below. Average data above the t_min time value and below the 
%t_max value determined below. Average over an integer number of cycle
%periods to avoid biasing the result up or down. 
t_min = 1;

%temporal cycle period
lambda_td = 1./omega_td;

%find the maximum time, as a multiple of the cycle period, with as many
%periods as will enter in t_end minus t_min to make sure we don't go over
%t_end
t_max =  floor((t_end-t_min)*omega_td).*lambda_td + t_min;

%i = 1; 
for i = 1:length(omega_td)
    %input luminance patterns
    lum_vect0 = alpha*cos(2*pi*(omega_x*x0-omega_td(i)*t_vect+phi));
    lum_vect1 = alpha*cos(2*pi*(omega_x*x1-omega_td(i)*t_vect+phi));
    
    %compute responses of the 1st order low-pass filters
    lv0_s_num = lp1filt_fn(lum_vect0,dt,tau_s);
    lv0_f_num = lp1filt_fn(lum_vect0,dt,tau_f);
    lv1_s_num = lp1filt_fn(lum_vect1,dt,tau_s);
    lv1_f_num = lp1filt_fn(lum_vect1,dt,tau_f);
    
    % This implements B .~ A_lp
    ang1 = (1+tanh(heavi_scale*((lv1_f_num-lv0_s_num)-0.5)))/2;
    
    % This implements A .~ B_lp
    ang2 = (1+tanh(heavi_scale*((lv0_f_num-lv1_s_num)-0.5)))/2;
        
    %This implements B - (B .~ A_lp)
    r_bl1_an = lv1_f_num - ang1;
    
    %This implements A - (A .~ B_lp)
    r_bl2_an = lv0_f_num - ang2;
    
    %full opponent response
    r_bl_an = r_bl1_an - r_bl2_an;
    
    %compute steady-state signal amplitude over an integer number of cycles
    ss_amp(i) = mean(r_bl_an( (t_vect>=t_min)&(t_vect<t_max(i)) ));
    
    %uncomment to compute oscillation amplitude and SNR
    %steady-state oscillation amplitude
    %sso_amp_an(i) = max(r_bl_an(t_vect>=t_min))-min(r_bl_an(t_vect>=t_min));
    
    %compute signal to noise ratio
    %snr_bl_an(i) = sss_amp_an(i)/sso_amp_an(i);
    
    if ( plot_debug == 1 )
        %theoretical steady-state response of correlator model
        r_ss_th = sin(2*pi*omega_x*dx)*alpha^2*2*pi*omega_td(i)*(tau_s-tau_f)/ ...
            ((1 +(2*pi*omega_td(i)*tau_s)^2)*(1+(2*pi*omega_td(i)*tau_f)^2));
        
        %plot AN output
        figure;
        subplot(2,1,1);
        hp1 = plot(t_vect,r_bl1_an,'r:');
        hold on;
        hp2 = plot(t_vect,r_bl2_an,'g:');
        legend([hp1 hp2],{'half-bl pref dir','half-bl null dir'});
        title('single AN BL correlator responses');
        subplot(2,1,2);
        hp3 = plot(t_vect,r_bl_an);
        hold on;
        hp4 = plot([t_start t_end],r_ss_th*[1 1],'r--');
        hp5 = plot(t_vect,ones(size(t_vect))*ss_amp,'g:');
        legend([hp3 hp4 hp5],{'full bl','ss theor','AN time averaged'});
        xlabel('time (s)');
        ylabel('bl and correlator responses');
    end
    
end

if ( plot_debug == 1 )
    %debug
    figure;
    ha = axes;
    semilogx(omega_td,ss_amp,'x');
    set(ha,'XLim',[1e-1 1e3]);
end

%% Compare the tuning curves of HR and extended BL
max_R_ss = max(R_ss);
R_ss_n = R_ss/max_R_ss;

h_f = figure; 
h_ax = axes;
h_p1 = semilogx(omega_tv,R_ss_n,'k');
hold on;

rssth_v_n = rssth_v/max_R_ss;
%plot(omega_td,r_num_ss_n,'ok');

[max_val, ind_max] = max(rssth_v_n);
ss_amp_n = (max_val/ss_amp(ind_max))*ss_amp;
h_p2 = plot(omega_td,ss_amp_n,'or');
xlabel('temporal frequency \omega_t (c/s)');
ylabel('normalized response');
legend([h_p1 h_p2],{'HR','eBL'});
set(h_ax,'TickDir','out');

%print(h_f,'figures/bl_response2.eps','-depsc');

h_f = figure; 
h_ax = axes;
semilogx(omega_td,ss_amp./rssth_v,'ko');
set(h_ax,'XLim',[1e-1 1e3],'TickDir','out');
xlabel('temporal frequency \omega_t (c/s)');
ylabel('eBL gain relative to HR');

%print(h_f,'figures/bl_response3.eps','-depsc');

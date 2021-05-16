%% Standard correlator response
%plot the output of a low pass filter with time constants of 20 and 2 ms
%as well as the output of a high pass filter with a time constant of 50 ms
%Plot the response of the two correlators at 3 delays of 12.5, 15 and 17.5
%ms. The input signal is a 10 ms long step in contrast. 

%luminance step after 5 ms
t_facet1 = 5e-3; 
%duration of the pulse (in s)
duration = 10e-3; 

%time step in s 
dt = 0.05e-3;
t_max = 100e-3;

tv = 0:dt:t_max-dt;
n_t = length(tv);

lum1 = (tv > t_facet1) & (tv < t_facet1 + duration);

h_f = figure;
h_a(1) = subplot(2,1,2);
h_p1 = plot(tv,lum1);
%hold on;
set(h_a(1),'YLim',[-0.1 1.1]);
%legend(h_p1,{'pulse'},'Box','off','Location','northwest');
xlabel('time (s)');
set(h_a(1),'TickDir','out');

%low-pass filter time constant in ms 
tau_lp = 20e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% low pass and high pass filter luminance pulse
lp_lum1 = lp1filt_fn(lum1,dt,tau_lp);

%figure; 
h_a(2) = subplot(2,1,1);
plot(tv,lp_lum1,'b');
hold on; %other stuff to plot below

debug = 0;
if ( debug == 1 )
    %compute and plot exact solution
    lp_lum1_th = zeros(size(tv));
    inds = find( (tv>t_facet1) & (tv < t_facet1 + duration) );
    lp_lum1_th(inds) = 1-exp( (t_facet1-tv(inds))/tau_lp );
    inds = find( tv >= t_facet1+duration);
    lp_lum1_th(inds) = exp( (t_facet1+duration-tv(inds))/tau_lp ) - exp( (t_facet1-tv(inds))/tau_lp );
    plot(tv,lp_lum1_th,'r:');
end

tau_lpf = 2e-3; 
lpf_lum1 = lp1filt_fn(lum1,dt,tau_lpf);

plot(tv,lpf_lum1,'k');

tau_hp = 50e-3; 

hp_lum1 = lum1 - lp1filt_fn(lum1,dt,tau_hp);

plot(tv,hp_lum1,'r');
set(h_a(2),'YLim',[-0.2 1.1]);
set(h_a(2),'TickDir','out');
ylabel('filter response');

%print(h_f,'figures/stdmod_resp1.eps','-depsc');

debug = 0;
if ( debug == 1 )
    %compute and plot exact solution
    hp_lum1_th = zeros(size(tv));
    inds = find( (tv>t_facet1) & (tv < t_facet1 + duration) );
    hp_lum1_th(inds) = exp( (t_facet1-tv(inds))/tau_hp );
    inds = find( tv >= t_facet1+duration);
    hp_lum1_th(inds) = exp( (t_facet1-tv(inds))/tau_hp ) - exp( (t_facet1+duration-tv(inds))/tau_hp );
    plot(tv,hp_lum1_th,'r:');
end

%% Response at 3 delays of the two models

pulse_delay = [12.5e-3 15e-3 17.5e-3]; %s

h_f = figure; 
for i = 1:length(pulse_delay)
    t_facet2 = t_facet1 + pulse_delay(i);
    lum2 = (tv > t_facet2) & (tv < t_facet2 + duration);
    
    lp_lum2 = lp1filt_fn(lum2,dt,tau_lp);
    lpf_lum2 = lp1filt_fn(lum2,dt,tau_lpf);
    hp_lum2 = lum2 - lp1filt_fn(lum2,dt,tau_hp);
    
    %correlator responses
    r1 = lp_lum1.*lpf_lum2 - lp_lum2.*lpf_lum1;
    r2 = lp_lum1.*hp_lum2 - lp_lum2.*hp_lum1;
    
    h_a = subplot(3,1,i);
    h_p1 = plot(tv,r1);
    hold on;
    plot(tv,r2,'r');
    legend(h_p1,{['delay: ' num2str(pulse_delay(i)) ' ms']});
    set(h_a,'YLim',[0 0.4]);
    set(h_a,'TickDir','out');
end
xlabel('time (s)');
ylabel('correlator response');
%print(h_f,'figures/stdmod_resp2.eps','-depsc');

%% Response to various delays

pulse_delay = [-50:0.1:50]*1e-3; %in s
int_resp1 = zeros(size(pulse_delay));
int_resp2 = zeros(size(pulse_delay));

for i = 1:length(pulse_delay);
    %note: the first pulse and its correlation delayed version do not
    %change; we need only to update the second pulse and its delayed
    %version
    t_facet2 = t_facet1 + abs(pulse_delay(i));
    lum2 = (tv > t_facet2) & (tv < t_facet2 + duration);
    
    lp_lum2 = lp1filt_fn(lum2,dt,tau_lp);
    lpf_lum2 = lp1filt_fn(lum2,dt,tau_lpf);
    hp_lum2 = lum2 - lp1filt_fn(lum2,dt,tau_hp);
    
    %correlator responses
    r1 = lp_lum1.*lpf_lum2 - lp_lum2.*lpf_lum1;
    r2 = lp_lum1.*hp_lum2 - lp_lum2.*hp_lum1;
    
    if (pulse_delay(i) < 0)
        %effectively swap the order of the pulses
        r1 = -r1;
        r2 = -r2;
    end

    int_resp1(i) = sum(r1)*dt;
    int_resp2(i) = sum(r2)*dt;
    
end

%convert to ms
int_resp1 = int_resp1*1e3;
int_resp2 = int_resp2*1e3;

h_f = figure;
h_a = axes;
plot(pulse_delay,int_resp1);
hold on;
plot(pulse_delay,int_resp2);
set(h_a,'TickDir','out');
xlabel('delay between pulses (s)');
ylabel('integrated response (ms)');
%print(h_f,'figures/stdmod_resp3.eps','-depsc');

%% Basic correlation model 
% plot two 10 ms pulses separated by 20 ms onset times
% plot the two pulses shifted by 20 ms
% plot output of correlation detector for 3 delays of 5, 10, 15 ms
% plot integrated response as a function of pulse delay

%define the pulses 
dt = 0.05; %ms
t_max = 60; %ms
t_vect = 0:dt:t_max;

%start of first pulse
t1 = 10;
%width of the pulse
wd = 10;

pulse1 = (t_vect>=t1)&(t_vect<(t1+wd));

%delay between the two pulses
pulse_delay = 15; %ms
t2 = t1 + pulse_delay;
pulse2 = (t_vect>=t2)&(t_vect<(t2+wd));

h_f = figure; 
h_a(1) = subplot(2,1,1);
h_p1 = plot(t_vect,pulse1,'r');
hold on;
h_a(2) = subplot(2,1,2);
h_p2 = plot(t_vect,pulse2,'k');
hold on;

%correlator delay
corr_delay = 20; %ms
%pulse1 onset after delay stage
t1_cd = t1+corr_delay; 
pulse1_cd = (t_vect>=t1_cd)&(t_vect<(t1_cd+wd));

%pulse2 onset after delay stage
t2_cd = t2+corr_delay;
pulse2_cd = (t_vect>=t2_cd)&(t_vect<(t2_cd+wd));

axes(h_a(1));
h_p1d = plot(t_vect,pulse1_cd,'r--');
legend([h_p1 h_p1d],{'pulse at facet 1','after delay stage'});
set(h_a(1),'TickDir','out');
axes(h_a(2));
h_p2d = plot(t_vect,pulse2_cd,'k--');
legend([h_p2 h_p2d],{'pulse at facet 2', 'after delay stage'});
xlabel('time (ms)');
ylabel('stimulus contrast');
set(h_a(2),'TickDir','out');

%print(h_f,'figures/pulses.eps','-depsc');

%% Compute ouptput of correlator for 3 delays
pulse_delay = [12.5 15 17.5]; %ms

h_f = figure;
for i = 1:length(pulse_delay)
    %note: the first pulse and its correlation delayed version do not
    %change; we need only to update the second pulse and its delayed
    %version
    t2 = t1 + pulse_delay(i);
    pulse2 = (t_vect>=t2)&(t_vect<(t2+wd));
    
    t2_cd = t2+corr_delay;
    pulse2_cd = (t_vect>=t2_cd)&(t_vect<(t2_cd+wd));
    
    r = pulse1_cd.*pulse2 - pulse2_cd.*pulse1;
    
    h_a = subplot(3,1,i);
    h_p(i) = plot(t_vect,r,'k');
    legend(h_p(i),{['delay: ' num2str(pulse_delay(i)) ' ms']});
    set(h_a,'TickDir','out');
end

xlabel('time (ms)');
ylabel('model output');

%print(h_f,'figures/basic_hrc_exresp.eps','-depsc');

%% Plot integrated response as a function of delay
pulse_delay = -40:0.1:40;
int_resp = zeros(size(pulse_delay));

for i = 1:length(pulse_delay);
    %note: the first pulse and its correlation delayed version do not
    %change; we need only to update the second pulse and its delayed
    %version
    t2 = t1 + abs(pulse_delay(i));
    pulse2 = (t_vect>=t2)&(t_vect<(t2+wd));
    
    t2_cd = t2+corr_delay;
    pulse2_cd = (t_vect>=t2_cd)&(t_vect<(t2_cd+wd));
    
    r = pulse1_cd.*pulse2 - pulse2_cd.*pulse1;

    if (pulse_delay(i) < 0)
        %effectively swap the order of the pulses
        r = -r;
    end

    int_resp(i) = sum(r)*dt;
end

h_f = figure;
h_a = axes;
plot(pulse_delay,int_resp,'k');
xlabel('delay between pulses (ms)');
ylabel('integrated response (ms)');
set(h_a,'TickDir','out');

%print(h_f,'figures/deltuning_basic.eps','-depsc');
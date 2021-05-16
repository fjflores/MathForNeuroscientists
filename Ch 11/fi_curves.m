%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

function fi_curves
%plot the theoretical fI curve of the LIF
[fr, i_inj] = lif_fr(1,1,20,0);

h_f1 = figure; 
h_a1 = axes;
line('Parent',h_a1,'XData',i_inj,'YData',fr);

%perfect integrator
rin = 20; % in MOhms
tau = 30; % in msec
tref = 1; % in msec
vth = 16; % in mV

%tref = 0;
fr_int = 1./(tref + vth*tau./(i_inj*rin));%in spk/ms
fr_int = 1000*fr_int; %in spk/s

line('Parent',h_a1,'XData',i_inj,'YData',fr_int,'Color','r');

%LIF with reset voltage set at 8 mV
[fr_res, i_inj] = lif_fr(1,1,20,8);
line('Parent',h_a1,'XData',i_inj,'YData',fr_res,'LineStyle','--');

set(gca,'YLim',[0 600]);
xlabel('current (nA)');
ylabel('firing frequency (spk/s)');

%print(h_f1,'-depsc','fi_lif.eps');

function [fr,i_inj] = lif_fr (i_min, i_step, i_max,vreset)
%
% [fr,i_inj] = lif_fr(i_min, i_step, i_max)
%
%lif_fr returns the f(i_inj) curve within the parameter range of currents.
%If the minimum current is negative, returns a sensible range of ...
%current values


%Biophysical parameters of the LIF model
rin = 20; % in MOhms
tau = 30; % in msec
tref = 1; % in msec
vth = 16; % in mV
%vreset = 0;

%Compute minimum current needed to get spikes (nA)
i_rheo = vth/rin; 

%Saturation current needed to get within 5% firing rate
%of the maximal firing rate set by the refractory period
i_sat = vth/(rin*(1-exp(-(0.05/(1-0.05))*(tref/tau))));

if ( i_min < 0 )
    %display rheobase and saturation current
    rheo_str = sprintf('Rheobase current: %f',i_rheo);
    disp(rheo_str);
    sat_str = sprintf('Saturation current: %f',i_sat);
    disp(sat_str);
    fr = [];
    i_inj = [];
    return;
end;

%A couple of obvious checks
if ( i_min <= i_rheo )
    warn_str1 = sprintf('warning i_min = %f smaller than rheobase: %f',...
	i_min,i_rheo);
    warn_str2 = sprintf('setting i_min to rheobase');
    i_min = i_rheo + i_step;
end;
  
if ( i_min > i_max )
    disp('fatal error: i_min larger than i_max')
    return;
end;


i_inj = i_min:i_step:i_max;
fr = 1./(tref - tau*log(1 - (vth-vreset)./(rin*i_inj))); %in kspk/sec

fr = 1000*fr; %in spk/sec


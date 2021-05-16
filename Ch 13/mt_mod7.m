%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%
%  called by mt_mod6.m

function [tv, Pv, Vv, Qv, ind_v] = mt_mod7(ap_freq, t_max)
%
% ap_freq: stimulation frequency (Hz)
% t_max: end of simulation (in ms)
% 
%
% returns the time vector corresponding to the simulations (tv)
% the resource vector Pv, with maximum value of 1
% the membrane potential vector Vv (in mV)
% the charge vector Qv 
% a vector of indices of the time of action potentials (ind_v)
%
%basic setup for the Markram-Tsodyks model
%This file generates a stochastic MT model
%

t_rec = 800; %ms
dt = 0.5; %ms

t_frac = dt/t_rec;
den_t = (1 + t_frac);

Pvmax = 1;
Pv0 = Pvmax;
p = 0.55;

if ( isempty(t_max) )
    t_max = 1000;
end;

tv = 0:dt:t_max;
n_tv = length(tv);

Pv = zeros(size(tv));
Pv(1) = Pv0;

if ( isempty(ap_freq) )
    ap_freq = 20; % 1/s
end;

%convert to time interval in ms and round to the closest bin
dt_ap = round((1000/ap_freq)/dt);

%time of first pulse 
t_aps = 20; %in ms
t_ind = round(t_aps/dt);
%indices of pulses
ind_v = t_ind:dt_ap:n_tv;

%leaky IF membrane potential
Vv = zeros(size(tv));
V0 = 0; %rest is at zero
Vv(1) = V0;

tm = 50; %ms
tm_frac = dt/tm;
den_v = (1 + tm_frac);

%synaptic factor
Rin = 250; %MOhm
qmax = 1.5; %pC
syn = p*Rin*qmax/(Pvmax*tm);
ap = (p*qmax/Pvmax); %in pC
Qv = zeros(size(tv));

for i = 2:n_tv
    Pv(i) = (Pv(i-1) + t_frac * Pvmax)/den_t;

    %backward Euler on membrane potential
    Vv(i) = Vv(i-1)/den_v;
    
    if ( ~isempty(find(ind_v == i)) )
        %time for a pulse
        
        %first update V and Q by taking the fraction p of the 
        %available resources
        Vv(i) = Vv(i) + syn*Pv(i)/den_v;
        Qv(i) = ap*Pv(i);
        
        %then update the available resources as 1-p those immediately
        %before the AP
        Pv(i) = (1-p)*Pv(i);

    end;
end;

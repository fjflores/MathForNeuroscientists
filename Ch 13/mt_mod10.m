%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  called by mt_mod12.m

function [tv, Pv, Uv, Vv, Qv, ind_v] = mt_mod10(ap_freq, t_max)
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
%Implements facilitation as well as depression
%using the second MT model

t_rec = 130; %ms
dt = 0.5; %ms

t_frac = dt/t_rec;
den_t = (1 + t_frac);

Pvmax = 1;
Pv0 = Pvmax;

u0 = 0.03;

if ( isempty(t_max) )
    t_max = 1000;
end;

tv = 0:dt:t_max;
n_tv = length(tv);

Pv = zeros(size(tv));
Pv(1) = Pv0;
t_facil = 530;
t_fracu = dt/t_facil;
den_u = 1 + t_fracu;

Uv = zeros(size(tv));
Uv(1) = u0;

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

tm = 60; %ms
tm_frac = dt/tm;
den_v = (1 + tm_frac);

%synaptic factor
Rin = 1000; %MOhm
qmax = 2e-3*1540; %pC because 3ms * 1540 pA
syn = Rin*qmax/(Pvmax*tm);
ap = (qmax/Pvmax); %in pC
Qv = zeros(size(tv));

for i = 2:n_tv
    Pv(i) = (Pv(i-1) + t_frac * Pvmax)/den_t;

    %backward Euler on membrane potential
    Vv(i) = Vv(i-1)/den_v;
    
    Uv(i) = (Uv(i-1) + t_fracu * u0)/den_u; 
    
    if ( ~isempty(find(ind_v == i)) )
        %time for a pulse
        
        %first update V and Q by taking the fraction p of the 
        %available resources
        Vv(i) = Vv(i) + syn*Uv(i)*Pv(i)/den_v;
        Qv(i) = ap*Uv(i)*Pv(i);
        
        %then update the available resources as 1-p those immediately
        %before the AP
        Pv(i) = (1-Uv(i))*Pv(i);
        Uv(i) = Uv(i) + u0*(1-Uv(i));
        
    end;
end;

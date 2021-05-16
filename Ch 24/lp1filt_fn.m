function lp1vect = lp1filt_fn(lum_vect,dt, tau_lp)
% Numerical implementation of first order low-pass filter for correlator models
%
% Usage:   vect = lp1filt_fn(lum_vec,dt,tau)
%
% where lum_vect is to be filtered,
% tau_lp is the time constant in s,
% dt is time step in s. 

if ( nargin ~= 3 )
    disp('Unexpected number of inputs. Aborting...');
    return;
end

%%%%%%%%%%%%%%%%%%
% Model parameters

if ( isempty(tau_lp) )
    %time constants of 1st order low pass (ms)
    tau_lp = 20e-3; %s
end

if ( isempty(dt) )
    dt = 0.05e-3; %s
end

%setup constants for backward euler
lp_const1 =  dt/(tau_lp + dt);
lp_const2 = 1 - lp_const1;

lp1vect = zeros(size(lum_vect));

%initial condition 
lp1vect(1) = lum_vect(1);

%marching algorithm
for i = 2:length(lum_vect)
    lp1vect(i) = lp_const1*lum_vect(i)+lp_const2*lp1vect(i-1);
end

end
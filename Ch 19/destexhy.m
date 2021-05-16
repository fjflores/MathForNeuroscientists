%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  called by destex_f1.m

%
% destexhy.m
%
% set up the destexhe equations on a patch
% and solve via hybrid euler
%
% usage		v = destexhy(dt,t1,t2,T,Istim)
%
% where     dt = timestep (ms)
%	        t1 = start time of current pulse (ms)
%	        t2 = stop time of current pulse (ms)
%	        T = final time (ms)
%	        Istim = size of current pulse (micro A/cm^2)
%
%		v = voltage (mV)
% example
%		v = destexhy(0.01,10,100,200,-1)
%		v = destexhy(0.01,10,100,200,3)
%
% original	Destexhe et al. Neuroscience 2001 
%           Destexhe and Pare J. Neurophysiol 1999
%           modelDB

function v = destexhy(dt,t1,t2,T,Istim)

% maximal conductances, in mS / (cm)^2
%gK = 0;
%gM = 0;
%gNa = 0;
gK = 10; %delayed rectifier            
gM = 0.5; %M current
gNa = 51.6; %HH-type sodium channel        
gl = 0.045; %leak           

% reversal potentials, in mV		
VK = -90;             
VNa = 50;           
Vl = -80;         

%capacitance
C = 1;                  % micro F / (cm)^2

N = floor((T-dt)/dt);
v = zeros(N,1);		% allocate space for v
%t = (1:N)'*dt;

% initial conditions
v(1) = -80.3935;

%steady-state values
m = am(v(1))/(am(v(1)) + bm(v(1)));
h = ah(v(1))/(ah(v(1)) + bh(v(1)));
n = an(v(1))/(an(v(1)) + bn(v(1)));
p = ap(v(1))/(ap(v(1)) + bp(v(1)));

j = 2;			% time counter

while j*dt < T,		% run for T ms

t = j*dt;

    Iapp = 0;		% apply Istim between t1 and t2 
    if (t>t1 & t<t2)
       Iapp = Istim;
    end

    % advance the gating variables, m, h, n, and p using hybrid Euler

    tmp = am(v(j-1));
    m = (m+dt*tmp)/(1+dt*(tmp+bm(v(j-1))));

    tmp = ah(v(j-1));
    h = (h+dt*tmp)/(1+dt*(tmp+bh(v(j-1))));

    tmp = an(v(j-1));
    n = (n+dt*tmp)/(1+dt*(tmp+bn(v(j-1))));

    tmp = ap(v(j-1));
    p = (p+dt*tmp)/(1+dt*(tmp+bp(v(j-1))));

    % compute the conductances of the K, Na, NaP, h & Ca currents

    cIK = gK*(n^4);
    cINa = gNa*m^3*h;
    
    %M current 
    cIM = gM*p;

    % advance the voltage, v

    tmp1 = cINa*VNa+(cIK+cIM)*VK+gl*Vl+Iapp;
    tmp2 = cINa+cIK+cIM+gl;

    %the correct formula is 
    %v(j) = (C*v(j-1) + dt*tmp1)/(C + dt*tmp2)
    %but with C = 1 the following will do.
    v(j) = (v(j-1) + dt*tmp1)/(1 + dt*tmp2);
   
    j = j + 1;

end

return

% steady-state activation, time constant and opening and closing
% functionals from Destexhe and  Pare


%forward rate of delayed rectifier
function val = an(V)
VT = -58;
tmp = (V-VT-15);
val = -0.032*tmp./(exp(-tmp/5)-1);

%backward rate of delayed rectifier
function val = bn(V)
VT = -58;
val = 0.5*exp(-(V-VT-10)/40);

%forward rate of fast sodium current  activation
function val = am(V)
VT = -58; 
tmp = (V-VT-13);
val = -0.32*tmp/(exp(-tmp/4)-1);

%backward rate of fast sodium current activation
function val = bm(V)
VT = -58;
tmp = (V-VT-40);
val = 0.28*tmp/(exp(tmp/5)-1);

%forward rate of fast sodium current inactivation
function val = ah(V)
VT = -58;
VS = -10;
tmp = (V-VT-VS-17);
val = 0.128*exp(-tmp/18);

%backward rate of fast sodium current inactivation
function val = bh(V)
VT = -58;
VS = -10;
tmp = (V-VT-VS-40);
val = 4/(1+exp(-tmp/5));

%forward rate of non-inactivating K current
function val = ap(V)
tmp = (V+30);
val = 0.0001*tmp./(1-exp(-tmp/9));

%backward rate of non-inactivating K current
function val = bp(V)
tmp = (V+30);
val = -0.0001*tmp./(1-exp(tmp/9));








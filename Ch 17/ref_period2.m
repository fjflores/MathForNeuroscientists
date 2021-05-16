%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

function ref_period2
%initial amplitude of second pulse (in pA)
I1_min = 60;
I1_max = 10*I1_min;

%simulation endpoint
Tfin = 100;

%for debugging purposes
%t2_v = 20:10:80;

t2_v = 15:1:80;
ith_v = zeros(size(t2_v));

disp('running, this may take some time...');
for i = 1:length(t2_v)
    t2 = t2_v(i);

    I1 = I1_min;
    spcnt = getRefrac(t2,Tfin,I1);

    while ( (spcnt < 2) & (I1 < I1_max) )
    
        I1 = I1 + 1;
    
        spcnt = getRefrac(t2,Tfin,I1);
    
    end

    if ( I1 < I1_max )
        ith_v(i) = I1;
        fprintf('Current threshold is %0.4g pA at t2 = %.2f \n',I1,t2);
    else
        ith_v(i) = -1;
        err_msg = sprintf('max. current %.1f pA exceeded',I1_max);
        disp(err_msg);
    end;
end;

h_f1 = figure; 
h_a1 = axes;
I0 = 67;
ith_vn = ith_v/I0;
line('Parent',h_a1,'XData',t2_v,'YData',ith_vn);
line('Parent',h_a1,'XData',[0 80],'YData',[1 1],'LineStyle','--');
xlabel('Interpulse interval (ms)');
ylabel('I_1/I_0');

h_f2 = figure;
h_a2 = subplot(2,1,1);
[spcnt1,v,t] = getRefrac(t2_v(1),Tfin,ith_v(1));
line('Parent',h_a2,'XData',t,'YData',v); %,'LineWidth',1);
set(h_a2,'XLim',[0 100],'YLim',[-10 120]);
title(h_a2,sprintf('t2 = %.2f Ith = %0.05g pA',t2_v(1),ith_v(1)));

%for debugging
%i_rep2 = 3;
i_rep2 = 26;

h_a3 = subplot(2,1,2);
[spcnt1,v,t] = getRefrac(t2_v(i_rep2),Tfin,ith_v(i_rep2));
line('Parent',h_a3,'XData',t,'YData',v); %,'LineWidth',1);
set(h_a3,'XLim',[0 100],'YLim',[-10 120]);
title(h_a3,sprintf('t2 = %.2f Ith = %0.05g pA',t2_v(i_rep2),ith_v(i_rep2)));
xlabel(h_a3,'Time (ms)');
ylabel(h_a3,'membrane potential (mV)');

%print(handles.figure1,'-depsc2','ref_period2.eps');

function [spcnt,V,T] = getRefrac(t2,Tfin,I1)
%returns the number of spikes and the membrane potential and time vector
%of the corresponding simulation. t2 is the onset time of the second, 1ms
%long pulse and I1 its intensity. Tfin is the last simulation time. 

%debug
t2 = 40.0;
Tfin = 100.0;
I1 = 60.0;

%this is the smallest current yielding an action potential at pA resolution
I0 = 67; %pA = 10e-12A; ultimately the injected current will be in micro[mu] A/cm^2
vth = 20;
dt = 0.01;

A = 4*pi*(1e-6);% cm^2
vK = -6;	% mV
vNa = 127;	% mV
vCl = 2.8417; 	% mV

Cm = 1;		% micro F/cm^2
GK =  36;	% mS/(cm)^2
GNa = 120;	% mS/(cm)^2
GCl = 0.3;	% mS/(cm)^2

N = ceil(Tfin/dt);
T = zeros(1,N);   % allocate space for long vectors
V = T;

t = 0;
v = 0;
n = an(0)/(an(0)+bn(0));  
m = am(0)/(am(0)+bm(0)); 
h = ah(0)/(ah(0)+bh(0));

up = 0;
spcnt = 0;
Iapp0 = I0*(1e-6)/A; %1 pA = 1e-6 muA
Iapp1 = I1*(1e-6)/A;

%simulate the double current pulse
j = 1;

while t < Tfin+2

    t = t + dt;

    Istim = Iapp0*(t>2)*(t<3) + Iapp1*(t>t2)*(t<t2+1);

    %Istim = Iapp0*(t>2)*(t<3);
    %one pulse of 1 ms. 
    
    ant = an(v);
    n = ( n + dt*ant )/(1 + dt*(ant+bn(v)) );
    amt = am(v);
    m = ( m + dt*amt )/(1 + dt*(amt+bm(v)) );
    aht = ah(v);
    h = ( h + dt*aht )/(1 + dt*(aht+bh(v)) );

    iNa = GNa*m^3*h;
    iK = GK*n^4;
    
    top = Cm*v+dt*(iNa*vNa + iK*vK + GCl*vCl + Istim);

    bot = Cm + dt*(iNa + iK + GCl);

    v = top/bot; 

    if v > vth
       up = 1;  	     % mark upward crossing
    else
       spcnt = spcnt + up;   % increment spcnt on downward crossing
       up = 0;
    end
    
    T(j) = t; V(j) = v;
    j = j + 1;

end

return

function val = an(v)
val = .01*(10-v)./(exp(1-v/10)-1);

function val = bn(v)
val = .125*exp(-v/80);

function val = am(v)
val = .1*(25-v)./(exp(2.5-v/10)-1);

function val = bm(v)
val = 4*exp(-v/18);

function val = ah(v)
val = 0.07*exp(-v/20);

function val = bh(v)
val = 1./(exp(3-v/10)+1);




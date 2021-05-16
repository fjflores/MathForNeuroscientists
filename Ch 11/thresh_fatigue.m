%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%


%
%constant current injection
%

%integration time step in ms (0.1 ms)
dt = 0.1; 

%capacitance, nF
Cm = 2;

%membrane time constant
tau_v = 20; %ms

%threshold in mV
vm_thres0 = 8;

%increment in threshold for each spike in mV
thres_incr = 4;

%threshold recovery time constant ms
tau_t = 80; 

%1 second
t_v = 0:dt:350-dt;

%current vector nA
i_v = [zeros(size([0:dt:50-dt])) 2*ones(size([50:dt:300-dt])) zeros(size([300:dt:350-dt]))];

%membrane potential vector mV
vm_v = zeros(size(t_v));

%threshold values during simulations
thres_v = zeros(size(t_v));

%initializes to baseline threshold 
thres_v(1) = vm_thres0;

%spike vector
s_v = zeros(size(t_v));

%scale by dt to speed up computations
i_v_n = i_v*dt/Cm;

%membrane time constant: precompute values used in loop to save time
dtau_m1 = dt/tau_v;
m_fact = (1-dtau_m1);

%threshold time constant: calculations to speed up the simulations
dtau_t1 = dt/tau_t;
t_fact = 1- dtau_t1;
t0_fact = vm_thres0*dtau_t1;

%integrate with forward euler
for k = 2:length(t_v)
    %membrane potential evolution
    vm_v(k) = vm_v(k-1)*m_fact + i_v_n(k-1);
    
    %threshold evolution
    thres_v(k) = thres_v(k-1)*t_fact + t0_fact;
    
    if (vm_v(k) > thres_v(k))
        vm_v(k) = 0; %reset
        s_v(k) = 1; %register spike
        thres_v(k) = thres_v(k) + thres_incr; %increment threshold
    end;
    
end;

ind = find(s_v > 0.5);
t_spks = t_v(ind);

h_f1 = figure;
h_a1 = subplot(3,1,1);
h_a2 = subplot(3,1,2);
h_a3 = subplot(3,1,3);

line('Parent',h_a1,'XData',t_v,'YData',vm_v);

for i =1:length(t_spks)
    line('Parent',h_a1,'XData',[t_v(ind(i)-1) t_v(ind(i)-1)],'YData',[vm_v(ind(i)-1) 50]);
end;

set(h_a1,'XLim',[0 350],'YLim',[0 50]);

%without threshold fatigue
vm_v2 = zeros(size(t_v));
s_v2 = zeros(size(t_v));

%integrate with forward euler
for k = 2:length(t_v)
    %membrane potential evolution
    vm_v2(k) = vm_v2(k-1)*m_fact + i_v_n(k-1);
        
    if (vm_v2(k) > vm_thres0)
        vm_v2(k) = 0; %reset
        s_v2(k) = 1; %register spike
    end;
    
end;

ind2 = find(s_v2 > 0.5);
t_spks2 = t_v(ind2);

line('Parent',h_a2,'XData',t_v,'YData',vm_v2);
ylabel(h_a2,'membrane potential (mV)');

for i =1:length(t_spks2)
    line('Parent',h_a2,'XData',[t_v(ind2(i)-1) t_v(ind2(i)-1)],'YData',[vm_v2(ind2(i)-1) 50]);
end;

set(h_a2,'XLim',[0 350],'YLim',[0 50]);

line('Parent',h_a3,'XData',t_v,'YData',i_v);
set(h_a3,'XLim',[0 350]);
xlabel(h_a3,'time (ms)');

%print(handles.figure1,'-depsc','thresh_fatigue.eps');

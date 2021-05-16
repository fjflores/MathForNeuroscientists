%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%

%% 
%simulate an integrate and fire neuron with a time-varying
%threshold that is incremented after an action potential and
%relaxes exponentially to its resting value,between action potentials

%Does the same thing for a static threshold and compares with the above
%model. Thiss is a modified version of lif1 and lif2

%integration time step in ms (0.1 ms)
dt = 0.1; 

%time vector; 100000 ms = 100 s
t_v = 0:dt:100000-dt;

%capacitance, nF
Cm = 2;

%current vector nA
i_v = ones(size(t_v));

%white noise vector nA
i_rnd = normrnd(0,1.5,1,length(t_v));
i_v = i_v + i_rnd;

%membrane potential vector mV
vm_v = zeros(size(t_v));

%threshold in mV
vm_thres0 = 8;

%increment in threshold for each spike in mV
thres_incr = 4;

%threshold recovery time constant ms
tau_t = 80; 

%threshold values during simulations
thres_v = zeros(size(t_v));
%initializes to baseline threshold 
thres_v(1) = vm_thres0;

%calculations to speed up the simulations
dtau_t1 = dt/tau_t;
t_fact = 1- dtau_t1;
t0_fact = vm_thres0*dtau_t1;

%spike vector
s_v = zeros(size(t_v));

%membrane time constant
tau_v = 20;
dtau_m1 = dt/tau_v;
m_fact = (1-dtau_m1);

%scale by dt to speed up computations
i_v_n = i_v*dt/Cm;

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
    
%compute ISI serial correlation coefficient
spk_times = t_v(find(s_v>0.5));

isi_v = diff(spk_times);

mean_isi = mean(isi_v);

v1 = isi_v(2:end) - mean_isi;
v2 = isi_v(1:end-1) - mean_isi;

%compute the serial correlation coefficient from the 
%formulas in the lecture notes
scc1 = mean(v1.*v2)/(std(v1)*std(v2));

%%%%%%%%%%%%%%%%%%%%%%%
%
%static threshold model
%
%%%%%%%%%%%%%%%%%%%%%%%

%membrane potential vector mV
vm_s = zeros(size(t_v));

%threshold in mV; adjusted manually to 
vm_thress = 10.4;

%spike vector
s_s = zeros(size(t_v));

%integrate with forward euler
for k = 2:length(t_v)
    %membrane potential evolution
    vm_s(k) = vm_s(k-1)*m_fact + i_v_n(k-1);
    
    if (vm_s(k) > vm_thress)
        vm_s(k) = 0; %reset
        s_s(k) = 1; %register spike
    end;
    
end;
    
%compute ISI serial correlation coefficient
spk_times_s = t_v(find(s_s>0.5));

isi_s = diff(spk_times_s);

mean_isi_s = mean(isi_s);

v1_s = isi_s(2:end) - mean_isi_s;
v2_s = isi_s(1:end-1) - mean_isi_s;

%compute the serial correlation coefficient from the 
%formulas in the lecture notes
scc1_s = mean(v1_s.*v2_s)/(std(v1_s)*std(v2_s));

%%

h_f1 = figure; 
h_a1 = subplot(2,2,1);
[n,xout] = hist(isi_v,20);
h = bar(h_a1,xout,n);
set(h,'FaceColor','r');
set(h,'EdgeColor','r');
set(h_a1,'XLim',[0 300]);

h_a2 = subplot(2,2,2);
[n,xout] = hist(isi_s,20);
h = bar(h_a2,xout,n);
set(h,'FaceColor','b');
set(h,'EdgeColor','b');
set(h_a2,'XLim',[0 300]);
xlabel(h_a2,'ISI (ms)');

h_a3 = subplot(2,2,3);
line('Parent',h_a3,'XData',v2+mean_isi,'YData',v1+mean_isi,...
    'Marker','x','Color','r','LineStyle','none');
set(h_a3,'XLim',[0 200],'YLim',[0 200]);
xlabel(h_a3,'ISI_n (ms)');
ylabel(h_a3,'ISI_n_+_1 (ms)');
title_str = sprintf('r = %.2g',scc1);
title(h_a3,title_str);

h_a4 = subplot(2,2,4);
line('Parent',h_a4,'XData',v2_s+mean_isi_s,'YData',v1_s+mean_isi_s,...
    'Marker','x','Color','b','LineStyle','none');
set(h_a4,'XLim',[0 200],'YLim',[0 200]);
title_str = sprintf('r = %.2g',scc1_s);
title(h_a4,title_str);

%%

h_f2 = figure;
h_a5 = subplot(3,1,1);
t_pl = spk_times(find(spk_times<2000));
stem(h_a5,t_pl,ones(size(t_pl)),'Marker','none','Color','r');
set(h_a5,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);

h_a6 = subplot(3,1,2);
t_pl_s = spk_times_s(find(spk_times_s<2000));
stem(h_a6,t_pl_s,ones(size(t_pl_s)),'Marker','none','Color','k');
set(h_a6,'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);

h_a7 = subplot(3,1,3);
tv_pl = t_v(find(t_v < 2000));
wn_pl = i_v(find(t_v < 2000));
line('Parent',h_a7,'XData',tv_pl,'YData',wn_pl);
xlabel(h_a7,'time (ms)');
ylabel(h_a7,'current (nA)');

%print(handles.figure1,'-depsc','neg_corr.eps'); 

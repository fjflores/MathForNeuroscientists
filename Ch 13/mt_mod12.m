%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%
%  requires mt_mod10.m

h_f1 = figure; 
h_a1 = subplot(5,1,1);
h_a2 = subplot(5,1,2);
h_a3 = subplot(5,1,3);
h_a4 = subplot(5,1,4);
h_a5 = subplot(5,1,5);

plot_intermediate = 1;

%stimulation frequencies
stim_freq = [1 2:2:100];
nstim_freq = length(stim_freq);

freq_q = zeros(1,nstim_freq);
freq_v = zeros(1,nstim_freq);
freq_m = zeros(1,nstim_freq);

for j = 1:nstim_freq
%for j = 6:6
    [tv, Pv, Uv, Vv, Qv, ind_v] = mt_mod10(stim_freq(j),6000);

    ind_max = length(tv);

    ind_ss = find(ind_v > 3000);
    nind_ss = length(ind_ss);

    Q_vals = zeros(1,nind_ss);
    V_vals = zeros(1,nind_ss);

    for i = 1:nind_ss-1
        Q_vals(i) = max(Qv(ind_v(ind_ss(i)):ind_v(ind_ss(i+1))));
        V_vals(i) = max(Vv(ind_v(ind_ss(i)):ind_v(ind_ss(i+1))));
    end;

    Q_vals(nind_ss) = max(Qv(ind_v(ind_ss(nind_ss)):ind_max));
    V_vals(nind_ss) = max(Vv(ind_v(ind_ss(nind_ss)):ind_max));

    freq_q(j) = mean(Q_vals);
    freq_v(j) = mean(V_vals);

    freq_m(j) = mean(Vv(ind_ss(1):ind_max));

    if ( (j == 6) & (plot_intermediate == 1) )
        xlim_val = 2000;
        line('Parent',h_a1,'XData',tv,'YData',Pv);
        set(h_a1,'XLim',[0 xlim_val]);
        
        line('Parent',h_a2,'XData',tv,'YData',Uv);
        set(h_a2,'XLim',[0 xlim_val]);
        
        for i = 1:length(ind_v)
            line('Parent',h_a3,'XData',[tv(ind_v(i)) tv(ind_v(i))],'YData',[0 1]);
        end;
        set(h_a3,'XLim',[0 xlim_val]);
        
        line('Parent',h_a4,'XDAta',tv,'YData',Vv);
        set(h_a4,'XLim',[0 xlim_val]);
        ylabel(h_a4,'Vm (mV)');
        
        line('Parent',h_a5,'XDAta',tv,'YData',Qv);
        set(h_a5,'XLim',[0 xlim_val]);
        xlabel(h_a5,'time (ms)');
        ylabel(h_a5,'charge (pC)')
    end;
end;

h_f2 = figure; 
h_a6 = subplot(3,1,1);
h_a7 = subplot(3,1,2);
h_a8 = subplot(3,1,3);

line('Parent',h_a6,'XData',stim_freq,'YData',freq_q);
ylabel(h_a6,'charge (pC)');
line('Parent',h_a7,'XData',stim_freq,'YData',freq_v);
ylabel(h_a7,'peak Vm (mV)');
line('Parent',h_a8,'XData',stim_freq,'YData',freq_m);
xlabel(h_a8,'stimulation frequency (Hz)');
ylabel(h_a8,'average Vm (mV)');


%theoretical curves
qmax = 2e-3*1540; %pC
u0 = 0.03;
tau_facil = 530; % ms
tau_facil_sec = tau_facil * 1e-3; %in seconds
v_u = exp(-1./(stim_freq*tau_facil_sec));

uss = u0./(1-(1-u0)*v_u);

tau_rec = 130; %ms
tau_rec_sec = tau_rec *1e-3; %in seconds
 
v_e = exp(-1./(stim_freq*tau_rec_sec));

q_ss = (qmax*uss).*(1-v_e)./(1-(1-uss).*v_e);

line('Parent',h_a6,'XData',stim_freq(1:2:end),'YData',q_ss(1:2:end),...
    'LineStyle','none','Marker','o','MarkerFaceColor','none',...
    'MarkerEdgeColor','r','MarkerSize',4); %'LineStyle',':','Color','r');

Rin = 1000; %MOhm
tm = 60; %ms
tm_sec = tm*1e-3; %in seconds

v_e2 = 1./(1-exp(-1./(stim_freq*tm_sec)));

v_ss = (q_ss*Rin/tm).*v_e2;

line('Parent',h_a7,'XData',stim_freq(1:2:end),'YData',v_ss(1:2:end),...
    'LineStyle','none','Marker','o','MarkerFaceColor','none',...
    'MarkerEdgeColor','r','MarkerSize',4); %,':','Color','r');

v_av = q_ss*Rin .* stim_freq*1e-3;
line('Parent',h_a8,'XData',stim_freq(1:2:end),'YData',v_av(1:2:end),...
    'LineStyle','none','Marker','o','MarkerFaceColor','none',...
    'MarkerEdgeColor','r','MarkerSize',4);

lim_val = qmax*Rin/tau_rec;
 
line('Parent',h_a7,'XData',[0 100],'YData',[lim_val lim_val],...
    'LineStyle','--');

line('Parent',h_a8,'XData',[0 100],'YData',[lim_val lim_val],...
    'LineStyle','--');

%print(h_f1,'-depsc2','mt_mod12.eps');
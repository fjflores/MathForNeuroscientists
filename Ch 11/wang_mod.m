%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  requires: wanghy.m
%

%time in ms
dt = 0.01;
tcurr_start = 10;
tcurr_end = 100;
tsim_end = 200;
curr_val = -1; 

[v, h, H2] = wanghy(dt,tcurr_start,tcurr_end,tsim_end,curr_val);
N = length(v);
t = (1:N)*dt;
curr_v = curr_val*ones(1,N).*(t>tcurr_start).*(t<tcurr_end);

h_f1 = figure; 
h_a1 = subplot(4,2,1);
h_a2 = subplot(4,2,2);
h_a3 = subplot(4,2,3);
h_a4 = subplot(4,2,4);
h_a5 = subplot(4,2,5);
h_a6 = subplot(4,2,6);
h_a7 = subplot(4,2,7);
h_a8 = subplot(4,2,8);

line('Parent',h_a1,'XData',t,'YData',v);
line('Parent',h_a3,'XData',t,'YData',h);
line('Parent',h_a5,'XData',t,'YData',H2);
line('Parent',h_a7,'XData',t,'YData',curr_v);

curr_val = 3;

[v, h, H2] = wanghy(dt,tcurr_start,tcurr_end,tsim_end,curr_val);
curr_v = curr_val*ones(1,N).*(t>tcurr_start).*(t<tcurr_end);

line('Parent',h_a2,'XData',t,'YData',v);
line('Parent',h_a4,'XData',t,'YData',h);
line('Parent',h_a6,'XData',t,'YData',H2);
line('Parent',h_a8,'XData',t,'YData',curr_v);

set(h_a1,'YLim',[-90 20]);
set(h_a2,'YLim',[-90 20]);
set(h_a4,'YLim',[0 0.5]);
set(h_a6,'YLim',[0 0.04]);
set(h_a7,'YLim',[-1 3]);
set(h_a8,'YLim',[-1 3]);
xlabel(h_a7,'time (ms)');
xlabel(h_a8,'time (ms)');
ylabel(h_a1,'membrane potential (mV)');

%print(handles.figure1,'-depsc','wang_mod.eps');

%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  requires prsolve.m
%

[t,y] = prsolve(1500,[0.75 0.0 2.1]);

h_f1 = figure; 
h_a1 = axes;
line('Parent',h_a1,'XData',t,'YData',y(:,1));
set(h_a1,'XLim',[0 1200],'YLim',[-20 120]);
%xlabel(h_a1,'time (ms)');
%ylabel(h_a1,'membrane potential (mV)');

h_f3 = figure; 
h_a3 = axes;

ca_sf = 1e-3;
line('Parent',h_a3,'XData',t,'YData',y(:,7));
line('Parent',h_a3,'XData',t,'YData',ca_sf*y(:,8),'Color','r');
set(h_a3,'XLim',[0 1200],'YLim',[0 0.4]);

[t,y] = prsolve(500,[2.5 0.0 2.1]);

h_f2 = figure; 
h_a2 = axes;

line('Parent',h_a2,'XData',t,'YData',y(:,1));
set(h_a2,'XLim',[0 500],'YLim',[-20 120]);
xlabel(h_a2,'time (ms)');
ylabel(h_a2,'membrane potential (mV)');

h_f4 = figure; 
h_a4 = axes;

line('Parent',h_a4,'XData',t,'YData',y(:,7));
line('Parent',h_a4,'XData',t,'YData',ca_sf*y(:,8),'Color','r');
set(h_a4,'XLim',[0 500],'YLim',[0 0.4]);
xlabel(h_a4,'time (ms)');
ylabel(h_a4,'IK-AHP activation and Ca');

%print(handles.figure1,'-depsc2','pr_modes.eps');

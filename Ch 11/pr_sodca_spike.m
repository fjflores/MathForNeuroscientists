%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%
%  requires: prsolve.m
%

[t,y] = prsolve(100,[0 0.68 0]);

h_f1 = figure; 
h_a1 = axes;
line('Parent',h_a1,'XData',t,'YData',y(:,1));
line('Parent',h_a1,'XData',t,'YData',y(:,2),'Color','r');
set(h_a1,'YLim',[-20 120]);
xlabel(h_a1,'time (ms)');
ylabel(h_a1,'membrane potential (mV)');

%print(handles.figure1,'-depsc2','pr_sodca_spike.eps');

%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  requires destexhy2.m dest_f2.m dest_f2b.m dest_f2c.m

 dt = 0.05;
 v = destexhy2(dt,10,1010,1020,5);
 N = length(v);
 t = (1:N)*dt;
 
 h_f1 = figure; 
 h_a1 = subplot(2,2,1);
 line('Parent',h_a1,'XData',t,'YData',v)
 set(h_a1,'XLim',[200 800],'YLim',[-80 50]);
 xlabel(h_a1,'time (ms)');
 ylabel(h_a1,'Vm (mV)');
 
 disp('computing 1st panel...');
 [n_mean s_mean m_isis c_isis] = dest_f2;
 h_a2 = subplot(2,2,2);
 line('Parent',h_a2,'XData',m_isis,'YData',c_isis,'Marker','o',...
    'LineStyle','none','MarkerFaceColor','k','MarkerSize',4);
 set(h_a2,'XLim',[0 1500],'YLim', [0 1.1]);
 xlabel(h_a2,'mean ISI (ms)');
 ylabel(h_a2,'CV');
 

 disp('computing 2nd panel...');
 [currents n_mean n_mean2 n_mean3] = dest_f2b;
 h_a3 = subplot(2,2,3);
 line('Parent',h_a3,'XData',currents,'YData',n_mean,...
     'Marker','o','MarkerFaceColor','k','MarkerSize',4);
 line('Parent',h_a3,'XData',currents,'YData',n_mean2,...
     'Marker','o','MarkerFaceColor','r','MarkerSize',4);
 line('Parent',h_a3,'XData',currents,'YData',n_mean3,...
     'Marker','o','MarkerFaceColor','w','MarkerSize',4);
 set(h_a3,'YLim',[0 300]);
 xlabel(h_a3,'current (muA/cm2)');
 ylabel(h_a3,'firing frequency (spk/s)');
 
 disp('computing 3rd panel...');
 [currents n_mean n_mean2 n_mean3] = dest_f2c;
 h_a4 = subplot(2,2,4);
 line('Parent',h_a4,'XData',currents,'YData',n_mean,...
     'Marker','o','MarkerFaceColor','k','MarkerSize',4);
 line('Parent',h_a4,'XData',currents,'YData',n_mean2,...
     'Marker','o','MarkerFaceColor','r','MarkerSize',4);
 line('Parent',h_a4,'XData',currents,'YData',n_mean3,...
     'Marker','o','MarkerFaceColor','w','MarkerSize',4);
 set(h_a4,'YLim',[0 300]);
 xlabel(h_a4,'current (muA/cm2)');
 
%print(handles.figure1,'-depsc2','destex_f2.eps');

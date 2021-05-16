%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

stdobsy2 = load('stdcol_obs_y2');
h_f1 = figure; 
h_a1 = axes;
line('Parent',h_a1,'XData',stdobsy2(:,1),'YData',stdobsy2(:,2),...
    'Marker','o','MarkerFaceColor','k','MarkerSize',3,'LineStyle','-','Color','k');
set(h_a1,'XLim',[360 800]);
set(h_a1,'XTick',[400:100:800]);
xlabel('wavelength (nm)');
ylabel('Spectral luminous efficiency, v(\lambda)');

%print(handles.figure1,'-depsc2','stdobs_plot.eps');

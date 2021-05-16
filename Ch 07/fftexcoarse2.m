%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

%
% fftexcoarse
%

N = 10;
T = 10;
dt = T/N;
t = 0:dt:T-dt;
u = t.*(4-t).*(10-t);

h_f1 = figure;
h_a1 = axes;
line('Parent',h_a1,'XData',t,'YData',u,...
    'Color','k','Marker','x','LineStyle','-');
axis(h_a1,'tight'); 
box(h_a1,'off');
xlabel(h_a1,'t  (s)','fontsize',10)
ylabel('u','fontsize',10)

h_f2 = figure;
h_a2 = axes;
uhat = fft(u);
line('Parent',h_a2,'XData',[0:N-1],'YData',real(uhat),...
    'Color','k','Marker','x','LineStyle','-');
line('Parent',h_a2,'XData',[0:N-1],'YData',imag(uhat),...
    'Color','r','Marker','o','LineStyle','-');
legend(h_a2,'real part','imaginary part','location','best')
xlabel(h_a2,'index','fontsize',10)
ylabel(h_a2,'uhat','fontsize',10)
axis(h_a2,'tight');
box(h_a2,'off');

h_f3 = figure;
h_a3 = axes;
f = [0:floor(N/2) -(ceil(N/2)-1:-1:1)]/T;
f_shift = circshift(f,[0 (N-2)/2]);
ruhat_shift = circshift(real(uhat),[0 (N-2)/2]);
iuhat_shift = circshift(imag(uhat),[0 (N-2)/2]);
line('Parent',h_a3,'XData',f_shift,'YData',ruhat_shift,...
    'Color','k','Marker','x','LineStyle','-');
line('Parent',h_a3,'XData',f_shift,'YData',iuhat_shift,...
    'Color','r','Marker','o','LineStyle','-');

%legend(handles.axes3,'real part','imaginary part','location','best')
xlabel(h_a3,'\omega  (Hz)','fontsize',10)
%ylabel('uhat','fontsize',10)
axis(h_a3,'tight');
box(h_a3,'off');

%print(h_f1,'-depsc2','fftexcoarse2.eps');

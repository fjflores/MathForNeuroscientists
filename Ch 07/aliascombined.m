%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

x = 0:0.01:1;
y = sin(2*pi*x);
z = sin(2*pi*3*x);

h_f1 = figure;
h_a1 = axes;
line('Parent',h_a1,'XData',x,'YData',y,'Color','k');
line('Parent',h_a1,'XData',x,'YData',z,'Color','r');
line('Parent',h_a1,'XData',[0 0.5 1],'YData',[0 0 0],'Marker','x','MarkerSize',8,'Color','k');

set(h_a1,'TickDir','out');
xlabel('time (s)');
ylabel('response (arbitrary units)');

h_f2 = figure;
h_a2 = axes;
line('Parent',h_a2,'XData',x,'YData',y+z,'Color','k');
line('Parent',h_a2,'XData',x,'YData',2*y,'Color','r');
line('Parent',h_a2,'XData',[0 0.5 1],'YData',[0 0 0],'Marker','x','MarkerSize',8,'Color','k');

set(h_a2,'TickDir','out');

%%%

h_f3 = figure;
h_a3 = axes;
N = 16;
T = 10;
dt = T/N;
t = 0:dt:T-dt;
u = sin(2*pi*t);
line('Parent',h_a3,'XData',t,'YData',u,'Color','k');

h_f4 = figure;
h_a4 = axes;
f = (0:N/2)/T;
c = fft(u)/N;
line('Parent',h_a4,'XData',f,'YData',abs(c(1:1+N/2)),'Color','k');

N =32;

dt = T/N;
t = 0:dt:T-dt;
u = sin(2*pi*t);
line('Parent',h_a3,'XData',t,'YData',u,'Color','r');
axis(h_a3,'tight'); 
box(h_a3,'off');
legend(h_a3,'N=16','N=32','location','best');
legend('boxoff');
xlabel(h_a3,'t  (s)','fontsize',10)
ylabel(h_a3,'u','fontsize',10)

f = (0:N/2)/T;
c = fft(u)/N;
line('Parent',h_a4,'XData',f,'YData',abs(c(1:1+N/2)),'Color','r');
axis(h_a4,'tight');
box(h_a4,'off');
legend(h_a4,'N=16','N=32','location','best');
legend('boxoff');
xlabel(h_a4,'\omega (Hz)','fontsize',10)
ylabel(h_a4,'|c|','fontsize',10)

%print(h_f1,'-depsc2','aliascombined.eps');

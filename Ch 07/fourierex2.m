%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

% fourierex

x = -0.5:.005:0.5;

h_f =figure; 
h_a1 = subplot(1,2,1);
h_a2 = subplot(1,2,2);

f = ones(size(x))/pi;
for n=1:10,
    f = f + 2*sin(n)*cos(n*2*pi*x)/pi/n;
end
line('Parent',h_a1,'XData',x,'YData',f,'Color','r');
f = ones(size(x))/pi;
for n=1:100,
    f = f + 2*sin(n)*cos(n*2*pi*x)/pi/n;
end
line('Parent',h_a1,'XData',x,'YData',f,'Color','k');
axis(h_a1,'tight');
box(h_a1,'off');
legend(h_a1,'N=10','N=100')
xlabel(h_a1,'x','fontsize',14)

f = ones(size(x));
for n=1:10,
    f = f + 2*cos(n*2*pi*x);
end
line('Parent',h_a2,'XData',x,'YData',f,'Color','r');
f = ones(size(x));
for n=1:100,
    f = f + 2*cos(n*2*pi*x);
end
line('Parent',h_a2,'XData',x,'YData',f,'Color','k');
axis(h_a2,'tight');
box(h_a2,'off');
legend(h_a2,'N=10','N=100');
xlabel(h_a2,'x','fontsize',14)


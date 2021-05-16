%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

h_f = figure;
h_a1 = subplot(2,3,1);
h_a2 = subplot(2,3,2);
h_a3 = subplot(2,3,3);
h_a4 = subplot(2,3,4);
h_a5 = subplot(2,3,5);
h_a6 = subplot(2,3,6);

%normal distribution
x = -3.5:0.05:3.5;
y = normpdf(x,0,1);
line('Parent',h_a5,'XData',x,'YData',y);
title(h_a5,'normal');

%exponential distribution
x = 0:0.05:3;
y = exppdf(x,1);
line('Parent',h_a3,'XData',x,'YData',y);
title(h_a3,'exponential');

%binomial 1 distribution
n = 100;
p = 0.1;

x = 0:1:20;
y = binopdf(x,n,p);
line('Parent',h_a1,'XData',x,'YData',y,'LineStyle','-','Marker','x');
title(h_a1,'binomial');

%binomial 2 distribution
n = 12;
p = 0.83;

x = 0:1:20;
y = binopdf(x,n,p);
line('Parent',h_a2,'XData',x,'YData',y,'LineStyle','-','Marker','x');
title(h_a2,'binomial');

%poisson distribution
lambda = 100*0.1;

x = 0:1:20;
y = poisspdf(x,lambda);
line('Parent',h_a4,'XData',x,'YData',y,'LineStyle','-','Marker','x');
title(h_a4,'poisson');

%gamma distribution
x = 0:0.05:10;
y = gampdf(x,3,1);
line('Parent',h_a6,'XData',x,'YData',y,'LineStyle','-');
title(h_a6,'gamma');


%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

m0 = 2;
m1 = 4;

x0 = 0:6;
y0 = poisspdf(x0,m0);
h_f1 = figure; 
h_a1 = subplot(2,2,1);
h_a2 = subplot(2,2,2);
h_a3 = subplot(2,2,3);
h_a4 = subplot(2,2,4);

line('Parent',h_a1,'XData',x0,'YData',y0,'LineStyle','--','Marker','o');

x1 = 0:9;
y1 = poisspdf(x1,m1);
line('Parent',h_a1,'XData',x1,'YData',y1,'LineStyle','-','Marker','s');
xlabel(h_a3,'number of spikes');
ylabel(h_a1,'probability of occurence');

%for thresholds k_0 of zero to six
p_fa = [1 1 - poisscdf(0:5,m0)];
p_d  = [1 1 - poisscdf(0:5,m1)];

line('Parent',h_a2,'XData',p_fa,'YData',p_d,'Marker','^');

m0 = 4;
m1 = 10;

%for thresholds k_0 of zero to six
p_fa = [1 1 - poisscdf(0:8,m0)];
p_d  = [1 1 - poisscdf(0:8,m1)];

line('Parent',h_a2,'XData',p_fa,'YData',p_d,'Marker','*');

line('Parent',h_a2,'XData',[0 1],'Ydata',[0 1],'Linestyle','--');
set(h_a2,'XLim',[0 1],'YLim',[0 1]);
xlabel(h_a2,'probability of false-alarm');
ylabel(h_a2,'probability of detection');

m0 = 30;
m1 = 50;
sigma = 10;

x0 = 0:70;
y0 = normpdf(x0,m0,sigma);

line('Parent',h_a3,'XData',x0,'YData',y0);

x1 = 10:90;
y1 = normpdf(x1,m1,sigma);

line('Parent',h_a3,'XData',x1,'YData',y1);
ylabel(h_a3,'probability of occurence');

%compute the 
d = (m1-m0)/sigma;

xi = -3:0.1:3;
p_fa = 1 - normcdf(xi);
p_d =  1 - normcdf(xi - d);

line('Parent',h_a4,'XData',p_fa,'YData',p_d);

%compute the case d = 1
d = 1;
p_d = 1 - normcdf(xi - d);

line('Parent',h_a4,'XData',p_fa,'YData',p_d,'Color','r');

line('Parent',h_a4,'XData',[0 1],'Ydata',[0 1],'Linestyle','--');
set(h_a4,'XLim',[0 1],'YLim',[0 1]);
xlabel(h_a4,'probability of false-alarm');
ylabel(h_a4,'probability of detection');

%print(handles.figure1,'-depsc2','poiss1.eps');

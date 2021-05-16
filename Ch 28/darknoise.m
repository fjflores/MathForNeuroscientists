%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

%this m-file reproduces the Barlow Fig. 2

function darknoise

%second data set
ss_dat2 = [ ...
    23.5 0.0;...
    37.1 0.0;...
    58.5 12.0;...
    92.9 44.0;...
    148.6 94.0;...
    239.3 100.0];

%log base 10 of the mean number of photons
x_hsp2s = log10(ss_dat2(:,1));

%convert to percents
y_hsp2s = ss_dat2(:,2)/100;

h_f1 = figure; 
h_a1 = subplot(1,2,1);
line('Parent',h_a1,'XData',x_hsp2s,'YData',y_hsp2s,...
    'Marker','x', 'LineStyle','none');

%compare with the Barlow model
x_hb = 1.3:0.01:2.7;
y_hb = p_c(x_hb,21,0.13,8.9);
 
line('Parent',h_a1,'XData',x_hb,'YData',y_hb,'Color','r');


%Barlow model and fig. 1 
x_hb1 = 1.3:0.01:2.3;
y_hb1 = p_c(x_hb1,19,0.13,9.8);
 
%note: we use 9.8 instead of 8.9 for the dark noise since it obviously 
%works much better with the data. 

h_a2 = subplot(1,2,2);
line('Parent',h_a2,'XData',x_hb1,'YData',y_hb1,'Color','r');

%Barlow model and fig. 1 
x_hb2 = 1.3:0.01:2.3;
y_hb2 = p_c(x_hb2,17,0.13,9.8);
 
line('Parent',h_a2,'XData',x_hb2,'YData',y_hb2,'Color','k');

%coordinates measured in points, relative to (1.5 0.2)
%0.2 in x direction is 28.88 pts
%0.2 in y direction is 28.58 pts
cf_x = 0.2/28.88;
cf_y = 0.2/28.58;

ptsx = [ 18.16; 38.1; 61.62; 82.16; 105.38 ];
ptsy = [ -7.74; 9.82; 54.48; 97.64; 107.17 ];

ptsx_c = 1.5 + ptsx*cf_x;
ptsy_c = 0.2 + ptsy*cf_y;

line('Parent',h_a2,'XData',ptsx_c,'YData',ptsy_c,...
    'LineStyle','none','Marker','x');
set(h_a2,'XLim',[1.3 2.3]);

ptsy2 = [ 12.21; 44.06; 83.65; 103.89; 111.04];
ptsy2_c = 0.2 + ptsy2*cf_y;
line('Parent',h_a2,'XData',ptsx_c,'YData',ptsy2_c,...
    'LineStyle','none','Marker','+');

xlabel(h_a1,'log10 nphotons');
ylabel(h_a1,'detection probability');

%print(handles.figure1,'-depsc2','darknoise.eps');

function y = p_c(x,c,alpha,d_p)
%
% function y = p_c(x, c, alpha, d_p)
%  
% cumulative poisson function
% for photon counts
%  
% x is the log10 photon count at the retina
% alpha is the absorption yield factor
% c is the threshold for seeing 
% d_p is the dark photon number
  
%convert back x to linear space
xl = 10.^x;

%corresponding poisson parameter

p_p = alpha*xl + d_p;

y = 1 - poisscdf(c,p_p);






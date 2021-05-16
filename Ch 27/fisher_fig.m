%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

%%
%rotation angle in radians
alpha = 20* (pi/180);
%rotation matrix
R = [cos(alpha) -sin(alpha);sin(alpha) cos(alpha)];

x = -4:0.2:4;
y = -2:0.2:4;
n_x = length(x);
n_y = length(y);
C0 = R'*[1 -0.3;-0.3 1]*R;
C0i = inv(C0);
dC0 = det(C0);
c_z0 = 1/(2*pi*sqrt(dC0));

n0 = [-1; 1];
z = zeros(n_y,n_x);
for i = 1:n_x
    for j = 1:n_y
        v = [x(i); y(j)];
        vz = v - n0;
        z(j,i) = c_z0*exp(-0.5*vz'*C0i*vz);
    end;
end;

h_f1 = figure; 
h_a1 = axes;
contour(h_a1,x,y,z,'EdgeColor','black');
hold on;

%%%%%%%%%%%%

x = -4:0.2:4;
y = -4:0.2:4;
n_x = length(x);
n_y = length(y);
C1 = R'*[1 0.8;0.8 1]*R;
C1i = inv(C1);
dC1 = det(C1);
c_z1 = 1/(2*pi*sqrt(dC1));

n1 = [1; 1];
z = zeros(n_x,n_x);
for i = 1:n_x
    for j = 1:n_y
        v = [x(i); y(j)];
        vz = v - n1;
        z(j,i) = c_z1*exp(-0.5*vz'*C1i*vz);
    end;
end;

%mesh(h_a2,x,x,z,'EdgeColor','black');
contour(h_a1,x,x,z,'EdgeColor','red');
set(h_a1,'XLim',[-4 4],'YLim',[-4 4],'TickDir','out');


%compute the best direction according to the Fisher method
wopt = inv(C1+C0)*(n1-n0);
wopt = wopt/norm(wopt);

line('Parent',h_a1,'XData',[0 wopt(1)],'YData',[0 wopt(2)],'Color','b');
line('Parent',h_a1,'XData',[0 wopt(1)],'YData',[0 wopt(2)],'Color','b');
xlabel(h_a1,'Data dimension 1');
ylabel(h_a1,'Data dimension 2');

%%
mu0 = wopt'*n0;
sig0 = wopt'*C0*wopt;

dx0 = 6*sig0/100;
x0 = mu0-3*sig0:dx0:mu0+3*sig0;
y0 = normpdf(x0,mu0,sig0);

h_f2 = figure;
h_a2 = subplot(2,1,1);
line('Parent',h_a2,'XData',x0,'YData',y0);

mu1 = wopt'*n1;
sig1 = wopt'*C1*wopt;

dx1 = 6*sig1/100;
x1 = mu1-3*sig1:dx1:mu1+3*sig1;
y1 = normpdf(x1,mu1,sig1);
line('Parent',h_a2,'XData',x1,'YData',y1,'Color','r');
set(h_a2,'XLim',[-4 4]);
xlabel(h_a2,'Projection on optimal discrimination direction');
ylabel(h_a2,'Probability density');

wper = [-wopt(2); wopt(1)];

mu0p = wper'*n0;
sig0p = wper'*C0*wper;

dx0p = 6*sig0p/100;
x0p = mu0p-3*sig0p:dx0p:mu0p+3*sig0p;
y0p = normpdf(x0p,mu0p,sig0p);
h_a3 = subplot(2,1,2);
line('Parent',h_a3,'XData',x0p,'YData',y0p);

mu1p = wper'*n1;
sig1p = wper'*C1*wper;

dx1p = 6*sig1p/100;
x1p = mu1p-3*sig1p:dx1p:mu1p+3*sig1p;
y1p = normpdf(x1p,mu1p,sig1p);
line('Parent',h_a3,'XData',x1p,'YData',y1p,'Color','r');
set(h_a3,'XLim',[-4 4]);
xlabel(h_a3,'Orthogonal projection');
ylabel(h_a3,'Probability density');

%print(handles.figure1,'-depsc2','fisher_fig.eps');



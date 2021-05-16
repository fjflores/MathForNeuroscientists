%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%

%plot the properties of the ISI distribution and forward recurrence time
%of the equilibrium gamma order 2 process

%mean isi
misi = 25; %ms

%corresponding gamma rate
rho = 2/misi; %kHz

%b parameter of matlab
b = 1/rho;

dx = 0.5;
x = 0:dx:100;
y1 = rho*(rho*x).*exp(-rho*x);
h = figure;
plot(x,y1,'k');

z1 = (rho/2)*(rho*x+1).*exp(-rho*x);
hold on;
plot(x,z1,'r');
xlabel('time (ms)');
ylabel('probability density');
set(gca,'TickDir','out');

%uncomment to print figure
%print(h,'-depsc2','frecfig_ex.eps');

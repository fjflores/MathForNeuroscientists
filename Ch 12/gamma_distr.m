%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

h_f = figure; 
h_a = axes;

m_inter = 20e-3; %mean isi in sec
t_vect = (0:0.1:100)*1e-3; %in sec

%Poisson
rho = 1/m_inter; %firing rate in Hz = 1/sec

p_dens = rho*exp(-rho*t_vect);

line('Parent',h_a,'XData',t_vect*1e3,'YData',p_dens);

%gamma of order 2
n = 2;
rho = n/m_inter;

p_dens = rho*(rho*t_vect).^(n-1).*exp(-rho*t_vect)/factorial(n-1);

line('Parent',h_a,'XData',t_vect*1e3,'YData',p_dens);


%gamma of order 5
n = 5;
rho = n/m_inter;

p_dens = rho*(rho*t_vect).^(n-1).*exp(-rho*t_vect)/factorial(n-1);

line('Parent',h_a,'XData',t_vect*1e3,'YData',p_dens);

%gamma of order 10
n = 10;
rho = n/m_inter;

p_dens = rho*(rho*t_vect).^(n-1).*exp(-rho*t_vect)/factorial(n-1);

line('Parent',h_a,'XData',t_vect*1e3,'YData',p_dens);

axes(h_a);
xlabel('ISI (msec)');
ylabel('Probability density');


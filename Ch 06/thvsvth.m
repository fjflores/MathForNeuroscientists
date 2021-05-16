%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% th vs vth
%
% contrast the spectra of the continuous and N-compartment cables
%
%  e.g.,  thvsvth(100)
%

function thvsvth(N)
ell = 1000;
n = 0:N-1;
vth = -(n*pi/ell).^2;
plot(n,vth,'k+')
th = -4*(N*sin(n*pi/2/N)/ell).^2;
hold on
plot(n,th,'ro')
hold on
legend('\vartheta_n','\theta_n')
hold off
box off
xlabel('index n','fontsize',14)
ylabel('eigenvalue  ( \mum^{-2})','fontsize',14)

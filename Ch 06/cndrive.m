%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% contrast accuracy of trap and be on passive fiber
%
% cndrive.m
%
function cndrive

dt = [1e-1 5e-2 1e-2 5e-3 1e-3];
for j=1:5,
    [verrbe(j),verr(j)] = cnpfib(.001,1,100,dt(j),10,0);
end
figure(2)
loglog(dt,verrbe,'kx-','linewidth',2)
hold on
loglog(dt,verr,'ro-','linewidth',2)
box off
legend('Backward Euler','Trapezoid','location','best')
legend('boxoff')
xlabel('dt  (ms)','fontsize',14)
ylabel('maximum absolute error  (mV)','fontsize',14)
set(gca,'tickdir','out')

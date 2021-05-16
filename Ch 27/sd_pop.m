%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

%average SNR for single neurons
d_single = 1;

%no correlation curve is diagonal
x = 1:50;
y_noc = 1:50;

h_f1 = figure;
h_a1 = axes;
line('Parent',h_a1,'XData',x,'Ydata',y_noc,'LineStyle','--');

%correlation among neurons
c = 0.1;

rho = [0 0.5 1];
col = ['r' 'g' 'b'];
for i =1:length(rho)
    d_pop2 = d_single^2 * (x/(1+rho(i))) .* ( (1./(x*c - c +1)) + (1/(1-c))*rho(i));
    h_rho(i) = line('Parent',h_a1,'XData',x,'YData',d_pop2,'Color',col(i));
end;

set(h_a1,'XLim',[0 50]);
xlabel(h_a1,'number of neurons');
ylabel(h_a1,'population SNR/mean single neuron SNR');
legend(h_rho,'r = 0','r = 0.5','r = 1');

%print(handles.figure1,'-depsc2','sd_pop.eps');


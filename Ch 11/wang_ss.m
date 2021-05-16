%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

v = -120:1:0;

%steady-state activation curve of low-threshold Ca current
sinf = 1./(1+exp(-(v+65)/7.8));

%steady-state inactivation of low-threshold Ca current
thh = -81;
kh = 6.25;
hinf = 1./(1+exp((v-thh)/kh));

%steady-state activation of H current
Hinf = 1./(1+exp((v+69)/7.1));

h_f1 = figure;
h_a1 = axes;
%plot results
line('Parent',h_a1,'XData',v,'YData',sinf);
line('Parent',h_a1,'XData',v,'YData',hinf,'LineStyle',':');
line('Parent',h_a1,'XData',v,'YData',Hinf,'Color','r');
set(h_a1,'XLim',[-120 -10]);
xlabel(h_a1,'membrane potential (mV)');

%time constant of (de-)inactivation of T-type calcium current
tauh = hinf.*exp((v+162.3)/17.8) + 20;
tauh = tauh/2; %temperature scaling factor

%time constant of H current
tauH = 1000./(exp((v+66.4)/9.3)+exp(-(v+81.6)/13));

h_f2 = figure; 
h_a2 = axes;
line('Parent',h_a2,'XData',v,'YData',tauh);
set(h_a2,'XLim',[-120 -10]);
ylabel(h_a2,'I_T inactivation time constant (ms)');

yyaxis right
line('Parent',h_a2,'XData',v,'YData',tauH,'Color','r');
set(h_a2,'XLim',[-120 -10],'Color','none');
set(h_a2,'YLim',[0 1500]);
xlabel(h_a2,'membrane potential (mV)');
ylabel(h_a2,'I_H activation time constant (ms)');
h_a2.YColor = 'r';

%print(handles.figure1,'-depsc2','wang_ss.eps');


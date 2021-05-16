%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%
%  requires prsolve.m 
%

function pr_complex

[t,y] = prsolve(1500,[0.75 0.0 2.1]);

h_f1 = figure; 
h_a1 = axes;

line('Parent',h_a1,'XData',t,'YData',y(:,1));
line('Parent',h_a1,'XData',t,'YData',y(:,2),'Color','r');
set(h_a1,'XLim',[970 1000],'YLim',[-10 90]);
ylabel(h_a1,'membrane potential (mV)');

%compute the corresponding currents
z = currents_from(y,2.1);

h_f2 = figure; 
h_a2 = axes;

line('Parent',h_a2,'XData',t,'YData',z(:,1));
line('Parent',h_a2,'XData',t,'YData',z(:,2));
line('Parent',h_a2,'XData',t,'YData',z(:,3));
set(h_a2,'XLim',[970 1000],'YLim',[-1800 500]);
xlabel(h_a2,'time (ms)');
ylabel(h_a2,'current (mA/cm2)');

h_f3 = figure; 
h_a3 = axes;

line('Parent',h_a3,'XData',t,'YData',z(:,4));
line('Parent',h_a3,'XData',t,'YData',z(:,5));
line('Parent',h_a3,'XData',t,'YData',z(:,6));
line('Parent',h_a3,'XData',t,'YData',50*z(:,7),'Color','r');
set(h_a3,'XLim',[970 1000],'YLim',[-800 800]);

h_f4 = figure; 
h_a4 = axes;

line('Parent',h_a4,'XData',t,'YData',z(:,8));
set(h_a4,'XLim',[970 1000]); %,'YLim',[-20 120]);
xlabel(h_a4,'time (ms)');
ylabel(h_a4,'current (mA/cm2)');

%print(handles.figure1,'-depsc2','pr_complex.eps');
return;

function z = currents_from(y,gc)

%maximal conductances in mS/cm^2
gL = 0.1; %leak
gNa = 30; %fast sodium
gKDR = 15; %delayed rectifier
gCa = 10; %fast calcium current
gKAHP = 0.8; %calcium dependent potassium current (slow)
gKC = 15; %voltage and calcium dependent potassium current (fast)

%reversal potentials (in mV)
VNa = 120; VCa = 140; VK = -15; VL = 0;

%fraction of cable length assigned to soma (1-p for dendrite)
p = 0.5;

%somatic leak current
Ils = gL*(y(:,1)-VL);

%steady-state sodium activation (instantaneous)
minf = am(y(:,1))./(am(y(:,1))+bm(y(:,1)));

%sodium current (y(3) is h, inactivation of sodium current)
INa = gNa*minf.^2.*y(:,3).*(y(:,1)-VNa);

%delayed rectifier current (y(4) is n, activation of DR)
IKDR = gKDR*y(:,4).*(y(:,1)-VK);

%dendritic leak current
Ild = gL*(y(:,2)-VL);

%dendritic calcium current (y(5) is s activation variable)
ICa = gCa*y(:,5).^2.*(y(:,2)-VCa);

%voltage and calcium dependent K current (y(6) is c activation variable,
%y(8) is Ca)
IKC = gKC*y(:,6).*min(y(:,8)/250,1).*(y(:,2)-VK);

%calcium dependent K current (y(7) is q activation variable)
IKAHP = gKAHP*y(:,7).*(y(:,2)-VK);

Idtos = (gc/p)*(y(:,2)-y(:,1));

z = [Ils INa IKDR Ild ICa IKC IKAHP Idtos];

return;

%
%For following rate constants, see eq. 6 of [PR94] and erratum
%

%forward rate constant for fast sodium
function val = am(v)
val = 0.32*(13.1-v)./(exp((13.1-v)/4)-1);
return;

%backward rate constant for fast sodium
function val = bm(v)
val = 0.28*(v-40.1)./(exp((v-40.1)/5)-1);
return;
%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

%%
%mean excitatory conductance
ge0 = 0.012; %mu S 

%relaxation time constant and steady-state SD
tau_e = 2.7; %ms
sig_e = 0.003; %mu S

%time step and time vector
dt = 0.5; %ms

%this gives 16 segments to average over for the 
%power spectrum
n_seg = 2;
seg_len = 2048;
n_v = 0:(n_seg*seg_len-1);
t_v = n_v*dt;
n_t_v = length(t_v);

%numerical simulation factors
ef_1 = exp(-dt/tau_e);
ef_2 = sqrt(1-exp(-2*dt/tau_e));

%initial value of the random conductance variable
x0 = 0;

%random white noise increment
w_ve = normrnd(0,1,1,n_t_v);
u_ve = sig_e*ef_2*w_ve;

%save space for random vector
x_ve = zeros(1,n_t_v);

%initial condition
x_ve(1) = x0;

for i =2:n_t_v
    %stochastic update
    x_ve(i) = x_ve(i-1)*ef_1 + u_ve(i);
end;

%naive method, periodogram
%[] = rectangular window
%[] = standard length, equal to vector length
[p_perio,f_perio] = periodogram(x_ve,[],[],2000,'twosided');
n_max = (n_seg*seg_len)/2 + 1;

h_f1 = figure;
h_a1 = subplot(3,1,1);
line('Parent',h_a1,'XData',f_perio(1:n_max),'YData',p_perio(1:n_max));

%%
%2048 = window length (hamming)
%1024 = overlap (50 percent)
%2048 = number of samples for fft
%2000 = sampling frequency (0.5 ms = 2kHz)
[p_welch,f_welch] = pwelch(x_ve,seg_len,seg_len/2,seg_len,2000,'twosided');
h_a2 = subplot(3,1,2);
line('Parent',h_a2,'XData',f_welch(1:seg_len/2 + 1),'YData',p_welch(1:seg_len/2+1));

%32 means 2*32 -1 data windows
%[] = use standard length, whole data
[p_mtm,p_c,f_mtm] = pmtm(x_ve,14,[],2000,'twosided');
h_a3 = subplot(3,1,3)
line('Parent',h_a3,'XData',f_mtm(1:n_max),'YData',p_mtm(1:n_max));
%line('Parent',h_a3,'XData',f_mtm(1:n_max),'YData',p_c(1:n_max,1),'LineStyle',':');
%line('Parent',h_a3,'XData',f_mtm(1:n_max),'YData',p_c(1:n_max,2),'LIneStyle',':');

%%
top_ylim = 10e-8;
%theoretical power spectrum
f = 0:0.1:200;
pfe = 2*sig_e^2*tau_e*1e-3./(1 + (2*pi*f*tau_e*1e-3).^2);
line('Parent',h_a3,'XData',f,'YData',pfe,'Color','r');
xlabel(h_a3,'frequency (Hz)');
ylabel(h_a3,'Power density (\mu S^2 / Hz)');
set(h_a3,'XLim',[0 200],'YLim',[0 top_ylim]);

line('Parent',h_a2,'XData',f,'YData',pfe,'Color','r');
xlabel(h_a2,'frequency (Hz)');
ylabel(h_a2,'Power density (\mu S^2 / Hz)');
set(h_a2,'XLim',[0 200],'YLim',[0 top_ylim]);

line('Parent',h_a1,'XData',f,'YData',pfe,'Color','r');
xlabel(h_a1,'frequency (Hz)');
ylabel(h_a1,'Power density (\mu S^2 / Hz)');
set(h_a1,'XLim',[0 200],'YLim',[0 top_ylim]);

%%
h_f2 = figure;
h_a4 = axes;
wind = hamming(2048);
line('Parent',h_a4,'XData',1:2048,'YData',wind);
set(h_a4,'XLim',[0 2049]);

%print(handles.figure1,'-depsc','power_spec.eps'); 

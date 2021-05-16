%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

%%
fc = 5.5; %Hz
fc1 = 2*pi*fc*1e-3; %converts to kHz and circular frequency

%sampling step in msec
dt = 1;
fs = 1/(dt*1e-3); %sampling frequency
fn = fs/2; %nyquist frequency
nlgn = 512;
t_lgn = (0:nlgn-1)*dt; %in msec
v_lgn1 = exp(-fc1*t_lgn);
v_lgn = t_lgn.*v_lgn1 - (fc1/2) * (t_lgn.^2) .* v_lgn1;

max_fr = 20;
v_lgn = max_fr*(v_lgn./max(v_lgn));

nfft = 1024;
nw = 32*nfft;
sigma_w = 0.25;
wh_n = sigma_w*randn(1,nw);

%the matlab definition of the convolution requires the lgn 
%filter to be first and the stimulus to be second. To get
%no phase delay we need to take the full wave form and truncate
%it, as seen by convolving v_lgn with a [1 0 0 ... 0] sequence
fr_v = conv(v_lgn,wh_n);
fr_v = fr_v(1:length(wh_n));

[gf_lgn,f_lgn] = tfestimate(wh_n,fr_v,nfft,nfft/2,nfft,fs,'twosided');

gt_lgn = ifft(gf_lgn);

h_f1 = figure; 
h_a1 = subplot(2,1,1);
line('Parent',h_a1,'XData',t_lgn,'YData',gt_lgn(1:nlgn));
line('Parent',h_a1,'XData',t_lgn,'YData',v_lgn,...
     'LineStyle','--','Color','r');
set(h_a1,'XLim',[0 t_lgn(256)],'YLim',[-20 30]);
xlabel(h_a1,'time (ms)');
ylabel(h_a1,'firing rate change (spk/s)');

%%
%threshold the firing rate to obtain a positive rate
fr_vth = max(fr_v,-100) + 100;

%generate Poisson spike train
vm_vect = zeros(size(fr_vth));
spk_vect = zeros(size(fr_vth));

vm_thres = exprnd(1);
vm_vect(1) = 0;
rho_vect = (fr_vth/1000); %in spk/msec

%save some time
dtrho_vect = dt*rho_vect;

for i = 2:length(fr_vth)
    vm_vect(i) = vm_vect(i-1) + dtrho_vect(i);
    if (vm_vect(i) > vm_thres)
      vm_vect(i) = 0;
      spk_vect(i) = 1;
      vm_thres = exprnd(1);
    end;
  end;

spk_vect = spk_vect*1000; %in spk/sec

h_f2 = figure;
h_a3 = subplot(2,1,1);
h_a4 = subplot(2,1,2);
nplot = 500;
line('Parent',h_a3,'XData',(1:nplot)*dt,'YData',wh_n(1:nplot));
line('Parent',h_a4,'XData',(1:nplot)*dt,'YData',fr_vth(1:nplot));

inds_spk = find(spk_vect(1:nplot) > 500);
for k = 1:length(inds_spk)
    line('Parent',h_a4,'XData',[inds_spk(k) inds_spk(k)]*dt,'YData',[0 30]);
end;
set(h_a4,'YLim',[0 200]);

[gf_spk,f_spk] = tfestimate(wh_n,spk_vect,nfft,nfft/2,nfft,fs,'twosided');

gt_spk = ifft(gf_spk);

figure(h_f1);
h_a2 = subplot(2,1,2);
line('Parent',h_a2,'XData',t_lgn,'YData',gt_spk(1:nlgn));
line('Parent',h_a2,'XData',t_lgn,'YData',v_lgn,...
     'LineStyle','--','Color','r');
set(h_a2,'XLim',[0 t_lgn(256)],'YLim',[-20 30]);
xlabel(h_a2,'time (ms)');
ylabel(h_a2,'firing rate change (spk/s)');

%%

%reverse correlation between stimulus and spike train obtained after 
%passing through a LGN filter. The white noise has a cut-off frequency
%of 50 Hz.

fc = 5.5; %Hz
fc1 = 2*pi*fc*1e-3; %converts to kHz and circular frequency

%sampling step in msec
dt = 1;
deltat = 1e-3; % in sec
fs = 1/(deltat); %sampling frequency
fn = fs/2; %nyquist frequency
nlgn = 512;
t_lgn = (0:nlgn-1)*dt; %in msec
v_lgn1 = exp(-fc1*t_lgn);
v_lgn = t_lgn.*v_lgn1 - (fc1/2) * (t_lgn.^2) .* v_lgn1;

max_fr = 20;
v_lgn = max_fr*(v_lgn./max(v_lgn));

%shift by 10 ms
%v_lgn = [ zeros(1,10) v_lgn(1:end-10)];

nfft = 1024;
nw = 32*nfft;
sigma_w = 0.25;

fc = 50; %cutoff frequency in Hz 

sigma = sqrt((sigma_w^2*nw)/(4*fc*deltat));
 
%typically, a random seed controls the state of the random number 
%generator. To make the following sequence reproducible it should
%be initialized here 

whitef = zeros(1,nw);
   
%zero frequency component
ai = randn(1) * sigma;
whitef(1) = complex(ai,0);

%all component except for nyquist
for i=2:nw/2
    if ( i/(nw*deltat) < fc )
        ai = randn(1) * sigma;
        bi = randn(1) * sigma;
        whitef(i) = complex(ai,bi);
        whitef(nw-i+2) = complex(ai,-bi);
    end;
end;

ai = randn(1) * sigma;
whitef(nw/2+1) = complex(ai,0);

wh_n = fft(whitef);
wh_n = wh_n/nw;
wh_n = real(wh_n);


%the matlab definition of the convolution requires the lgn 
%filter to be first and the stimulus to be second. To get
%no phase delay we need to take the full wave form and truncate
%it, as seen by convolving v_lgn with a [1 0 0 ... 0] sequence
fr_v = conv(v_lgn,wh_n);
fr_v = fr_v(1:length(wh_n));

[gf_lgn,f_lgn] = tfestimate(wh_n,fr_v,nfft,nfft/2,nfft,fs,'twosided');

gt_lgn = ifft(gf_lgn);

%figure; 
%plot(t_lgn,gt_lgn(1:nlgn)); hold on;
%plot(t_lgn,v_lgn,'r--');

%threshold the firing rate to obtain a positive rate
fr_vth = max(fr_v,-100) + 100;

%generate Poisson spike train
vm_vect = zeros(size(fr_vth));
spk_vect = zeros(size(fr_vth));

vm_thres = exprnd(1);
vm_vect(1) = 0;
rho_vect = (fr_vth/1000); %in spk/msec

%save some time
dtrho_vect = dt*rho_vect;

for i = 2:length(fr_vth)
    vm_vect(i) = vm_vect(i-1) + dtrho_vect(i);
    if (vm_vect(i) > vm_thres)
      vm_vect(i) = 0;
      spk_vect(i) = 1;
      vm_thres = exprnd(1);
    end;
  end;

spk_vect = spk_vect*1000; %in spk/sec

nplot = 500;
line('Parent',h_a3,'XData',(1:nplot)*dt,'YData',wh_n(1:nplot),'Color','r');
line('Parent',h_a4,'XData',(1:nplot)*dt,'YData',fr_vth(1:nplot),'Color','r');

inds_spk = find(spk_vect(1:nplot) > 500);
for k = 1:length(inds_spk)
    line('Parent',h_a4,'XData',[inds_spk(k) inds_spk(k)]*dt,'YData',[220 250],'Color','r');
end;

set(h_a3, 'XLim',[0 500],'YLim',[-1 1]);
xlabel(h_a3,'time (ms)');

set(h_a4, 'XLim',[0 500],'YLim',[0 250]);
xlabel(h_a4,'time (ms)');
ylabel(h_a4,'firing rate (spk/s)');

[gf_spk,f_spk] = tfestimate(wh_n,spk_vect,nfft,nfft/2,nfft,fs,'twosided');

inds_over = find( (f_spk > fc) & (f_spk < fs-fc) );
gf_spk(inds_over) = 0;

gt_spk = ifft(gf_spk);

%figure;
%plot(t_lgn,gt_spk(1:nlgn)); hold on;
%plot(t_lgn,v_lgn,'r--');

%%
h_f3 = figure;
h_a5 = axes;
line('Parent',h_a5,'XData',t_lgn,'YData',gt_spk(1:nlgn));
line('Parent',h_a5,'XData',t_lgn,'YData',v_lgn,...
     'LineStyle','--','Color','r');
set(h_a5,'XLim',[0 t_lgn(256)],'YLim',[-20 30]);
xlabel(h_a5,'time (ms)');
ylabel(h_a5,'firing rate change (spk/s)');

%print(handles.figure1,'-depsc2','lgn_revcor_wn3.eps');

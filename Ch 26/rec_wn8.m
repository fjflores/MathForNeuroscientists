%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

%%
%stimulus estimation

dt_ms = 0.5; %sampling step, milliseconds
dt_s = dt_ms*1e-3; %sampling step seconds
fs = 1/dt_s; %sampling frequency, Hz
ni = 64;
nfft = 8192;
nw = ni*nfft;
sigma_w = 1;

fc = 10; %cutoff frequency in Hz 

sigma = sqrt((sigma_w^2*nw)/(4*fc*dt_s));
 
%typically, a random seed controls the state of the random number 
%generator. To make the following sequence reproducible it should
%be initialized here 

whitef = zeros(1,nw);
   
%zero frequency component
ai = randn(1) * sigma;
whitef(1) = complex(ai,0);

%all component except for nyquist
for i=2:nw/2
    if ( i/(nw*dt_s) < fc )
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

%low pass filter the white noise
tau_ef = 20e-3; %in s
ms_max = 6*tau_ef;
npt1 = ceil(ms_max/dt_s);

%look for the closest power of two
p2 = ceil(log(npt1)/log(2));
npt2 = 2^p2;
t_ef = (0:npt2-1)*dt_s;
kt = exp(-t_ef/tau_ef);

%The function filter(kt,1,wh_n) yields the same value
wh_nf = conv(kt,wh_n)*dt_s;
wh_nf = wh_nf(1:length(wh_n));

mfr_v = [10:10:100]; %target firing rate
%debug
%mfr_v = [ 50 ];
mfr_v2 = zeros(size(mfr_v)); %true firing rate
nerr_v = zeros(size(mfr_v));

disp('processing 10 different firing rates. This may take a while...');
for mfr_ind = 1:length(mfr_v)
    info_str = sprintf('processing firing rate: %i spk/s',mfr_v(mfr_ind));
    disp(info_str);
    
    %scale to yield about mfr spk/s for each of the two spike trains
    mfr = mfr_v(mfr_ind); %in spk/s
    m_awh = 0.5*(mean(max(0,wh_nf)) - mean(min(0,wh_nf)));
    %the following scale factor will also affect the power spectrum of the
    %spike train
    scale_factor = mfr/m_awh;
    fr_vect = scale_factor*wh_nf;

    fr_p = max(fr_vect,0);
    fr_m = -min(fr_vect,0);

    %generate Poisson spike train
    %positive part of the stimulus
    vm_vect = zeros(size(fr_vect));
    spk_p = zeros(size(fr_vect));

    n_gamma = 1;
    vm_thres = gamrnd(n_gamma,1/n_gamma);
    vm_vect(1) = 0;
    rho_vect = (fr_p/1000); %in spk/msec

    %save some time
    dtrho_vect = dt_ms*rho_vect;

    for i = 2:length(fr_p)
        vm_vect(i) = vm_vect(i-1) + dtrho_vect(i);
        if (vm_vect(i) > vm_thres)
            vm_vect(i) = 0;
            spk_p(i) = 1;
            vm_thres = gamrnd(n_gamma,1/n_gamma);
        end;
    end;

    spk_p = spk_p*fs; %in spk/sec

    %generate Poisson spike train
    %negative part of the stimulus
    vm_vect = zeros(size(fr_vect));
    spk_m = zeros(size(fr_vect));
    
    vm_thres = gamrnd(n_gamma,1/n_gamma);
    vm_vect(1) = 0;
    rho_vect = (fr_m/1000); %in spk/msec

    %save some time
    dtrho_vect = dt_ms*rho_vect;

    for i = 2:length(fr_m)
        vm_vect(i) = vm_vect(i-1) + dtrho_vect(i);
        if (vm_vect(i) > vm_thres)
            vm_vect(i) = 0;
            spk_m(i) = 1;
            vm_thres = gamrnd(n_gamma,1/n_gamma);
        end;
    end;

    spk_m = spk_m*fs; %in spk/sec

    mfr_v2(mfr_ind) = 0.5*(mean(spk_p) + mean(spk_m));
    
    %combine the two spike trains
    spk_c = spk_p - spk_m;

    %compute the reconstruction filter
    [gspkf_lgn,f_lgn] = tfestimate(spk_c,wh_n,nfft,nfft/2,nfft,fs,'twosided');

    %set to zero irrelevant components, i.e., above the cutoff frequency 
    %of the stimulus
    inds_z = find( (f_lgn > fc) & (f_lgn < fs - fc) );
    gspkf_lgn(inds_z) = 0;

    gspkt_lgn = ifft(gspkf_lgn);

    nlgn = length(gspkt_lgn);

    %make sure we have a 1 x nlgn columns vector
    gspkt_lgn = gspkt_lgn(:)';

    %unwrap the filter
    gspkt_lgn2 = circshift(gspkt_lgn,[0  (nlgn/2-1)]);
    t_lgn2 = ((-nlgn/2+1):nlgn/2)*dt_s;

    %compute the reconstruction
    rec_v = fftfilt(gspkt_lgn2,spk_c);

    %compensate for the delay due to the negative components 
    rec_v2 = rec_v(nlgn/2:nw-nlgn/2);
    wh_n2 = wh_n(1:nw-nfft+1);

    %compute the normalized error
    wn_noise = wh_n2 - rec_v2;
    err_n = sqrt(mean((wn_noise).^2))/std(wh_n);
    nerr_v(mfr_ind) = err_n;
    
    if ( mfr == 50 )
        %plot stimulus and reconstruction
        n_plotpts = 2000;
        t_plot = (0:n_plotpts-1)*dt_ms;
        h_f1 = figure; 
        h_a1 = subplot(2,1,1);
        h_a2 = subplot(2,1,2);
        line('Parent',h_a1,'XData',t_plot,'YData',wh_n2(n_plotpts:2*n_plotpts-1),...
            'Color','k');
        line('Parent',h_a1,'XData',t_plot,'YData',rec_v2(n_plotpts:2*n_plotpts-1),...
            'Color','r');
        
        line('Parent',h_a2,'XData',t_plot,'YData',fr_p(n_plotpts:2*n_plotpts-1),...
            'Color','k');
        line('Parent',h_a2,'XData',t_plot,'YData',-fr_m(n_plotpts:2*n_plotpts-1),...
            'Color','r');
        
        line('Parent',h_a2,'XData',[t_plot(1) t_plot(n_plotpts)],'YData',[0 0],...
            'Color','k');
        xlabel(h_a2,'time (ms)');
        ylabel(h_a2,'firing rate (spk/s)');
        
        for i = n_plotpts:2*n_plotpts-1
            if ( spk_p(i) > fs/2 )
                t_c = t_plot(i-n_plotpts+1);
                line('Parent',h_a2,'XData',[t_c t_c],'YData',[400 500],'Color','k');
            end;
            if ( spk_m(i) > fs/2 )
                t_c = t_plot(i-n_plotpts+1);
                line('Parent',h_a2,'XData',[t_c t_c],'YData',[-500 -400],'Color','R');
            end;
        end;
        
        %compute and plot the SNR
        [pwh,fwh] = pwelch(wh_n,nfft,nfft/2,nfft,fs,'twosided');
        [pwn,fwn] = pwelch(wn_noise,nfft,nfft/2,nfft,fs,'twosided');

        %compute snr
        inds = find( (fwh < fc) | (fwh > fs-fc) );
        snr = zeros(size(pwh));
        snr(inds) = pwh(inds)./pwn(inds);
        h_f2 = figure; 
        h_a3 = axes;
        line('Parent',h_a3,'XData',fwh,'YData',snr);
        set(h_a3,'XLim',[0 15]);
        xlabel(h_a3,'frequency (Hz)');
        ylabel(h_a3,'SNR');
        
    end;
        
end;

%%
tau_omega = 2*pi*tau_ef*fc;
lambda_v = 2*mfr_v2; %spk/s

gamma_v = (pi^2/2)*(tau_ef*lambda_v/atan(tau_omega));

norm_err2 = (1 - (1/tau_omega)*(gamma_v./sqrt(1+gamma_v)).*atan(tau_omega./sqrt(1+gamma_v)));
norm_err = sqrt(norm_err2);

h_f3 = figure;
h_a4 = subplot(3,1,1);
line('Parent',h_a4,'XData',mfr_v2,'YData',nerr_v,...
    'Marker','s','Color','k','LineStyle','--');
line('Parent',h_a4,'XData',mfr_v2,'YData',norm_err,...
    'Marker','o','Color','r','LineStyle','-');
xlabel(h_a4,'firing rate (spk/s)');
ylabel(h_a4,'normalized error');

%simulate different gamma orders
gamman_v = [1 2 5 10 100]; %gamma order
nerr_v = zeros(size(gamman_v));
cv_v = zeros(size(gamman_v));

for g_ind = 1:length(gamman_v)
    %scale to yield about mfr spk/s for each of the two spike trains
    mfr = 50; %in spk/s
    m_awh = 0.5*(mean(max(0,wh_nf)) - mean(min(0,wh_nf)));
    %the following scale factor will also affect the power spectrum of the
    %spike train
    scale_factor = mfr/m_awh;
    fr_vect = scale_factor*wh_nf;

    fr_p = max(fr_vect,0);
    fr_m = -min(fr_vect,0);

    %generate Poisson spike train
    %positive part of the stimulus
    vm_vect = zeros(size(fr_vect));
    spk_p = zeros(size(fr_vect));

    n_gamma = gamman_v(g_ind);
    vm_thres = gamrnd(n_gamma,1/n_gamma);
    vm_vect(1) = 0;
    rho_vect = (fr_p/1000); %in spk/msec

    %save some time
    dtrho_vect = dt_ms*rho_vect;

    for i = 2:length(fr_p)
        vm_vect(i) = vm_vect(i-1) + dtrho_vect(i);
        if (vm_vect(i) > vm_thres)
            vm_vect(i) = 0;
            spk_p(i) = 1;
            vm_thres = gamrnd(n_gamma,1/n_gamma);
        end;
    end;

    spk_p = spk_p*fs; %in spk/sec

    %generate Poisson spike train
    %negative part of the stimulus
    vm_vect = zeros(size(fr_vect));
    spk_m = zeros(size(fr_vect));
    
    vm_thres = gamrnd(n_gamma,1/n_gamma);
    vm_vect(1) = 0;
    rho_vect = (fr_m/1000); %in spk/msec

    %save some time
    dtrho_vect = dt_ms*rho_vect;

    for i = 2:length(fr_m)
        vm_vect(i) = vm_vect(i-1) + dtrho_vect(i);
        if (vm_vect(i) > vm_thres)
            vm_vect(i) = 0;
            spk_m(i) = 1;
            vm_thres = gamrnd(n_gamma,1/n_gamma);
        end;
    end;

    spk_m = spk_m*fs; %in spk/sec
    
    %combine the two spike trains
    spk_c = spk_p - spk_m;

    %compute the reconstruction filter
    [gspkf_lgn,f_lgn] = tfestimate(spk_c,wh_n,nfft,nfft/2,nfft,fs,'twosided');

    %set to zero irrelevant components, i.e., above the cutoff frequency 
    %of the stimulus
    inds_z = find( (f_lgn > fc) & (f_lgn < fs - fc) );
    gspkf_lgn(inds_z) = 0;

    gspkt_lgn = ifft(gspkf_lgn);

    nlgn = length(gspkt_lgn);
   
    %make sure we have a 1 x nlgn columns vector
    gspkt_lgn = gspkt_lgn(:)';

    %unwrap the filter
    gspkt_lgn2 = circshift(gspkt_lgn,[0  (nlgn/2-1)]);
    t_lgn2 = ((-nlgn/2+1):nlgn/2)*dt_s;

    %compute the reconstruction
    rec_v = fftfilt(gspkt_lgn2,spk_c);

    %compensate for the delay due to the negative components 
    rec_v2 = rec_v(nlgn/2:nw-nlgn/2);
    wh_n2 = wh_n(1:nw-nfft+1);

    %compute the normalized error
    wn_noise = wh_n2 - rec_v2;
    err_n = sqrt(mean((wn_noise).^2))/std(wh_n);
    nerr_v(g_ind) = err_n;

    inds_p = find(spk_p > fs/2);
    isis_p = diff(inds_p)*dt_ms;
    inds_m = find(spk_m > fs/2);
    isis_m = diff(inds_m)*dt_ms;
    isis = [isis_p isis_m];
    cv_v(g_ind) = std(isis)/mean(isis);
end;

h_a5 = subplot(3,1,2);
line('Parent',h_a5,'XData',1:length(gamman_v),'YData',nerr_v,...
    'Marker','s','Color','k','LineStyle','--');
xlabel(h_a5,'gamma order');
set(h_a5,'XTick',1:length(gamman_v),'XTickLabel',gamman_v,'XLim',[1 length(gamman_v)]);
ylabel(h_a5,'normalized error');

h_a6 = subplot(3,1,3);
line('Parent',h_a6,'XData',1:length(gamman_v),'YData',cv_v,...
    'Marker','o','Color','r','LineStyle','--');
xlabel(h_a6,'gamma order');
set(h_a6,'XTick',1:length(gamman_v),'XTickLabel',gamman_v,'XLim',[1 length(gamman_v)]);
ylabel(h_a6,'coefficient of variation');

%print(handles.figure1,'-depsc2','rec_wn8.eps');

%% Start by plotting all the filters
% further clean up to get to the final figures for the paper
%
load('exp_filters.mat');


dgg_func3 = @(p,x) (p(1)*x.^p(2).*exp(-x.^p(4)/p(3)) + p(5)*x.^p(6).*exp(-x.^p(8)/p(7)))';

params_tm1 = [-2.4036e-1 1.24532e1 5.3642e-2 2.4366e-1 1.8255e-8 2.6814 2.036e2 9.9257e-1]; 
f1fit = dgg_func3(params_tm1,t);

h_f = figure; 
ha = subplot(1,2,1);
hp1d = plot(t,f_tm1,'Color',[0 0 0]);
hold on; 

hp1f = plot(t,f1fit,'Color',[1 0 0]);
legend([hp1d hp1f],{'Tm1 data', 'Tm1 fit'});
set(ha,'XLim',[0 2000],'TickDir','out');
xlabel('time (ms)');
ylabel('filter weights, per ms contrast'); 

params_tm2 = [-2.1473e3 15.4885 0.0364 0.2343 2.2721e-9 3.7628 31.9043 0.8539]; 
f2fit = dgg_func3(params_tm2,t);

ha = subplot(1,2,2);
hp2d = plot(t,f_tm2,'k');
hold on; 

hp2f = plot(t,f2fit,'r');
legend([hp2d hp2f],{'Tm2 data', 'Tm2 fit'});
set(ha,'XLim',[0 2000],'TickDir','out');
%print(h_f,'figures/tempfilters_off.eps','-depsc');

%%

params_tm3 = [2.2757e-8 7.8192 0.404 0.4516 -1.0482e-9 3.2796 179.2137 1.0173]; 
f3fit = dgg_func3(params_tm3,t);

h_f = figure; 
ha = subplot(1,2,1);
hp3d = plot(t,f_tm3,'k');
hold on; 

hp3f = plot(t,f3fit,'r');
legend([hp3d hp3f],{'Tm3 data', 'Tm3 fit'});
set(ha,'XLim',[0 2000],'TickDir','out');
xlabel('time (ms)');
ylabel('contrast weight per ms'); 

params_mi1 = [7.7293e-1 1.4552e1 4.0010e-2 2.1850e-1 -8.2252e-12 3.8217 291.4865 1.0578];
f4fit = dgg_func3(params_mi1,t);

ha = subplot(1,2,2);
hp4d = plot(t,f_mi1,'k');
hold on; 
hp4f = plot(t,f4fit,'r');
legend([hp4d hp4f],{'Mi1 data', 'Mi1 fit'});
set(ha,'XLim',[0 2000],'TickDir','out');
%print(h_f,'figures/tempfilters_on.eps','-depsc');

%% Plot pairs of filters together

% OFF pathway
h_f = figure; 
ha = axes;
hp1d = plot(t,f_tm1,'Color',[0 0 0]);
hold on; 
hp1f = plot(t,f1fit,'r');
hp2d = plot(t,f_tm2,'Color',[0.5 0.5 0.5]);
hold on; 
hp2f = plot(t,f2fit,'Color',[0.5 0 0]);
set(ha,'XLim',[0 500],'TickDir','out');
legend([hp1f hp2f],{'Tm1 fit', 'Tm2 fit'});
%print(h_f,'figures/compfilters_off.eps','-depsc');

%difference between peaks in experimental data
[tm1m tm1i] = min(f_tm1);
[tm2m tm2i] = min(f_tm2);
diff_exp_off = abs(t(tm1i)-t(tm2i));

%repeat for the fits
[ff1m ff1i] = min(f1fit);
[ff2m ff2i] = min(f2fit);
diff_fit_off = abs(t(ff1i)-t(ff2i));

disp(['experimental peak difference ' num2str(diff_exp_off) ' ms; fit peak difference: ' num2str(diff_fit_off)]);

%%
% ON pathway
h_f = figure; 
ha = axes;
hp3d = plot(t,f_tm3,'Color',[0 0 0]);
hold on; 
hp3f = plot(t,f3fit,'r');
hp4d = plot(t,f_mi1,'Color',[0.5 0.5 0.5]);
hold on; 
hp4f = plot(t,f4fit,'Color',[0.5 0 0]);
set(ha,'XLim',[0 500],'TickDir','out');
legend([hp3f hp4f],{'Tm3 fit', 'Mi1 fit'});
%print(h_f,'figures/compfilters_on.eps','-depsc');

%difference between peaks in experimental data
[tm3m tm3i] = max(f_tm3);
[mi1m mi1i] = max(f_mi1);
diff_exp_on = abs(t(tm3i)-t(mi1i));

%repeat for the fits
[ff3m ff3i] = max(f3fit);
[ff4m ff4i] = max(f4fit);
diff_fit_on = abs(t(ff3i)-t(ff4i));

disp(['experimental peak difference ' num2str(diff_exp_on) ' ms; fit peak difference: ' num2str(diff_fit_on)]);

%% fits of static non-linarities from nonlinear_fit.m
%

%typical input range
x = -30:0.1:30;

%average fit
p12 = [-10.46 38.20 7.773 6.8075];
nlf12_y = sig_func(p12,x);

h_f = figure; 
hp1 = plot(nl_tm1_x,nl_tm1_y,'k');
hold on;
hp2 = plot(nl_tm2_x,nl_tm2_y,'Color',[0.5 0.5 0.5]);
hp12 = plot(x,nlf12_y,'r');
legend([hp1 hp2 hp12],{'Tm1', 'Tm2', 'fit'});
set(gca,'XLim',[-30 30],'TickDir','out');
xlabel('linear filter output');
ylabel('non-linearity output');
%print(h_f,'figures/statnl_off.eps','-depsc');

%average fit
p34 = [-9.9460 26.505 3.1410 5.7985];
nlf34_y = sig_func(p34,x);

h_f = figure; 
hp3 = plot(nl_tm3_x,nl_tm3_y,'k');
hold on;
hp4 = plot(nl_mi1_x,nl_mi1_y,'Color',[0.5 0.5 0.5]);
hp34 = plot(x,nlf34_y,'r');
legend([hp3 hp4 hp34],{'Tm3', 'Mi1','fit'});
set(gca,'XLim',[-30 30],'TickDir','out');
xlabel('linear filter output');
ylabel('non-linearity output');
%print(h_f,'figures/statnl_on.eps','-depsc');

%% Compute tuning to sinewave stimuli
debug_flag = 0;

dt = 1e-3; %sampling step in s (1 ms)
t_sim = [0:dt:20-dt]; %10 s
t_pre = [-5:dt:-dt]; %5 s
n_tot = length(t_pre) + length(t_sim);

%spatial wavelength of sinusoidal stimulus
lambda_x = 30; % deg

%spacing between 2 sampling stations of the correlator
dphi = 5.1; %deg

%omega_t = [0.1 0.5 1 5 10]; c/s (debug)
omega_t = 10.^[-1:.05:1]; %c/s
amp = 0.11; 
for l = 1:length(omega_t)
    %stimulus at first facet
    stim_a = amp*sin(2*pi*omega_t(l) * t_sim);
    
    %stimulus at second facet, phase shifted by the ratio of the distance
    %between the facets and wavelength of the sinusoidal stimulus
    stim_b = amp*sin(2*pi*omega_t(l) * t_sim + 2*pi*dphi/lambda_x);
    
    %make sure we are at steady-state, compute response over 5 s before motion onset
    stim_a = [stim_a(1)*ones(1,5000),stim_a];
    stim_b = [stim_b(1)*ones(1,5000),stim_b];
    
    %output of Tm1 at position a
    h1a = filter(f1fit,1,stim_a);
    h1an = sig_func(p12,h1a);
    
    %output of Tm2 at position b
    h2b = filter(f2fit,1,stim_b);
    h2bn = sig_func(p12,h2b);

    %output of Tm1 at position b
    h1b = filter(f1fit,1,stim_b);
    h1bn = sig_func(p12,h1b);

    %output of Tm2 at position a
    h2a = filter(f2fit,1,stim_a);
    h2an = sig_func(p12,h2a);

    %compute full correlation model response
    hrc12n = (h1bn.*h2an - h1an.*h2bn);
    
    %average over the last 10 s
    mhrc12n(l) = mean(hrc12n(n_tot-10000:n_tot));
    
    if ( debug_flag == 1 )
        figure;
        hp2 = plot([t_pre t_sim],hrc12n);
        legend([hp2],{'non-linear HRC'});
        title([ 'Response at a temporal frequency of ' num2str(omega_t(l)) ' c/s']);
    end
    
    %output of Tm3 at position a
    h3a = filter(f3fit,1,stim_a);
    h3an = sig_func(p34,h3a);

    %output of Mi1 at position b
    h4b = filter(f4fit,1,stim_b);
    h4bn = sig_func(p34,h4b);

    %output of Tm3 at position b
    h3b = filter(f3fit,1,stim_b);
    h3bn = sig_func(p34,h3b);

    %output of Mi1 at position a
    h4a = filter(f4fit,1,stim_a);
    h4an = sig_func(p34,h4a);

    %compute full correlation model response
    hrc34n = -(h3bn.*h4an - h3an.*h4bn);

    %average over the last 10 s
    mhrc34n(l) = mean(hrc34n(n_tot-10000:n_tot));
end

mhrc12n_norm = mhrc12n/max(mhrc12n);
mhrc34n_norm = mhrc34n/max(mhrc34n);

%%
load('mHRC.mat');

h_f = figure;
ha = subplot(1,2,1);
hp4 = plot(omega_t,mhrc34n_norm,'r');
hold on;
hp4a = plot(omega,mHRC34_norm,'k');
set(ha,'xscale','log','YLim',[0 1],'TickDir','out');
%legend([hp4 hp4a],{'LN model','experimental'});
title('L1 (ON) pathway');
xlabel('temporal frequency (Hz)');
ylabel('normalized mean response amplitude');

ha = subplot(1,2,2);
hp2 = plot(omega_t,mhrc12n_norm,'r');
hold on;
hp2a = plot(omega,mHRC12_norm,'k');
set(ha,'xscale','log','YLim',[0 1],'TickDir','out');
legend([hp2 hp2a],{'LN model','experimental'});
title('L2 (OFF) pathway');
xlabel('temporal frequency (Hz)');
%print(h_f,'figures/drosofreq_tune.eps','-depsc');

%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

%%

fid = fopen('timeseries/ts001.bin','r');

%big endian format, ieee-be, 16 bits integer, unsigned, uint16
ts = fread(fid,inf,'uint16','ieee-be');

fclose(fid);

%in samples per sec
f_samp = 1200;

%in sec
dt = 1/f_samp;

%45 mins
n_ts = length(ts);

s_v = (1:n_ts);
tv = s_v*dt;

h_f1 = figure; 
h_a1 = subplot(3,1,1);
line('Parent',h_a1,'XData',tv,'YData',ts(s_v));
set(h_a1,'XLim',[0 2700]);

%45 sec
inds = find ( (tv > 960) & (tv <=1005) );
h_a2 = subplot(3,1,2);
line('Parent',h_a2,'XData',tv(inds),'YData',ts(inds));
ylabel(h_a2,'Luminance (arbitrary units)');

%4.5 sec
inds = find ( (tv > 960) & (tv <=964.5) );
h_a3 = subplot(3,1,3);
line('Parent',h_a3,'XData',tv(inds),'YData',ts(inds));
xlabel(h_a3,'time (s)');

%%
[n xout] = hist(ts,100);
h_f2 = figure;
h_a4 = subplot(2,2,1);
bar(h_a4,xout,n);
set(h_a4,'XLim',[-0.1e4 4e4],'TickDir','out');
xlabel(h_a4,'Luminance');
ylabel(h_a4,'Observations');

[nl xoutl] = hist(log10(ts),100);
h_a5 = subplot(2,2,2);
bar(h_a5,xoutl,nl);
set(h_a5,'TickDir','out','XLim',[1 5]);
xlabel(h_a5,'log_1_0 (Luminance)');

%the sampling rate is 1200 samples/s = 2*fNyquist
%we have 45 mins of data or 3,240,000 samples
%split in 4.5min data sets or 324,000 samples
%carry out fft on 4096 samples yielding a frequency resolution
%of 1200/4096 = 0.29 Hz. Since 324,000/4096 = 79.1, we use
%79*4096 = 323,584 samples with 416 (=324,000 - 323,584) unused
%or approx. a fraction of 416/324,000 = 0.0013 unused.

%centered at zero 
ts = ts - mean(ts);

l_ts = length(ts);
l_set = l_ts/10;

n_fft = 4096;
n_samp = floor(l_set/4096)*4096;

%positive part of the power spectrum
p_m = zeros(10,2049);

for i = 1:10
    
    %the parameters are:
    %1) the length of the window n_fft = 4096
    %overlap: n_fft/2 = 2048
    %fft size: n_ftt = 4096
    %maximal two-sided frequency 1200 samples/s
    samples = (i-1)*l_set + (1:n_samp);
    [p_v,w_v] = pwelch(ts(samples),n_fft,n_fft/2,n_fft,1200,'twosided');

    p_m(i,:) = p_v(1:2049);
    
    %take only positive frequencies
    flog_v = log10(w_v(2:2049));
    plog_v = log10(p_v(2:2049));

    %use the same df and xx vector for all the data
    %so compute the following only in the first pass
    if ( i == 1 )
        %linear fit to the data to avoid bias to high 
        %frequencies, we need to resample
        df = (flog_v(end) - flog_v(1))/100;
        xx = flog_v(1):df:flog_v(end);
    end;
    
    %resample log10 power spectrum at xx values
    yy = spline(flog_v,plog_v,xx);
    
    %linear fit
    %i_xx = [xx; ones(1,length(xx))]';
    %l_y = i_xx\yy';
    
    %uncomment to plot individual power spectra and fits
    %figure;
    %plot(flog_v,plog_v);
    %hold on;
    %plot(xx,i_xx*l_y,'g');
    
    if (i == 1)
        yy_m = zeros(10,length(yy));
    end;
    
    yy_m(i,:) = yy;
    
end;

p_mean = mean(p_m,1);
h_a6 = subplot(2,2,3);
line('Parent',h_a6,'XData',w_v(1:2049),'YData',p_mean,'Color','k');
set(h_a6,'XLim',[0 5],'YLim',[0 8e6]);
xlabel(h_a6,'Frequency (Hz)');
ylabel(h_a6,'Luminance^2/Hz');

yy_mean = mean(yy_m,1);
yy_std = std(yy_m,0,1);


%linear fit
i_xx = [xx; ones(1,length(xx))]';
l_y = i_xx\yy_mean';

%plot data
h_a7 = subplot(2,2,4);
line('Parent',h_a7,'XData',xx,'YData',yy_mean,'Color','k');
line('Parent',h_a7,'XData',xx,'YData',yy_mean+yy_std,'Color','k','LineStyle',':');
line('Parent',h_a7,'XData',xx,'YData',yy_mean-yy_std,'Color','k','LineStyle',':');

%plot fit
line('Parent',h_a7,'XData',xx,'YData',i_xx*l_y,'Color','r');
set(h_a7,'XLim',[-1 3],'YLim',[0 7]);
xlabel(h_a7,'log_1_0(Frequency)');
ylabel(h_a7,'log_1_0(Luminance^2/Hz)');

%print(handles.figure1,'-depsc2','disp_ts2.eps');

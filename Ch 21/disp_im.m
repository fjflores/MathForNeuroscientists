%% Reads in the data
%computes average luminance histograms of ten images from the 
%van Hateren data base in linear and log units 
i_max = 10;
w_im = 1536;
h_im = 1024;
im = zeros(i_max,w_im,h_im);

for i = 1:i_max
    if ( i ~= 10 )
        im_str = ['images/imk0000' num2str(i) '.iml'];
    else
        im_str = ['images/imk000' num2str(i) '.iml'];
    end;
    
    f1=fopen(im_str,'rb','ieee-be');
    
    im(i,:,:)=fread(f1,[w_im,h_im],'uint16');
    fclose(f1);
    
    

end;

%plot last image
figure; 
colormap(gray);
imagesc(squeeze(im(i,:,:))');

%% Power spectrum
%computes the power spectrum of ten images from the 
%van Hateren data base along the horizontal and vertical
%dimension. Averages the corresponding 20 power spectra
%and fits them to a power law in log-log space

%power spectra in linear units, 10 lines with horizontal spectra, 10 lines
%with vertical spectra (total 20 lines x 128 power spectrum points)
p_lin_h = zeros(i_max,128);
p_lin_v = zeros(i_max,128);

%power spectra in log10 units
p_log_h = zeros(i_max,128);
p_log_v = zeros(i_max,128);

for i = 1:i_max
    
    disp(['processing image ' num2str(i) ' ...']);
    %each image is 1536x1024 (width x height)
    %compute power spectrum along the y (vertical) axis

    %get rid of 3rd dimension
    a = squeeze(im(i,257:1536-256,:)); %lines are horizontal dimension and columns vertical dimension

    %reshape takes elements columnwise, so we first transpose to get each
    %vertical image stripe as a column
    b_v = reshape(a',prod(size(a)),1);
    c_v = b_v - mean(b_v);
    
    %the parameters are: 1) the length of the window 256
    %2) overlap: 128
    %3) fft size 256
    %maximal two-sided frequency 60 (each pixel is 1 min of arc)
    [p_v,w_v] = pwelch(c_v,256,128,256,60,'twosided');
    
    %take only positive frequencies
    flog_v = log10(w_v(2:129));
    plog_v = log10(p_v(2:129));

    %linear fit to the data
    %flog_ve = [flog_v'; ones(1,128)]';
    %x = flog_ve\plog_v;
    
    b_h = reshape(a,prod(size(a)),1);
    c_h = b_h - mean(b_h);
    [p_h,w_h] = pwelch(c_h,256,128,256,60,'twosided');
    flog_h = log10(w_h(2:129));
    plog_h = log10(p_h(2:129));

    %linear fit to the data
    %flog_he = [flog_h'; ones(1,128)]';
    %y = flog_he\plog_h;
    
    %plot individual power spectra and fits
    %figure;
    %plot(flog_v,plog_v);
    %hold on;
    %plot(flog_h,plog_h,'r');
    %plot(flog_v,flog_ve*x,'g');
    %plot(flog_h,flog_he*y,'y');
    
    %save data one line for the vertical and one line for the horizontal 
    %spectrum
    p_log_h(i,:) = plog_h';
    p_log_v(i,:) = plog_v';
    p_lin_h(i,:) = p_h(2:129)';
    p_lin_v(i,:) = p_v(2:129)';
    
end;

%% Plot results
 
p_log_mean = mean([p_log_h; p_log_v],1);
p_log_std = std([p_log_h; p_log_v],0,1);

p_lin_mean = mean([p_lin_h; p_lin_v],1);

h_f1 = figure; 
h_a1 = axes;
plot(flog_v,p_log_mean,'k');
hold on; 
plot(flog_v,p_log_mean-p_log_std,'k:');
plot(flog_v,p_log_mean+p_log_std,'k:');

%linear fit
p_log_meant = p_log_mean';
flog_ve = [flog_v'; ones(1,128)]';
x = flog_ve\p_log_meant;
plot(flog_v,flog_ve*x,'r');

set(h_a1,'TickDir','out');
xlabel('log_{10}(Frequency)');
ylabel('log_{10}(Luminance^2/Hz)');

h_f2 = figure; 
h_a2 = axes;
plot(w_v(2:129),p_lin_mean,'k');
set(h_a2,'TickDir','out');
xlabel('Frequency (Hz)');
ylabel('Luminance power spectral density (cd/m^2)^2/Hz)');

%% Print figures
%print(h_f1,'figures_ed2/nat_scene_psdlog.eps','-depsc');
%print(h_f2,'figures_ed2/nat_scene_psdlin10.eps','-depsc');

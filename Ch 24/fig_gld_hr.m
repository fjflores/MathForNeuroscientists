%% Plot glider stimuli and responses
%  generate glider stimuli with two checks of odd and even parity and 
%  compute standard correlation model response. Then generate diverging
%  3-point glider of odd and even parity and compute responses.

%position of each glider check with coordinates (x, y, t) in that order
glider_hr_lr = [0 0 0;1 0 1]; %from (0,0) to (1,1): translate from left to right
glider_d3 = [0 0 0;0 0 1;1 0 1]; %corresponds to |_
glider_c3 = [0 0 0;1 0 0;1 0 1]; %corresponds to -|
parity_e = 0; %even parity
parity_o = 1; %odd parity

x_max = 90;
dx = 2;
x_vect = dx:dx:x_max;
n_x = length(x_vect);

%set to 1 to see the output of a spatial array of correlation detectors
spatial_averaging = 0;

%time interval
dt = 25e-3; %s

%t_start = 0;
t_start = dt;
%t_end = 101*dt; %debug
%t_end = 1; %debug
t_end = 5; %s
t_vect = t_start:dt:t_end;
n_t = length(t_vect);

%For HR responses, resample at a resolution of 0.1 ms. Since each point is sampled at 25 ms
%resolution, we need 250 points
rs_rate = 250;
dt_rs = dt/rs_rate;
n_trs = n_t*rs_rate;

%HR low-pass filter time constants
tau_s = 20e-3; %20 ms
tau_f = 2e-3; %2 ms

%set to 1 for additional debug plots
plot_debug = 0; 

%% Compute a left-right 2-point glider with even parity
%returns a n_x x 1 x n_t logical array
s = testgen(n_x,1,n_t,glider_hr_lr,parity_e);

%get rid of the y dimension
s = squeeze(s);

%convert into contrast pattern
s2 = zeros(size(s));
s2(s==0) = -1;
s2(s==1) = 1;

%convert to n_t x n_x array, which is easier to show
stim = s2';

%clear high memory load data
clear s s2; 

h_f1 = figure; 
ha = axes;
imagesc(stim);
colormap(gray);
set(ha,'TickDir','out');

%these are the indices in the stim matrix
tick_x = get(ha,'XTick');
pos_x = x_vect(tick_x); %convert to position
%convert to cell array of strings, one per tick
label_x = cellstr(num2str(pos_x'));
set(ha,'XTickLabel',label_x);
xlabel('position (deg)');

tick_y = get(ha,'YTick');
pos_y = t_vect(tick_y); %convert to position
%convert to cell array of strings, one per tick
label_y = cellstr(num2str(pos_y'));
set(ha,'YTickLabel',label_y);
%get(ha,'XLim')
set(ha,'XLim',[1 25]);
ylabel('time (s)');
%colorbar(ha);
%title('Left-right even parity 2 point glider');

%print(h_f1,'figures/gld_stims1.eps','-depsc','-r300','-painters');
%print(h_f1,'figures/gld_stims1.tif','-dtiff','-r300','-painters');

%% Compute response to a left-right 2-point glider with even parity

%extract luminance patterns
lum_vect0 = stim(:,1); % at x = 0

%this first add rs_rate columns to create a n_t x rs_rate matrix. Then the
%matrix is transposed because reshape takes elements in one column before 
%proceeding to the next. Final result is a column vector
lv0_rs = reshape(repmat(lum_vect0,1,rs_rate)',[n_trs 1]);
t_vect_rs = (1:n_trs)*dt_rs;

if ( plot_debug == 1 )
    figure;
    hp1 = plot(t_vect,lum_vect0);
    hold on;
    hp2 = plot(t_vect_rs,lv0_rs,'r');
    xlabel('time (s)');
    legend([hp1 hp2],{'original','resampled'});
end

lum_vect1 = stim(:,2); % at x = 2
lv1_rs = reshape(repmat(lum_vect1,1,rs_rate)',[n_trs 1]);

%compute responses of the 1st order low-pass filters
lv0_s_num = lp1filt_fn(lv0_rs,dt_rs,tau_s);
lv0_f_num = lp1filt_fn(lv0_rs,dt_rs,tau_f);
lv1_s_num = lp1filt_fn(lv1_rs,dt_rs,tau_s);
lv1_f_num = lp1filt_fn(lv1_rs,dt_rs,tau_f);

%half correlator responses
r_hc1 = lv0_s_num.*lv1_f_num;
r_hc2 = lv1_s_num.*lv0_f_num;

%full opponent response
r_hr_single = r_hc1 - r_hc2;

%time average
av_r_hr_single = mean(r_hr_single);

%plot the correlator response
hf = figure; 
h_p1 = plot(t_vect_rs,r_hr_single,'k');
hold on;
h_p1a = plot([t_vect_rs(1) t_vect_rs(end)],av_r_hr_single*ones(1,2),'r');
xlabel('time (s)');
ylabel('correlation model response');
%legend([h_p1 h_p1a],{'correlation model','time-average'});

%title('single correlator response');

if ( spatial_averaging == 1 )
    %compute the correlator response averaged over space
    r_hr_mean = zeros(size(t_vect_rs));
    
    for j=1:n_x-1 %n_corr-1
        lum_vect0 = stim(:,j); %position1
        lum_vect1 = stim(:,j+1); % at position2 = position1 + dx
        
        %resample
        lv0_rs = reshape(repmat(lum_vect0,1,rs_rate)',[n_trs 1]);
        lv1_rs = reshape(repmat(lum_vect1,1,rs_rate)',[n_trs 1]);
        
        %compute responses of the 1st order low-pass filters
        lv0_s_num = lp1filt_fn(lv0_rs,dt_rs,tau_s);
        lv0_f_num = lp1filt_fn(lv0_rs,dt_rs,tau_f);
        lv1_s_num = lp1filt_fn(lv1_rs,dt_rs,tau_s);
        lv1_f_num = lp1filt_fn(lv1_rs,dt_rs,tau_f);
        
        %half correlator responses
        r_hc1 = lv0_s_num.*lv1_f_num;
        r_hc2 = lv1_s_num.*lv0_f_num;
        
        %full opponent response
        r_hr = r_hc1 - r_hc2;
        
        %add to the mean vector
        r_hr_mean = r_hr_mean + r_hr';
    end
    
    r_hr_mean = r_hr_mean/(n_x-1);
    
    %compute the mean response over the last 4 seconds of the response
    m_resp_1 = mean(r_hr_mean(t_vect_rs>1));
    
    %plot the averaged correlator response
    figure;
    h_p1 = plot(t_vect_rs,r_hr_single,'k');
    hold on;
    h_p2 = plot(t_vect_rs,r_hr_mean,'r');
    plot([t_vect_rs(1) t_vect_rs(end)],[m_resp_1 m_resp_1],'r--')
    xlabel('time (s)');
    ylabel('response');
    legend([h_p1 h_p2], {'single correlator','spatial average'});
    title('spatially averaged vs. single correlator response');
end

%% Create a left-right glider with odd parity
%returns a n_x x 1 x n_t logical array
s = testgen(n_x,1,n_t,glider_hr_lr,parity_o);

%get rid of the y dimension
s = squeeze(s);

%convert into contrast pattern
s2 = zeros(size(s));
s2(s==0) = -1;
s2(s==1) = 1;

%convert to n_t x n_x array, which is easier to show
stim = s2';

%clear high memory load data
clear s s2; 

h_f2 = figure; 
ha = axes;
imagesc(stim);
colormap(gray);
set(ha,'TickDir','out');
set(ha,'XLim',[1 25]);

%these are the indices in the stim matrix
tick_x = get(ha,'XTick');
pos_x = x_vect(tick_x); %convert to position
%convert to cell array of strings, one per tick
label_x = cellstr(num2str(pos_x'));
set(ha,'XTickLabel',label_x);
xlabel('position (deg)');

tick_y = get(ha,'YTick');
pos_y = t_vect(tick_y); %convert to position
%convert to cell array of strings, one per tick
label_y = cellstr(num2str(pos_y'));
set(ha,'YTickLabel',label_y);
ylabel('time (s)');
%colorbar(ha);
%title('Left-right odd parity 2 point glider')
%print(h_f2,'figures/gld_stims2.eps','-depsc','-r300','-painters');
%print(h_f2,'figures/gld_stims2.tif','-dtiff','-r300','-painters');

%% Response to a left-right glider with odd parity

%extract luminance patterns
lum_vect0 = stim(:,1); % at x = 0

%this first add rs_rate columns to create a n_t x rs_rate matrix. Then the
%matrix is transposed because reshape takes elements in one column before 
%proceeding to the next. Final result is a column vector
lv0_rs = reshape(repmat(lum_vect0,1,rs_rate)',[n_trs 1]);

lum_vect1 = stim(:,2); % at x = 2
lv1_rs = reshape(repmat(lum_vect1,1,rs_rate)',[n_trs 1]);

%compute responses of the 1st order low-pass filters
lv0_s_num = lp1filt_fn(lv0_rs,dt_rs,tau_s);
lv0_f_num = lp1filt_fn(lv0_rs,dt_rs,tau_f);
lv1_s_num = lp1filt_fn(lv1_rs,dt_rs,tau_s);
lv1_f_num = lp1filt_fn(lv1_rs,dt_rs,tau_f);

%half correlator responses
r_hc1 = lv0_s_num.*lv1_f_num;
r_hc2 = lv1_s_num.*lv0_f_num;

%full opponent response
r_hr_single = r_hc1 - r_hc2;

%time average
av_r_hr_single = mean(r_hr_single);

%plot the correlator response
%hf = figure;
figure(hf);
h_p2 = plot(t_vect_rs,r_hr_single,'k:');
hold on;
h_p2a = plot([t_vect_rs(1) t_vect_rs(end)],av_r_hr_single*ones(1,2),'r--');
xlabel('time (s)');
ylabel('correlation model response');
%legend([h_p2 h_p2a],{'correlation model','time-average'});
legend([h_p1 h_p1a h_p2 h_p2a],{'correlation model','time-average','odd parity','time-average'});
%ylims = get(gca,'YLim');
%ylims = [-ylims(2) ylims(2)];
%set(gca,'YLim',ylims);

%print(hf,'figures/gld_stims3.eps','-depsc','-r300','-painters');

%title('single correlator response');

if ( spatial_averaging == 1 )
    %compute the correlator response averaged over space
    r_hr_mean3 = zeros(size(t_vect_rs));
    
    for j=1:n_x-1 %n_corr-1
        lum_vect0 = stim(:,j); %position1
        lum_vect1 = stim(:,j+1); % at position2 = position1 + dx
        
        %resample
        lv0_rs = reshape(repmat(lum_vect0,1,rs_rate)',[n_trs 1]);
        lv1_rs = reshape(repmat(lum_vect1,1,rs_rate)',[n_trs 1]);
        
        %compute responses of the 1st order low-pass filters
        lv0_s_num = lp1filt_fn(lv0_rs,dt_rs,tau_s);
        lv0_f_num = lp1filt_fn(lv0_rs,dt_rs,tau_f);
        lv1_s_num = lp1filt_fn(lv1_rs,dt_rs,tau_s);
        lv1_f_num = lp1filt_fn(lv1_rs,dt_rs,tau_f);
        
        %half correlator responses
        r_hc1 = lv0_s_num.*lv1_f_num;
        r_hc2 = lv1_s_num.*lv0_f_num;
        
        %full opponent response
        r_hr = r_hc1 - r_hc2;
        
        %add to the mean vector
        r_hr_mean3 = r_hr_mean3 + r_hr';
    end
    
    r_hr_mean3 = r_hr_mean3/(n_x-1);
    
    %compute the mean response over the last 4 seconds of the response
    m_resp_3 = mean(r_hr_mean3(t_vect_rs>1));
    
    %plot the averaged correlator responses to even and odd gliders
    figure;
    h_p3 = plot(t_vect_rs,r_hr_mean);
    hold on;
    h_p4 = plot(t_vect_rs,r_hr_mean3,'r');
    plot([t_vect_rs(1) t_vect_rs(end)],[m_resp_1 m_resp_1],'b--')
    plot([t_vect_rs(1) t_vect_rs(end)],[m_resp_3 m_resp_3],'r--')
    xlabel('time (s)');
    ylabel('response');
    legend([h_p3 h_p4],{'even parity','odd parity'});
    title('spatially averaged response to left-right gliders even and odd parity');
end

%% Create a 3-point diverging glider even parity
%returns a n_x x 1 x n_t logical array
s = testgen(n_x,1,n_t,glider_d3,parity_e);

%get rid of the y dimension
s = squeeze(s);

%convert into contrast pattern
s2 = zeros(size(s));
s2(s==0) = -1;
s2(s==1) = 1;

%convert to n_t x n_x array, which is easier to show
stim = s2';

%clear high memory load data
clear s s2; 

h_f4 = figure; 
ha = axes;
imagesc(stim);
colormap(gray);
set(ha,'TickDir','out');
set(ha,'XLim',[1 25]);

%these are the indices in the stim matrix
tick_x = get(ha,'XTick');
pos_x = x_vect(tick_x); %convert to position
%convert to cell array of strings, one per tick
label_x = cellstr(num2str(pos_x'));
set(ha,'XTickLabel',label_x);
xlabel('position (deg)');

tick_y = get(ha,'YTick');
pos_y = t_vect(tick_y); %convert to position
%convert to cell array of strings, one per tick
label_y = cellstr(num2str(pos_y'));
set(ha,'YTickLabel',label_y);
ylabel('time (s)');
%colorbar(ha);
%title('Even parity diverging 3 point glider')
%print(h_f4,'figures/gld_stims4.eps','-depsc','-r300','-painters');
%print(h_f4,'figures/gld_stims4.tif','-dtiff','-r300','-painters');

%% Response to a 3-point diverging glider even parity
%extract luminance patterns
lum_vect0 = stim(:,1); % at x = 0

%this first add rs_rate columns to create a n_t x rs_rate matrix. Then the
%matrix is transposed because reshape takes elements in one column before 
%proceeding to the next. Final result is a column vector
lv0_rs = reshape(repmat(lum_vect0,1,rs_rate)',[n_trs 1]);

lum_vect1 = stim(:,2); % at x = 2
lv1_rs = reshape(repmat(lum_vect1,1,rs_rate)',[n_trs 1]);

%compute responses of the 1st order low-pass filters
lv0_s_num = lp1filt_fn(lv0_rs,dt_rs,tau_s);
lv0_f_num = lp1filt_fn(lv0_rs,dt_rs,tau_f);
lv1_s_num = lp1filt_fn(lv1_rs,dt_rs,tau_s);
lv1_f_num = lp1filt_fn(lv1_rs,dt_rs,tau_f);

%half correlator responses
r_hc1 = lv0_s_num.*lv1_f_num;
r_hc2 = lv1_s_num.*lv0_f_num;

%full opponent response
r_hr_single = r_hc1 - r_hc2;

%time average
av_r_hr_single = mean(r_hr_single);

%plot the correlator response
h_f5 = figure;
%figure(hf);
h_p2 = plot(t_vect_rs,r_hr_single,'k');
hold on;
h_p2a = plot([t_vect_rs(1) t_vect_rs(end)],av_r_hr_single*ones(1,2),'r--');
xlabel('time (s)');
ylabel('correlation model response');
legend([h_p2 h_p2a],{'correlation model','time-average'});
%legend([h_p1 h_p1a h_p2 h_p2a],{'correlation model','time-average','odd parity','time-average'});

%print(h_f5,'figures/gld_stims5.eps','-depsc','-r300','-painters');

if ( spatial_averaging == 1 )
    %compute the correlator response averaged over space
    r_hr_mean4 = zeros(size(t_vect_rs));
    
    
    for j=1:n_x-1 %n_corr-1
        lum_vect0 = stim(:,j); %position1
        lum_vect1 = stim(:,j+1); % at position2 = position1 + dx
        
        %resample
        lv0_rs = reshape(repmat(lum_vect0,1,rs_rate)',[n_trs 1]);
        lv1_rs = reshape(repmat(lum_vect1,1,rs_rate)',[n_trs 1]);
        
        %compute responses of the 1st order low-pass filters
        lv0_s_num = lp1filt_fn(lv0_rs,dt_rs,tau_s);
        lv0_f_num = lp1filt_fn(lv0_rs,dt_rs,tau_f);
        lv1_s_num = lp1filt_fn(lv1_rs,dt_rs,tau_s);
        lv1_f_num = lp1filt_fn(lv1_rs,dt_rs,tau_f);
        
        %half correlator responses
        r_hc1 = lv0_s_num.*lv1_f_num;
        r_hc2 = lv1_s_num.*lv0_f_num;
        
        %full opponent response
        r_hr = r_hc1 - r_hc2;
        
        %add to the mean vector
        r_hr_mean4 = r_hr_mean4 + r_hr';
    end
    
    r_hr_mean4 = r_hr_mean4/(n_x-1);
    
    %compute the mean response over the last 4 seconds of the response
    m_resp_4 = mean(r_hr_mean4(t_vect_rs>1));
    
    %plot the averaged correlator response
    figure;
    h_p5 = plot(t_vect_rs,r_hr_mean);
    hold on;
    h_p6 = plot(t_vect_rs,r_hr_mean4,'r');
    plot([t_vect_rs(1) t_vect_rs(end)],[m_resp_1 m_resp_1],'b--')
    plot([t_vect_rs(1) t_vect_rs(end)],[m_resp_4 m_resp_4],'r--')
    xlabel('time (s)');
    ylabel('response');
    legend([h_p5 h_p6],{'even spatial average','diverging-even'});
    title('spatially averaged response left-right even and diverging-even parity gliders');
end

%% Create a 3-point diverging glider odd parity
%returns a n_x x 1 x n_t logical array
s = testgen(n_x,1,n_t,glider_d3,parity_o);

%get rid of the y dimension
s = squeeze(s);

%convert into contrast pattern
s2 = zeros(size(s));
s2(s==0) = -1;
s2(s==1) = 1;

%convert to n_t x n_x array, which is easier to show
stim = s2';

%clear high memory load data
clear s s2; 

figure; 
ha = axes;
imagesc(stim);
colormap(gray);
set(ha,'TickDir','out');
set(ha,'XLim',[1 25]);

%these are the indices in the stim matrix
tick_x = get(ha,'XTick');
pos_x = x_vect(tick_x); %convert to position
%convert to cell array of strings, one per tick
label_x = cellstr(num2str(pos_x'));
set(ha,'XTickLabel',label_x);
xlabel('position (deg)');

tick_y = get(ha,'YTick');
pos_y = t_vect(tick_y); %convert to position
%convert to cell array of strings, one per tick
label_y = cellstr(num2str(pos_y'));
set(ha,'YTickLabel',label_y);
ylabel('time (s)');
colorbar(ha);
title('Odd parity diverging 3 point glider')

%% Response to a 3-point diverging glider odd parity
%extract luminance patterns
lum_vect0 = stim(:,1); % at x = 0

%this first add rs_rate columns to create a n_t x rs_rate matrix. Then the
%matrix is transposed because reshape takes elements in one column before 
%proceeding to the next. Final result is a column vector
lv0_rs = reshape(repmat(lum_vect0,1,rs_rate)',[n_trs 1]);

lum_vect1 = stim(:,2); % at x = 2
lv1_rs = reshape(repmat(lum_vect1,1,rs_rate)',[n_trs 1]);

%compute responses of the 1st order low-pass filters
lv0_s_num = lp1filt_fn(lv0_rs,dt_rs,tau_s);
lv0_f_num = lp1filt_fn(lv0_rs,dt_rs,tau_f);
lv1_s_num = lp1filt_fn(lv1_rs,dt_rs,tau_s);
lv1_f_num = lp1filt_fn(lv1_rs,dt_rs,tau_f);

%half correlator responses
r_hc1 = lv0_s_num.*lv1_f_num;
r_hc2 = lv1_s_num.*lv0_f_num;

%full opponent response
r_hr_single = r_hc1 - r_hc2;

%time average
av_r_hr_single = mean(r_hr_single);

%plot the correlator response
hf = figure;
%figure(hf);
h_p2 = plot(t_vect_rs,r_hr_single,'k');
hold on;
h_p2a = plot([t_vect_rs(1) t_vect_rs(end)],av_r_hr_single*ones(1,2),'r--');
xlabel('time (s)');
ylabel('correlation model response');
legend([h_p2 h_p2a],{'correlation model','time-average'});

if ( spatial_averaging == 1 )
    %compute the correlator response averaged over space
    r_hr_mean5 = zeros(size(t_vect_rs));
    
    for j=1:n_x-1 %n_corr-1
        lum_vect0 = stim(:,j); %position1
        lum_vect1 = stim(:,j+1); % at position2 = position1 + dx
        
        %resample
        lv0_rs = reshape(repmat(lum_vect0,1,rs_rate)',[n_trs 1]);
        lv1_rs = reshape(repmat(lum_vect1,1,rs_rate)',[n_trs 1]);
        
        %compute responses of the 1st order low-pass filters
        lv0_s_num = lp1filt_fn(lv0_rs,dt_rs,tau_s);
        lv0_f_num = lp1filt_fn(lv0_rs,dt_rs,tau_f);
        lv1_s_num = lp1filt_fn(lv1_rs,dt_rs,tau_s);
        lv1_f_num = lp1filt_fn(lv1_rs,dt_rs,tau_f);
        
        %half correlator responses
        r_hc1 = lv0_s_num.*lv1_f_num;
        r_hc2 = lv1_s_num.*lv0_f_num;
        
        %full opponent response
        r_hr = r_hc1 - r_hc2;
        
        %add to the mean vector
        r_hr_mean5 = r_hr_mean5 + r_hr';
    end
    
    r_hr_mean5 = r_hr_mean5/(n_x-1);
    
    %compute the mean response over the last 4 seconds of the response
    m_resp_5 = mean(r_hr_mean5(t_vect_rs>1));
    
    %plot the averaged correlator response
    figure;
    h_p7 = plot(t_vect_rs,r_hr_mean);
    hold on;
    h_p8 = plot(t_vect_rs,r_hr_mean5,'r');
    plot([t_vect_rs(1) t_vect_rs(end)],[m_resp_1 m_resp_1],'b--')
    plot([t_vect_rs(1) t_vect_rs(end)],[m_resp_5 m_resp_5],'r--')
    xlabel('time (s)');
    ylabel('response');
    legend([h_p7 h_p8],{'even spatial average','diverging-odd'});
    title('spatially averaged response lr-even and diverging-odd parity gliders');
end

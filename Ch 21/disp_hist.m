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

%% Compute histograms for each data set 
min_lum = min(im(:));
max_lum = max(im(:));

%this will give 100 bins when we do histograms
edges = linspace(min_lum,max_lum,101); %1x101
centers = (edges(2:end) + edges(1:end-1))/2; %1x100
hist_im = zeros(i_max,100);

for i = 1:i_max
    %compute and plots the histogram for each image in a separate figure
    h_f = figure; 
    x = histogram(im(i,:,:),edges); %alternatively, use hist
        
    %saves the data for averaging
    hist_im(i,:) = x.Values;
    
    %get rid of the intermediate figure
    delete(h_f);

end

%% Compute average histogram

hist_im_mean = mean(hist_im,1);
h_f1 = figure; 
h_a1 = axes;
bar(centers,hist_im_mean);
set(h_a1,'TickDir','out');
xlabel('luminance (cd/m2)');
ylabel('average number of observations');

%% Compute histograms for each data set in log coordinates

min_lum = min(log10(im(:)));
max_lum = max(log10(im(:)));

%this will give 100 bins when we do histograms
edges = linspace(min_lum,max_lum,101); %1x101
centers = (edges(2:end) + edges(1:end-1))/2; %1x100
hist_lim = zeros(i_max,100);

for i = 1:i_max
    %compute and plots the histogram for each image in a separate figure
    h_f = figure; 
    x = histogram(log10(im(i,:,:)),edges);
    
    %saves the data for averaging
    hist_lim(i,:) = x.Values;
    
    %get rid of the intermediate figure
    delete(h_f);
end

%% Compute average histogram

hist_lim_mean = mean(hist_lim,1);
h_f2 = figure; 
h_a2 = axes;
bar(centers,hist_lim_mean);
set(h_a2,'TickDir','out');
xlabel('log10 luminance (cd/m2)');
ylabel('average number of observations');

%% Print figures
%print(h_f1,'figures_ed2/nat_scene_histlin.eps','-depsc');
%print(h_f2,'figures_ed2/nat_scene_histlog.eps','-depsc');



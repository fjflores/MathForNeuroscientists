%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  called by destex_f2.m
%
%  requires destexhy3.m

function [currents n_mean n_mean2 n_mean3] = dest_f2b
%compute f-i curve of destexhy3 model
%changes peak conductance of inhibitory 
%conductance
dt = 0.05;

currents = [ 4:2:28 ];
n_currents = length(currents);

n_trials = 10;
n_spk = zeros(n_currents,n_trials);

for j = 1:n_currents
    %info_str = sprintf('processing current %i ...',j);
    %disp(info_str);
    isis = [];
    
    for k = 1:n_trials
        v = destexhy3(dt,10,1010,1020,currents(j),[]);

        %find spike times
        spk = [];
        inds = find(v > 0);
        if ( ~isempty(inds) )
            %at least one spike
            gaps = diff(inds);
            gi = find(gaps>1);
            if ( ~isempty(gi) )
                %at least two spikes
                [valm indm] = max(v(inds(1:gi(1))));
                spk(1) = indm + (inds(1)-1);

                for i =2:length(gi)
                    %if there are 3 or more spikes
                    [valm indm] = max(v(inds(gi(i-1)+1:gi(i))));
                    spk(i) = indm + (inds(gi(i-1)+1)-1);
                end;
        
                %last spike
                [valm indm] = max(v(inds(gi(end)+1:end)));
                spk(end+1) = indm + (inds(gi(end)+1)-1);
            else
                %get the only spike
                [valm indm] = max(v(inds(1:end)));
                spk(1) = indm + (inds(1) - 1);
            end;
        end;
        n_spk(j,k) = length(spk);
    end;
end;

n_mean = mean(n_spk,2);
s_mean = std(n_spk,0,2);

% figure; 
% h = axes; 
% line('Parent',h,'XData',currents,'YData',n_mean,'Marker','o','MarkerFaceColor','k');
% for i = 1:n_currents
%     line('Parent',gca,'XData',[currents(i) currents(i)],...
%         'YData',[n_mean(i)-s_mean(i) n_mean(i)+s_mean(i)]); 
% end;
% set(h,'XLim',[currents(1)-1 currents(n_currents)+1]);


n_spk2 = zeros(n_currents,n_trials);
syn_par = [0.012e-3 0.003e-3 2*0.057e-3 0.0066e-3];

for j = 1:n_currents
    %info_str = sprintf('processing current %i ...',j);
    %disp(info_str);
    isis = [];
    
    for k = 1:n_trials
        v = destexhy3(dt,10,1010,1020,currents(j),syn_par);

        %find spike times
        spk = [];
        inds = find(v > 0);
        if ( ~isempty(inds) )
            %at least one spike
            gaps = diff(inds);
            gi = find(gaps>1);
            if ( ~isempty(gi) )
                %at least two spikes
                [valm indm] = max(v(inds(1:gi(1))));
                spk(1) = indm + (inds(1)-1);

                for i =2:length(gi)
                    %if there are 3 or more spikes
                    [valm indm] = max(v(inds(gi(i-1)+1:gi(i))));
                    spk(i) = indm + (inds(gi(i-1)+1)-1);
                end;
        
                %last spike
                [valm indm] = max(v(inds(gi(end)+1:end)));
                spk(end+1) = indm + (inds(gi(end)+1)-1);
            else
                %get the only spike
                [valm indm] = max(v(inds(1:end)));
                spk(1) = indm + (inds(1) - 1);
            end;
        end;
        n_spk2(j,k) = length(spk);
    end;
end;

n_mean2 = mean(n_spk2,2);
s_mean2 = std(n_spk2,0,2);

% line('Parent',h,'XData',currents,'YData',n_mean2,'Marker','o','MarkerFaceColor','r');
% for i = 1:n_currents
%     line('Parent',gca,'XData',[currents(i) currents(i)],...
%         'YData',[n_mean2(i)-s_mean2(i) n_mean2(i)+s_mean2(i)],'Color','r'); 
% end;

n_spk3 = zeros(n_currents,n_trials);
syn_par = [0.012e-3 0.003e-3 0.5*0.057e-3 0.0066e-3];

for j = 1:n_currents
    %info_str = sprintf('processing current %i ...',j);
    %disp(info_str);
    isis = [];
    
    for k = 1:n_trials
        v = destexhy3(dt,10,1010,1020,currents(j),syn_par);

        %find spike times
        spk = [];
        inds = find(v > 0);
        if ( ~isempty(inds) )
            %at least one spike
            gaps = diff(inds);
            gi = find(gaps>1);
            if ( ~isempty(gi) )
                %at least two spikes
                [valm indm] = max(v(inds(1:gi(1))));
                spk(1) = indm + (inds(1)-1);

                for i =2:length(gi)
                    %if there are 3 or more spikes
                    [valm indm] = max(v(inds(gi(i-1)+1:gi(i))));
                    spk(i) = indm + (inds(gi(i-1)+1)-1);
                end;
        
                %last spike
                [valm indm] = max(v(inds(gi(end)+1:end)));
                spk(end+1) = indm + (inds(gi(end)+1)-1);
            else
                %get the only spike
                [valm indm] = max(v(inds(1:end)));
                spk(1) = indm + (inds(1) - 1);
            end;
        end;
        n_spk3(j,k) = length(spk);
    end;
end;

n_mean3 = mean(n_spk3,2);
s_mean3 = std(n_spk3,0,2);

% line('Parent',h,'XData',currents,'YData',n_mean3,'Marker','o','MarkerFaceColor','b');
% for i = 1:n_currents
%     line('Parent',gca,'XData',[currents(i) currents(i)],...
%         'YData',[n_mean3(i)-s_mean3(i) n_mean3(i)+s_mean3(i)],'Color','b'); 
% end;

%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  called by destex_f2.m

function [n_mean s_mean m_isis c_isis] = dest_f2

%compute f-i curve of destexhy2 model
%as well as cv 

dt = 0.05;

currents = [ 4.5:0.5:9 ];
n_currents = length(currents);

n_trials = 10;
n_spk = zeros(n_currents,n_trials);
m_isis = zeros(1,n_currents);
s_isis = zeros(1,n_currents);
n_isis = zeros(1,n_currents);

for j = 1:n_currents
    %info_str = sprintf('processing current %i ...',j);
    %disp(info_str);
    isis = [];
    
    for k = 1:n_trials
        v = destexhy2(dt,10,1010,1020,currents(j));

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
        isis = [isis diff(spk)];
        n_spk(j,k) = length(spk);
    end;
    m_isis(j) = mean(isis);
    s_isis(j) = std(isis);
    n_isis(j) = length(isis);
end;

n_mean = mean(n_spk,2);
s_mean = std(n_spk,0,2);
c_isis = s_isis./m_isis;

%figure; 
%h = axes; 
%line('Parent',h,'XData',currents,'YData',n_mean,'Marker','o','MarkerFaceColor','k');
%for i = 1:n_currents
%    line('Parent',gca,'XData',[currents(i) currents(i)],...
%        'YData',[n_mean(i)-s_mean(i) n_mean(i)+s_mean(i)]); 
%end;
%set(h,'XLim',[currents(1)-1 currents(n_currents)+1]);

%figure; 
%plot(m_isis,s_isis./m_isis,'o');



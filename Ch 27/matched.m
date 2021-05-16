%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

t = 0:0.01:0.25;

alpha = 2*pi*5.5;

f = t.*(1 - (alpha/2)*t).*exp(-alpha*t);
scale_factor = 5/max(f);
f = scale_factor*f;

baseline = 25;
f = f + baseline;
l_f = length(f);

%for display purposes
t_fine = 0:0.001:0.25;
f_fine = t_fine.*(1 - (alpha/2)*t_fine).*exp(-alpha*t_fine);
f_fine = scale_factor*f_fine;
f_fine = f_fine + baseline;

h_f1 = figure; 
h_a1 = subplot(2,2,1);
line('Parent',h_a1,'XData',t_fine,'YData',f_fine,'Color','r');
line('Parent',h_a1,'XData',[0 0.25],'YData',baseline*[1 1],'Color','k','LineStyle','--');

n_trials = 1000;
p_val = n_trials/10;
ms_fr50 = zeros(1,n_trials);
ms_fr250 = zeros(1,n_trials);
ms_match = zeros(1,n_trials);

variability = 5;
for i = 1:n_trials
    r = normrnd(f,variability);
    ms_fr50(1,i) = mean(r(1:5));
    ms_fr250(1,i) = mean(r);
    ms_match(1,i) = (r - baseline)*(f - baseline)';
    if ( mod(i,p_val) == 0 )
        line('Parent',h_a1,'XData',t-0.002,'YData',r,...
         'Marker','o','MarkerSize',1,'MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none');
    end;
end;

mn_fr50 = zeros(1,n_trials);
mn_fr250 = zeros(1,n_trials);
mn_match = zeros(1,n_trials);

for i = 1:n_trials
    r = normrnd(baseline*ones(size(f)),variability);
    mn_fr50(1,i) = mean(r(1:5));
    mn_fr250(1,i) = mean(r);
    mn_match(1,i) = (r - baseline)*(f - baseline)';
    if ( mod(i,p_val) == 0 )
        line('Parent',h_a1,'XData',t+0.002,'YData',r,...
         'Marker','o','MarkerSize',1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineStyle','none');
    end;
end;

set(h_a1,'XLim',[0 0.25],'YLim',[0 45]);
xlabel(h_a1,'time (s)');
ylabel(h_a1,'firing rate (spk/s)');

[h1, x1] = hist([mn_fr50; ms_fr50]',   20); 
h1 = h1/n_trials; 
h_a2 = subplot(2,2,2);
bar(h_a2,x1,h1);
set(h_a2,'XLim',[15 35],'TickDir','out');
xlabel(h_a2,'mean firing rate 0-50 ms (spk/s)');
ylabel(h_a2,'probability');

h2 = hist([mn_fr250; ms_fr250]', 20); 
h2 = h2/n_trials;

[h3, x3] = hist([mn_match; ms_match]', 20); 
h3 = h3/n_trials; 
h_a3 = subplot(2,2,3);
bar(h_a3,x3,h3);
set(h_a3,'XLim',[-180 220],'TickDir','out');
xlabel(h_a3,'log-likelihood (matched filter)')
ylabel(h_a3,'probability'); 

%yields P_FA and P_D
ch1 = 1 - cumsum(h1,1);
ch2 = 1 - cumsum(h2,1);
ch3 = 1 - cumsum(h3,1);

h_a4 = subplot(2,2,4);
line('Parent',h_a4,'XData',ch1(:,1),'YData',ch1(:,2));
line('Parent',h_a4,'XData',ch2(:,1),'YData',ch2(:,2),'LineStyle','--');
line('Parent',h_a4,'XData',ch3(:,1),'YData',ch3(:,2),'Color','r');
xlabel(h_a4,'probability of false-alarm');
ylabel(h_a4,'probability of detection');

%compute the minimum error
err1 = 0.5*ch1(:,1) + 0.5*(1-ch1(:,2));
%figure; plot(err1); %for debugging
info_str = sprintf('error 1: %.2f',min(err1));
disp(info_str);

err2 = 0.5*ch2(:,1) + 0.5*(1-ch2(:,2));
info_str = sprintf('error 2: %.2f',min(err2));
disp(info_str);

err3 = 0.5*ch3(:,1) + 0.5*(1-ch3(:,2));
info_str = sprintf('error 3: %.2f',min(err3));
disp(info_str);

%print(handles.figure1,'-depsc2','matched.eps');


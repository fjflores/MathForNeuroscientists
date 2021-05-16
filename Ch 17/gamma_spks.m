%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

%%
n_samples = 10000;

m_isi = 10; %in ms (without ref period)
tref = 2; %in ms

max_plot = 400;

%10 ms isi, n_samples samples
v_1 = exprnd(m_isi,n_samples,1);

%add a 2ms refractory period
v_1r = v_1 + tref;

%compute spike times
spk_1r = cumsum(v_1r);

inds = find(spk_1r <= max_plot);

h_f1 = figure; 
h_a1 = subplot(5,2,1);

for i = 1:length(inds);
    line('Parent',h_a1,'XData',[spk_1r(inds(i)) spk_1r(inds(i))],'YData',[0 1]);
end;
set(h_a1,'XLim',[0 max_plot]);

h_a2 = subplot(5,2,2);
[n, xout] = hist(v_1r,[0:60]);
bar(h_a2,xout,n/n_samples);
set(h_a2,'XLim',[0 60],'TickDir','out');
m_v1r = mean(v_1r);
s_v1r = std(v_1r);
cv_1 = s_v1r/m_v1r;
info_str = sprintf('n = 1 mean ISI: %.2g',mean(v_1r));
disp(info_str);
title(h_a2,sprintf('n = 1 CV = %.2g',cv_1));

y = exppdf(xout-2,m_isi);
y = y/sum(y);
line('Parent',h_a2,'XData',xout,'YData',y','Color','r');

%%
%gamma of order 2

a = 2; 

%10 ms isi, n_samples samples
v_2 = gamrnd(a,m_isi/a,n_samples,1);

%add a 2ms refractory period
v_2r = v_2 + tref;

%compute spike times
spk_2r = cumsum(v_2r);

inds = find(spk_2r <= max_plot);

h_a3 = subplot(5,2,3);
for i = 1:length(inds);
    line('Parent',h_a3,'XData',[spk_2r(inds(i)) spk_2r(inds(i))],'YData',[0 1]);
end;
set(h_a3,'XLim',[0 max_plot]);

h_a4 = subplot(5,2,4);
[n, xout] = hist(v_2r,[0:60]);
bar(h_a4,xout,n/n_samples);
set(h_a4,'XLim',[0 60],'TickDir','out');
m_v2r = mean(v_2r);
s_v2r = std(v_2r);
cv_2 = s_v2r/m_v2r;
info_str = sprintf('n = 2 mean ISI: %.2g',mean(v_2r));
disp(info_str);
title(h_a4,sprintf('n = 2 CV = %.2g',cv_2));

y = gampdf(xout-2,a,m_isi/a);
y = y/sum(y);
line('Parent',h_a4,'XData',xout,'YData',y','Color','r');

%%
%gamma of order 5

a = 5;

%10 ms isi, n_samples samples
v_5 = gamrnd(a,m_isi/a,n_samples,1);

%add a 2ms refractory period
v_5r = v_5 + tref;

%compute spike times
spk_5r = cumsum(v_5r);

inds = find(spk_5r <= max_plot);

h_a5 = subplot(5,2,5);
for i = 1:length(inds);
    line('Parent',h_a5,'XData',[spk_5r(inds(i)) spk_5r(inds(i))],'YData',[0 1]);
end;
set(h_a5,'XLim',[0 max_plot]);

h_a6 = subplot(5,2,6);
[n, xout] = hist(v_5r,[0:60]);
bar(h_a6,xout,n/n_samples);
set(h_a6,'XLim',[0 60],'TickDir','out');
m_v5r = mean(v_5r);
s_v5r = std(v_5r);
cv_5 = s_v5r/m_v5r;
info_str = sprintf('n = 5 mean ISI: %.2g',mean(v_5r));
disp(info_str);
title(h_a6,sprintf('n = 5 CV = %.2g',cv_5));

y = gampdf(xout-2,a,m_isi/a);
y = y/sum(y);
line('Parent',h_a6,'XData',xout,'YData',y','Color','r');

%%
%gamma of order 10

a = 10;

%10 ms isi, n_samples samples
v_10 = gamrnd(a,m_isi/a,n_samples,1);

%add a 2ms refractory period
v_10r = v_10 + tref;

%compute spike times
spk_10r = cumsum(v_10r);

inds = find(spk_10r <= max_plot);

h_a7 = subplot(5,2,7);
for i = 1:length(inds);
    line('Parent',h_a7,'XData',[spk_10r(inds(i)) spk_10r(inds(i))],'YData',[0 1]);
end;
set(h_a7,'XLim',[0 max_plot]);

h_a8 = subplot(5,2,8);
[n, xout] = hist(v_10r,[0:60]);
bar(h_a8,xout,n/n_samples);
set(h_a8,'XLim',[0 60],'TickDir','out');
m_v10r = mean(v_10r);
s_v10r = std(v_10r);
cv_10 = s_v10r/m_v10r;
info_str = sprintf('n = 10 mean ISI: %.2g',mean(v_10r));
disp(info_str);
title(h_a8,sprintf('n = 10 CV = %.2g',cv_10));

y = gampdf(xout-2,a,m_isi/a);
y = y/sum(y);
line('Parent',h_a8,'XData',xout,'YData',y','Color','r');

%%
%perfect integrator

n_spk = ceil(max_plot/(m_isi + tref));
spk_infr = (m_isi + tref)*[0:1:n_spk+1];
inds = find(spk_infr <= max_plot);

h_a9 = subplot(5,2,9);
for i = 1:length(inds);
    line('Parent',h_a9,'XData',[spk_infr(inds(i)) spk_infr(inds(i))],'YData',[0 1]);
end;
set(h_a9,'XLim',[0 max_plot]);

h_a10 = subplot(5,2,10);
line('Parent',h_a10,'XData',[m_isi+tref m_isi+tref],'YData',[0 1]);
set(h_a10,'XLim',[0 60],'TickDir','out');
info_str = sprintf('n = infinity mean ISI: 12');
disp(info_str);
title(h_a10,sprintf('n = \\infty CV = 0'));

xlabel(h_a9,'time (ms)');
xlabel(h_a10,'ISI (ms)');

%print(handles.figure1,'-depsc2','gamma_spks.eps');

%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

%%
%peak firing rate 
pfr = 30; %in spk/s
vfr = 5; %variability in spk/s

%compute the average mean square error between stimulus and reconstruction
n_rand = 5000; %number of random samples

%space dimension
d_v = 2;


%number of even spaced vectors phi{1}, ... phi{n} (in addition this 
%implies phi{-1}, ... phi{-n}, where phi{-k} = - phi{k}
n_basis_v = 3:2:15;

%normalized root mean squared error
rmse_n = zeros(1,length(n_basis_v));
%corresponding standard deviation
rmss_n = zeros(1,length(n_basis_v));

for j = 1:length(n_basis_v)
    %generate basis vectors evenly distributed on the
    %unit circle
    n_basis = n_basis_v(j);
    ang_basis = 2*pi*[0:n_basis-1]/n_basis;
    e = [cos(ang_basis); sin(ang_basis)];

    %normalized square error for each sample 
    err2_n = zeros(1,n_rand);
    
    fr = zeros(n_basis,1);
    rfr = zeros(n_basis,1);

    for i=1:n_rand
        %select a random vector
        r_num = rand(1,2);
        ang = 2*pi * r_num(1); %evenly distributed between 0 and 2pi
        len = 3*pfr * r_num(2); %evenly distributed between 0 and 3 times the peak firing rate
        v = len*[cos(ang); sin(ang)];
        
        %corresponding firing rates along the basis 
        %vectors
        fr = e'*v;
    
        %generate a random response
        rfr = randn(n_basis,1)*vfr + fr;
    
        %compute the corresponding vector estimates
        fr_mat = (d_v/n_basis) * e*rfr;
        err_coord = fr_mat - v;
        err2_tot = sum(err_coord.^2,1);%squared distance between random sample and stimulus
        
        %normalize by expected variance
        err2_n(i) = err2_tot/vfr^2;
    end;
        
    rmse_n(j) = sqrt(mean(err2_n));
    rmss_n(j) = sqrt(std(err2_n));
end;

h_f1 = figure;
h_a1 = subplot(2,2,1);
line('Parent',h_a1,'XData',n_basis_v,'YData',rmse_n,...
    'Marker','o','MarkerFaceColor','k','MarkerSize',4);
y2 = d_v./sqrt(n_basis_v);
line('Parent',h_a1,'XData',n_basis_v,'YData',y2,'Color','r');
for i = 1:length(n_basis_v)
    line('Parent',h_a1,'XData', [n_basis_v(i) n_basis_v(i)],...
        'YData',[rmse_n(i)-rmss_n(i) rmse_n(i)+rmss_n(i)]);
end;
set(h_a1,'XLim',[2 15]);
set(h_a1,'YLim',[-1 3]);
xlabel(h_a1,'neuron pairs');
ylabel(h_a1,'normalized RMSE');

%%
%use the latest base to generate a sample firing rate 
%distribute the angle randomly between pi/2 (90deg) and 3pi/2 (270 deg)
ang = (pi/2) + pi*rand(1);
ang_d = (180/pi) * ang;
v = pfr*[cos(ang); sin(ang)];
fr = e'*v;
rfr = randn(n_basis,1)*vfr + fr;
ang_basis_d = (180/pi) * ang_basis;
ang_basis_m = rem(ang_basis_d + 180,360);
fr_val = zeros(1,2*n_basis);
for i = 1:n_basis
    if rfr(i) > 0
        fr_val(i) = rfr(i);
    else
        fr_val(n_basis + i) = -rfr(i);
    end;
end;
    
ang_basis_b = [ ang_basis_d ang_basis_m ];
[ang_sort, inds] = sort(ang_basis_b);
fr_sort = fr_val(inds);

h_a2 = subplot(2,2,2);
line('Parent',h_a2,'XData',ang_sort,'YData',fr_sort,...
    'Marker','o','MarkerFaceColor','k','MarkerSize',4,'LineStyle',':');
line('Parent',h_a2,'XData',[ang_d ang_d],'YData',[0 30],'Color','r');
set(h_a2,'XLim',[0 360]);
xlabel(h_a2,'preferred direction (deg)');
ylabel(h_a2,'firing rate (spk/s)');

h_a3 = subplot(2,2,3);
%generate a plot of the corresponding directions
for i = 1:15
line('Parent',h_a3,'XData',[0 e(1,i)],'YData',[0 e(2,i)]);
line('Parent',h_a3,'XData',[0 -e(1,i)],'YData',[0 -e(2,i)],'Color','r');
end

%print(handles.figure1,'-depsc2','over_rep.eps');

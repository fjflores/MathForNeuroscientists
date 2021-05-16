%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%


%first plot the correlation function 
c = 0.38;
rho = 1;

x = -180:1:180;
x_r = (pi/180)*x;
y_rho = c*exp(-abs(x_r)/rho);

h_f1 = figure; 
h_a1 = subplot(2,2,1);
line('Parent',h_a1,'XData',x,'YData',y_rho);

%change correlation length
nu = 0.125;
y_nu = c*exp(-abs(x_r)/nu);

line('Parent',h_a1,'XData',x,'YData',y_nu,'Color','r');
xmin = -45; xmax = 45;
set(h_a1,'XLim',[xmin xmax],'YLim',[-0.02 0.4]);

%plot the position of the basis vectors
n_basis = 95;

%in radians, from 0 to 2pi 
ang_basis = 2*pi*[0:n_basis-1]/n_basis;

%in degrees, from 0 to 360
ang_basis_d = (180/pi) * ang_basis;

%in degrees, from -180 to 180
ang_basis_f = ang_basis_d - 180; 

%same thing for the negative vectors
ang_basis_m = mod(ang_basis_d + 180,360) - 180;

for i = 1:length(ang_basis_f)
    c_ang = ang_basis_f(i);
    if ( (c_ang <= xmax) & (c_ang >= xmin) )
        line('Parent',h_a1,'XData',ang_basis_f(i),'YData',-0.015,...
            'Marker','o','MarkerSize',2,'Color','r');
    end;
end;

for i = 1:length(ang_basis_m)
    c_ang = ang_basis_m(i);
    if ( (c_ang <= xmax) & (c_ang >= xmin) )
        line('Parent',h_a1,'XData',ang_basis_m(i),'YData',-0.015,...
            'Marker','o','MarkerSize',2,'Color','k');
    end;
end;
xlabel(h_a1,'angle (deg)');
ylabel(h_a1,'correlation coeff.');

%space dimension
d_v = 2;

%basis vectors
e = [cos(ang_basis); sin(ang_basis)];

%peak firing rate 
pfr = 30; %in spk/s
vfr = 5; %variability in spk/s
n_rand = 5000; %number of random samples

%do reconstructions for a range of correlation coefficients
c_v = [0:0.1:0.5];

%root mean square error obtained by the normal (non-optimal)
%method
norm_rmse_n = zeros(size(c_v));
%corresponding sd
norm_rmss_n = zeros(size(c_v));
%optimal rmse    
opt_rmse_n = zeros(size(c_v));
%corresponding sd
opt_rmss_n = zeros(size(c_v));

%rho %debugging only
for k = 1:length(c_v)
    %set up the covariance matrix
    cmat = zeros(n_basis,n_basis);
    %c = 0; %for debugging gets rid of all correlations
    for i = 1:n_basis
        for j = 1:n_basis
            if ( i == j )
                cmat(i,j) = vfr^2;
            else
                cmat(i,j) = vfr^2*c_v(k)*exp(-abs(ang_basis(i) - ang_basis(j))/rho);
            end;
        end;
    end;

    %compute the Choleski factorization
    rmat = chol(cmat);

    %compute the optimal reconstruction matrix
    %in our notation e = U^t and the inverse is given by 
    %(U^t C^(-1) U)^(-1) U^t C^(-1) 
    U = e';
    cmati = inv(cmat);
    U_d = inv(U'*cmati*U)*U'*cmati;

    %normalized root mean square error for each reconstruction
    %using the uncorrelated optimal algorithm
    err2_n = zeros(1,n_rand);
    
    %normalized root mean square error for each reconstruction 
    %using the optimal algorithm for correlated firing rates
    cerr2_n = zeros(1,n_rand);
    
    %mean firing rate and random component
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
        rfr = rmat'*randn(n_basis,1) + fr;
    
        %compute the corresponding vector estimates
        fr_mat = (d_v/n_basis) * e*rfr;
        err_coord = fr_mat - v;
        err2_tot = sum(err_coord.^2,1); %distance between random sample and stimulus
    
        %normalize by vector length
        err2_n(i) = err2_tot/vfr^2;
    
        %do the same thing with the optimal estimator
        fr_mat_o = U_d * rfr;
    
        err_coord_o = fr_mat_o - v;
        err2_tot_o = sum(err_coord_o.^2,1); %distance between random sample and stimulus
    
        %normalize by vector length
        cerr2_n(i) = err2_tot_o/vfr^2;
    end;

    %root mean squared error
    norm_rmse_n(k) = sqrt(mean(err2_n));
    %root sd
    norm_rmss_n(k) = sqrt(std(err2_n));
   
    opt_rmse_n(k) = sqrt(mean(cerr2_n));
    opt_rmss_n(k) = sqrt(std(cerr2_n));
end;

h_a2 = subplot(2,2,2);
line('Parent',h_a2,'XData',c_v,'YData',norm_rmse_n,'Marker','x');
line('Parent',h_a2,'XData',c_v,'YData',opt_rmse_n,...
    'Marker','o','MarkerSize',2,'MarkerFaceColor','k');
set(h_a2,'XLim',[0 0.55],'YLim',[0.15 0.6]);
set(h_a2,'XTick',[0 0.1 0.2 0.3 0.4 0.5]);
xlabel(h_a2,'correlation coeff.');
ylabel(h_a2,'normalized RMSE');

%do reconstructions for a range of correlation lengths
%rho_v = [1:0.5:5];
rho_v = [0.125 0.25 0.5 1 2 4 8 16];

%root mean square error obtained by the normal (non-optimal)
%method
norm_rmse_n = zeros(size(rho_v));
%corresponding sd
norm_rmss_n = zeros(size(rho_v));
%optimal rmse    
opt_rmse_n = zeros(size(rho_v));
%corresponding sd
opt_rmss_n = zeros(size(rho_v));
   
%c %debugging only
for k = 1:length(rho_v)
    %set up the covariance matrix
    cmat = zeros(n_basis,n_basis);
    %c = 0; %for debugging gets rid of all correlations
    for i = 1:n_basis
        for j = 1:n_basis
            if ( i == j )
                cmat(i,j) = vfr^2;
            else
                cmat(i,j) = vfr^2*c*exp(-abs(ang_basis(i) - ang_basis(j))/rho_v(k));
            end;
        end;
    end;

    %compute the Choleski factorization
    rmat = chol(cmat);

    %compute the optimal reconstruction matrix
    %in our notation e = U^t and the inverse is given by 
    %(U^t C^(-1) U)^(-1) U^t C^(-1) 
    U = e';
    cmati = inv(cmat);
    U_d = inv(U'*cmati*U)*U'*cmati;

    %normalized root mean square error for each reconstruction
    %using the uncorrelated optimal algorithm
    err2_n = zeros(1,n_rand);
    
    %normalized root mean square error for each reconstruction 
    %using the optimal algorithm for correlated firing rates
    cerr2_n = zeros(1,n_rand);
    
    %mean firing rate and random component
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
        rfr = rmat'*randn(n_basis,1) + fr;
    
        %compute the corresponding vector estimates
        fr_mat = (d_v/n_basis) * e*rfr;
        err_coord = fr_mat - v;
        err2_tot = sum(err_coord.^2,1); %distance between random sample and stimulus
    
        %normalize by vector length
        err2_n(i) = err2_tot/vfr^2;
    
        %do the same thing with the optimal estimator
        fr_mat_o = U_d * rfr;
    
        err_coord_o = fr_mat_o - v;
        err2_tot_o = sum(err_coord_o.^2,1); %distance between random sample and stimulus
    
        %normalize by vector length
        cerr2_n(i) = err2_tot_o/vfr^2;
    end;

    %root mean squared error
    norm_rmse_n(k) = sqrt(mean(err2_n));
    %root sd
    norm_rmss_n(k) = sqrt(std(err2_n));
   
    opt_rmse_n(k) = sqrt(mean(cerr2_n));
    opt_rmss_n(k) = sqrt(std(cerr2_n));
end;

h_a3 = subplot(2,2,3);
line('Parent',h_a3,'XData',1:length(rho_v),'YData',norm_rmse_n,'Marker','x');
line('Parent',h_a3,'XData',1:length(rho_v),'YData',opt_rmse_n,...
    'Marker','o','MarkerSize',2,'MarkerFaceColor','k');
 
set(h_a3,'XTick',[1 3 5 7],'XTickLabel',{'0.125', '0.5',  '2', '8'});
set(h_a3,'XLim',[0 9],'YLim',[0.2 0.55]);
xlabel(h_a3,'correlation length (rad)');
ylabel(h_a3,'normalized RMSE');

%generate example firing rates in the case of correlated and uncorrelated
%neurons. 
ang = pi;
ang_d = (180/pi) * ang;
v = pfr*[cos(ang); sin(ang)];

%uncorrelated case
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

h_a4 = subplot(2,2,4);
line('Parent',h_a4,'XData',ang_sort,'YData',fr_sort);
line('Parent',h_a4,'XData',[ang_d ang_d],'YData',[0 30],'Color','r');

%correlated case

%First, set up the covariance matrix
cmat = zeros(n_basis,n_basis);
c = 0.38; 
rho = 1;
for i = 1:n_basis
    for j = 1:n_basis
    	if ( i == j )
        	cmat(i,j) = vfr^2;
        else
        	cmat(i,j) = vfr^2*c*exp(-abs(ang_basis(i) - ang_basis(j))/rho);
        end;
    end;
end;

%compute the Choleski factorization
rmat = chol(cmat);
    
fr = e'*v;
rfr = rmat'*randn(n_basis,1) + fr;
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

line('Parent',h_a4,'XData',ang_sort,'YData',fr_sort,'Color',[.5 0.5 0.5]);
set(h_a4,'XLim',[0 360],'YLim',[0 50]);
set(h_a4,'XTick',[0 90 180 270 360]);
xlabel(h_a4,'angle (deg)');
ylabel(h_a4,'firing rate (spk/s)');

%print(handles.figure1,'-depsc2','pop_corr.eps');

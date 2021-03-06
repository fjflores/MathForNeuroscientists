%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%

%plot the autocorrelation and chi_c function of a gamma order
%two and 10 process.

clear all; close all;

%paramters of gamma process
rho = 160; %in spk/s
n_g = 2;

%mean rate in spk/s
chi = rho/n_g; 

tv_ms = [0:0.1:30];
tv_s = tv_ms*1e-3;
 
%compute chi_c
k_v = 1;
n_k = length(k_v);
zk = exp(2*pi*i*k_v/n_g);
t_mat = repmat(tv_s,n_k,1);
t_mat2 = zeros(size(t_mat));
for j = 1:n_k
    t_mat2(j,:) = rho*(zk(j)-1)*t_mat(j,:);
end;
rhoc = chi*(1+ zk*exp(t_mat2));

%get rid of potential imaginary part due 
%to numerical errors
rhoc_r = real(rhoc);
figure(1);
plot(tv_ms,rhoc_r,'k');
set(gca,'TickDir','out');
xlabel('time (ms)');
ylabel('xc(t) ((spk/s)^2)');

%uncomment to print figure
%print(1,'-depsc2','gamma_corr_p1.eps');

%autocorrelation
tv_ms2 = [-tv_ms(end:-1:2) tv_ms];
rhoc_s = [rhoc_r(end:-1:2) rhoc_r];
c_g = chi*(rhoc_s-chi);


figure(2); 
plot(tv_ms2,c_g,'k');
hold on;
plot([-30 30],[0 0],'k--');
plot([0 0],[-7000 1000],'k--');
set(gca,'TickDir','out');
xlabel('time (ms)');
ylabel('C(t) ((spk/s)^2)');

%uncomment to print figure
%print(2,'-depsc2','gamma_corr_p2.eps');

%set(gca,'YLim',[-4000 4000]);

%clear all;
rho = 800; %in spk/s
n_g = 10;

%mean rate in spk/sec
chi = rho/n_g;

tv_ms = [0:0.1:40];
tv_s = tv_ms*1e-3;
k_v = [1:1:9];
n_k = length(k_v);
zk = exp(2*pi*i*k_v/n_g);
t_mat = repmat(tv_s,n_k,1);
t_mat2 = zeros(size(t_mat));
for j = 1:n_k
    t_mat2(j,:) = rho*(zk(j)-1)*t_mat(j,:);
end;
rhoc = chi*(1+ zk*exp(t_mat2));
rhoc_r = real(rhoc);

tv_ms2 = [-tv_ms(end:-1:2) tv_ms];
rhoc_s = [rhoc_r(end:-1:2) rhoc_r];
c_g = chi*(rhoc_s-chi);
figure(3);
plot(tv_ms2,c_g,'k');
hold on;
plot([-40 40],[0 0],'k--');
plot([0 0],[-7000 3000],'k--');
set(gca,'TickDir','out');
xlabel('time (ms)');
ylabel('C(t) ((spk/s)^2)');

%uncomment to print figure
%print(3,'-depsc2','gamma_corr_p3.eps');
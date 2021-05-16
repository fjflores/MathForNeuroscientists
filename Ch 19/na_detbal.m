%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%

%target membrane potential
v_v = -37.94;

am_v = .1*(25-(v_v+71))./(exp(2.5-(v_v+71)/10)-1);
bm_v = 4*exp(-(v_v+71)/18);
taum_v = 1./(am_v+bm_v);
minf_v = am_v.*taum_v;

ah_v  = 0.07*exp(-(v_v+71)/20);
bh_v = 1./(exp(3-(v_v+71)/10)+1);
tauh_v = 1./(ah_v+bh_v);
hinf_v = ah_v.*tauh_v;

%rates, converted in 1/s = Hz from 1/ms = kHz
am = am_v*1e3;
bm = bm_v*1e3;
ah = ah_v*1e3;
bh = bh_v*1e3;

%mean closed and open state durations for m and h
m_am = 1/am;
m_bm = 1/bm;
m_ah = 1/ah;
m_bh = 1/bh;

%average transition rate
tm = 0.5*(am + am);
th = 0.5*(ah + ah);

%average time between transitions
m_tm = 1/tm;
m_th = 1/th;


info_str = sprintf('m gate: mean time in closed state = %.2g  open state: %.2g (in s)',m_am,m_bm);
disp(info_str);
info_str = sprintf('        mean time between transitions = %.2g (in s)\n',m_tm);
disp(info_str);

info_str = sprintf('h gate: mean time in closed state = %.2g  open state: %.2g (in s)',m_ah,m_bh);
disp(info_str);
info_str = sprintf('        mean time between transitions = %.2g (in s)\n\n',m_th);
disp(info_str);


%define transition matrices
Qm = [-am am; bm -bm];
Qh = [-ah ah; bh -bh];

%compute EV for Qm
[vm, dm] = eig(Qm');

%find the smallest eigenvalue
[c_min, i_min] = min(abs(diag(dm)));
info_str = sprintf('minimum Qm eigenvalue: %.2g, index: %i',c_min,i_min);
disp(info_str);

[c_max, i_max] = max(abs(diag(dm)));
info_str = sprintf('maximum Qm eigenvalue: %.2g, index: %i',c_max,i_max);
disp(info_str);

info_str = sprintf('ratio of min to max: %.2g\n',abs(c_min/c_max));
disp(info_str);

%compute equilibrium probabilities
vm_ss = vm(:,i_min)/sum(vm(:,i_min));

%check detailed balance
delta_dbm = vm_ss(1)*Qm(1,2) - vm_ss(2)*Qm(2,1);
info_str = sprintf('Qm deviation from detailed balance: %.2g',delta_dbm);
disp(info_str);

mean_dbm = 0.5*(vm_ss(1)*Qm(1,2) + vm_ss(2)*Qm(2,1));
info_str = sprintf('Qm deviation relative to mean: %.2g\n\n',delta_dbm/mean_dbm);
disp(info_str);

%Repeat for Qh
[vh, dh] = eig(Qh');

%find the smallest eigenvalue
[c_min, i_min] = min(abs(diag(dh)));
info_str = sprintf('minimum eigenvalue: %.2g, index: %i',c_min,i_min);
disp(info_str);

[c_max, i_max] = max(abs(diag(dh)));
info_str = sprintf('maximum eigenvalue: %.2g, index: %i',c_max,i_max);
disp(info_str);

info_str = sprintf('ratio of min to max: %.2g\n',abs(c_min/c_max));
disp(info_str);

%compute equilibrium probabilities
vh_ss = vh(:,i_min)/sum(vh(:,i_min));

%check detailed balance
delta_dbh = vh_ss(1)*Qh(1,2) - vh_ss(2)*Qh(2,1);
info_str = sprintf('Qh deviation from detailed balance: %.2g',delta_dbh);
disp(info_str);

mean_dbh = 0.5*(vh_ss(1)*Qh(1,2) + vh_ss(2)*Qh(2,1));
info_str = sprintf('Qh deviation relative to mean: %.2g\n\n',delta_dbh/mean_dbh);
disp(info_str);

%build the transition matrix for the entire channel
n_st = 8; %number of states
Qna = zeros(n_st,n_st);

%fill in the non-zero diagonal elements
Qna(1,2) = 3*am;
Qna(1,5) = ah;

Qna(2,1) = bm; 
Qna(2,3) = 2*am;
Qna(2,6) = ah;

Qna(3,2) = 2*bm;
Qna(3,4) = am;
Qna(3,7) = ah;

Qna(4,3) = 3*bm;
Qna(4,8) = ah;

Qna(5,1) = bh;
Qna(5,6) = 3*am;

Qna(6,2) = bh;
Qna(6,5) = bm;
Qna(6,7) = 2*am;

Qna(7,3) = bh;
Qna(7,6) = 2*bm;
Qna(7,8) = am;

Qna(8,4) = bh;
Qna(8,7) = 3*bm;

%sets the diagonal element values
for i = 1:n_st
    Qna(i,i) = -sum(Qna(i,:));
end;

%find the EV of Qna
[vna, dna] = eig(Qna');

%find the equilibrium vector
[cna ina] = min(abs(diag(dna)));
vna_ss = vna(:,ina)/sum(vna(:,ina));

%compare to predicted values
vna_pred(1,1) = (vm_ss(1))^3*vh_ss(1);
vna_pred(2,1) = 3*vm_ss(2)*(vm_ss(1))^2*vh_ss(1);
vna_pred(3,1) = 3*(vm_ss(2))^2*(vm_ss(1))*vh_ss(1);
vna_pred(4,1) = (vm_ss(2))^3*vh_ss(1);
vna_pred(5,1) = (vm_ss(1))^3*vh_ss(2);
vna_pred(6,1) = 3*vm_ss(2)*(vm_ss(1))^2*vh_ss(2);
vna_pred(7,1) = 3*(vm_ss(2))^2*(vm_ss(1))*vh_ss(2);
vna_pred(8,1) = (vm_ss(2))^3*vh_ss(2);

vna_diff = vna_ss - vna_pred;
vna_maxdiffrel = max(abs(vna_diff))/min(abs(vna_ss));

str1 = sprintf('%.4g ',vna_ss);
str2 = sprintf('calculated vna_ss: %s',str1);
disp(str2);

str1 = sprintf('%.4g ',vna_pred);
str2 = sprintf('predicted  vna_ss: %s',str1);
disp(str2);

str1 = sprintf('maximal relative difference: %.2g',vna_maxdiffrel);
disp(str1);





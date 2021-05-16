%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%
%  requires destexhy.m destexhy2.m

%compute the distribution of the membrane potential
dt = 0.05;
v_stat = destexhy2(dt,10,100,2000,0);

%skip the first 50 ms
n1 = ceil(50/dt);
mv = mean(v_stat(n1:end));
sv = std(v_stat(n1:end));

text_str = sprintf('mean vm: %.3g std: %.3g',mv,sv);
disp(text_str);

%histogram of values
[nv, xoutv] = hist(v_stat(n1:end),50);

h_f1 = figure; 
h_a2 = subplot(2,2,2);
bar(h_a2,xoutv,nv,'r');
set(h_a2,'TickDir','out','XLim',[-80 -60]);
xlabel(h_a2,'Vm (mV)');
ylabel(h_a2,'number of occurrences');

%plot sample noise trace
t = (1:length(v_stat))*dt;
inds = find( (t>500)&(t<1000));
h_a1 = subplot(2,2,1);
line('Parent',h_a1,'XData',t(inds),'YData',v_stat(inds));
set(h_a1,'XLim',[500 1000]);

vn = destexhy(dt,150,350,550,-0.25);
N = length(vn);
t = (1:N)*dt;
h_a3 = subplot(2,2,3);
line('Parent',h_a3,'XData',t,'YData',vn);
xlabel(h_a3,'time (ms)');
ylabel(h_a3,'Vm (mv)');

reps = 1000;
vb = zeros(N,reps);

disp('iterating 1000 times. This may take a while...');
for i = 1:reps
    if ( mod(i,100) == 0 )
        info_str = sprintf('at iteration %i',i);
        disp(info_str);
    end;
    
    vb(:,i) = destexhy2(dt,150,350,550,-0.25);
end;
vba = mean(vb,2);

h_a4 = subplot(2,2,4);
line('Parent',h_a4,'XData',t,'YData',vba);
xlabel(h_a4,'time (ms)');

inds1 = find( (t>40) & (t<140) );
inds2 = find( (t>200) & (t<300) );

vpre = mean(vba(inds1));
vpul = mean(vba(inds2));
vdiff = vpre-vpul;

text_str2 = sprintf('mean vm prepulse: %.2g during: %.2g difference: %.2g',vpre,vpul,vdiff);
disp(text_str2);

vpren = mean(vn(inds1));
vpuln = mean(vn(inds2));
vdiffn = vpren-vpuln;

text_str3 = sprintf('mean vm prepulse: %.2g during: %.2g difference: %.2g',vpren,vpuln,vdiffn);
disp(text_str3);

%print(handles.figure1,'-depsc2','destex_f1.eps');

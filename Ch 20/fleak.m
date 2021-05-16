%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

%%
N = 32;

x = 0:N-1;
x_fine = 0:0.01:N-1;

f_osc = 1/4;
y = cos( 2*pi*x*f_osc );
y_fine = cos( 2*pi*x_fine*f_osc );

msize = 3;

h_f1 = figure;
h_a1 = subplot(3,1,1);
line('Parent',h_a1,'XData',x,'YData',y,'LineStyle','none','Marker','o',...
    'MarkerFaceColor','k','MarkerSize',msize);
line('Parent',h_a1,'XData',x_fine,'YData',y_fine,'LineStyle','-','Marker','none');
set(h_a1,'XLim',[0 N]);
xlabel(h_a1,'time (normalized)');

%%
z = fft(y);
mz2 = (1/N)*(z.*conj(z));
mz2_c = circshift(mz2,[0 15]); %line vector

%debug: normalized such that sum(y.^2) = sum(mz2)
%sy2 = sum(y.^2)
%smz2 = sum(mz2)

%compute discrete spheroidal sequence of order zero
[e,v] = dpss(32,1);
dpss_win = e(:,1);

[pxx, w] = periodogram(y,dpss_win,N,1,'twosided');
pxxc = circshift(pxx,[15 0]); %column vector

%debug: same normalization as above sum(pxx) = sum(y.^2)
%spxx = sum(pxx)

df = 1/N;
f = (0:N-1)*df;
f_corr = (-N/2+1:N/2)*df;
i_res1 = find(f == f_osc);
i_res2 = find(f == f_osc + 0.5);
inds_lin = [i_res1-1:i_res1+1 i_res2-1:i_res2+1];

i_res3 = find(f_corr == f_osc);
i_res4 = find(f_corr == -f_osc);
inds_linc = [i_res4-1:i_res4+1 i_res3-1:i_res3+1];

%linear scale
h_a2 = subplot(3,1,2);
line('Parent',h_a2,'XData',f_corr,'YData',mz2_c,'LineStyle','none','Marker','o',...
    'MarkerFaceColor','k','MarkerSize',msize);
line('Parent',h_a2,'XData',f_corr(inds_linc),'YData',pxxc(inds_linc),'LineStyle','none','Marker','x',...
    'MarkerEdgeColor','r','MarkerSize',msize);
set(h_a2,'XLim',[-0.5 0.5],'YLim',[0 10]);
ylabel(h_a2,'power (normalized)');


%log scale
h_a3 = subplot(3,1,3);
line('Parent',h_a3,'XData',f_corr,'YData',pow2db(mz2_c),'LineStyle','none','Marker','o',...
    'MarkerFaceColor','k','MarkerSize',msize);
line('Parent',h_a3,'XData',f_corr,'YData',pow2db(pxxc),'LineStyle','none','Marker','x',...
    'MarkerEdgeColor','r','MarkerSize',msize);
set(h_a3,'XLim',[-0.5 0.5],'YLim',[-50 25]);
xlabel(h_a3,'frequency (normalized)','Color','k');
ylabel(h_a3,'dB(power)');

%debug
%figure; subplot(2,1,1); plot(f,mz2); hold on; plot(w,pxx,'xr');
%subplot(2,1,2); plot(f,pow2db(mz2)); hold on; plot(w,pow2db(pxx),'xr');

%line('Parent',h_a3,'XData',f,'YData',mz2,'LineStyle','none','Marker','o',...
%    'MarkerFaceColor','k','MarkerSize',msize);
%set(h_a3,'XLim',[0 1],'YLim',[200 300]);

%%
f1 = f_osc + 0.5*df;
y1 = cos( 2*pi*f1*x );
y1_fine = cos( 2*pi*f1*x_fine );

h_f2 = figure; 
h_a4 = subplot(3,1,1);
line('Parent',h_a4,'XData',x,'YData',y1,'LineStyle','none','Marker','o', ...
    'MarkerFaceColor','k','MarkerSize',msize);
line('Parent',h_a4,'XData',x_fine,'YData',y1_fine,'LineStyle','-','Marker','none');
set(h_a4,'XLim',[0 N]);
xlabel(h_a4,'time (normalized)');


z1 = fft(y1);
mz12 = (1/N)*(z1.*conj(z1));
mz12_c = circshift(mz12,[0 15]); %line vector

%debug: normalized so that sum(y1.^2) = sum(mz12)
%sy12 = sum(y1.^2)
%smz12 = sum(mz12)

[pxx2, w2] = periodogram(y1,dpss_win,N,1,'twosided');
pxx2c = circshift(pxx2,[15 0]); %column vector
%debug: same normalization as above sum(pxx) = sum(y.^2)
%spxx2 = sum(pxx2)

%linear scale
h_a5 = subplot(3,1,2);
line('Parent',h_a5,'XData',f_corr,'YData',mz12_c,'LineStyle','none','Marker','o',...
    'MarkerFaceColor','k','MarkerSize',msize);
line('Parent',h_a5,'XData',f_corr(inds_linc),'YData',pxx2c(inds_linc),'LineStyle','none','Marker','x',...
    'MarkerEdgeColor','r','MarkerSize',msize);
set(h_a5,'XLim',[-0.5 0.5],'YLim',[0 10]);

%log scale
h_a6 = subplot(3,1,3);
line('Parent',h_a6,'XData',f_corr,'YData',pow2db(mz12_c),'LineStyle','none','Marker','o',...
    'MarkerFaceColor','k','MarkerSize',msize);
line('Parent',h_a6,'XData',f_corr,'YData',pow2db(pxx2c),'LineStyle','none','Marker','x',...
    'MarkerEdgeColor','r','MarkerSize',msize);
set(h_a6,'XLim',[-0.5 0.5],'YLim',[-50 25]);
xlabel(h_a6,'frequency (normalized)');

%debug
%figure; subplot(2,1,1); plot(f,mz12); hold on; plot(w,pxx2,'xr');
%subplot(2,1,2); plot(f,pow2db(mz12)); hold on; plot(w,pow2db(pxx2),'xr');

%line('Parent',h_a6,'XData',f,'YData',mz12,'LineStyle','none','Marker','o', ...
%    'MarkerFaceColor','k','MarkerSize',msize);
%set(h_a6,'XLim',[0 1],'YLim',[200 300]);

%%

%time and frequency domain kernels
%dt = 1 so T = N = 32
T = 32;
omega = -0.5:0.001:0.5;
ind0 = (length(omega)-1)/2 + 1;

%continuous (Fejer) kernel
ghat = sin(pi*omega*T).^2./(T*pi^2*omega.^2);
ghat(ind0) = T;

%discrete (Dirichlet) kernel
%ghat2 = sin(pi*omega*T).^2./(T*sin(pi*omega).^2);
%ghat2(ind0) = T;

[e,v] = dpss(32,1);
e_win = e(:,1);

%column vector
ehat = zeros(size(omega));
k_vect = 0:31;

for i = 1:length(ehat)
    arg_exp = (-1)*k_vect*2*pi*omega(i);
    arg_comp = complex(cos(arg_exp),sin(arg_exp));
    ehat(i) = arg_comp*e_win;
end;

ehat = ehat.*conj(ehat);

h_f3 = figure; 
h_a7 = subplot(2,1,1);
line('Parent',h_a7,'XData',1:32,'YData',e_win,'Color','r');
line('Parent',h_a7,'XData',[1 1 32 32],'YData',[0 1/sqrt(32) 1/sqrt(32) 0]);
set(h_a7,'XLim',[0 32]);
xlabel(h_a7,'time (normalized)');

h_a8 = subplot(2,1,2);
line('Parent',h_a8,'XData',omega,'YData',pow2db(ghat)); 
line('Parent',h_a8,'XData',omega,'YData',pow2db(ehat),'Color','r');
line('Parent',h_a8,'XData',[0.5/32 0.5/32],'YData',[-25 10],'LineStyle','--');
line('Parent',h_a8,'XData',[-0.5/32 -0.5/32],'YData',[-25 10],'LineStyle','--');
set(h_a8,'XLim',[-0.15 0.15],'YLim',[-25 20],'XTick',[-0.15 -0.1 -0.05 0 0.05 0.1 0.15]);
xlabel(h_a8,'frequency (normalized)');
ylabel(h_a8,'dB');

%figure; plot(omega,pow2db(ghat));
%hold on; plot(omega,pow2db(ghat2),'r');
%plot(omega,pow2db(ehat),'g');
%set(gca,'YLim',[-25 25]);

%print(handles.figure1,'-depsc2','fleak.eps');

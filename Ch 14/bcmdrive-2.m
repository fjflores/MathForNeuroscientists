%
% Gabbiani and Cox, Mathematics for Neuroscientists
%
%  bcmdrive
%
w = zeros(100,6);
c = zeros(100,6);
outrate = zeros(1,6);

Tminutes = 166; %simulation time in bcmsim, in minutes
mc = zeros(Tminutes,6);
mw = mc;

for r = 1:6,
    rate = r*10;
    [c(:,r), w(:,r), outrate(r), mc(:,r), mw(:,r)] = bcmsim(rate);
end
%save bcmrun2

h_f1 = figure(1);
c = 0:.01:1;
Omega = 1./(1+exp(40*(0.4-c))) - 0.5./(1+exp(60*(0.25-c)));
plot(c,Omega,'k','linewidth',1.5)
xlabel('calcium concentration, c (uM)','fontsize',14)
ylabel('\Omega','fontsize',14)
text(0.1,0.3,'(A)','fontsize',20)
box off
set(gca,'TickDir','out');
%print(h_f1,'figs_fg/omega_func.eps','-depsc');

h_f2 = figure(2);
[nb,xb]=hist(w(:,3));
bh = bar(xb,nb);
set(bh,'facecolor',[1 1 1],'barwidth',1,'linewidth',1);
%xlim([4.6 5])
xlabel('synaptic weight, w','fontsize',14)
ylabel('count','fontsize',14)
text(4.65,20,'(D)','fontsize',20)
box off
set(gca,'TickDir','out');
%print(h_f2,'figs_fg/weight_hist.eps','-depsc');

h_f3 = figure(3);     
f = 10:10:60;
mw_syn = mean(w,1);
[ax,h1,h2] = plotyy(f,mw_syn,f,outrate);
set(ax(1),'ycolor',[1 0 0])
set(h1,'color','red','linewidth',1.5,'marker','o')
set(ax(2),'ycolor',[0 0 0])
set(h2,'color','black','linewidth',1.5,'marker','o')
xlabel('presynaptic spike rate   (Hz)','fontsize',14)
ylabel(ax(1),'average weight','fontsize',14,'color','red')
ylabel(ax(2),'output spike rate (Hz)','fontsize',14)
text(50,4.5,'(C)','fontsize',20)
box off
set(gca,'TickDir','out');
%print(h_f3,'figs_fg/prerate_weightoutrate.eps','-depsc');


h_f4 = figure(4);     
T = 1:Tminutes;
[ax,h1,h2] = plotyy(T,mw(:,3),T,mc(:,3));
set(ax(1),'ycolor',[1 0 0])
set(h1,'color','red','linewidth',1.5); %,'marker','o')
set(ax(2),'ycolor',[0 0 0])
set(h2,'color','black','linewidth',1.5); %,'marker','o')
xlabel('time (minutes)','fontsize',14)
ylabel(ax(1),'average weight','fontsize',14,'color','red')
ylabel(ax(2),'average calcium concentration (uM)','fontsize',14)
text(50,4.5,'(B)','fontsize',20)
box off
set(gca,'TickDir','out');
%print(h_f4,'figs_fg/time_weightca.eps','-depsc');

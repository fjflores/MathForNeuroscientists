%
% Gabbiani and Cox, Mathematics for Neuroscientists
%
% ampareg1.m    Plot activities
%
%  From Castellani et al 2001, PNAS 98:12772-12777
%
%  usage:    ampareg1
%

close all

h_f0 = figure(1);   % first draw the reaction diagram

text(0,5-.2,'a_0','fontsize',18)   % upper branch
hold on
plot([0.5 4.85],[5 5],'k')
plot([4.5 4.85],[5.1 5],'k')
text(2.4,5.35,'k_1(c)','fontsize',14)
plot([0.5 4.85],[4.75 4.75],'k')
plot([0.5 0.85],[4.75 4.65],'k')
text(2.4,4.5-.1,'k_{-1}(c)','fontsize',14)
text(5,5-.2,'a_1','fontsize',18)

text(0,5-.2-5,'a_2','fontsize',18)   % lower branch
plot([0.5 4.85],[5 5]-5,'k')
plot([4.5 4.85],[5.1 5]-5,'k')
text(2.4,5.35-5,'k_1(c)','fontsize',14)
plot([0.5 4.85],[4.75 4.75]-5,'k')
plot([0.5 0.85],[4.75 4.65]-5,'k')
text(2.4,4.5-.1-5,'k_{-1}(c)','fontsize',14)
text(5,5-.2-5,'a_3','fontsize',18)

plot([0 0],[0.1 4.4],'k')     % left branch
plot([-.1 0],[4.05 4.4],'k')
h = text(-.4,2.1,'k_{-2}(c)','fontsize',14);
set(h,'rotation',90)
plot([0.25 0.25],[0.1 4.4],'k')
plot([0.35 0.25],[0.45 0.1],'k')
h = text(0.55,2.1,'k_2(c)','fontsize',14);
set(h,'rotation',90)

plot([0 0]+5,[0.1 4.4],'k')     % right branch
plot([-.1 0]+5,[4.05 4.4],'k')
h = text(-.4+5,2.1,'k_{-2}(c)','fontsize',14);
set(h,'rotation',90)
plot([0.25 0.25]+5,[0.1 4.4],'k')
plot([0.35 0.25]+5,[0.45 0.1],'k')
h = text(0.55+5,2.1,'k_2(c)','fontsize',14);
set(h,'rotation',90)

ylim([-1 6])
xlim([-1 6])

text(-1,5.8,'(A)','fontsize',20)

hold off
axis off
%print(h_f0,'diagram.eps','-depsc');

% now draw the LTD/LTP curve

c = 0:0.02:5; %calcium concentration, uM

km1 = 1+30*(6*c).^2./(1+(6*c).^2);
km2 = km1;
k1 = 1+100*(6*c).^2./(64+(6*c).^2);
k2 = k1;

h_f1 = figure; 
h_p1 = plot(c, km1,'k','linewidth',1);
hold on;
h_p2 = plot(c,k1,'r','linewidth',1);
xlabel('c (uM)','fontsize',14)
ylabel('rate (1/s)','fontsize',14)
set(gca,'TickDir','out');
legend([h_p1 h_p2],{'k_{-1}','k_1'})
box off
%print(h_f1,'phosphokin_rates.eps','-depsc');

den = km1.*km2 + k1.*km2 + k2.*km1 + k1.*k2;

a0 = km1.*km2./den;
a1 = k1.*km2./den;
a2 = k2.*km1./den;
a3 = k1.*k2./den;

h_f2 = figure; 
h_p1 = plot(c,a0,'k');
hold on;
h_p2 = plot(c,a1,'r--');
%plot(c,a2,'g');
h_p3 = plot(c,a3,'r');
xlabel('c (uM)','fontsize',14)
ylabel('steady-state probability','fontsize',14)
set(gca,'TickDir','out');
legend([h_p1 h_p2 h_p3],{'a_0','a_1, a_2', 'a_3'})
box off
%print(h_f2,'phosphokin_states.eps','-depsc');

g = a0 + 2*(a1+a2) + 4*a3;

h_f3 = figure;
plot(c,g/g(1),'k','linewidth',1)
hold on
plot(c,ones(size(c)),'r--')
xlabel('c (uM)','fontsize',14)
ylabel('g_{AMPA}(c)/g_{AMPA}(0)','fontsize',14)
set(gca,'TickDir','out');
box off
%print(h_f3,'gampa_rel.eps','-depsc');

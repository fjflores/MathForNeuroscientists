
%
% gKCa.m
%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% plot the K,Ca functionals
%

function gKCa

v = -120:.01:20;
c = .05;
lint = {'-','--',':'};

figure(3)
for j=1:3,
    c = 10*c;
    a = amKCa(v,c); b = bmKCa(v,c);
    KCatau = 1./(a+b);
    KCainf = a.*KCatau;
    subplot(2,1,1)
    plot(v,KCainf,lint{j},'color','k')
    set(gca,'xticklabel',[])
    hold on
    box off
    subplot(2,1,2)
    plot(v,KCatau,lint{j},'color','k')
    hold on
    box off
end
axis tight
xlabel('V  (mV)','fontsize',14)
ylabel('\tau_{m,KCa}  (ms)','fontsize',14)
legend('c = 0.05 \muM','c = 0.5 \muM', 'c = 5.0 \muM')
legend boxoff
hold off
box off
subplot(2,1,1)
ylabel('m_{\infty,KCa}','fontsize',14)
box off
return

function val = amKCa(v,c)
VT = 25.8;
val = 0.28*c./(c+.48*exp(-1.7*v/VT));

function val = bmKCa(v,c)
VT = 25.8;
val = 0.48./(1+(c/.13e-3).*exp(2*v/VT));


  

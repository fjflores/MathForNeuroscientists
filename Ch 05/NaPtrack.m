
%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% Rest potential and eigenvalues of the quasi-active HH+I_NaP system
% over a range of channel conductances
%
% usage: NaPtrack
%

function NaPtrack

close all               % closes all open figures

for g = 0:.01:1,

    [Vr,B] = hhsymNaP(g);     % find rest potential and quasiactive matrix

    figure(1)
    plot(g,Vr,'.','color',[g 0 0])
    hold on

    figure(2)
    e = eig(B);
    plot(real(e),imag(e),'.','color',[g 0 0])
    hold on

end

figure(1)
xlabel('I_{NaP} (mS/cm^2)','fontsize',14)
ylabel('V_r (mV)','fontsize',14)
box off

figure(2)
xlabel('real(z) (1/ms)','fontsize',14)
ylabel('imag(z) (1/ms)','fontsize',14)
box off

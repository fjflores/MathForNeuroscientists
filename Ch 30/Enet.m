%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  Enet.m
%
% simulate a random (or saved) net of LIF E cells 
%
%  usage:    Enet(N,dt,Tfin,P,den,mk,fname)
%
%  where:    N = number of cells
%           dt = timestep (ms)
%         Tfin = final time (ms)
%            P = period of stimulus (ms)
%          den = connection density (between 0 and 1)
%           mk = spike marker string
%
%  e.g.:   Enet(20,.05,200,30,0.15,'k+')
%          Enet(40,.05,200,50,.07,'k+')
%     Enet(20,.05,200,50,.07,'k+','EnetW20')  % to reproduce Fig. 27.11
%     Enet(40,.05,200,50,.07,'k+','EnetW40')  % to reproduce Fig. 27.12
%

function W = Enet(N,dt,Tfin,P,den,mk,fname)

close all

Nt = ceil(Tfin/dt);

tauE = 2; 

gL = 0.3; VL = -68; Cm = 1;

if nargin == 7
  load(fname)
  N = size(W,1);
else
  W = sprand(N,N,den);
  W = W - spdiags(diag(W),0,N,N);
end

Win = zeros(N,1);
Win(1:ceil(N*.2)) = 1;

figure(1)
imagesc(W)
colormap(1-gray)
colorbar
xlim([0.5 N+0.5])
ylim([0.5 N+0.5])
axis equal
set(gca,'xtick',get(gca,'ytick'))
set(gca,'xaxislocation','top','tickdir','out')
grid

aE = (2*tauE-dt)/(2*tauE+dt);  bE = 2/(2*tauE+dt);

tref = 3; Vthr = -50; Vres = -70;

gE = zeros(N,1); V = VL*ones(N,1); t = zeros(Nt,1);

sp = zeros(N,1); T = -10*ones(N,1); ref = zeros(N,1);

a = 2*Cm/dt; b = 2*gL*VL; c = a + gL;

figure(2)

for j=1:Nt-1,

    t = j*dt;

    spin = (t/P==round(t/P));

    gEo = gE;

    gE = aE*gE + bE*(W*sp + spin*Win);

    V = ( (a - (gL + gEo)).*V + b ) ./  (c + gE);

    ref = find( t-T < tref );
    V(ref) = Vres;

    sp = (V>Vthr);

    if sum(sp)>0
        fsp = find(sp>0);
        T(fsp) = t;
        V(fsp) = Vres;
        plot(t,fsp,mk)
        hold on
    end

end  % for j

hold off
ylim([0.5 N+0.5]);
set(gca,'ytick',[2:2:N],'tickdir','out')
box off
xlabel('t  (ms)','fontsize',14)
ylabel('cell','fontsize',14)

%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% threecellI.m
%
% simulate the three cell integrate and fire net
%
%  usage:    threecellI(dt,Tfin,P)
%
%  where:    dt = timestep (ms)
%          Tfin = final time (ms)
%             P = period of stimulus (ms)
%
%  e.g.:   threecellI(0.01,50,2)
%

function threecellI(dt,Tfin,P)

Nt = ceil(Tfin/dt);

tauE = 2; VsynE = 0;
tauI = 2; VsynI = -70;

gL = 0.3; VL = -68; Cm = 1;

Win = [0.5 0 0]';

W = [0 0 0; 0.5 0 0; 0.5 0.5 0];

Winh = 3;

aE = (2*tauE-dt)/(2*tauE+dt);  bE = 2/(2*tauE+dt);

aI = (2*tauI-dt)/(2*tauI+dt);  bI = 2/(2*tauI+dt);

tref = 3;
Vthr = -50;
Vres = -70;

gE = zeros(3,Nt); 
gI = zeros(1,Nt); 
V = VL*ones(3,Nt);
t = zeros(Nt,1);

sp = zeros(3,1); T = -10*ones(3,1); ref = zeros(3,1);

for j=1:Nt-1,

    t(j+1) = j*dt;

    spin = (t(j+1)/P==round(t(j+1)/P));

    gE(:,j+1) = aE*gE(:,j) + bE*(W*sp + spin*Win);

    gI(j+1) = aI*gI(j) + bI*Winh*sp(3);

    top = (2*Cm/dt - (gL+gE(:,j))).*V(:,j) + 2*gL*VL;
    top(1) = top(1) - gI(j)*V(1,j) + (gI(j)+gI(j+1))*VsynI;
    bot = 2*Cm/dt + gL + gE(:,j+1);
    bot(1) = bot(1) + gI(j+1);
    V(:,j+1) = top./bot;

    ref = find( t(j+1)-T < tref );
    V(ref,j+1) = Vres;

    sp = (V(:,j+1)>Vthr);

    if sum(sp)>0
        fsp = find(sp>0);
        T(fsp) = t(j+1);
        V(fsp,j+1) = Vres;
    end


end  % for j

for n=1:3
    subplot(4,1,n)
    [ax,h1,h2] = plotyy(t,gE(n,:),t,V(n,:));
    set(ax,'xtick',[],'tickdir','out')
    set(h1,'color','k')
    set(h2,'linestyle','-','color','r')
    set(ax(1),'ytick',[0 .2 .4])
    set(ax(2),'ytick',[-70:10:-50])
    set(ax(1),'ycolor','k')
    set(ax(2),'ycolor','r')
    box off
    ylabel(ax(1),['g_{E,' num2str(n) '}'],'fontsize',14)
    ylabel(ax(2),['V_' num2str(n)],'fontsize',14)
end

subplot(4,1,4)
plot(t,gI,'k')
set(gca,'tickdir','out')
box off
ylabel('g_{I,1}','fontsize',14)

%set(ax(1),'ytick',[0:5]/10)
%set(ax(2),'ytick',[-70:5:-50])
xlabel('t (ms)','fontsize',14)

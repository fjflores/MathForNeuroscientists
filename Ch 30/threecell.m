%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% threecell.m
%
% simulate the three cell feed forward integrate and fire net
%
%  usage:    threecell(dt,Tfin,P)
%
%  where:    dt = timestep (ms)
%          Tfin = final time (ms)
%             P = period of stimulus (ms)
%
%  e.g.:   threecell(0.01,50,5)
%

function threecell(dt,Tfin,P)

Nt = ceil(Tfin/dt);

tauE = 2; VE = 0;

gL = 0.3; VL = -68; Cm = 1;

Win = [0.5 0 0]';

w = 0.5;   %0.26;
W = [0 0 0; 0.5 0 0; w w 0];

aE = (2*tauE-dt)/(2*tauE+dt);  bE = 2/(2*tauE+dt);

tref = 3;
Vthr = -50;
Vres = -70;

gE = zeros(3,Nt); 
V = VL*ones(3,Nt);
t = zeros(Nt,1);

sp = zeros(3,1); T = -10*ones(3,1); ref = zeros(3,1);

for j=1:Nt-1,

    t(j+1) = j*dt;

    spin = (t(j+1)/P==round(t(j+1)/P));

    gE(:,j+1) = aE*gE(:,j) + bE*(W*sp + spin*Win);

    top = (2*Cm/dt - (gL+gE(:,j))).*V(:,j) + 2*gL*VL;
    bot = 2*Cm/dt + gL + gE(:,j+1);
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
    subplot(3,1,n)
    [ax,h1,h2] = plotyy(t,gE(n,:),t,V(n,:));
    if n<3
       set(ax,'xtick',[])
    end
    set(h2,'linestyle','-','color','r')
    set(h1,'color','k')
    ylim(ax(1),[0 0.4])
    set(ax(1),'ytick',[0:4]/10)
    set(ax(2),'ytick',[-70:5:-50])
    set(ax(1),'ycolor','k')
    set(ax(2),'ycolor','r')
    set(ax,'tickdir','out')
    box off
    ylabel(ax(1),['g_{E,' num2str(n) '}'],'fontsize',14)
    ylabel(ax(2),['V_' num2str(n)],'fontsize',14)
end

%set(ax(1),'ytick',[0:5]/10)
%set(ax(2),'ytick',[-70:5:-50])
xlabel('t (ms)','fontsize',14)

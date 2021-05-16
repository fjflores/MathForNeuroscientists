%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% threecellIrast.m
%
% simulate the three cell integrate and fire net
%
%  usage:    threecellI(dt,Tfin,P)
%
%  where:    dt = timestep (ms)
%          Tfin = final time (ms)
%             P = period of stimulus (ms)
%
%  e.g.:   threecellIrast(0.01,75,2,'ro')
%

function threecellIrast(dt,Tfin,P,mk)

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

gE = zeros(3,1); 
gI = 0;
V = VL*ones(3,1);

sp = zeros(3,1); T = -10*ones(3,1); ref = zeros(3,1);

a = 2*Cm/dt; b = 2*gL*VL; c = a + gL;

for j=1:Nt-1,

    t = j*dt;

    spin = (t/P==round(t/P));

    gEo = gE;
    gE = aE*gE + bE*(W*sp + spin*Win);

    gIo = gI;
    gI = aI*gI + bI*Winh*sp(3);

    top = (a - (gL + gEo)).*V + b;
    top(1) = top(1) - gIo*V(1) + (gIo+gI)*VsynI;
    bot = c + gE;
    bot(1) = bot(1) + gI;
    V = top./bot;

    ref = find( t-T < tref );
    V(ref) = Vres;

    sp = (V>Vthr);

    if sum(sp)>0
        fsp = find(sp>0);
        T(fsp) = t;
        V(fsp) = Vres;
        plot(t,fsp*Tfin/10,mk)
        hold on
    end


end  % for j

hold off
ylim([0.5 3*Tfin/10 + 5]);
set(gca,'ytick',[1:3]*Tfin/10)
set(gca,'yticklabel',[1:3])
set(gca,'tickdir','out')
box off
xlabel('t  (ms)','fontsize',14)
ylabel('cell','fontsize',14)


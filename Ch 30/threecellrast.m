%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% threecellrast.m
%
% simulate the three cell feed forward integrate and fire net
%
%  usage:    threecellrast(dt,Tfin,P,mk)
%
%  where:    dt = timestep (ms)
%          Tfin = final time (ms)
%             P = period of stimulus (ms)
%            mk = plot marker
%
%  e.g.:   threecellrast(0.01,50,2,'k+')
%

function threecellrast(dt,Tfin,P,mk)

Nt = ceil(Tfin/dt);

tauE = 2; VE = 0;

gL = 0.3; VL = -68; Cm = 1;

Win = [0.5 0 0]';

W = [0 0 0; 0.5 0 0; 0.5 0.5 0];

aE = (2*tauE-dt)/(2*tauE+dt); bE = 2/(2*tauE+dt);

tref = 3; Vthr = -50; Vres = -70;

gE = zeros(3,1); V = VL*ones(3,1);

sp = zeros(3,1); T = -10*ones(3,1); ref = zeros(3,1);

a = 2*Cm/dt; b = 2*gL*VL; c = a + gL;

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
        plot(t,fsp*Tfin/10,mk)
        hold on
    end

end  % for j

hold off
%axis equal
%axis tight
ylim([0.5 3*Tfin/10 + 5]);
set(gca,'ytick',[1:3]*Tfin/10)
set(gca,'yticklabel',[1:3])
set(gca,'tickdir','out')
box off
xlabel('t  (ms)','fontsize',14)
ylabel('cell','fontsize',14)

%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% EInetH.m
%
% simulate a net of E and I cells and Hebbian Plasticity on WEE 
%
%  usage:    EInet(N,dt,Tfin,P,den)
%
%  where:    N = struct('E',# of E cells,'I',# of I cells)
%           dt = j (ms)
%         Tfin = final time (ms)
%            P = period of stimulus (ms)
%          den = struct('EE',den of EE cnxs,'EI',lw,'II',lw,'IE') 
%  e.g.:   
%          N = struct('E',80,'I',20);
%          den = struct('EE',.2,'EI',.25,'IE',.25,'II',.05);
%          EInetH(N,.01,100,30,den)
%
%  	EInetH(N,.01,2000,50,den,'EInet8020');
%

function [W,Waft] = EInetH(N,dt,Tfin,P,den,fname)

close all

Nt = ceil(Tfin/dt);

tauE = 2; VsynE = 0;
tauI = 1; VsynI = -70;

gL = 0.3; VL = -68; Cm = 1;

if nargin == 6

   load(fname)
   WEE = W(1:N.E,1:N.E);
   WII = W(N.E+1:end,N.E+1:end);
   WIE = W(1:N.E,N.E+1:end);
   WEI = W(N.E+1:end,1:N.E);

else

   WEE = sprand(N.E,N.E,den.EE);
   WEE = (WEE - spdiags(diag(WEE),0,N.E,N.E))/5;

   WEI = sprand(N.I,N.E,den.EI)/4;

   WIE = sprand(N.E,N.I,den.IE)/4;

   WII = sprand(N.I,N.I,den.II);
   WII = (WII - spdiags(diag(WII),0,N.I,N.I))/4;

   W = [WEE WIE; WEI WII];

end

WinEE = zeros(N.E,1);                        % drive 1st 20% of E cells
WinEE(1:ceil(N.E*.2)) = 1;

Wmax = 0.2;

tauP = 5; tauD =  5;
AP =  0.1; AD = 0.3;

pre = cell(N.E,1); post = cell(N.E,1);

for i=1:N.E,
    pre{i} = find(WEE(i,:));
    post{i} = find(WEE(:,i));
end

%W = [WEE WIE; WEI WII];
%Wbef = W;
Ntot = N.E+N.I;
figure(1)
imagesc(W)
colormap(1-gray)
colorbar
hold on
plot([0 Ntot+0.5],[N.E+0.5 N.E+0.5],'r')
plot([N.E+0.5 N.E+0.5],[0 Ntot+0.5],'r')
axis square
set(gca,'xtick',get(gca,'ytick'))
set(gca,'xaxislocation','top')
grid
set(gca,'tickdir','out')
hold off

aE = (2*tauE-dt)/(2*tauE+dt); bE = 2/(2*tauE+dt);

aI = (2*tauI-dt)/(2*tauI+dt); bI = 2/(2*tauI+dt);

tref = 3; Vthr = -50; Vres = -70;

gEE = zeros(N.E,1); VE = VL*ones(N.E,1); t = zeros(Nt,1);

gII = zeros(N.I,1); VI = VL*ones(N.I,1); 

gEI = zeros(N.I,1); 

gIE = zeros(N.E,1); 

spE = zeros(N.E,1); TE = -10*ones(N.E,1); 

spI = zeros(N.I,1); TI = -10*ones(N.I,1); 

a = 2*Cm/dt; b = 2*gL*VL; c = a + gL;

for j=1:Nt-1,

    t = j*dt;

    spin = (t/P==round(t/P));

    gEEo = gEE;
    gEE = aE*gEE + bE*(WEE*spE + spin*WinEE);

    gEIo = gEI;
    gEI = aE*gEI + bE*WEI*spE;

    gIEo = gIE;
    gIE = aI*gIE + bI*WIE*spI; 

    gIIo = gII;
    gII = aI*gII + bI*WII*spI;

    VE = ( (a - (gL + gEEo + gIEo)).*VE + b + (gIE+gIEo)*VsynI ) ./  (c + gEE + gIE);

    VI = ( (a - (gL + gIIo + gEIo)).*VI + b + (gII+gIIo)*VsynI ) ./  (c + gII + gEI);

    refE = find( t-TE < tref );
    VE(refE) = Vres;

    spE = (VE>Vthr);

    if sum(spE)>0
        fspE = find(spE>0);
        TE(fspE) = t;
        VE(fspE) = Vres;
        figure(2)
        plot(t,fspE,'k.')
        hold on

        nEs = length(fspE);

        for n = 1:nEs,

            i = fspE(n);     % cell i has spiked

            ipre = pre{i};

            if ~isempty(ipre)    % potentiate presyn weights

                Deltat = TE(ipre) - t;
                dWP = AP*exp(Deltat/tauP)';
                WEE(i,ipre) = Wmax*dWP + WEE(i,ipre).*(1 - dWP);

            end

            ipost = post{i};

            if isempty(ipost)
               continue
            end

            Deltat = TE(ipost)-t;       % depress postsyn weights
            dWD = AD*exp(Deltat/tauD);
            WEE(ipost,i) = WEE(ipost,i).*(1 - dWD);

        end

    end

    winc = 600*10;

    if mod(j,winc) == 0

       figure(3)
       [iw,jw,sw] = find(WEE);
       nbin = 20;
       hw = hist(sw,nbin);
       wax = linspace(0,Wmax,nbin);
       plot3(wax,t*ones(nbin,1),hw/length(sw),'k')
       hold on

    end

    refI = find( t-TI < tref );
    VI(refI) = Vres;

    spI = (VI>Vthr);

    if sum(spI)>0
        fspI = find(spI>0);
        TI(fspI) = t;
        VI(fspI) = Vres;
        figure(2)
        plot(t,N.E+fspI,'r.')
        hold on
    end

end % for j

Waft = [WEE WIE; WEI WII];

figure(2)
hold off
box off
xlim([10 590])
ylim([0.5 Ntot+0.5]);
set(gca,'tickdir','out')
xlabel('t  (ms)','fontsize',14)
ylabel('cell','fontsize',14)

figure(3)
hold off
view(27.5,42)
xlabel('weight','fontsize',14)
ylabel('time','fontsize',14)
zlabel('fraction','fontsize',14)

figure(4)
imagesc(Waft)
colormap(1-gray)
colorbar
hold on
plot([0 Ntot+0.5],[N.E+0.5 N.E+0.5],'r')
plot([N.E+0.5 N.E+0.5],[0 Ntot+0.5],'r')
axis square
set(gca,'xtick',get(gca,'ytick'))
set(gca,'xaxislocation','top')
grid
set(gca,'tickdir','out')
hold off


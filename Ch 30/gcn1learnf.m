%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% gcn1learnf.m
%
% simulate the learning phase of a grid cell net of E and I cells 
% following Widloski and Fiete
%
% usage:    gcn1learnf
%

function [WII,WIER,WIEL,WERI,WELI,xE,xI] = gcn1learnf

close all

T_long = 14400;     %total sim time (s)
T = 10;             %sim time blocks (s)
dt = 1/2000;       %step size of numerical integration (s)

N = struct('E',200,'I',80);   % numbers of cells

w0 = 0.001;
WII = w0*rand(N.I); %(2*rand(N.I)-1);      % initialize weights
WERI = w0*rand(N.E,N.I); % (2*rand(N.E,N.I)-1);
WELI = w0*rand(N.E,N.I); % (2*rand(N.E,N.I)-1);
WIER = w0*rand(N.I,N.E); % (2*rand(N.I,N.E)-1);
WIEL = w0*rand(N.I,N.E); %(2*rand(N.I,N.E)-1);

eta = 0.015; % +0.005;   % learning rate

TER = -100*ones(N.E,1);    % most recent spike time
TEL = -100*ones(N.E,1); 
TI = -100*ones(N.I,1); 

xE = (1:N.E)'/N.E-1/2/N.E;  %  m      
xI = (1:N.I)'/N.I-1/2/N.I;  %  m

AE = ones(size(xE));
ind = find(abs(xE-.50)>.28);
AE(ind) = exp(-60*((abs(xE(ind)-.50)-.28)/.72).^2);

AI = ones(size(xI));
ind = find(abs(xI-.50)>.28);
AI(ind) = exp(-60*((abs(xI(ind)-.50)-.28)/.72).^2);

FER = zeros(N.E,1);
FEL = zeros(N.E,1);
FI = zeros(N.I,1);

thrER = exprnd(ones(N.E,1));
thrEL = exprnd(ones(N.E,1));
thrI = exprnd(ones(N.I,1));

%STDP parameters
tau_STDP = 0.006;   %STDP time constant

%Sinusoidal animal trajectory for training
x = (sin((dt:dt:T)*2*pi/(T/2)-pi)+1)/2;
v = diff(x)/dt; v = [v(1),v];
s2 = (0.02)^2;

tim = dt:dt:T;
NT = length(tim);
% tim1 = tim(1:NT/4);
% xcv(1:NT/4) = abs(tim1-1.25)/2/1.25;
% tim2 = tim(NT/4+1:NT/2);
% xcv(NT/4+1:NT/2) = 1 - abs(tim2-3.75)/2/1.25;
% xcv(NT/2+1:NT) = xcv(1:NT/2);
% x = xcv;
% v = diff(x)/dt; v = [v(1),v];

FEf = zeros(N.E,length(x));    % precompute exponentials
for j=1:N.E,
    FEf(j,:) = 10*AE(j)*exp(-(x-xE(j)).^2/s2);
end
FIf = zeros(N.I,length(x));
for j=1:N.I,
    FIf(j,:) = 50*AI(j)*exp(-(x-xI(j)).^2/s2);
end

for tt = 1:T_long/T %OUTER time loop (T_long sec)
    
    for t = 1:T/dt %INNER time loop (T sec)
        
        realt = T*(tt-1) + t*dt;
        
        FERnf = (1+0.9*v(t))*FEf(:,t);
        FELnf = (1-0.9*v(t))*FEf(:,t);
        
        FER = FER + dt*FERnf;              % cumm spike rates
        FEL = FEL + dt*FELnf;     
        FI = FI + dt*FIf(:,t);
        
        spER = find(FER > thrER);    % spiking ER cells
        TER(spER) = realt;
        spEL = find(FEL > thrEL);
        TEL(spEL) = realt;
        spI = find(FI > thrI);
        TI(spI) = realt;
        
        if sum(spER)>0    % process ER spikes
         
            FER(spER) = 0;    % reset cumm spike rate
            nERs = length(spER);
            thrER(spER) = exprnd(ones(nERs,1));   % reset thresholds
            
            for n = 1:nERs,
                
                k = spER(n);     % cell k has spiked
                Deltat = TI - realt;
                dWP = 2*eta*exp(Deltat/(2*tau_STDP));
                WERI(k,:) = max(WERI(k,:) - dWP',0);  % potentiate presyn weights
                dWD = eta*exp(Deltat/(3*tau_STDP));
                WIER(:,k) = max(WIER(:,k) - dWD,0);  % depress postsyn weights
                
            end
            
        end
        
        if sum(spEL)>0    % adjust weights of spiking EL cells
            
            FEL(spEL) = 0;
            nELs = length(spEL);
            thrEL(spEL) = exprnd(ones(nELs,1));
            
            for n = 1:nELs,
                
                k = spEL(n);     % cell k has spiked
                Deltat = TI - realt;
                dWP = 2*eta*exp(Deltat/(2*tau_STDP));
                WELI(k,:) = max(WELI(k,:) - dWP',0);  % potentiate presyn weights
                dWD = eta*exp(Deltat/(3*tau_STDP));
                WIEL(:,k) = max(WIEL(:,k) - dWD,0);  % depress postsyn weights
                
            end
            
        end
        
        if sum(spI)>0
            
            FI(spI) = 0;
            nIs = length(spI);
            thrI(spI) = exprnd(ones(nIs,1));

            for n = 1:nIs,
                
                k = spI(n);     % cell k has spiked
                
                Deltat = TER - realt;
                ed4 = eta*exp(Deltat/(4*tau_STDP));
                dWP = 1.2*ed4; % 1.2*eta*exp(Deltat/(4*tau_STDP));
                WIER(k,:) = WIER(k,:) + dWP';  % potentiate presyn weights
                dWD = ed4; %eta*exp(Deltat/(4*tau_STDP));
                WERI(:,k) = WERI(:,k) + dWD;  % depress postsyn weights
                
                Deltat = TEL - realt;
                ed4 = eta*exp(Deltat/(4*tau_STDP));
                dWP = 1.2*ed4; % 1.2*eta*exp(Deltat/(4*tau_STDP));
                WIEL(k,:) = WIEL(k,:) + dWP';  % potentiate presyn weights
                dWD = ed4; % eta*exp(Deltat/(4*tau_STDP));
                WELI(:,k) = WELI(:,k) + dWD;  % depress postsyn weights
                
                Deltat = TI - realt;
                dWP = 3.5*eta*exp(Deltat/(4*tau_STDP)); % swapped 2 and 4 on May 9 and 7 and 3.5
                WII(k,:) = WII(k,:) + dWP';  % potentiate presyn weights
                dWD = 7*eta*exp(Deltat/(2*tau_STDP));
                WII(:,k) = max(WII(:,k) - dWD,0);  % depress postsyn weight
                
            end
            
        end
        
    end  % t
    
    disp(['elapsed time: ', num2str(realt),' sec'])
    
    %if mod(tt,60) == 0
    if mod(realt/60-30,70)==0
        figure(1)
        plot3((1:N.I)',(tt/6)*ones(N.I,1),WIEL(:,100),'k','linewidth',1)
        hold on
        plot3((1:N.I)',(tt/6)*ones(N.I,1),WIER(:,100),'r','linewidth',1)
 
        figure(2)
        plot3((1:N.E)',(tt/6)*ones(N.E,1),WELI(:,40),'k','linewidth',1)
        hold on
        plot3((1:N.E)',(tt/6)*ones(N.E,1),WERI(:,40),'r','linewidth',1)
       
        figure(3)
        plot3((1:N.I)',(tt/6)*ones(N.I,1),WII(:,N.I/2),'k','linewidth',1)
        hold on
        
        drawnow
    end
    
end % tt

figure(1)
legend('W_{IE_L}','W_{IE_R}')
ylabel('time (minutes)','fontsize',14)
xlabel('Inhibitory neuron','fontsize',14)
zlabel('synaptic weight','fontsize',14)
set(gca,'ytick',[30 100 170 240])
ylim([30 240])
hold off

figure(2)
legend('W_{E_LI}','W_{E_RI}')
ylabel('time (minutes)','fontsize',14)
xlabel('Excitatory neuron','fontsize',14)
zlabel('synaptic weight','fontsize',14)
set(gca,'ytick',[30 100 170 240])
ylim([30 240])
hold off

figure(3)
ylabel('time (minutes)','fontsize',14)
xlabel('Inhibitory neuron','fontsize',14)
zlabel('synaptic weight, W_{II}','fontsize',14)
set(gca,'ytick',[30 100 170 240])
ylim([30 240])
hold off

return




%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% gcn1test.m
%
% simulate the test phase of a grid cell net of E and I cells 
% following Widloski and Fiete
%
%  usage:    gcn1test
%
%  requires:   gcnrun.mat    good memories
%              fastsmooth.m  from 
%  https://www.mathworks.com/matlabcentral/fileexchange/19998-fast-smoothing-function
%

%function [FERn,FELn,FIn] = gcn1test(WII,WIER,WIEL,WERI,WELI)

load gcnrun 

N = struct('E',200,'I',80);
T_long = 60*60*1; %14400;     %total sim time (s)
T = 10;             %sim time blocks (s)
dt = 1/2000;       %step size of numerical integration (s)

xE = (1:N.E)'/N.E-1/2/N.E;  %  m      
xI = (1:N.I)'/N.I-1/2/N.I;  %  m

taus = 30e-3;  % s,       terms for updating synaptic "activations"
a = (2*taus - dt)/(2*taus + dt);
b = 2*taus/(2*taus + dt);

FER = zeros(N.E,1);  %gER = FER;     % spike rates, synaptic activations
FEL = zeros(N.E,1);  %gEL = FEL;
FI = zeros(N.I,1);   %gI = FI;

yERI = zeros(N.E,1);  yELI = yERI;       % synaptic currents
yIER = zeros(N.I,1);  yIEL = yIER; yII = yIER;   

thrER = exprnd(ones(N.E,1));    % thresholds for inhomogeneous Poisson
thrEL = exprnd(ones(N.E,1));
thrI = exprnd(ones(N.I,1));

%Sinusoidal animal trajectory for training
tim = dt:dt:T;
NT = length(tim);
x = (sin((dt:dt:T)*2*pi/(T/2)-pi)+1)/2;
v = diff(x)/dt; v = [v(1),v];
s2 = 0.02^2;

xpath = x(2500:7500);
nx = hist(xpath,100);
timeperbin = nx*dt;
xEhist = 0.005:.01:0.995;
ssEav = zeros(1,100);
ssIav = ssEav;

FEx = zeros(N.E,length(x));
for n=1:N.E,
    FEx(n,:) = 10*exp(-(x-xE(n)).^2/s2);
end 

FIx = zeros(N.I,length(x));
for n=1:N.I,
    FIx(n,:) = 50*exp(-(x-xI(n)).^2/s2);
end

gridE = zeros(T_long/T/6,100);
gridI = gridE;

for tt = 1:T_long/T %OUTER time loop (T_long sec)
    
   spkrecER = zeros(N.E,length(x));
   spkrecEL = zeros(N.E,length(x));
   spkrecI = zeros(N.I,length(x));
    
    for t = 1:T/dt %INNER time loop (T sec)
             
        realt = T*(tt-1) + t*dt;
        
        FERn = max(15 + (1+0.9*v(t))*(50 + FEx(:,t) - yERI),0); % instantaneous spike rate
        FER = FER + dt*FERn;           % cumm spike rate
        
        FELn = max(15 + (1-0.9*v(t))*(50 + FEx(:,t) - yELI),0);
        FEL = FEL + dt*FELn;
        
        FIn = max(15 + 50 + FIx(:,t) + yIER + yIEL - yII,0);
        FI = FI + dt*FIn;
        
        spER = find(FER > thrER);   % check for ER spikes
        nERs = length(spER);        % number of spiking ER cells
        yIER = a*yIER;
        
        if nERs > 0    % spiking ER cells
            
            if abs(t-5000) <= 2500 %|| abs(t-15000) < 2500,
               spkrecER(spER,t) = x(t);
%                if mod(t,50)==0
%                    figure(1)
%                    plot3(((1/200:1/200:1) -1/400)',(t-2500)*ones(200,1),FERn)
%                    hold on
%                end
            end

            FER(spER) = 0;                        % reset cumm spike rate
            thrER(spER) = exprnd(ones(nERs,1));   % reset threshold
            yIER = yIER + b*sum(WIER(:,spER),2);

        end
        
        spEL = find(FEL > thrEL);
        nELs = length(spEL);
        yIEL = a*yIEL;
        
        if nELs > 0    %  spiking EL cells
            
            spkrecEL(spEL,t) = x(t);
            
            FEL(spEL) = 0;
            thrEL(spEL) = exprnd(ones(nELs,1));
            yIEL = yIEL + b*sum(WIEL(:,spEL),2);
            
        end
        
        spI = find(FI > thrI);
        nIs = length(spI);
        yERI = a*yERI;
        yELI = a*yELI;
        yII = a*yII;
        
        if nIs > 0    %  spiking I cells
            
            if abs(t-5000)<2500 %|| abs(t-15000)<2500
              spkrecI(spI,t) = x(t);
            end
            
            FI(spI) = 0;
            thrI(spI) = exprnd(ones(nIs,1));
            yERI = yERI + b*sum(WERI(:,spI),2);
            yELI = yELI + b*sum(WELI(:,spI),2);
            yII = yII + b*sum(WII(:,spI),2);
            
        end
        
    end % t
    
%     drawnow
%     view(15.5,84)
%     %keyboard
%     hold off
    
    disp(['elapsed time: ', num2str(realt),' sec'])
    
    if mod(tt,6) == 0 
        figure(1)
        plot(FELn,'k','linewidth',1.5)
        hold on
        plot(FERn,'r','linewidth',1.5)
        hold off
        legend('FELn','FERn')
        xlabel('neuron','fontsize',14)
        ylabel('spks/s','fontsize',14)
        
        figure(2)
        plot(FIn,'k','linewidth',1.5)
        legend('FIn')
        xlabel('neuron','fontsize',14)
        ylabel('spk/s','fontsize',14)
        
        figure(3)
        ssE = spkrecER(80,:);
        ssE2 = ssE(ssE>0);
        hssE2 = hist(ssE2,100);
        nhssE2 = hssE2./timeperbin;
        nhssEs = fastsmooth(nhssE2,5,1,1);
        gridE(tt/6,:) = nhssEs;
%         ssEav = ssEav + nhssE2s;
        plot(xEhist,nhssEs)
        hold on
        
        ssE = spkrecI(40,:);
        ssE2 = ssE(ssE>0);
        hssE2 = hist(ssE2,100);
        nhssE2 = hssE2./timeperbin;
        nhssIs = fastsmooth(nhssE2,5,1,1);
        gridI(tt/6,:) = nhssIs;
%         ssIav = ssIav + nhssE2s;
        plot(xEhist,nhssIs,'g')
        
        hold off
        
        drawnow
   end
    
end % tt


figure(4)
imagesc(gridE)
h=colorbar;
colormap(1-gray)
xlabel('x  (cm)','fontsize',14)
ylabel('t  (minute)','fontsize',14)
ylabel(h,'spike rate density  (Hz/cm)','fontsize',14)
title('Cell E_R 80','fontsize',14)

figure(5)
imagesc(gridI)
h=colorbar;
colormap(1-gray)
xlabel('x  (cm)','fontsize',14)
ylabel('t  (minute)','fontsize',14)
ylabel(h,'spike rate density  (Hz/cm)','fontsize',14)
title('Cell I 40','fontsize',14)

figure(6)
plot(1:100,gridE(45,:),'r','linewidth',1)
hold on
plot(1:100,gridI(45,:),'k','linewidth',1)
hold off
xlabel('x  (cm)','fontsize',14)
ylabel('spike rate density  (Hz/cm)','fontsize',14)
legend('E_R 80','I 40')
box off

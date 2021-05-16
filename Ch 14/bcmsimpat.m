%
% Gabbiani and Cox
%
% bcmsimpat.m    from Yeung, Shouval, Blais and Cooper, PNAS 101(41)
%
%  presentation of patterned input 
%  the two parameters are firing rates in spikes/second
%
%   usage:   bigw = bcmsimpat(ratelo,ratehi);
%
%   e.g.  bigw = bcmsimpat(10,40);
%

function bigw = bcmsimpat(ratelo,ratehi)

N = struct('E',100,'I',20);   % number of cells in each population
dt = 1;   % timestep, ms
Tminutes = 150;   % duration of simulation

T = -100;           % most recent output spike time
TE = -100*ones(N.E,1);    % most recent E-presynaptic spike time

FE = zeros(N.E,1);   % cumulative E cell spike rate
gE = 0*rand(N.E,1); 
GE = sum(gE);    % total excitatory conductance
w = 0.5*ones(N.E,1);    % synaptic E weights
bigw = zeros(N.E,Tminutes);   % for reporting purposes
c = 0.2*ones(N.E,1);  % calcium levels at E synapses
FI = zeros(N.I,1); % cummulative I cell spike rate
gI = 0*rand(N.I,1); 
GI = sum(gI);   % total inhibitory conductance

thrE = exprnd(ones(N.E,1));   % Poisson E thresholds
thrI = exprnd(ones(N.I,1));   % Poisson I thresholds

Erate(:,1) = [ratehi*ones(25,1); ratelo*ones(75,1)];
Erate(:,2) = circshift(Erate(:,1),25);   % spike rate of E cells, Hz
Erate(:,3) = circshift(Erate(:,1),50); 
Erate(:,4) = circshift(Erate(:,1),75);

Irate = 10;   % spike rate of I cells, Hz
out = 0;   % number of output spikes

Vreset = -65;   % Integrate and Fire reset potential
Vm = -65;       % membrane potential of output cell
gNMDA = 1.5e-3;  % gt
t = 0;

%stim = ceil(rand*4);  % choose a stimulus pattern
stim = 2;

figure(1)
plot([1:100]',Erate(:,2),'k','linewidth',1.5)
xlabel('E synapse','fontsize',14)
ylabel('presynaptic spike rate  (Hz)','fontsize',14)
ylim([0 50])
box off
text(10,40,'(A)','fontsize',20)

figure(2)
plot3([1:100]',0*ones(100,1),w,'k')
hold on

for m=1:Tminutes,
    for j=1:60e3,     % 60 thousand time steps per minute
        
        %if mod(j,1000) == 0    % switch stimulus every 500 ms
        %    stim = ceil(rand*4);
        %end
        stim = 2;
        
        t = t+dt;
        
        FE = FE + 1e-3*dt*Erate(:,stim);
        GE = GE - dt*GE/5;
        spE = find(FE > thrE);    % presyn E spikes
        
        if sum(spE) > 0    % process E spikes
            
            TE(spE) = t;    % record time
            FE(spE) = 0;    % reset cumulative spike rate
            nEs = length(spE);
            thrE(spE) = exprnd(ones(nEs,1));   % reset thresholds
            GE = GE + 0.03*sum(w(spE));   % augment E conductance
            
        end
        
        tiltE = t - TE;
        A = 0.7*exp(-tiltE/50) + 0.3*exp(-tiltE/200);
        
        FI = FI + 1e-3*dt*Irate;     
        GI = GI - dt*GI/5;
        spI = find(FI > thrI);    % presyn I spikes
        
        if sum(spI) > 0     % process I spikes
            
            FI(spI) = 0;    % reset cumm spike rate
            nIs = length(spI);
            thrI(spI) = exprnd(ones(nIs,1));   % reset thresholds
            GI = GI + 0.1*nIs;
            
        end
        
        tiltm = t - T;
        B = 42*(0.75*exp(-tiltm/3) + 0.25*exp(-tiltm/35));
        
        gNMDA = gNMDA + dt*( (8e-7)*((5e-3)-gNMDA) - (8e-9)*B^2*gNMDA );
        
        %FG 05/03/16
        %JNMDA = gNMDA*H(Vm+B)*A;
        JNMDA = gNMDA*H(Vreset+B)*A;
        
        c = c + dt*( JNMDA - c/20);
        
        w = max(w + dt*( (2e-5)*c.*(Omega(c)-w/200) ),0);
        
        Vreset = Vreset - dt*(Vreset+65)/100;
        
        %Vm = Vm + dt*( Vreset-Vm - Vm*sum(gE) - (Vm+65)*sum(gI) )/20;
        num = Vm + dt*(Vreset-65*GI)/20;
        den = 1 + dt*(1+GE+GI)/20;
        Vm = num/den;   % via backward Euler
        
        if Vm > -55
            Vreset = Vreset - 2;
            Vm = Vreset;
            T = t;
            out = out + 1;
        end
        
    end % j loop (one minute)
    
    if mod(m,5) == 0
        figure(2);
        plot3([1:100]',m*ones(100,1),w,'k')
        drawnow
    end
    bigw(:,m) = w;
    
end   % m loop (over minutes)

xlabel('E synapse','fontsize',14)
ylabel('t  (minutes)','fontsize',14)
zlabel('weight, w','fontsize',14)
text(80,10,8,'(B)')
hold off

return

function val = Omega(c)  % "U" shaped calcium function
val = 1./(1+exp(40*(0.4-c))) - 0.5./(1+exp(60*(0.25-c))); % code

function val = H(V)   % Jahr Steven Function
val = (130-V)/(1+exp(-0.062*V)/3.57);


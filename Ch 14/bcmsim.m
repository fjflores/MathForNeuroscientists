%
% Gabbiani and Cox, Mathematics for Neuroscientists
%
% bcmsim.m    from Yeung, Shouval, Blais and Cooper, PNAS 101(41)
%
%  called by bcmdrive.m
%

function [c, w, outrate, mc, mw] = bcmsim(rate)

N = struct('E',100,'I',20);   % number of cells in each population
dt = 1;   % timestep, ms
Tminutes = 166;   % duration of simulation
mc = zeros(Tminutes,1);
mw = mc;

T = -100;           % most recent output spike time
TE = -100*ones(N.E,1);    % most recent E-presynaptic spike time

FE = zeros(N.E,1);   % cumulative E cell spike rate
gE = 0*rand(N.E,1); %excitatory synaptic conductance vector, initialized at 0
GE = sum(gE);    % total excitatory conductance
w = 0.5*ones(N.E,1);    % synaptic E weights
c = 0.2*ones(N.E,1);  % calcium levels at E synapses
FI = zeros(N.I,1); % cumulative I cell spike rate
gI = 0*rand(N.I,1); %inhibitory synaptic conductance vector, initialized at 0
GI = sum(gI);   % total inhibitory conductance

thrE = exprnd(ones(N.E,1));   % Poisson E thresholds
thrI = exprnd(ones(N.I,1));   % Poisson I thresholds

Erate = rate;   % spike rate of E cells, Hz
Irate = rate;   % spike rate of I cells, Hz
out = 0;   % number of output spikes

Vreset = -65;   % Integrate and Fire reset potential
Vm = -65;       % membrane potential of output cell
gt = 5e-3;
gNMDA = 1.5e-3;  % gt/2; % 4.5e-3;  % gt
t = 0;

for m=1:Tminutes,
    for j=1:60e3,     % 60 thousand time steps per minute since dt = 1 ms
        
        t = t+dt;
        
        FE = FE + 1e-3*dt*Erate;
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
        %A = 0.75*exp(-tiltE/50) + 0.25*exp(-tiltE/200);
        
        FI = FI + 1e-3*dt*Irate;     
        GI = GI - dt*GI/5;
        spI = find(FI > thrI);    % presyn I spikes
        
        if sum(spI) > 0     % process I spikes
            
            FI(spI) = 0;    % reset cumulative spike rate
            nIs = length(spI);
            thrI(spI) = exprnd(ones(nIs,1));   % reset thresholds
            GI = GI + 0.1*nIs;
            
        end
        
        tiltm = t - T;
        B = 42*(0.75*exp(-tiltm/3) + 0.25*exp(-tiltm/35));
        %B = 60*(0.75*exp(-tiltm/3) + 0.25*exp(-tiltm/30));
        
        gNMDA = gNMDA + dt*( (8e-7)*((4.5e-3)-gNMDA) - (8e-9)*B^2*gNMDA ); % gNMDA(0) = 1.5e-3
       
        %V = Vm + B;
        %FG 04/15/2016
        V = Vreset + B;
        
        JNMDA = gNMDA*H(V)*A;
        
        c = c + dt*( JNMDA - c/20);
        
        %w = w + dt*( (2e-5)*c.*(Omega(c)-w/200) );
        %w = max(w + dt*( (2e-5)*c.*(Omega(c)) ),0);
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
    
%     subplot(1,2,1)    % show c and w every minute
%     hist(c)
%     xlabel('calcium')
%     subplot(1,2,2)
%     hist(w)
%     xlabel('weight')
%     drawnow 
     
    mc(m) = mean(c);
    mw(m) = mean(w);
    
    disp(['rate = ' num2str(rate) ', minutes = ' num2str(m) ...
         ', mean(w) = ' num2str(mean(w)) ', gNMDA = ' num2str(gNMDA) ...
         ', mean(c) = ' num2str(mean(c))])
    
end   % m loop (over minutes)

outrate = (out/(Tminutes*60));

return

function val = Omega(c)  % "U" shaped calcium function
val = 1./(1+exp(40*(0.4-c))) - 0.5./(1+exp(60*(0.25-c))); % code

function val = H(V)   % Jahr Stevens Function
val = (130-V)/(1+exp(-0.062*V)/3.57);


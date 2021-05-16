%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% stEqatrain.m
%
% staggered Euler on the quasi active sphere
%
% usage:  stEqatrain(dt,Tfin,G,stim)
%     G = struct('K', 36, 'Na', 120, 'Cl', 0.3);
%     stim = struct('amp',1,'per',50,'tau1',0.6,'tau2',0.5);
%

function [t,v] = stEqatrain(dt,Tfin,G,stim,ind)

A = 4*pi*(1e-6); % (cm)^2

[Vr,B] = hhsym(G,ind);

Nt = ceil(Tfin/dt);

t = zeros(Nt,1); 
y = zeros(4,Nt);
f = zeros(4,1);

td = 1/dt;
B1 = 2*td*eye(4) - B;
B2 = 2*td*eye(4) + B;

Istim0 = 0;
pstim = t;
toff = 2;

for j=2:Nt,

      t(j) = (j-1)*dt;
      
      tch = mod(t(j)-toff,stim.per(ind));

      Istim1 = stim.amp(ind)*(exp(-tch/stim.tau1(ind))-exp(-tch/stim.tau2(ind)))*(1e-6);
      pstim(j) = Istim1;

      f(1) = (Istim0 + Istim1)/A;

      y(:,j) = B1 \ ( B2*y(:,j-1) + f );

      Istim0 = Istim1;

end

v = y(1,:)';
% figure(1)
% plot(t,v,'r')
% figure(2)
% plot(t,pstim)

return

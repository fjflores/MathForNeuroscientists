%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% stEforkthreshloc.m
%

function Itheta = stEforkthreshloc

fork = struct('rad',1e-4*[1 1 1],'ell',250e-4*[1 1 1],'dx',5e-4,'dt',0.05);
stim = struct('t1',1,'t2',2,'amp',4e-4,'loc',5e-4,'Tfin',10);
g = struct('K', 36, 'Na', 120, 'Cl', 1/15);
Vr = forkrest(fork,g);

N = fork.ell/fork.dx;
Hlen = sum(N)+1;
eloc = [1:1:N(1) N(1)+N(2)+1:1:Hlen];
deloc = [0:1:N(1)+N(3)]*fork.dx;
eloc = fliplr(eloc);
nloc = length(eloc);

for j=1:nloc,

   stim.loc = eloc(j);
   
   xloc(j) = deloc(j);

   Itheta(j) = getItheta(fork,g,stim,Vr);
   
   %disp('paused')
   %pause
    
   j

end

plot(xloc,1e6*Itheta,'k')
box off
xlabel('x_s  (cm)','fontsize',14)
ylabel('I_\theta  (pA)','fontsize',14)

%
% find the threshold current via bisection on vmax
%

function Ith = getItheta(fork,g,stim,Vr)
a = 10e-6;
b = 250e-6;
stim.amp = a;
[t,v] = stEforkwVr(fork,g,stim,0,Vr);
va = max(v);
stim.amp = b;
[t,v] = stEforkwVr(fork,g,stim,0,Vr);
vb = max(v);

Ith = (a+b)/2;

while  b-a > .25e-6   % solve to .25 pA resolution

       stim.amp = Ith;
       [t,v] = stEforkwVr(fork,g,stim,0,Vr);
       vmid = max(v);
 
       if va*vmid < 0
          b = Ith;
          vb = vmid;
       else
          a = Ith;
          va = vmid;
       end

       Ith = (a+b)/2;

end

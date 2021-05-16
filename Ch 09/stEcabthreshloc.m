%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% stEcabthreshloc.m
%
%  Compute the current amplitude, of a 1ms pulse injected at site x, 
%  needed to elicit an action potential, for each x along the cable.
%
%  For the uniform cable type     stEcabthreshloc(1);
%
%  For the nonuniform cable type  stEcabthreshloc(0);
%

function Itheta = stEcabthreshloc(uni)

cab = struct('rad',1e-4,'ell',1e-1,'dx',1e-4,'dt',0.02);
stim = struct('t1',1,'t2',2,'amp',4e-4,'loc',0.03,'Tfin',20);
if uni
   g = struct('K', 36, 'Na', 120, 'Cl', 1/15);
else
   x = cab.dx/2:cab.dx:cab.ell-cab.dx/2; hot = ((x>0.005).*(x<0.01))';
   g = struct('K', 40-20*hot, 'Na', 44+560*hot, 'Cl', 1/15);
end

Vr = cabrest(cab,g);

for j=1:99,

   stim.loc = j*0.001;
   
   xloc(j) = stim.loc;

   Itheta(j) = getItheta(cab,g,stim,Vr);
   
   j

end

plot(xloc,1e6*Itheta,'k')
box off
xlabel('x_s  (cm)','fontsize',14)
ylabel('I_\theta  (pA)','fontsize',14)

%
% find the threshold current via bisection on vmax
%
function Ith = getItheta(cab,g,stim,Vr)
a = 10e-6;
b = 500e-6;
stim.amp = a;
[t,v] = stEcabwVr(cab,g,stim,0,Vr);
va = max(v);
stim.amp = b;
[t,v] = stEcabwVr(cab,g,stim,0,Vr);
vb = max(v);

Ith = (a+b)/2;

while  b-a > .25e-6   % solve to .25 pA resolution

       stim.amp = Ith;
       [t,v] = stEcabwVr(cab,g,stim,0,Vr);
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

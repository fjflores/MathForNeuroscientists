%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% trapforksyngain
%
function trapforksyngain

stim = struct('t1',1,'tau',1,'Gsyn',5e-7,'vsyn',70,'loc',100,'Tfin',15);

ell = [2.5 2.5 2.5]/100;
dx = .0001;
N = ell/dx;
Hlen = sum(N)+1;

eloc = [10:10:N(1) N(1)+N(2)+10:10:Hlen];
deloc = [0:10:N(1)+N(3)-1]*dx;
eloc = fliplr(eloc);
nloc = length(eloc);
vsynmax = zeros(nloc,1);
vsomamax = zeros(nloc,1);

for jloc = 1:nloc

    jloc

    stim.loc = eloc(jloc);

    [t,vrec] = trapforksyn(stim,0);

    vsomamax(jloc) = max(vrec(2,:));

    vsynmax(jloc) = max(vrec(1,:));

end

deloc = deloc*10000; % cm to microns
plot(deloc,vsomamax,'k')
hold on
plot(deloc,vsynmax,'r');
hold off
box off
xlabel('distance from soma  (\mum)','fontsize',14)
ylabel('(mV)','fontsize',14)
legend('Peak Soma Potential','Peak Synaptic Potential','location','best')

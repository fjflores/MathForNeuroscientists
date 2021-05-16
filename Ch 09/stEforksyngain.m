%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
%  stEforksyngain
%
% Contrast peak soma and synaptic potentials as the synapse travels
% away from the soma
%
function stEforksyngain

stim = struct('t1',1,'tau',1,'Gsyn',5e-7,'Vsyn',0,'loc',100,'Tfin',15);

ell = [2.5 2.5 2.5]/100;
dx = .0001;
N = ell/dx;
Hlen = sum(N)+1;

eloc = [10:10:N(1) N(1)+N(2)+10:10:Hlen];
deloc = [0:10:N(1)+N(3)-1]*dx;
eloc = fliplr(eloc);
nloc = length(eloc);

g = struct('K', 36, 'Na', 120, 'Cl', 1/15);
vsomamax = zeros(nloc,1);
vsynmax = zeros(nloc,1);

for jloc = 1:nloc

    stim.loc = eloc(jloc);
    [t,vrec] = stEforksyn(stim,g,0);
    plot(t,vrec(1,:),'r',t,vrec(2,:),'k')
    drawnow
    vsomamax(jloc) = max(vrec(2,:)) - vrec(2,1);
    vsynmax(jloc) = max(vrec(1,:))-vrec(1,1);

end

deloc = deloc*10000; % cm to microns
plot(deloc,vsomamax,'k');
hold on
plot(deloc,vsynmax,'r');
hold off
legend('Peak Soma Potential','Peak Synaptic Potential','location','best')
xlabel('distance from soma (\mum)','fontsize',14)
ylabel('(mV)','fontsize',14)
box off
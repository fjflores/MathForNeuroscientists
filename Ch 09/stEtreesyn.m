
%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  stEtreesyn.m
%
%  Compute the response of an active tree to random synaptic input
%
%  usage:    [t,Vrec] = stEtreesyn(filename,stim)
%
%  e.g.: 
%       filename = 'sep12a.asc'; 'AF-6-01-04-B.asc'
%       stim = struct('num',20,'inum',5,'tau',1,'Gsyn',.75e-6,'Tfin',50,...
%              'dx',5,'dt',0.02);
%
%  Note: filename must be of either asc or swc type
%
%  Required m-files:  asc_converter.m
%                     swc_converter.m
%                     makemd.m
%                     treeplot.m
%                     makeH.m
%

function [t,Vrec,H] = stEtreesyn(filename,stim)

if filename(end-3:end) == '.asc'
    nc = asc_converter(filename);
else 
    nc = swc_converter(filename,0);
end

md = makemd(nc, stim.dx);  %neur_loadmorph_direct(nc,stim.dx);
md.filename = filename;

treeplot(md,nc,0)    % plot the tree

Ra = 0.3;                 % k Ohm cm
Cm = 1;                    % uF / cm^2
As = nc.somadata.As*1e-8;

% build H
H = makeH(md,Ra,As);
%keyboard
Nx = length(H);
dH = diag(H);

dt = stim.dt;
Nt = ceil(stim.Tfin/dt);
e = ones(Nx,1);

g = struct('K', 40*e, 'Na', 40*e, 'Cl', 1/15*e);
%g.K(end) = 20;
%g.Na(end) = 600;
%g = struct('K', 36*e, 'Na', 120*e, 'Cl', 1/15*e);
E = struct('K', -77, 'Na', 56, 'Cl', -68);
options = optimset('Jacobian','on','maxiter',10000);
V = fsolve(@(V) Iss(V,E,g,H),-70*e,options);    % initial conditions

nstim = stim.num;
Iapp = zeros(nstim,1);
eloc = Iapp;
stim_ind = ceil(rand(nstim,1)*md.count.L);      % random indices
stim.t1 = 1 + rand(nstim,1)*(90/100)*stim.Tfin;  % random times
Vexc = 0;
Vinh = -80;
stim.Vsyn = Vexc*ones(nstim,1);
stim.Vsyn(1:stim.inum) = Vinh;

cum_ind = cumsum(md.interp.size) - [1:md.count.N];
% translate indicies into branch and comp num for right dx and rad
for i=1:nstim   
    tmp = stim_ind(i) > cum_ind;
    br = 1 + sum(tmp);
    mc = length(md.interp.radii{br});
    brcomp = min(max(cum_ind(br) - stim_ind(i),1),mc);
    eloc(i) = stim_ind(i);   
    rad = md.interp.radii{br}(brcomp);
    dx = md.grid.hstep(br);
    Iapp(i) = stim.Gsyn/rad/2/pi/dx;
    x = nc.celldata{br}(:,1);
    s = brcomp/md.interp.size(br);
    xloc = (1-s)*x(1) + s*x(end);
    y = nc.celldata{br}(:,2);
    yloc = (1-s)*y(1) + s*y(end);
    z = nc.celldata{br}(:,3);
    zloc = (1-s)*z(1) + s*z(end);
    h = text(xloc,yloc,zloc,num2str(round(stim.t1(i))));
    if i<= stim.inum
        set(h,'Color','k','FontSize',8);    % label inh synapses
    else
        set(h,'Color','r','FontSize',8);    % label exc synapses
    end
end

axis off
hold off

Vrec = zeros(Nt,nstim+1);
Vrec(1,:) = V([eloc; Nx]);

n = an(V)./(an(V)+bn(V)); 
m = am(V)./(am(V)+bm(V)); 
h = ah(V)./(ah(V)+bh(V)); 

gsyn = zeros(Nx,1);

for j=2:Nt,

      t = (j-1)*dt;
      t2 = t-dt/2;
      
      gsyn(eloc) = Iapp.*((t2-stim.t1)./stim.tau).*exp(1-(t2-stim.t1)./stim.tau).*(t2>stim.t1);

      a = an(V);  c = (a+bn(V))/2;
      n = ( (1/dt-c).*n + a) ./ (1/dt + c); n4 = n.^4;

      a = am(V);  c = (a+bm(V))/2;
      m = ( (1/dt-c).*m + a) ./ (1/dt + c);

      a = ah(V);  c = (a+bh(V))/2;
      h = ( (1/dt-c).*h + a) ./ (1/dt + c); m3h = m.^3.*h;

      d = g.Na.*m3h + g.K.*n4 + g.Cl + gsyn;

      f = g.Na.*m3h*E.Na + g.K.*n4*E.K + g.Cl.*E.Cl;
      f(eloc) = f(eloc) + gsyn(eloc).*stim.Vsyn;

      H(1:Nx+1:end) = dH + d + 2*Cm/dt;         % update the diagonal

      Vmid = H\(2*Cm*V/dt + f);
      
      V = 2*Vmid - V;
  
      Vrec(j,:) = V([eloc; Nx]);

end

t = linspace(0,stim.Tfin,Nt)';

figure(2)                   % plot soma response
plot(t,Vrec(:,end),'k')
box off
xlabel('t  (ms)','fontsize',14)
ylabel('V_{soma}  (mV)','fontsize',14)

figure(3)                   % plot responses at all synapses
[ts,tsind] = sort(stim.t1);

for j=1:nstim
    plot3(t,j*ones(Nt,1),Vrec(:,tsind(j)),'k')
    hold on
end
plot3(t,0*stim.Tfin*ones(Nt,1),Vrec(:,end),'r')
hold off

return

function [val, jac] = Iss(V,E,g,B)
Nx = length(V);
a = am(V); b = bm(V);
m = a./(a+b);
dm = (dam(V).*b - a.*dbm(V))./(a+b).^2;

a = ah(V); b = bh(V);
h = a./(a+b);
dh = (dah(V).*b - a.*dbh(V))./(a+b).^2;

a = an(V); b = bn(V);
n = a./(a+b);
dn = (dan(V).*b - a.*dbn(V))./(a+b).^2;

m3h = m.^3.*h;
n4 = n.^4;

val = B*V + g.Na.*m3h.*(V-E.Na) + g.K.*n4.*(V-E.K) + g.Cl.*(V-E.Cl);

dj = g.Na.*((3*dm.*m.^2.*h + m.^3.*dh).*(V-E.Na) + m3h) + ...
     g.K.*(4*dn.*n.^3.*(V-E.K) + n4) + g.Cl;
jac = B + spdiags(dj,0,Nx,Nx);

function val = an(v)
val = .01*(10-(v+71))./(exp(1-(v+71)/10)-1);

function val = dan(v)
tmp = exp(-(61+v)/10);
val = -( tmp.*(71+v) - 10 )./(tmp-1).^2/1000;

function val = bn(v)
val = .125*exp(-(v+71)/80);

function val = dbn(v)
val = -exp(-(v+71)/80)/640;

function val = am(v)
val = .1*(25-(v+71))./(exp(2.5-(v+71)/10)-1);

function val = dam(v)
tmp = exp(-(46+v)/10);
val = -( tmp.*(56+v) - 10 )./(tmp-1).^2/100;

function val = bm(v)
val = 4*exp(-(v+71)/18);

function val = dbm(v)
val = -(2/9)*exp(-(v+71)/18);

function val = ah(v)
val = 0.07*exp(-(v+71)/20);

function val = dah(v)
val = -(7/2000)*exp(-(v+71)/20);

function val = bh(v)
val = 1./(exp(3-(v+71)/10)+1);

function val = dbh(v)
tmp = exp(-(v+41)/10);
val = tmp./(tmp+1).^2/10;
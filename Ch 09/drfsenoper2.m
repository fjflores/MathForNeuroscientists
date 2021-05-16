%
%  drfsenoper.m
%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  Compute the response of an active tree to random synaptic input
%
%  usage:    [t,Vrec] = stEtreesyn(filename,stim)
%
%  e.g.: 
%       filename = 'sep12a.asc'; 'AF-6-01-04-B.asc'; 'hefork.swc'
%
%  stim = struct('Vc',-65,'Tfin',200,'dx',2,'dt',0.005)
%
%  Note: filename must be of either asc or swc type
%
%  Required m-files:  asc_converter.m
%                     swc_converter.m
%                     makemd.m
%                     treeplot.m
%                     makeH.m
%

function [Vc, Ic, t] = drfsenoper(filename,stim)

if filename(end-3:end) == '.asc'
    nc = asc_converter(filename);
else 
    nc = swc_converter(filename,0);
end

md = makemd(nc,stim.dx);
md.filename = filename;

close all

treeplot(md,nc,0)    % plot the tree
set(gca,'tickdir','out')
box off
xlabel('\mum','fontsize',14)
ylabel('\mum','fontsize',14)

Ra = 0.3;                 % k Ohm cm
Cm = 1;                    % uF / cm^2
As = nc.somadata.As*1e-8;

% build H
[H,SA,psa] = makeH(md,Ra,As);
%keyboard
Nx = length(H);

% modify for clamp

HLC = H(1:Nx-1,Nx);   % compartments adjacent to soma
roots = find(HLC);
HLCR = -HLC(roots).*SA;

H(Nx,:) = [];
H(:,Nx) = [];
Nx = Nx - 1;
dH = diag(H);

dt = stim.dt;
Nt = ceil(stim.Tfin/dt);
e = ones(Nx,1);

g = struct('K', 40*e, 'Na', 40*e, 'Cl', 1/15*e);
g.K(end) = 20;
g.Na(end) = 600;
E = struct('K', -77, 'Na', 56, 'Cl', -68);

u = rand(Nx,1);
exc = find(u>=.2);
inh = find(u<.2);
Ne = length(exc);
Ni = length(inh);
gbarE = 1e-6;
gbarI = 2e-6;
gE = gbarE*(1 + max(randn(Ne,1)/10,-.9)); 
gI = gbarI*(1 + max(randn(Ni,1)/10,-.9)); 
tauE = 0.5;
tauI = 1;

figure
t = 0:.01:10;
plot(t,1e6*gbarE*(t/tauE).*exp(1-t/tauE),'k')
hold on
plot(t,1e6*gbarI*(t/tauI).*exp(1-t/tauI),'r')
hold off
legend('g_E','g_I')
set(gca,'tickdir','out','fontsize',14)
box off
xlabel('t  (ms)')
ylabel('nS')

Vexc = 0;
Vinh = -80;

Vc = [-55 -65];
Ic = zeros(2,Nt);
Vr = -75*e;

for cexp = 1:2,

    VcVec = Vc(cexp)*HLC;

    options = optimset('Jacobian','on');
    H(1:Nx+1:end) = dH;       % restore diagonal of H
    Vr = fsolve(@(V) Iss(V,E,g,H,VcVec),Vr,options);    % initial conditions
    V = Vr;

    n = an(V)./(an(V)+bn(V));
    m = am(V)./(am(V)+bm(V));
    h = ah(V)./(ah(V)+bh(V));

    gsyn = zeros(Nx,1);

    Ic(cexp,1) = g.Cl(end)*(As*(Vc(cexp)-E.Cl) - HLCR'*(V(roots)-Vc(cexp)));

    Nesp = 800;
    spetime = sort(stim.Tfin*rand(1,Nesp));
    specomp = 1 + round((Ne-1)*rand(1,Nesp));

    dispecomp = diff(specomp);   % remove doubles
    dbls = find(dispecomp==0);
    specomp(dbls) = [];
    spetime(dbls) = [];
    Nesp = length(specomp);

    te1 = stim.Tfin*ones(Ne,1);
    for i=1:Ne,                 % initialize 1st spk times for each comp
        ehit = find(specomp==i);
        if length(ehit) > 0
            te1(i) = min(spetime(ehit));
        end
    end

    specnt = 1;

    Nisp = 200;
    spitime = sort(stim.Tfin*rand(1,Nisp));
    spicomp = 1 + round((Ni-1)*rand(1,Nisp));

    dispicomp = diff(spicomp);   % remove doubles
    dbls = find(dispicomp==0);
    spicomp(dbls) = [];
    spitime(dbls) = [];
    Nisp = length(spicomp);
    
    figure
    plot(spetime,exc(specomp),'k+')
    hold on
    plot(spitime,inh(spicomp),'rx')
    hold off
    set(gca,'tickdir','out','fontsize',14)
    box off
    xlabel('t  (ms)')
    ylabel('compartment number')
   
    ti1 = stim.Tfin*ones(Ni,1);
    for i=1:Ni,                 % initialize 1st spk times for each comp
        ihit = find(spicomp==i);
        if length(ihit) > 0
            ti1(i) = min(spitime(ihit));
        end
    end

    spicnt = 1;

    for j=2:Nt,

        t = (j-1)*dt;
        t2 = t-dt/2;

        if t2 > spetime(specnt)
            specnt = min(specnt + 1,Nesp);
            te1(specomp(specnt)) = spetime(specnt);
        end

        if t2 > spitime(spicnt)
            spicnt = min(spicnt + 1,Nisp);
            ti1(spicomp(spicnt)) = spitime(spicnt);
        end

        gsyn(exc) = gE.*((t2-te1)/tauE).*exp(1-(t2-te1)/tauE).*(t2>te1)./psa(exc);
        gsyn(inh) = gI.*((t2-ti1)/tauI).*exp(1-(t2-ti1)/tauI).*(t2>ti1)./psa(inh);

        a = an(V);  c = (a+bn(V))/2;
        n = ( (1/dt-c).*n + a) ./ (1/dt + c); n4 = n.^4;

        a = am(V);  c = (a+bm(V))/2;
        m = ( (1/dt-c).*m + a) ./ (1/dt + c);

        a = ah(V);  c = (a+bh(V))/2;
        h = ( (1/dt-c).*h + a) ./ (1/dt + c); m3h = m.^3.*h;

        d = g.Na.*m3h + g.K.*n4 + g.Cl + gsyn;

        f = g.Na.*m3h*E.Na + g.K.*n4*E.K + g.Cl.*E.Cl;
        f(exc) = f(exc) + gsyn(exc)*Vexc;
        f(inh) = f(inh) + gsyn(inh)*Vinh;
        f(roots) = f(roots) - HLC(roots)*Vc(cexp);

        H(1:Nx+1:end) = dH + d + 2*Cm/dt;         % update the diagonal

        Vmid = H\(2*Cm*V/dt + f);

        V = 2*Vmid - V;

        Ic(cexp,j) = g.Cl(end)*(As*(Vc(cexp)-E.Cl) - HLCR'*(V(roots)-Vc(cexp)));

    end  % j

end % cexp

t = linspace(0,stim.Tfin,Nt);

a11 = Vc(1) - Vexc;
a12 = Vc(1) - Vinh;

a21 = Vc(2) - Vexc;
a22 = Vc(2) - Vinh;
A = [a22 -a12; -a21 a11]/(a11*a22-a12*a21);

%I = [i65'-i65(1); i75'-i75(1)];
%I = [i65'; i75'];
g = A*Ic;

figure
hist(1e6*g(1,:),20)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','k','EdgeColor','w')
set(gca,'tickdir','out')
set(gca,'fontsize',14)
box off
xlabel('g_e  (nS)')

figure
hist(1e6*g(2,:),20)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w')
set(gca,'tickdir','out')
set(gca,'fontsize',14)
box off
xlabel('g_i  (nS)')

figure
plot(t,1e6*Ic(1,:),'k',t,1e6*Ic(2,:),'r')
set(gca,'tickdir','out')
set(gca,'fontsize',14)
box off
xlabel('t (ms)')
ylabel('I_c  (pA)')
legend(['V_c = ' num2str(Vc(1))],['V_c = ' num2str(Vc(2))],'location','best')

%plot time course of conductances
figure
plot(t,1e6*g(1,:),'k',t,1e6*g(2,:),'r')
set(gca,'tickdir','out')
set(gca,'fontsize',14)
box off
xlabel('t (ms)')
ylabel('g_e, g_i  (nS)')
legend(['g_e = '],['g_i = '],'location','best')

figure; 
%%%%%
%we subsample at half the step size to compute the power spectra
dt_eff = 2*dt*1e-3; %converts to sec
fs = 1/dt_eff; %compute the cut-off frequency associated with subsampling step
g_es = g(1,1:2:end);
g_es = g_es - mean(g_es); %subtract the mean
[p_e,p_c,f] = pmtm(g_es,14,[],fs,'twosided');

g_is = g(2,1:2:end);
g_is = g_is - mean(g_is);
[p_i,p_c,f] = pmtm(g_is,14,[],fs,'twosided');

inds_plot = find(f < 500);

plot(f(inds_plot),p_e(inds_plot)./max(p_e),'k')
hold on
plot(f(inds_plot),p_i(inds_plot)./max(p_i),'r')

set(gca,'tickdir','out')
set(gca,'fontsize',14)
box off
hold off
xlim([0 500])
xlabel('frequency')
ylabel('power')

return

function [val, jac] = Iss(V,E,g,B,VcVec)
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

val = B*V + g.Na.*m3h.*(V-E.Na) + g.K.*n4.*(V-E.K) + g.Cl.*(V-E.Cl) + VcVec;

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












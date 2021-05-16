%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% getVrJacfull.m
%
% find the rest potential of the full HH cell 
% by providing the Jacobian
%
% To precompute the derivatives of the gating functionals I recommend
%    >> syms V
%    >> bn = .125*exp(-(V+71)/80);
%    >> dbn = simple(diff(bn,V))
%

function Vr = getVrJacfull

E = struct('K', -77, 'Na', 56, 'Cl', -68); % reversal potentials, mV

G = struct('K', 36, 'Na', 120, 'Cl', 0.3);  % channel conductances, mS/cm^2

options = optimset('display','iter','jacobian','on');

Vr = fsolve(@(V) IssJac(V,E,G),-71,options);

function [val,jac] = IssJac(V,E,G)

an = .01*(10-(V+71))./(exp(1-(V+71)/10)-1);
bn = .125*exp(-(V+71)/80);
ninf = an./(an+bn);

dan = -1/1000*(71*exp(-61/10-1/10*V)-10+exp(-61/10-1/10*V).*V)/(exp(-61/10-1/10*V)-1).^2;
dbn = -1/640*exp(-1/80*V-71/80);
dninf = (bn.*dan-an.*dbn)./(an+bn).^2;

am = .1*(25-(V+71))./(exp(2.5-(V+71)/10)-1);
bm = 4*exp(-(V+71)/18);
minf = am./(am+bm);

dam = -1/100*(56*exp(-23/5-1/10*V)-10+exp(-23/5-1/10*V).*V)./(exp(-23/5-1/10*V)-1).^2;
dbm = -2/9*exp(-1/18*V-71/18);
dminf = (bm.*dam-am.*dbm)./(am+bm).^2;

ah = 0.07*exp(-(V+71)/20);
bh = 1./(exp(3-(V+71)/10)+1);
hinf = ah./(ah+bh);

dah = -7/2000*exp(-1/20*V-71/20);
dbh = 1/10./(exp(-41/10-1/10*V)+1).^2.*exp(-41/10-1/10*V);
dhinf = (bh.*dah-ah.*dbh)./(ah+bh).^2;

val = G.Na*minf.^3.*hinf.*(V-E.Na) + G.K*ninf.^4.*(V-E.K) + G.Cl*(V-E.Cl);

jac = G.Na*((3*minf.^2.*dminf.*hinf+minf.^3*dhinf).*(V-E.Na) + minf.^3*hinf) + ...
      G.K*(4*ninf.^3.*dninf.*(V-E.K) + ninf.^4) + G.Cl;


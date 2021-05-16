%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% getVrJac.m
%
% find the rest potential of teh cell with leak and delayed rectifier K current,
% by providing the Jacobian
%
% To precompute the derivatives of the gating functionals I recommend
%    >> syms V
%    >> bn = .125*exp(-(V+71)/80);
%    >> dbn = simple(diff(bn,V))
%

function Vr = getVrJac

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

val = G.K*ninf^4.*(V-E.K) + G.Cl*(V-E.Cl);

jac = G.K*(4*ninf.^3.*dninf.*(V-E.K) + ninf.^4) + G.Cl;


%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  hhsymh.m
%
%  set up the Quasi Active HH+I_h system
%

function [Vr,B] = hhsymh(g)

syms V m n h q

VK = -77;               % mV
VNa = 56;              % mV
VL = -68;            % mV
Vh = -40; %-55;

an = .01*(61+V)/(1-exp(-(V+61)/10));
bn = .125*exp(-(V+71)/80);
ninf = an/(an+bn);

am = .1*(46+V)/(1-exp(-(46+V)/10));
bm = 4*exp(-(V+71)/18);
minf = am/(am+bm);

ah = 0.07*exp(-(V+71)/20);
bh = 1/(exp(-(41+V)/10)+1);
hinf = ah/(ah+bh);

qinf = 1/(1+exp((V+69)/7.1));
tauq = 100*10/(exp((V+66.4)/9.3)+exp(-(V+81.6)/13));

I = g.K*(ninf^4)*(V-VK) + g.Na*(minf^3)*hinf*(V-VNa) + g.Cl*(V-VL) + g.h*(qinf^2)*(V-Vh);
i = inline(char(I));
Vr = fsolve(@(V) i(V), -30);

nr = subs(ninf,V,Vr);
mr = subs(minf,V,Vr);
hr = subs(hinf,V,Vr);
qr = subs(qinf,V,Vr);
 
F(1) = -g.K*(n^4)*(V-VK) - g.Na*(m^3)*h*(V-VNa) - g.Cl*(V-VL) - g.h*(q^2)*(V-Vh);

F(2) = an*(1-n) - bn*n;

F(3) = am*(1-m) - bm*m;

F(4) = ah*(1-h) - bh*h;

F(5) = (qinf - q)/tauq;

J = jacobian(F,[V n m h q]);

B = subs(J,{V,n,m,h,q},{Vr,nr,mr,hr,qr});


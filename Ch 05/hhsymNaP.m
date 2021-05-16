%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  hhsymNaP.m
%
%  set up the Quasi Active Hodgkin Huxley + I_NaP system
%

function [Vr,B] = hhsymNaP(gNaP)

syms V m n h

gK = 36;             % mS / (cm)^2
gNa = 120;           % mS / (cm)^2
gl = 0.3;            % mS / (cm)^2

VK = -77;               % mV
VNa = 56;              % mV
VL = -68;            % mV

an = .01*(61+V)/(1-exp(-(V+61)/10));
bn = .125*exp(-(V+71)/80);
ninf = an/(an+bn);

am = .1*(46+V)/(1-exp(-(46+V)/10));
bm = 4*exp(-(V+71)/18);
minf = am/(am+bm);

ah = 0.07*exp(-(V+71)/20);
bh = 1/(exp(-(41+V)/10)+1);
hinf = ah/(ah+bh);

pinf = 1/(1+exp(-(V+49)/5));

I = gK*(ninf^4)*(V-VK) + (gNa*(minf^3)*hinf+gNaP*pinf)*(V-VNa) + gl*(V-VL);
i = inline(char(I));
Vr = fsolve(@(V) i(V), -30);

nr = subs(ninf,V,Vr);
mr = subs(minf,V,Vr);
hr = subs(hinf,V,Vr);
 
F(1) = -gK*(n^4)*(V-VK) - (gNa*(m^3)*h+gNaP*pinf)*(V-VNa) - gl*(V-VL);

F(2) = an*(1-n) - bn*n;

F(3) = am*(1-m) - bm*m;

F(4) = ah*(1-h) - bh*h;

J = jacobian(F,[V n m h]);

B = subs(J,{V,n,m,h},{Vr,nr,mr,hr});


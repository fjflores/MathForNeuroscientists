%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% exocytosis.m
%
%  Exercise the 2-pool model of Sorensen et al
%  (with af as the mutant switch) of the role of
%  synaptotagmin in exocytosis
%
as = 0.8;
af = 4.9;
A = zeros(11);
A(1:4,1:3) = as*[-3 0 0; 3 -2 0; 0 2 -1; 0 0 1];
A(5:9,5:8) = af*[-4 0 0 0;4 -3 0 0;0 3 -2 0; 0 0 2 -1; 0 0 0 1];

k1 = 0.12; km1 = 0.1;
km2 = 0.05;
gs = 20; gf = 1450;
bs = 4; bf = 56; b = 0.55;
B = zeros(11);
B(1:4,2:4) = bs*[1 0 0;-1 2 0;0 -2 3;0 0 -3];
B(5:9,6:9) = bf*[1 0 0 0;-1 2*b 0 0; 0 -2*b 3*b^2 0; 
                0 0 -3*b^2 4*b^3; 0 0 0 -4*b^3];
B(4,4) = B(4,4) - gs;
B(1,1) = -k1-km2;
B(1,5) = km1;
B(5,1) = k1;
B(5,5) = -km1;
B(9,9) = B(9,9) - gf;
B(10,4) = gs;
B(11,9) = gf;

rmax = 55;
KD = 2.3;
dt = .05;
T = 8;
Nt = T/dt;

x = zeros(11,1);
e1 = x;
x(1) = 200;
x(5) = 300;
e1(1) = 1;

ES = zeros(Nt,1); EF = ES;
c = ES;
t = 0;
a2 = log(10)/2.5;
a1 = exp(-1.5*a2);
c(1) = min(a1*exp(t*a2),20);


for j=2:Nt
    
    t = (j-1)*dt;
    
    Y = 2*eye(11) + dt*(c(j-1)*A+B);
    c(j) = min(a1*exp(t*a2),20);
    f = dt*rmax*(c(j-1)/(c(j-1)+KD) + c(j)/(c(j)+KD))*e1;
    W = 2*eye(11) - dt*(c(j)*A+B);
    x = W\(Y*x+f);
    ES(j) = x(10); 
    EF(j) = x(11);
       
end

t = linspace(0,T,Nt);
subplot(3,1,1)
semilogy(t,c,'k')
box off
%xlabel('t (s)','fontsize',14)
ylabel('[Ca^{2+}]  (\muM)','fontsize',14)

subplot(3,1,2)
plot(t,ES,'r--')
hold on
plot(t,EF,'r')
plot(t,EF+ES,'k')
hold off
box off
ylim([0 800])
legend('E_{SRP}','E_{RRP}','E_{SRP}+E_{RRP}','location','best')
legend boxoff
ylabel('\Delta C_m  (fF)','fontsize',14)
text(3,700,'WT')

af = 2.4;   % run the mutant
A(5:9,5:8) = af*[-4 0 0 0;4 -3 0 0;0 3 -2 0; 0 0 2 -1; 0 0 0 1];
x = zeros(11,1);
e1 = x;
x(1) = 200;
x(5) = 300;
e1(1) = 1;

ES = zeros(Nt,1); EF = ES;
c = ES;
t = 0;
a2 = log(10)/2.5;
a1 = exp(-1.5*a2);
c(1) = min(a1*exp(t*a2),20);


for j=2:Nt

    t = (j-1)*dt;

    Y = 2*eye(11) + dt*(c(j-1)*A+B);
    c(j) = min(a1*exp(t*a2),20);
    f = dt*rmax*(c(j-1)/(c(j-1)+KD) + c(j)/(c(j)+KD))*e1;
    W = 2*eye(11) - dt*(c(j)*A+B);
    x = W\(Y*x+f);
    ES(j) = x(10);
    EF(j) = x(11);

end

t = linspace(0,T,Nt);
subplot(3,1,3)
plot(t,ES,'r--')
hold on
plot(t,EF,'r')
plot(t,EF+ES,'k')
hold off
box off
ylim([0 800])
legend('E_{SRP}','E_{RRP}','E_{SRP}+E_{RRP}','location','best')
legend boxoff
text(3,700,'R233Q')
xlabel('t (s)','fontsize',14)
ylabel('\Delta C_m  (fF)','fontsize',14)
%axis tight
    
    
    
    

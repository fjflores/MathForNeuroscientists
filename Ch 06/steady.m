%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% steady.m
%
% compute the steady state response with with
% 	N compartments   and a stimulus at compartment Nloc
%
% e.g.,   steady(41,[1 21])
%

function steady(N,Nloc)

a = 1e-4;
ell = .1;
I0 = 1e-5;
C_m = 1;		% micro F / cm^2
G_L = 1/15; %0.3;     		% mS / cm^2
R_2 = 0.3; %0.034;		% k Ohm cm

dx = ell/N;		% patch length
x = dx/2:dx:ell-dx/2;	% vector of patch midpoints

S = (2*eye(N)-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1))/dx/dx;
S(1,1) = 1/dx/dx;
S(N,N) = 1/dx/dx;

tau = C_m/G_L;
lambda = sqrt(a/(2*R_2*G_L))

[Q,th] = eig(S);

th = diag(th);

z = (1+lambda^2*th)/tau;

col = {'k' 'r'};

for k=1:2

    P = Q*diag(Q(Nloc(k),:)./z');
    w = I0*N/(2*pi*a*ell)*sum(P,2);
    plot(x,w,col{k})
    hold on

end

hold off
box off
xlabel('x  (cm)','fontsize',14)
ylabel('v_\infty  (mV)','fontsize',14)
%legend('s=1','s=21','location','NE')


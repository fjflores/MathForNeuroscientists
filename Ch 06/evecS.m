
%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% evecS.m
%
% compute and display the first 4 eigenvectors of the
% the N-by-N second difference matrix for a cell of length ell
%
%    evecS(ell,N)    e.g.,   evecS(0.1,40)
%

function evecS(ell,N)

dx = ell/N;	
x = dx/2:dx:ell-dx/2;

S = (2*eye(N)-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1))/dx/dx;
S(1,1) = 1/dx/dx;
S(N,N) = 1/dx/dx;

[Q,Z] = eig(S);

z = diag(Z);

z(1:4)

plot(x,Q(:,1),'k',x,Q(:,2),'k--',x,Q(:,3),'r',x,Q(:,4),'r--')
box off
legend('q_0','q_1','q_2','q_3','location','SE')
xlabel('x  (cm)','fontsize',14)
ylabel('eigenvector','fontsize',14)


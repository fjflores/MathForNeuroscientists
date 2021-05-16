%
% Gabbiani and Cox, Mathematics for Neuroscientists
%
% test the xor net trained by mlptrain
%
% usage,   xordrive
%

function [V,W] = xordrive

close all

P = [0 1 0 1
     0 0 1 1];    % the 4 input patterns of xor

t = [0 1 1 0]';   % their associated xor values

Nhid = 2;
maxiter = 2e5;
rate = 0.05;

[V,W] = mlptrain(P,t,Nhid,maxiter,rate);  % train up the system

out = zeros(1,4);

for j=1:4,

    p = P(:,j);	% grab the first input
    h = s(V*p);		% process through level 1
    out(j) = s(W*h)';	% process through level 2

end

out

return

%
% mlptrain
%
% train a 2 level feedforward net to associate states and targets
%
% usage [V,W] = ocrtrain(P,t,Nhid,maxiter,rate)
%
% where		state is the Nin-by-Npat  array of input patterns
%
%	        maxiter is the number of training iterations
%		rate is the learning rate
%

function [V,W] = mlptrain(state,target,Nhid,maxiter,rate)

[Nin Npat] = size(state);   
Nout = size(target,2);    
V = randn(Nhid,Nin);  % initialize V,
W = randn(Nout,Nhid);  % initialize W,
Vtr = zeros(4,maxiter);
Wtr = zeros(2,maxiter);
%V = [6 6;1 1];
%W = [12 -15];

for iter=1:maxiter,

      j = ceil(rand*Npat);		% index of input state
      i = ceil(rand*Nout);	% index of output bit

      p = state(:,j);   % input pattern
      h = s(V*p);		% hidden layet pattern
      o = s(W(i,:)*h);	% output pattern

      tmp = (o-target(j,i))*o*(1-o);    % the common part of the update
  
      V = V - rate*tmp*(W(i,:)'.*h.*(1-h))*p';  % update V
      
      W(i,:) = W(i,:) - rate*tmp*h';	% update row i of W
      
      Vtr(:,iter) = V(:);
      Wtr(:,iter) = W';
      
end

figure(1)
plot(Wtr(1,:),'r')
hold on
plot(Wtr(2,:),'k')
box off
legend('W_1','W_2','location','best')
xlabel('iteration','fontsize',14)
ylabel('synaptic weight','fontsize',14)
hold off

figure(2)
plot(Vtr(1,:),'k')
hold on
plot(Vtr(3,:),'k--')
plot(Vtr(2,:),'r')
plot(Vtr(4,:),'r--')
box off
legend('V(1,1)','V(2,1)','V(1,2)','V(2,2)','location','best')
xlabel('iteration','fontsize',14)
ylabel('synaptic weight','fontsize',14)
hold off

return

function val = s(x)   % the soft theshold
val = 1./(1+exp(-x+.5));

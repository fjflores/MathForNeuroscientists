%
% Gabbiani and Cox
%
% Test ocrtrain on a set of hand written numbers from
% http://yann.lecun.com/exdb/mnist/
%
%  using readMNIST.m available at
%  https://www.mathworks.com/matlabcentral/fileexchange/27675-read-digits-and-labels-from-mnist-database/content/readMNIST.m
%
% usage,   ocrdrive
%

function [V,W] = ocrdrive

close all

Npat = 200;

[imgs labels]=readMNIST('train-images.idx3-ubyte','train-labels.idx1-ubyte',Npat,0);
[c1, ia1] = unique(labels);

state = zeros(400,8);

figure(1)
for i=1:8, 
    
    subplot(2,4,i) 
    imagesc(round(imgs(:,:,ia1(i))));
    colormap('gray')
    axis equal
    axis off
    %title(num2str(labels(ia1(i+5))))
    
    state(:,i) = round(reshape(imgs(:,:,ia1(i)),400,1));
    
end

%target = [0 0;0 1;1 0;1 1];
target = [0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1];

Nhid = 400;
maxiter = 5e3;
rate = 0.3;

[V,W] = ocrtrain(state,target,Nhid,maxiter,rate);  % train up the system

for j=1:8,    % test for errors

    p = state(:,j);	% grab the first encoded letter

    h = s(V*p);		% process through level 1
    
    out = s(W*h)'

end

return

%
% ocrtrain
%
% train a 2 level feedforward net to associate states and targets
%
% usage [V,W] = ocrtrain(state,target,Nhid,maxiter,rate)
%
% where		state is the Nin-by-Npat  array of input patterns
%
%	        maxiter is the number of training iterations
%		rate is the learning rate
%

function [V,W] = ocrtrain(state,target,Nhid,maxiter,rate)

[Nin Npat] = size(state);   % imgs are 20-by-20
Nout = size(target,2);    % we need only classify numbers 0-9
V = (2*rand(Nhid,Nin)-1)/100; %randn(Nhid,Nin)/100;  % initialize V,
W = (2*rand(Nout,Nhid)-1)/100; %randn(Nout,Nhid)/100;  % initialize W,

Vtr = zeros(maxiter,1);
Wtr = Vtr;
Ztr = Vtr;
Ztrs = Vtr;

for iter=1:maxiter,

      j = ceil(rand*Npat);		% index of input state
      i = ceil(rand*Nout);		% index of output bit

      p = state(:,j);   % input pattern
      h = s(V*p);		% hidden layet pattern
      o = s(W(i,:)*h);	% output pattern

      tmp = (o-target(j,i))*o*(1-o);    % the common part of the update
  
      gv = tmp*(W(i,:)'.*h.*(1-h))*p';
      V = V - rate*gv;  % update V
      
      gw = tmp*h';
      W(i,:) = W(i,:) - rate*gw;	% update row i of W
      
      if mod(iter,1)==0

          Vtr(iter) = mean(W(:));
          Wtr(iter) = std(W(:));
          Ztr(iter) = mean(V(:));
          Ztrs(iter) = std(V(:));
          iter
      end

end

figure(2)
plot(Vtr,'r')
hold on
plot(Vtr+Wtr,'r--')
plot(Ztr,'k')
plot(Ztr+Ztrs,'k--')
legend('mean(W)','mean(W) \pm std(W)','mean(V)','mean(V) \pm std(V)')
plot(Ztr-Ztrs,'k--')
plot(Vtr-Wtr,'r--')
hold off
xlabel('iteration','fontsize',14)
ylabel('synaptic weight','fontsize',14)
ylim([-.15 .15])
box off
return

function val = s(x)   % the soft theshold
val = 1./(1+exp(-x+.5));
return


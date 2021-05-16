%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
%  spikepca.m
%
%  Generate a noisy trace corresponding to the extracellular potential
%  in the neighborhood of 3 cells. Align the active events into a data
%  matrix, compute the associated singular value decomposition and 
%  finally plot the triplet of scores for each event. You should be
%  able to dynamically rotate the axes of figure 4 and arrive at 3
%  relevant clusters (as in Fig 14.2B).
%
%  This routine requires:  stEqatrain.m    and   hhsym.m
%

function spikepca

close all

G = struct('K', [36 38 40], 'Na', [140 120 100], 'Cl', [2 1 0.3]);
stim = struct('amp',[.95 0.9 0.7]*1e3,'per',[50 90 70],...
               'tau1',[0.6 0.5 0.7],'tau2',[0.5 0.4 0.6]);
                                             
dt = 0.02;
Tfin = 3000;

[t,v1] = stEqatrain(dt,Tfin,G,stim,1);   % get train from cell 1
figure(1)
subplot(4,1,1)
plot(t/1000,v1,'k')
set(gca,'xticklabel',[])
ylim([-1 3.5])
%axis tight
box off
ylabel('cell 1','fontsize',14)

vtot = v1;
    
[t,v2] = stEqatrain(dt,Tfin,G,stim,2);   % get train from cell 2
subplot(4,1,2)
plot(t/1000,v2,'k')
set(gca,'xticklabel',[])
ylim([-1 3.5])
%axis tight
box off
ylabel('cell 2','fontsize',14)
vtot = vtot + v2;

[t,v3] = stEqatrain(dt,Tfin,G,stim,3);  % get train from cell 3
subplot(4,1,3)
plot(t/1000,v3,'k')
set(gca,'xticklabel',[])
ylim([-1 3.5])
%axis tight
box off
ylabel('cell 3','fontsize',14)
vtot = vtot + v3;

vtot = vtot + randn(size(v3))/10;    % noisy sum of trains

subplot(4,1,4)
plot(t/1000,vtot,'k')
ylim([-1 6])
%axis tight
box off
ylabel('cell sum','fontsize',14)
xlabel('t  (s)','fontsize',14)

% gather spikes

vtot = [vtot; vtot(end)*ones(600,1)];
t = [t; t(end)*ones(600,1)];

hot = find(vtot>1.5 & t>50)';     % find large parts of vtot

dh = diff(hot);
dhi = find(dh>50);
lo = [hot(1) hot(dhi+1)];
spN = length(lo);

figure(2)

row = 0;
st = linspace(0,12,601)';
for j=1:spN
   
    vrec = vtot(lo(j)-150:lo(j)+450);   % back 3, ahead 12
    
    if max(vrec) < 4
       plot(st(1:4:end),vrec(1:4:end),'k');  % plot every 4th val because too dense
       hold on
       row = row + 1;
       spdat(row,:) = vrec;
       v = [v1(lo(j)) v2(lo(j)) v3(lo(j))];
       [jnk,clab(row)] = max(v);		% label for later confirmation
    end
    
end

box off
xlim([0 12])
xlabel('t  (ms)','fontsize',14)
ylabel('mV','fontsize',14)
hold off

spdat = spdat';

[spm,spn] = size(spdat);
spdat0 = spdat;
spdat = spdat - repmat(mean(spdat,2),1,spn);
[Y,Sig,X] = svd(spdat);
sig = diag(Sig);
figure(3)
semilogy(sig(sig>1),'kx-')
box off
xlabel('index','fontsize',14)
ylabel('singular value','fontsize',14)

sc1 = spdat'*Y(:,1);
sc2 = spdat'*Y(:,2);
sc3 = spdat'*Y(:,3);
figure(4)

sc11 = spdat(:,find(clab==1))'*Y(:,1);
sc21 = spdat(:,find(clab==1))'*Y(:,2);
sc31 = spdat(:,find(clab==1))'*Y(:,3);
plot3(sc11,sc21,sc31,'ro','markersize',14)
hold on
sc12 = spdat(:,find(clab==2))'*Y(:,1);
sc22 = spdat(:,find(clab==2))'*Y(:,2);
sc32 = spdat(:,find(clab==2))'*Y(:,3);
plot3(sc12,sc22,sc32,'rs','markersize',14)
sc13 = spdat(:,find(clab==3))'*Y(:,1);
sc23 = spdat(:,find(clab==3))'*Y(:,2);
sc33 = spdat(:,find(clab==3))'*Y(:,3);
plot3(sc13,sc23,sc33,'rd','markersize',14)

plot3(sc1,sc2,sc3,'k+','markersize',10) 	% the scores
grid
xlabel('s_1','fontsize',14)
ylabel('s_2','fontsize',14)
zlabel('s_3','fontsize',14)
hold off


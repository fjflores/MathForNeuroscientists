%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
%  treeplot  (see stEtreesyn for usage)
%

function treeplot(md,nc,pbr)

close all
figure(1)
hold on

terms = repmat(1,1,md.count.N);
terms(md.par) = 0;
termdex = find(terms == 1);

for jj = 1:length(md.par)

    xyz = nc.celldata{md.par(jj)};
    plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k');

end

for jj = 1:length(termdex)

    xyz = nc.celldata{termdex(jj)};
    plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k');

end

for jj = 1:length(nc.roots)

    xyz = nc.celldata{nc.roots(jj)};
    plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k');

end

sd = nc.somadata;
xs = sd.x; xs(end+1) = sd.x(1);  
ys = sd.y; ys(end+1) = sd.y(1);  
zs = sd.z; zs(end+1) = sd.z(1);
plot3(xs,ys,zs,'k','linewidth',2)

% xlim = get(gca,'XLim');
% ylim = get(gca,'YLim');
% 
% mlen = round(1/8*(xlim(2)-xlim(1))*1e-1)*1e1;
% 
% xp = [xlim(1) xlim(1)+mlen];
% yp1 = [ylim(2)-1/8*(ylim(2)-ylim(1)) ylim(2)-1/8*(ylim(2)-ylim(1))];
% yp2 = [ylim(2)-3/16*(ylim(2)-ylim(1)) ylim(2)-3/16*(ylim(2)-ylim(1))];
% 
% plot(xp,yp1,'k','linewidth',2)
% mscale = text(xp(1),yp2(1),[num2str(mlen) '\mum']);

% ax = axis;
% mscale = text(ax(1),ax(3),md.filename);

if pbr  % display branch numbers

    for i = 1:length(nc.celldata);
        hold on
        xyz = nc.celldata{i};
        h = text(xyz(1,1),xyz(1,2),xyz(1,3),num2str(i));
        set(h,'Color','k','FontSize',8);
    end

end

return

function allPolys(polys)
% polys = resRoCoS.runs{3,1}.polys2D;
% opt=res.origVal;

% a=ones(length(polys),1)*opt;
% b=ones(length(polys),1)*opt;

for k =1:length(polys)
    b = cell2mat(polys{k}(2:end, 4));
    L(k)= min(b);
    U(k)= max(b);
end

n = round(sqrt(length(polys)));

for k =1:length(polys)
    
    x = cell2mat(polys{k}(2:end, 2));
    y = cell2mat(polys{k}(2:end, 3));
    z = cell2mat(polys{k}(2:end, 4));

%     sx = model.S(groups{gp(k,1)}.row, groups{gp(k,1)}.column);
%     sy = model.S(groups{gp(k,2)}.row, groups{gp(k,2)}.column);

%     ax(k)=subplot(n+1,n,k);    
    ax(k)=subplot(n,n,k);
    scatter(x, y,35,z,'filled','s');
    hold on;
%     title(k)
     bla=polys{k};
      xlabel(bla{1,2});
      ylabel(bla{1,3});
    caxis([min(L) max(U)])
%     xlim([-50 50]);
%     ylim([-50 50]);
end
% caxis([min(L) max(U)])
% h=colorbar;
% set(h, 'Position', [.95 .12 .02 .8])
% 
% for i=1:length(polys)
%     pos=get(ax(i), 'Position');
%     set(ax(i), 'Position', [pos(1) pos(2) pos(3) pos(4)]);
% end
%  hold off;
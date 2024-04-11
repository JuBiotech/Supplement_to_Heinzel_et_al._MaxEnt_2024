function PolygoneSubplot(polys, opt, w)
% PolygoneSubplot(resP.runs{3,1}.polys2D, res.origVal, [7,8,64])

% gp = nchoosek(1:length(groups), 2);

% a=ones(length(polys),1)*opt;
% b=ones(length(polys),1)*opt;

for k =1:length(polys)
    b = cell2mat(polys{k}(2:end, 4));
    L(k)= min(b);
    U(k)= max(b);
end

for k =1:length(w)
    
    x = cell2mat(polys{w(k)}(2:end, 2));
    y = cell2mat(polys{w(k)}(2:end, 3));
    z = cell2mat(polys{w(k)}(2:end, 4));

%     sx = model.S(groups{gp(w(k),1)}.row, groups{gp(w(k),1)}.column);
%     sy = model.S(groups{gp(w(k),2)}.row, groups{gp(w(k),2)}.column);

    grid on;
    ax(k)=subplot(1,length(w),k);
    scatter(x*100, y*100,180,z,'filled','s');
    set(gca,'FontSize',30)
      bla=polys{w(k)};
      xlabel(bla{1,2});
      ylabel(bla{1,3});
      
      caxis([min(L) max(U)])
end
hold on;
grid on;
caxis([min(L) max(U)])
h=colorbar;
set(h, 'Position', [.95 .12 .02 .8],'FontSize',30)

for i=1:length(w)
    pos=get(ax(i), 'Position');
    set(ax(i), 'Position', [pos(1) pos(2) pos(3) pos(4)]);
end
 hold off;
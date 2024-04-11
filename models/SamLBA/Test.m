% PolygoneSubplot(resMin18RoCoS.runs{3,1}.polys2D, resMin18.origVal, [7,8,64])
% PolygoneSubplot(resMaxBmRoCoS.runs{3,1}.polys2D, resMaxBm.origVal, [7,8,64])
% PolygoneSubplot(resMax26RoCoS.runs{3,1}.polys2D, resMax26.origVal, [7,8,64])

w=[7,8,64];

polys1=resMin18RoCoS.runs{3,1}.polys2D;

polys2=resMaxBmRoCoS.runs{3,1}.polys2D;

polys3=resMax26RoCoS.runs{3,1}.polys2D;

for k =1:length(polys1)
    b1 = cell2mat(polys1{k}(2:end, 4));
    L1(k)= min(b1);
    U1(k)= max(b1);
    b2 = cell2mat(polys2{k}(2:end, 4));
    L2(k)= min(b2);
    U2(k)= max(b2);
    b3 = cell2mat(polys3{k}(2:end, 4));
    L3(k)= min(b3);
    U3(k)= max(b3);
end

for k =1:length(w)
    
    x = cell2mat(polys1{w(k)}(2:end, 2));
    y = cell2mat(polys1{w(k)}(2:end, 3));
    z = cell2mat(polys1{w(k)}(2:end, 4));
    grid on;
    ax(k)=subplot(3,length(w),k);
    scatter(x*100, y*100,180,z,'filled','s');
    set(gca,'FontSize',30)
      bla=polys1{w(k)};
      xlabel(bla{1,2});
      ylabel(bla{1,3});
      caxis([min(L1) max(U1)])
end
hold on;
grid on;
caxis([min(L1) max(U1)])
h=colorbar;
set(h, 'Position', [.93 .7 .02 .25],'FontSize',30)

for i=1:length(w)
    pos=get(ax(i), 'Position');
    set(ax(i), 'Position', [pos(1) pos(2) pos(3) pos(4)]);
end
 hold off;
 
 for k =1:length(w)
    
    x = cell2mat(polys2{w(k)}(2:end, 2));
    y = cell2mat(polys2{w(k)}(2:end, 3));
    z = cell2mat(polys2{w(k)}(2:end, 4));

    grid on;
    ax(k)=subplot(3,length(w),3+k);
    scatter(x*100, y*100,180,z,'filled','s');
    set(gca,'FontSize',30)
      bla=polys2{w(k)};
      xlabel(bla{1,2});
      ylabel(bla{1,3});
      caxis([min(L2) max(U2)])
end
hold on;
grid on;
caxis([min(L2) max(U2)])
h=colorbar;
set(h, 'Position', [.93 .39 .02 .25],'FontSize',30)

for i=1:length(w)
    pos=get(ax(i), 'Position');
    set(ax(i), 'Position', [pos(1) pos(2) pos(3) pos(4)]);
end
 hold off;
 
 for k =1:length(w)
    
    x = cell2mat(polys3{w(k)}(2:end, 2));
    y = cell2mat(polys3{w(k)}(2:end, 3));
    z = cell2mat(polys3{w(k)}(2:end, 4));
    
    grid on;
    ax(k)=subplot(3,length(w),6+k);
    scatter(x*100, y*100,180,z,'filled','s');
    set(gca,'FontSize',30)
      bla=polys3{w(k)};
      xlabel(bla{1,2});
      ylabel(bla{1,3});
      caxis([min(L3) max(U3)])
end
hold on;
grid on;
caxis([min(L3) max(U3)])
h=colorbar;
set(h, 'Position', [.93 .09 .02 .25],'FontSize',30)

for i=1:length(w)
    pos=get(ax(i), 'Position');
    set(ax(i), 'Position', [pos(1) pos(2) pos(3) pos(4)]);
end
 hold off;
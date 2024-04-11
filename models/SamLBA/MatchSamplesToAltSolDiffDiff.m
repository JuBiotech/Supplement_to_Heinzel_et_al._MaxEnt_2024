% dpc=[0.1,0.5,1,2,5,8,10,12,15,17,20,23,25,28,30,50];
dpc=25;
% dpc=[1,5,10,25];

for i=1:length(dpc)
    load(['Lee_max26_mcBM_' num2str(dpc(i)) 'pc_100kSamp_truncnormal.mat']);
    [data,dists] = matchSamplesToAltSolutions( res, altOptmax26 );
    ind1=find(data(:,1)==1);
    ind2=find(data(:,1)==2);
    ind3=find(data(:,1)==3);
    ind4=find(data(:,1)==4);
%     subplot(2,2,i);
%     plot(res.diffs(ind1),data(ind1,2),'.y',res.diffs(ind2),data(ind2,2),'.k',res.diffs(ind3),data(ind3,2),'.r',res.diffs(ind4),data(ind4,2),'.b')
%     scatter(res.diffs,data(:,2),100,data(:,1))
%     subplot(4,4,i);plot(res.diffs,data(:,2),'.r')
%     title(dpc(i))
end

subplot(2,1,1)
plot(repmat(diffmax26(1),1,8),0:10:70,'-r',repmat(diffmax26(2),1,8),0:10:70,'-c',repmat(diffmax26(3),1,8),0:10:70,'-b',repmat(diffmax26(4),1,8),0:10:70,'-k','LineWidth',4);
hold on;
plot(res.diffs,dists(:,1),'.r',res.diffs,dists(:,2),'.c',res.diffs,dists(:,3),'.b',res.diffs,dists(:,4),'.k')
subplot(2,1,2)
plot(repmat(diffmax26(1),1,8),0:10:70,'-r',repmat(diffmax26(2),1,8),0:10:70,'-c',repmat(diffmax26(3),1,8),0:10:70,'-b',repmat(diffmax26(4),1,8),0:10:70,'-k','LineWidth',4);
hold on;
plot(res.diffs(ind1),data(ind1,2),'.r',res.diffs(ind2),data(ind2,2),'.c',res.diffs(ind3),data(ind3,2),'.b',res.diffs(ind4),data(ind4,2),'.k')

%% On the Road Plots:
% ind=1;
% figure;
% for k=0:1.1:39;
%     ZwSp=find(dists(:,1)>k & dists(:,1)<k+1);
% %     Scheibe3=dists(ZwSp,:);
%     
% %     figure; 
%     for j=1:size(dists(ZwSp,:),1) 
%         subplot(4,9,ind)
%         plot(1:4,dists(ZwSp(j),:),'-k');hold on; 
%     end
%     ind=ind+1;
% end
% figure;
% subplot(2,2,1); 
% stairs(0:0.5:65,histc(dists(:,1),0:0.5:65),'Color',[0.5 0.5 0.5],'LineWidth',2);hold on;
% stairs(0:0.5:65,histc(dists(ind1,1),0:0.5:65),'-r','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind2,1),0:0.5:65),'-c','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind3,1),0:0.5:65),'-b','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind4,1),0:0.5:65),'-k','LineWidth',2);
% subplot(2,2,2); 
% stairs(0:0.5:65,histc(dists(:,2),0:0.5:65),'Color',[0.5 0.5 0.5],'LineWidth',2);hold on;
% stairs(0:0.5:65,histc(dists(ind1,2),0:0.5:65),'-r','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind2,2),0:0.5:65),'-c','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind3,2),0:0.5:65),'-b','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind4,2),0:0.5:65),'-k','LineWidth',2);
% subplot(2,2,3); 
% stairs(0:0.5:65,histc(dists(:,3),0:0.5:65),'Color',[0.5 0.5 0.5],'LineWidth',2);hold on;
% stairs(0:0.5:65,histc(dists(ind1,3),0:0.5:65),'-r','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind2,3),0:0.5:65),'-c','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind3,3),0:0.5:65),'-b','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind4,3),0:0.5:65),'-k','LineWidth',2);
% subplot(2,2,4); 
% stairs(0:0.5:65,histc(dists(:,4),0:0.5:65),'Color',[0.5 0.5 0.5],'LineWidth',2);hold on;
% stairs(0:0.5:65,histc(dists(ind1,4),0:0.5:65),'-r','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind2,4),0:0.5:65),'-c','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind3,4),0:0.5:65),'-b','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind4,4),0:0.5:65),'-k','LineWidth',2);
% figure;
% subplot(2,2,1); hold on;
% stairs(0:0.5:65,histc(dists(ind1,1),0:0.5:65),'-r','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind1,2),0:0.5:65),'-r','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind1,3),0:0.5:65),'-r','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind1,4),0:0.5:65),'-r','LineWidth',2);
% subplot(2,2,2); hold on;
% stairs(0:0.5:65,histc(dists(ind2,1),0:0.5:65),'-c','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind2,2),0:0.5:65),'-c','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind2,3),0:0.5:65),'-c','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind2,4),0:0.5:65),'-c','LineWidth',2);
% subplot(2,2,3); hold on;
% stairs(0:0.5:65,histc(dists(ind3,1),0:0.5:65),'-b','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind3,2),0:0.5:65),'-b','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind3,3),0:0.5:65),'-b','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind3,4),0:0.5:65),'-b','LineWidth',2);
% subplot(2,2,4); hold on;
% stairs(0:0.5:65,histc(dists(ind4,1),0:0.5:65),'-k','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind4,2),0:0.5:65),'-k','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind4,3),0:0.5:65),'-k','LineWidth',2);
% stairs(0:0.5:65,histc(dists(ind4,4),0:0.5:65),'-k','LineWidth',2);
% figure;
% subplot(2,2,1); stairs(0:0.5:65,histc(dists(:,1),0:0.5:65),'-k')
% subplot(2,2,2); stairs(0:0.5:65,histc(dists(:,2),0:0.5:65),'-k')
% subplot(2,2,3); stairs(0:0.5:65,histc(dists(:,3),0:0.5:65),'-k')
% subplot(2,2,4); stairs(0:0.5:65,histc(dists(:,4),0:0.5:65),'-k')
% figure;
% subplot(2,2,1); stairs(0:0.5:65,histc(dists(ind1,1),0:0.5:65),'-r')
% subplot(2,2,2); stairs(0:0.5:65,histc(dists(ind1,2),0:0.5:65),'-r')
% subplot(2,2,3); stairs(0:0.5:65,histc(dists(ind1,3),0:0.5:65),'-r')
% subplot(2,2,4); stairs(0:0.5:65,histc(dists(ind1,4),0:0.5:65),'-r')
% figure;
% subplot(2,2,1); stairs(0:0.5:65,histc(dists(ind2,1),0:0.5:65),'-b')
% subplot(2,2,2); stairs(0:0.5:65,histc(dists(ind2,2),0:0.5:65),'-b')
% subplot(2,2,3); stairs(0:0.5:65,histc(dists(ind2,3),0:0.5:65),'-b')
% subplot(2,2,4); stairs(0:0.5:65,histc(dists(ind2,4),0:0.5:65),'-b')
% figure;
% subplot(2,2,1); stairs(0:0.5:65,histc(dists(ind3,1),0:0.5:65),'-c')
% subplot(2,2,2); stairs(0:0.5:65,histc(dists(ind3,2),0:0.5:65),'-c')
% subplot(2,2,3); stairs(0:0.5:65,histc(dists(ind3,3),0:0.5:65),'-c')
% subplot(2,2,4); stairs(0:0.5:65,histc(dists(ind4,4),0:0.5:65),'-c')
% figure;
% subplot(2,2,1); stairs(0:0.5:65,histc(dists(ind4,1),0:0.5:65),'-g')
% subplot(2,2,2); stairs(0:0.5:65,histc(dists(ind4,2),0:0.5:65),'-g')
% subplot(2,2,3); stairs(0:0.5:65,histc(dists(ind4,3),0:0.5:65),'-g')
% subplot(2,2,4); stairs(0:0.5:65,histc(dists(ind4,4),0:0.5:65),'-g')
% figure;
% subplot(2,2,1); hist(dists(:,1),50)
% subplot(2,2,2); hist(dists(:,2),50)
% subplot(2,2,3); hist(dists(:,3),50)
% subplot(2,2,4); hist(dists(:,4),50)
% figure;
% subplot(2,2,1); hist(dists(ind1,1),50);
% subplot(2,2,2); hist(dists(ind1,2),50);
% subplot(2,2,3); hist(dists(ind1,3),50);
% subplot(2,2,4); hist(dists(ind1,4),50);
% figure;
% subplot(2,2,1); hist(dists(ind2,1),50);
% subplot(2,2,2); hist(dists(ind2,2),50);
% subplot(2,2,3); hist(dists(ind2,3),50);
% subplot(2,2,4); hist(dists(ind2,4),50);
% figure;
% subplot(2,2,1); hist(dists(ind3,1),50);
% subplot(2,2,2); hist(dists(ind3,2),50);
% subplot(2,2,3); hist(dists(ind3,3),50);
% subplot(2,2,4); hist(dists(ind3,4),50);
% figure;
% subplot(2,2,1); hist(dists(ind4,1),50);
% subplot(2,2,2); hist(dists(ind4,2),50);
% subplot(2,2,3); hist(dists(ind4,3),50);
% subplot(2,2,4); hist(dists(ind4,4),50);
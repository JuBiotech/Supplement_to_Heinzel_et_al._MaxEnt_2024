clc;
figure;
hold on;
pic1 = subplot(2,2,1);
plot(metaTNcycle(:,1),metaTNcycle(:,2:4),'--ks','LineWidth',2,...
                'MarkerEdgeColor','r','MarkerSize',10);
datacursormode on
title('with Cycles');
pic2 = subplot(2,2,2); 
plot(metaTNnoCycle(:,1),metaTNnoCycle(:,2:4),'--ks','LineWidth',2,...
                'MarkerEdgeColor','b','MarkerSize',10);
title('without Cycles');
pic3 = subplot(2,2,3); 
hist(resCycle.diffs, 50);
ylim([0 1000]);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r')
set(pic3,'XMinorGrid','on')
pic4 = subplot(2,2,4); 
hist(resNoCycle.diffs, 50);
ylim([0 1000]);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b')
set(pic4,'XMinorGrid','on')
hold off;

% figure;
% hist([resCycle.diffs,resNoCycle.diffs], 50);
% set(get(gca,'child'),'facecolor','r');
% legend('Cycles','NoCycles')

% with (A) and without (B) allowing cycles
boundA=[1.05 17.6 23.4 29.2 40.0; 1.88 19.3 25.1 30.0 42];
boundB=[1.10 16.9 23.1 29.2 39.8; 1.98 18.7 24.8 30.1 45];

% lsgA(1,:) = resCycle.origSol;
% lsgA(2,:) = extractAverageSolution(resCycle, boundA(1,1), boundA(2,1));
% lsgA(3,:) = extractAverageSolution(resCycle, boundA(1,2), boundA(2,2));
% lsgA(4,:) = extractAverageSolution(resCycle, boundA(1,3), boundA(2,3));
% lsgA(5,:) = extractAverageSolution(resCycle, boundA(1,4), boundA(2,4));
% lsgA(6,:) = extractAverageSolution(resCycle, boundA(1,5), boundA(2,5));
% 
% lsgB(1,:) = resNoCycle.origSol;
% lsgB(2,:) = extractAverageSolution(resNoCycle, boundB(1,1), boundB(2,1));
% lsgB(3,:) = extractAverageSolution(resNoCycle, boundB(1,1), boundB(2,1));
% lsgB(4,:) = extractAverageSolution(resNoCycle, boundB(1,1), boundB(2,1));
% lsgB(5,:) = extractAverageSolution(resNoCycle, boundB(1,1), boundB(2,1));
% lsgB(6,:) = extractAverageSolution(resNoCycle, boundB(1,1), boundB(2,1));
% 
% % xlswrite('ComparisonFluxMaps.xlsx', lsgB, 'test','B8')
%

figure;
hold on;
for k=1:5
    disp(['k = ', num2str(k)]);
    indsA = (resCycle.diffs >= boundA(1,k)) & (resCycle.diffs <= boundA(2,k));
    vA = resCycle.sols(:,indsA);
    indsB = (resNoCycle.diffs >= boundB(1,k)) & (resNoCycle.diffs <= boundB(2,k));
    vB = resNoCycle.sols(:,indsB);

    varA=zeros(1,535);
    for i = 1:535, varA(i) = var(vA(i, :)); end
    m1 = max(varA);
    varB=zeros(1,535);
    for i = 1:535, varB(i) = var(vB(i, :)); end
    m2 = max(varB);
    x = round(10*max(m1,m2))/10;
    mw=repmat(x/10,10,length(10));
    
    disp(['   x = ', num2str(x/10)]);
    pic1=subplot(5,2,2*k-1);
    hist(varA,50); hold on;
    plot(mw,-1:11:100,'Color','green','LineWidth',2)
    datacursormode on;
    title('with Cycles');
    xlim([-.1 x]);
    ylim([0 46]);
    h1 = findobj(gca,'Type','patch');
    set(h1,'FaceColor','r')
    set(pic1,'XMinorGrid','on')
    pic2=subplot(5,2,2*k); 
    hist(varB,50);hold on;
    plot(mw,-1:11:100,'Color','green','LineWidth',2)
    title('without Cycles');
    ylim([0 46]);
    xlim([-.1 x]);
    h2 = findobj(gca,'Type','patch');
    set(h2,'FaceColor','b')
    set(pic2,'XMinorGrid','on')

%     pic1=subplot(5,1,k);
%     hist([varA;varB],50);
%     legend('varA','varB')
%     xlim([-.1 x]);
%     ylim([0 46]);
%     set(pic1,'XMinorGrid','on')
    
    indA=find(varA > x/10);
    indB=find(varB > x/10);
    L(k,1)=length(indA);
    L(k,2)=length(indB);
end
hold off;
L
% find(varA > 9.01 & varA < 9.91)
% find(varA > 9.91 & varA < 10.8)
% find(varA > 11.7 & varA < 12.6)
% find(varA > 44.1 & varA < 50.0)
% find(varB > 8.56 & varB < 9.42)
% find(varB > 10.3 & varB < 11.1)
% find(varB > 16.3 & varB < 17.1)
% find(varB > 41.9 & varB < 50.0)

% figure;
% subplot(2,1,1);
% hist(vA(44,:),50);
% subplot(2,1,2);
% hist(vB(44,:),50);
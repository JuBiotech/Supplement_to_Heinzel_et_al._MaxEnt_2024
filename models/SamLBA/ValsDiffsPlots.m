F4 = 'Lee_minATP_mcBM_^pc_100kSamp_truncnormal';
% F3 = 'Lee_max26_mcBM_^pc_100kSamp_truncnormal';
% F1 = 'Lee_min18_mcBM_^pc_100kSamp_truncnormal';
F2 = 'Lee_maxBM_mcBM_^pc_100kSamp_truncnormal';

files=regexp([F2 '°' F4], '\°','split');
% files=regexp([F1 '°' F2 '°' F3 '°' F4], '\°','split');
% pc = [0.1, 0.5, 1, 2, 5, 8, 10, 12, 15, 17, 20, 23, 25, 28, 30, 32, 35, 38, 40, 45, 50, 55, 65, 75, 85, 95];
% pc = [0.1, 0.5, 1, 2, 5, 10, 15, 20, 25, 30, 35, 40];
pc = [25, 35, 38, 40, 45, 50];%[1, 10, 25, 50, 75];
%% plotten
for p=1:length(files)
    A(p,:)=regexp(files(p), '\^','split');
end

folder='Lee2000/';
ending = '.mat';
n=length(pc);
m=length(files);
hb = 50;

figure;
k=1;
for i = 1:n
    for j = 1:m
        subplot(n,m,k);
        datacursormode on;
        pcS = num2str(pc(i));
        load([folder,A{j,1}{1,1},pcS,A{j,1}{1,2},ending]);
        t=find(abs(res.vals)>=0);
%         if length(find(abs(res.vals)==0))>0
%         t(k,1) = length(find(abs(res.vals)==0));
%         t(k,2) = pc(i);
%         end
        plot(res.vals(t),res.diffs(t),'.','Color',[0/255 91/255 130/255]);
        title([A{j,1}{1,1},pcS],'Interpreter','None')
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',[0,91/255,130/255],'EdgeColor',[128/255,128/255,128/255])
%         T(i,j)=max(res.vals(t));
%         T(i,j+m)=min(res.vals(t));
%         T(i,j+2*m)=max(res.vals(t))-min(res.vals(t));
%         T(i+n,j)=max(res.diffs(t));
%         T(i+n,j+m)=min(res.diffs(t));
        k = k+1;
    end
end
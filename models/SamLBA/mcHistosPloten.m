% F1 = 'Lee_maxBM_mcBM_^pc_100kSamp_truncnormal';
% F2 = 'Lee1k_maxBM_mcBM_^pc_100kSamp_truncnormal';
% F3 = 'Lee26_max26_mcBM_^pc_100kSamp_truncnormal';
% files=regexp([F1 '°' F2 '°' F3], '\°','split');

F4 = 'Lee_minATP_mcBM_^pc_100kSamp_truncnormal';
F3 = 'Lee_max26_mcBM_^pc_100kSamp_truncnormal';
F1 = 'Lee_min18_mcBM_^pc_100kSamp_truncnormal';
F2 = 'Lee_maxBM_mcBM_^pc_100kSamp_truncnormal';

files=regexp([F2 '°' F4], '\°','split');
% files=regexp([F1 '°' F2 '°' F3], '\°','split');
% files=regexp([F1 '°' F2 '°' F3 '°' F4], '\°','split');
pc = [1, 10, 25, 50, 75];%[38, 40];%

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
        hist(abs(res.vals(t)),hb);
        title([A{j,1}{1,1},pcS],'Interpreter','None')
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',[0,91/255,130/255],'EdgeColor',[128/255,128/255,128/255])
        k = k+1;
    end
end

figure;
k=1;
for i = 1:n
    for j = 1:m
        subplot(n,m,k);
        datacursormode on;
        pcS = num2str(pc(i));
        load([folder,A{j,1}{1,1},pcS,A{j,1}{1,2},ending]);
        t=find(abs(res.vals)>=0);
        hist(res.diffs(t),hb);
        title([A{j,1}{1,1},pcS],'Interpreter','None')
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',[0,91/255,130/255],'EdgeColor',[128/255,128/255,128/255])
        k = k+1;
    end
end

clear all
% a='testtest';
% b='^test°bla';
% c='^testblablub';
% z=[a,b,c];
% v=regexp(z, '\°','split')
% for p=1:length(v)
%     A(p,:)=regexp(v(p), '\^','split')
% end
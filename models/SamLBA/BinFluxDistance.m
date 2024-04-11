% function BinFluxDistance(res,center,edgeB,edgeA)
% res=resB;
center=2.28;
edgeB=1.6;
edgeA=3;

AltSolCenter=extractSolutionWithDiff(res,center);

i=find(res.diffs <= edgeA & res.diffs >= edgeB);
altSolBin=res.sols(:,i);

for k=1:length(altSolBin) 
    diffAltSol(k) = norm(AltSolCenter - altSolBin(:, k)); 
end

figure;
hist(diffAltSol)
title('max26 (40%) diff=2.28','Interpreter','None')
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0,91/255,130/255],'EdgeColor',[128/255,128/255,128/255])

% figure;
% subplot(1,3,1)
% plotPolyIsoLines(res10G6PKG.runs(1:9))
% title('r1=10-G6P-KG')
% subplot(1,3,2)
% plotPolyIsoLines(res10G6PR5P.runs(1:9))
% title('r1=10-G6P-R5P')
% subplot(1,3,3)
% plotPolyIsoLines(res10KGR5P.runs(1:9))
% title('r1=10-KG-R5P')
function [ result, resave ] = ppPlots( fileName, res, objValFormula )
%PPPLOTS Summary of this function goes here
%   Detailed explanation goes here

    resave = false;
    
    result = {};
    
    objVals = objValFormula(res.vals);

    clf;
    plotHistPercent(objVals, 30);
    saveas(gcf, [fileName(1:end-4) '_objVal.fig']);
    saveas(gcf, [fileName(1:end-4) '_objVal.png']);
end


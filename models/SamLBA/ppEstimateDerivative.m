function [ result, resave ] = ppEstimateDerivative( fileName, res )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    resave = false;
    p = polyfit(((res.samples./res.origData)' - 1)*100, res.vals,1);
    result = [res.origData, p];
end


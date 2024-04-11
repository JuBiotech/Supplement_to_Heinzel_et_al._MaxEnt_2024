clear all;clc;
changeCobraSolver('glpk','LP');changeCobraSolver('glpk','MILP');
clear ans

%% setup Model
M = struct('S', [], 'b', [], 'c', [], 'lb', [], 'ub', []);
M.S=[1 1 -2 0 1 0;
    1 -1 0 1 0 0;
    0 0 1 -2 0 -1.77;
    0 -1 0 1 -1 1.38];
[m,r]=size(M.S);
M.lb=[0; -10; -10; -10; 0; 0];
M.ub=ones(r,1)*10;M.ub(6)=20;
% upt=8.4; M.lb(1)=upt;M.ub(1)=upt;
M.rxns={'v1' 'v2' 'v3' 'v4' 'v5' 'v6'};
M.mets={'A' 'B' 'C' 'D'};        
M.c=zeros(r,1);M.c(6)=1;
M.b=zeros(m,1);

%% classical FBA
[nomSol,~, ~]=quickSolveFBA(M, true);

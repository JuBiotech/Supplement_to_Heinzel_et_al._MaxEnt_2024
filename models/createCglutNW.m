function [CoryneModel]=createCglutNW(filename)
RefModel=readCbModel('Ec_iJR904_GlcMM',1000,'SBML');

CoryneModel = removeRxns(RefModel,char(RefModel.rxns(1:end-1)));
CoryneModel.genes(1:end-1)='';
CoryneModel.rxnGeneMat(1:end-1)=[];
CoryneModel.subSystems='';

z='.xls';

[addNum, addTxt]=xlsread([filename z], 'used');

addID=addTxt(:,1);
addForm=addTxt(:,3);
addGRA=addTxt(:,4);
addSub=addTxt(:,6);
addRev=addNum(:,1);
addLB=addNum(:,2);
addUB=addNum(:,3);
addObj=addNum(:,4);

for i=1:length(addID)
    CoryneModel = addReaction(CoryneModel,addID{i},addForm{i},'',addRev(i),addLB(i),addUB(i),addObj(i),addSub{i},addGRA{i});
end

CoryneModel.metNames=CoryneModel.mets;
CoryneModel.metFormulas=CoryneModel.mets;
CoryneModel = removeRxns(CoryneModel,char(CoryneModel.rxns(1)));
CoryneModel.genes(1)='';
CoryneModel.rxns=CoryneModel.rxns';
% CoryneModel.lb=CoryneModel.lb;
% CoryneModel.ub=CoryneModel.ub;
% CoryneModel.c=CoryneModel.c;

cobraSolverName = 'lp_solve';
cobraSolverName = 'glpk';
LPsolverOK = changeCobraSolver(cobraSolverName, 'LP');
clc;
end


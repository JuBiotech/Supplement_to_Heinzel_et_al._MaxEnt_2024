changeCobraSolver('glpk','LP');changeCobraSolver('glpk','MILP');

[sols]=compareLooplessFBA(modelmax26,nomSol(26,1));

for j=1:size(sols,2)
    for k=j:size(sols,2)
        solsLee(j,k) = norm(sols(2:end,j)-sols(2:end,k)); 
        if solsLee(j,k)<5e-04
            solsLee(j,k)=0; 
        end
    end
end

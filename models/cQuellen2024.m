clear all; clc; 
% Solver laden
changeCobraSolver('glpk','LP');changeCobraSolver('glpk','MILP'); clear ans

% Netzwerk laden
M=readCbModel('iEZ481_2024');clc
% M=createCglutNW('CoryneNetwork_20150803');clc
M.rxns=M.rxns';
M.description='Cg_iEZ478_iEZ480_unchecked';
[m,r]=size(M.S);

% Standard Einstellungen
M=changeRxnBounds(M,'EX_o2_e',-30,'l');
M=changeRxnBounds(M,'PI_t_NA',0,'b');
M=changeRxnBounds(M,'PI_t_ATP',0,'b');
M=changeRxnBounds(M,'aceA',0,'b');
M=changeRxnBounds(M,'aceB',0,'b');
M=changeRxnBounds(M,'NH3_t_H',0,'b');
M=changeRxnBounds(M,'UREA_d',0,'b');
M=changeRxnBounds(M,'UREA_t_H',0,'b');
M=changeRxnBounds(M,'SUC_t_H',0,'b');
M=changeRxnBounds(M,'FUM_t_H',0,'b');
M=changeRxnBounds(M,'MAL_t_H',0,'b');
M=changeRxnBounds(M,'NA_t_H',0,'b');
M=changeRxnBounds(M,'LAC_L_t_H',0,'b');
M=changeRxnBounds(M,'LYS_t_H',0,'b');
M=changeRxnBounds(M,'AC_t_H',0,'b');
M=changeRxnBounds(M,'NO3_t_H',0,'b');
M=changeRxnBounds(M,'GLU_t_H',0,'b');
M=changeRxnBounds(M,'THR_t_H',0,'b');
M=changeRxnBounds(M,'nagB',0,'l');
M=changeRxnBounds(M,'PPI_d',0,'b');
M=changeRxnBounds(M,'tkt_1',0,'l');
M=changeRxnBounds(M,'tkt_2',0,'l');
M=changeRxnBounds(M,'tal',0,'l');
M=changeRxnBounds(M,'glgC',0,'b');
M=changeRxnBounds(M,'CO2_d',0,'l');
M=changeRxnBounds(M,'aspA',0,'l');
M=changeRxnBounds(M,'aspB',0,'l');
M=changeRxnBounds(M,'leakage',-50,'l');
M=changeRxnBounds(M,'leakage',50,'u');
M=changeRxnBounds(M,'folD2',-1,'l');
M=changeRxnBounds(M,'metF',-1,'l');
M=changeRxnBounds(M,'folD2',1,'u');
M=changeRxnBounds(M,'metF',1,'u');
M=changeRxnBounds(M,'pps',0,'b');
M=changeRxnBounds(M,'glpD',0,'b');
M=changeRxnBounds(M,'putA',0,'b');
M=changeRxnBounds(M,'tpiA',-100,'l');
M=changeRxnBounds(M,'aceA',-100,'l');
M=changeRxnBounds(M,'aceA',100,'u');
M=changeRxnBounds(M,'aceB',100,'u');
M=changeRxnBounds(M,'quiC',0,'b');
M=changeRxnBounds(M,'PROTEIN_b',0,'b');
M=changeRxnBounds(M,'DNA_a',0,'b');

for i=1:r
    if(i<=m)
        if(strcmp(M.mets{i},'h[e]')==1),           ind(1)=i;
        elseif(strcmp(M.mets{i},'h[c]')==1),       ind(2)=i;
        elseif(strcmp(M.mets{i},'pi[c]')==1),      ind(3)=i;
        elseif(strcmp(M.mets{i},'atp[c]')==1),     ind(4)=i;
        elseif(strcmp(M.mets{i},'adp[c]')==1),     ind(5)=i;           
        end 
    end
   
    if(strcmp(M.rxns{i},'biomass_a')==1), Rkt(1)=i; end
end
M.S(ind(4),Rkt(1))=-(29.2+6.4);
M.S(ind(5),Rkt(1))=29.2+6.4;
M.S(ind(3),Rkt(1))=29.2+6.4;
M=changeRxnBounds(M,'mt',0,'b'); % [1.08,1.8] nomineller Wert 1.4

%% Ende Vorabeinstellungen
%% Start C-Quellen Varianten:
%% PCA
M0=changeObjective(M,'biomass_a',0);% Biomasse als Optimierungsziel aus
M0=changeObjective(M0,'EX_pca',1); % PCA Uptake als Opt-Ziel an
M0=changeRxnBounds(M0,'EX_glc_e',0,'b'); %Glucose Uptake =0
M0=changeRxnBounds(M0,'EX_pca',-10,'l'); % PCA Uptake ermöglichen
M0=changeRxnBounds(M0,'biomass_a',0.1,'b'); % Wachstumsrate festsetzen auf 0.1 (b=both: lower =upper bound)
% Messwerte Wachstumsraten auf PCA als einziger C-Quelle
% M0=changeRxnBounds(M0,'biomass_a',0.099089,'b');
% M0=changeRxnBounds(M0,'biomass_a',0.096591,'b');
% M0=changeRxnBounds(M0,'biomass_a',0.096,'l');
% M0=changeRxnBounds(M0,'biomass_a',0.0991,'u');

[fluxMapM0,SolM0, ~]=quickSolveFBA(M0, false); % FBA
for i=1:r, if(abs(fluxMapM0(i))<0.00001),fluxMapM0(i)=0;end, end 
length(find(fluxMapM0~=0)) % Kleine Flüsse auf Null runden

% Indizes für FLüsse über Systemgrenze raus suchen 
n=1;
for i=1:r
    if(sum(abs(M0.S(:,i)))==1)
        ex(n)=i;n=n+1;
    end
end
% Flüsse über Systemgrenze (Bezeichnung + Wert) in Variable ausgeben
n=1;
for i=1:length(ex)
    if(fluxMapM0(ex(i))~=0)
        EX0{n,1}=M0.rxns(ex(i));
        EX0{n,2}=fluxMapM0(ex(i));
        n=n+1;
    end
end
clear ex
%Netzwerk in SBML überführen
% writeCbToSBML(M0,'Cglut20150803_PCA');

%% Glucose with PCA
M1=changeObjective(M,'biomass_a',0);% Biomasse als Optimierungsziel aus
M1=changeObjective(M1,'EX_glc_e',1);% Glucose Uptake als Opt-Ziel an
M1=changeRxnBounds(M1,'EX_glc_e',-30,'l');% Glucose Uptake ermöglichen Lower bound =30
M1=changeRxnBounds(M1,'EX_glc_e',30,'u');% Glucose Uptake ermöglichen Upper bound =30

% Messwerte Wachstumsraten auf Glucose+PCA als C-Quelle
% M1=changeRxnBounds(M1,'biomass_a',0.6,'b');
% M1=changeRxnBounds(M1,'biomass_a',0.493689,'b');
% M1=changeRxnBounds(M1,'biomass_a',0.522519,'b');

% Wachstumsrate festsetzen auf [0.493;0.523]
M1=changeRxnBounds(M1,'biomass_a',0.493,'l');
M1=changeRxnBounds(M1,'biomass_a',0.523,'u');
M1=changeRxnBounds(M1,'EX_pca',-1.1,'b');% PCA Uptake festsetzen auf 1.1
% FBA
[fluxMapM1,SolM1, ~]=quickSolveFBA(M1, false); 
% Kleine Flüsse auf Null runden
for i=1:r, if(abs(fluxMapM1(i))<0.00001),fluxMapM1(i)=0;end, end
length(find(fluxMapM1~=0))
% Indizes für FLüsse über Systemgrenze raus suchen 
n=1;
for i=1:r
    if(sum(abs(M0.S(:,i)))==1)
        ex(n)=i;n=n+1;
    end
end
% Flüsse über Systemgrenze (Bezeichnung + Wert) in Variable ausgeben
n=1;
for i=1:length(ex)
    if(fluxMapM1(ex(i))~=0)
        EX1{n,1}=M1.rxns(ex(i));EX1{n,2}=fluxMapM1(ex(i));
        n=n+1;
    end
end
clear ex
%Netzwerk in SBML überführen
% writeCbToSBML(M0,'Cglut20150803_GLC-PCA');

%% Plots
% pcaUpt=0:0.1:2;
% for i=1:length(pcaUpt)
%     ZWS=changeRxnBounds(M1,'EX_pca',-pcaUpt(i),'b');
%     [~,optS, ~]=quickSolveFBA(ZWS, true);    
%     objFluxS0(i)=-1*optS;
% end
% plot(pcaUpt,objFluxS0, 'sb',pcaUpt,abs(SolM1)*ones(1,length(objFluxS0)), '-r' )
% xlabel('PCA uptake','FontSize',20);
% ylabel('Glc uptake','FontSize',20);

%% BHI 
% Optimierungziel per Voreinstellung Wachstumsrate
% M2=changeRxnBounds(M,'EX_glc_e',0,'b');% Glucose Uptake =0
M2=changeRxnBounds(M,'EX_glc_e',-30,'l');% Glucose Uptake ermöglichen
M2=changeRxnBounds(M2,'EX_glc_e',30,'u');% Glucose Uptake ermöglichen
% Wachstumsrate festsetzen auf [0.1;0.89]
M2=changeRxnBounds(M2,'biomass_a',0.1,'l'); % 0.79,'l');
M2=changeRxnBounds(M2,'biomass_a',0.89,'u');

% (Aminosäure/Fettsäure) Uptake ermöglichen
M2=changeRxnBounds(M2,'EX_ala_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_arg_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_asn_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_asp_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_cys_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_gln_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_glu_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_gly_e',-10,'l');
M2=changeRxnBounds(M2,'EX_his_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_ile_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_leu_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_lys_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_met_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_phe_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_pro_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_ser_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_thr_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_trp_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_tyr_L_e',-10,'l');
M2=changeRxnBounds(M2,'EX_val_L_e',-10,'l');

M2=changeRxnBounds(M2,'EX_suc_e',0,'b');
M2=changeRxnBounds(M2,'EX_cit_e',0,'b');
M2=changeRxnBounds(M2,'EX_for_e',0,'b');
M2=changeRxnBounds(M2,'EX_ac_e',0,'b');

M2=changeRxnBounds(M2,'EX_C140',-10,'l');
M2=changeRxnBounds(M2,'EX_C150',-10,'l');
M2=changeRxnBounds(M2,'EX_C160',-10,'l');
M2=changeRxnBounds(M2,'EX_C180',-10,'l');
M2=changeRxnBounds(M2,'EX_C161',-10,'l');
M2=changeRxnBounds(M2,'EX_C181',-10,'l');

% FBA
[fluxMapM2,SolM2, ~]=quickSolveFBA(M2, true); 
% Kleine Flüsse auf Null runden
for i=1:r, if(abs(fluxMapM2(i))<0.0001),fluxMapM2(i)=0;end, end
% Indizes für FLüsse über Systemgrenze raus suchen 
n=1;
for i=1:r
    if(sum(abs(M2.S(:,i)))==1)
        ex(n)=i;n=n+1;
    end
end
% Flüsse über Systemgrenze (Bezeichnung + Wert) in Variable ausgeben
n=1;
for i=1:length(ex)
    if(fluxMapM2(ex(i))~=0)
        EX2{n,1}=M2.rxns(ex(i));
        EX2{n,2}=fluxMapM2(ex(i));
        n=n+1;
    end
end
%Netzwerk in SBML überführen

% writeCbToSBML(M0,'M0_biomass');
% writeCbToSBML(M1,'M1_glc');
% writeCbToSBML(M2,'M2_BHI');

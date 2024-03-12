%% FINAL PIPELINE - LOOPLESS+COMPACTED+ROOM test 4 & test 1 
%Test 4: Loopless+Compact+ParallelSL+ROOM_SoS
%Test 1: Loopless+Compact+ParallelSL
%Notes on the script: None of the comments from the PC3 pipeline have been removed as they may serve as reference.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUNE 2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Integrate Experimental Boundaries and Simplificate
%load('Recon1_tunning_original.mat')
%boundsFilename = 'FlxBound_initial_M_text.txt';
     %Call a function that updates the reaction names as they are different
     %between the experimental data and the model
%          MatchingNames = MatchRxnNamesPC3(boundsFilename);
%         for i = 1:size(MatchingNames,1)
%             model = changeRxnBounds(model,MatchingNames{i,1},...
%                 MatchingNames{i,2},"l");
%             model = changeRxnBounds(model,MatchingNames{i,1},...
%                 MatchingNames{i,3},"u");
%         end 
% [model_compacted,modelLoopless] = modelSimplification(model,MatchingNames(:,1));
%% Until Here
%LOAD working models with boundaries compact and loopless
%     load('RECON1CompactLoopless_S_FinalVersion.mat');
%     modelOriginal = model;


%diary('EcoTest_iJO1366.txt')
%% Load and constraint original EcoK12 model to match the BW25113 strain
pathEcoK12='iJO1366.mat'
EcoK12=load(pathEcoK12)
EcoK12=EcoK12.iJO1366
[EcoK12] = generateRules(EcoK12)

%%Load Trancriptomic data 
pathToTranscri='MetabolicEco_Average.csv'
MetbolicEco= readtable(pathToTranscri,'Format','%s%f','Delimiter',',','TextType','string')
%Parameters for model compaction
%Find ExcReactions (avoided in loopless analysis) 
ExcRxnsBool = findExcRxns(EcoK12);
ExcRxns= EcoK12.rxns(ExcRxnsBool);
%Biomass Rxns
biomass=EcoK12.rxns(find(EcoK12.c));
%Transport & exchange rxns (avoided in model compaction)
TransportRxns=findTransRxns(EcoK12,true) 
TransportRxns=vertcat(TransportRxns,biomass)


EcoK12=changeRxnBounds(EcoK12,'BIOMASS_Ec_iJO1366_core_53p95M',0.001,'l')
%Generate Compact Loopless model
[EcoCompacted,EcoLoopless,~,EcoRemovedRxns] = GenerateCompactLooplessTissueModel(EcoK12,biomass,ExcRxns,TransportRxns)

modelLoopless=EcoLoopless
model_compacted=EcoCompacted
modelOriginal=EcoK12
%%
EcoRemovedRxns(ismember(EcoRemovedRxns,EcoCompacted.rxns)) %There was an extra reaction that must be removed from the removed list
EcoRemovedRxns(23) = [];

%%
%%Characterize model (parameter for lethal)

optimalBM    = optimizeCbModel(model_compacted);
BM           = optimalBM.f;
lowerBoundBM = 0;
model_compacted.lb(find(model_compacted.c)) = lowerBoundBM;
model_compacted.ub(find(model_compacted.c)) = BM;
[EcoSumTable,EcoExpressionNumeric]=mapExpressionToRxn_Compacted(EcoK12,EcoRemovedRxns,MetbolicEco)


%%% PARAMETERS FOR CIMAT
expressionRxns=EcoExpressionNumeric
threshold_lb= prctile(expressionRxns(expressionRxns~=-1),40);
threshold_ub= prctile(expressionRxns(expressionRxns~=-1),60);
runtime = 300;
epsilon= 1;
tol = 1e-8;
core= {};
logfile = 'MILPlog';
 

 %Call iMAT and save MILPProblem (includes the A matrix)
  ExcRxns            = findExcRxns(model_compacted);
 IgnoreForLethality = model_compacted.rxns(ExcRxns);
 IgnoreForLethality = [IgnoreForLethality; model_compacted.rxns(find(model_compacted.c))];
 
 
  disp('Starting CiMAT...')
 [tissueModel,MILPproblem,~,~,solutionIMAT] = CiMAT(modelLoopless, expressionRxns, threshold_lb, threshold_ub, tol, core, logfile, runtime, epsilon,model_compacted);
                    MILPproblem.c(find(model_compacted.c)) = 1;
                    MILPproblem.S = full(MILPproblem.A);
                    model         = model_compacted;
                    BiomassiMAT   = solutionIMAT.cont(find(model_compacted.c));
                    MILPproblem.ub(find(model_compacted.c)) = BiomassiMAT;
  %Just in case: This will be used when identifying synthetic lethals
%if certain issues appear.
    disp('Starting iMAT_2...')
            [~,MILPproblemNormal,~,~] = iMAT_2(modelOriginal, expressionRxns, threshold_lb, threshold_ub, tol, core, logfile, runtime, epsilon);
                    MILPproblemNormal.c(find(model_compacted.c)) = 1;
                    MILPproblemNormal.S = full(MILPproblemNormal.A);
        
%clearvars -except  solutionIMAT model tolerance MILPproblem model_compacted modelOriginal tissueModel modelLoopless BiomassiMAT%% Characterize model - Sampling, Continuos and Discrete Distance, Find most dense space of
% solutions
 disp('Starting sampling...')
       samplingTime = 10;
      tic
        [minSolution,maxSolution] = SamplingWithNoKOsMaxMin(model,samplingTime,'Continuous',MILPproblem);
      toc
   clear MILPproblem
%%
cutoff = 0.05;
[Jsl,Jdl,Jtl] =  parallelSL_23(modelOriginal,cutoff,IgnoreForLethality);
%In case some error was returned after
%%
[JslCompactedFinal,JdlCompactedFinal] = Lethals4CompactedModelSL(modelLoopless,model_compacted,Jsl,Jdl);
JtlCompactedFinal={}
ArrayWithLethals = [{JslCompactedFinal};{JdlCompactedFinal};{JtlCompactedFinal}];
%%
disp('SLROOM') %NOTE! The fastSLROOM2 was modified to return only ERs and SLs - not triplets. 
[JslROOM,JdlROOM,BMReduction]= fastSLROOM2(model,cutoff,ArrayWithLethals,maxSolution,minSolution);%% It is probable that most lethals will return no answer in the ROOM, because the lethals in
%the simplified model could make the model infeasible. Therefore the following if
%statement is to perform again the fastSLROOM2 with the original model in
%case every Lethal returned NaN
%%
condition1 = sum(cellfun(@isnan,BMReduction{1,1}(:,3))) == size(BMReduction{1,1},1);
condition2 = sum(cellfun(@isnan,BMReduction{2,1}(:,3))) == size(BMReduction{2,1},1);
%condition3 = sum(cellfun(@isnan,BMReduction{3,1}(:,4))) == size(BMReduction{3,1},1);
            if  condition1 && condition2 && condition3
                disp('Performing FastSLROOM with original model, simplified model was too constrained')
                  samplingTime = 10;
                  tic
                    [minSolution,maxSolution] = SamplingWithNoKOsMaxMin(modelOriginal,samplingTime,'Continuous',MILPproblemNormal);
                  toc
                      %Eliminate possible Lethals that were in loops
                        JslLoopless      = Jsl(ismember(Jsl,modelLoopless.rxns),:);
                        tempVar(:,1)     = ismember(Jdl(:,1),modelLoopless.rxns);
                        tempVar(:,2)     = ismember(Jdl(:,2),modelLoopless.rxns);
                        JdlLoopless      = Jdl(sum(tempVar,2) == 2,:);
                        tempVar2(:,1) = ismember(Jtl(:,1),modelLoopless.rxns);
                        tempVar2(:,2) = ismember(Jtl(:,2),modelLoopless.rxns);
                        tempVar2(:,3) = ismember(Jtl(:,3),modelLoopless.rxns);
                        JtlLoopless   = Jtl(sum(tempVar2,2) == 3,:);
                        ArrayWithLethals                                        = [{JslLoopless};{JdlLoopless};{JtlLoopless}];
                        model = modelOriginal;
                        [JslROOM,JdlROOM,JtlROOM,BMReduction]                   = fastSLROOM2(model,cutoff,ArrayWithLethals,maxSolution,minSolution);
                        
            end
%%
%It is probable that most lethals will return no answer as the lethals tend
%to make the model infeasible
    [pointsCellArrayMaxSingles, pointsCellArrayMinSingles]  = SamplingROOMSynLethals(minSolution,maxSolution,model,JslROOM,samplingTime);
    %%
    [pointsCellArrayMaxPairs, pointsCellArrayMinPairs]      = SamplingROOMSynLethals(minSolution,maxSolution,model,JdlROOM,samplingTime);
    %[pointsCellArrayMaxTriplets, pointsCellArrayMinTriplets]= SamplingROOMSynLethals(minSolution,maxSolution,model,JtlROOM,samplingTime);
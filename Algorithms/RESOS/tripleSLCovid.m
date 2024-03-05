function [efficientSingleTargets,effectiveTargetPair,effectiveTargetTriplet] = tripleSLCovid(method,targets,modelDrug,modelPlacebo, placeboEuclMinSolution,placeboEuclMaxSolution,drugEuclMaxSolution,drugEuclMinSolution,threshold)
%%  [slist_id,dlist_id,tlist_id]=tripleSL(model,cutoff,eliList,atpm)
% INPUT
% method
    %1 = Original fastSL algorithm
    %2 =  ROOM (Requires two extra inputs: maxSolution minSolution)
% targets
    %FORMAT TO BE DEFINED YET
    %list of targets of a certain drug
% modelDrug
%   Specific model for patients who received the drug
%model Placebo
%   Specific model for patients who received a placebo
% placeboEuclMaxSolution = Using modelPlacebo, After a sampling and euclidean distance filtration, this
%   solution is the one that has the maximum euclidean distance value within the
%   valid solutions.
% placeboEuclMinSolution = Using modelPlacebo After a sampling and euclidean distance filtration, this
%   solution is the one that has minimum euclidean distance value within the
%   valid solutions.

% drugEuclMaxSolution = Using modelDrug, After a sampling and euclidean distance filtration, this
%   solution is the one that has the maximum euclidean distance value within the
%   valid solutions.
% drugEuclMinSolution = Using modelDrug After a sampling and euclidean distance filtration, this
%   solution is the one that has minimum euclidean distance value within the
%   valid solutions.
            
            
%OUTPUT
% EffectiveSingleTargets        Targets that manage to reduce the euclidean
            % distance between the drugWT and the original Placebo WT
% EffectivePairTargets        Indices of pair of reactions  that manage to reduce the euclidean
            % distance between the drugWT and the original Placebo WT
% EffectiveTripletsTargets        Indices of triplets of reactions  that manage to reduce the euclidean
            % distance between the drugWT and the original Placebo WT

%%

if iscell(targets)
    targetsIdx = find(ismember(modelPlacebo.rxns,targets));
end

if ~exist('method','var') ||  isempty(method)
    method = 1;
    disp('Variable "method" not set. Using traditional fastSL algorithm...')
end

if method == 1
    if or(~exist('modelPlacebo','var'),~exist('modelDrug','var'))
        error('For original fastSL algorithm, one of the following, or both, variables is missing: modelDrug modelPlacebo')
    end
else 
    if or(or(~exist('placeboEuclMinSolution','var'),~exist('placeboEuclMaxSolution','var')),or(~exist('drugEuclMaxSolution','var'),~exist('drugEuclMinSolution','var')))
        error('For ROOM fastSL algorithm, one of the following, variables is missing: placeboEuclMinSolution placeboEuclMaxSolution drugEuclMaxSolution drugEuclMinSolution')
    end
    if ~exist('threshold','var') ||isempty(threshold)
        disp('No threshold set for ROOM fastSL, default threshold will be used...threshold = .5')
        threshold = 0.8;
    end
end


%Wildtype FBA solution
%Step1 Identify Single Lethal Reactions...
switch method 
    
    case 1
        disp('Using FastSL original algorithm')
        fluxPlaceboWT    = optimizeCbModel(modelPlacebo);
        fluxPlaceboWT    = fluxPlaceboWT.x;
        fluxDrugWT       = optimizeCbModel(modelDrug);
        fluxDrugWT       = fluxDrugWT.x;
        originalDistance = norm(fluxDrugWT - fluxPlaceboWT);
        %Recovers efficient single targets to not use them again
        [efficientSingleTargets, ~] = singleSLCovid(method,targets, modelDrug,modelPlacebo,placeboEuclMinSolution,placeboEuclMaxSolution,drugEuclMaxSolution,drugEuclMinSolution,threshold);
        targets      = targets(~ismember(targets,efficientSingleTargets));
        
    case 2 
        origDistMaxSolutions = norm(placeboEuclMaxSolution - drugEuclMaxSolution);
        acceptedDistMax = origDistMaxSolutions * threshold;
        origDistMinSolutions = norm(placeboEuclMinSolution - drugEuclMinSolution);
        acceptedDistMin = origDistMinSolutions * threshold;
        %Recovers efficient single targets to not use them again
        [efficientSingleTargets, ~] = singleSLCovid(method,targets, modelDrug,modelPlacebo,placeboEuclMinSolution,placeboEuclMaxSolution,drugEuclMaxSolution,drugEuclMinSolution,threshold);
        efficientSingleTargets = find(ismember(modelDrug.rxns,efficientSingleTargets));
         %Eliminate those rxns which already have been classified as essential rxns
         %in singleSL
        %targets = targets(~ismember(targets,efficientSingleTargets));
        
end

                                                
effectiveTargetPair = [];
effectiveTargetTriplet = [];
%%
switch method
    case 1 
disp('Original fastSL algortihm')
h = waitbar(0,'0.00','Name','Identifying effective target pairs & effective target triplets - Part 1 of 1...');
modelSL = modelPlacebo;

for iRxn = 1:length(targets)
    currentSingleSL = targets(iRxn);
    modelSL.lb(currentSingleSL) = 0; 
    modelSL.ub(currentSingleSL) = 0;
    targetPairs = targets(~ismember(targets,currentSingleSL));
    for jRxn = 1:length(targetPairs)
            currentPotentialPair  = targetPairs(jRxn);
            if currentPotentialPair == currentSingleSL
                continue
            end
            idxsPair = [currentSingleSL currentPotentialPair];
             %Depending on the format of the variable 'targets' this could
            %change a little.
                modelSL.lb(currentPotentialPair) = 0;
                modelSL.ub(currentPotentialPair) = 0;
            testPairSL            = optimizeCbModel(modelSL);
            currentDist2WT        = norm(fluxDrugWT - testPairSL.x);
            
            if (originalDistance > currentDist2WT)    %  || isnan(solKO_i.f))
                effectiveTargetPair = [effectiveTargetPair;idxsPair]; %#ok<*AGROW>
            end
            %Cutoff of 10% for distance, not continue looking for triplets.
            %If the distance reduced is more than 90%, then not continue
            %looking for triplets
                if currentDist2WT <= originalDistance*.1
                    modelSL = modelPlacebo;
                    continue
                end
            %Finding triplets
                targetsTriplets = targetPairs(~ismember(targetPairs,idxsPair));
                
                for kRxn = 1:length(targetsTriplets)                   
                    currentPotentialTripletRxnidx = targetsTriplets(kRxn);
                    idxsTriplet                   = [currentSingleSL currentPotentialPair currentPotentialTripletRxnidx];
                    modelSL.lb(currentPotentialTripletRxnidx) = 0;
                    modelSL.ub(currentPotentialTripletRxnidx) = 0;
                    testTripletSL = optimizeCbModel(modelSL);
                    currentDist2WT        = norm(fluxDrugWT - testTripletSL.x);
                    if (originalDistance > currentDist2WT)    %  || isnan(solKO_i.f))
                        effectiveTargetTriplet = [effectiveTargetTriplet;idxsTriplet]; %#ok<*AGROW>
                    end
                    %Reset bounds on idx reactION
                    modelSL.lb(currentPotentialTripletRxnidx) = model.lb(currentPotentialTripletRxnidx);
                    modelSL.ub(currentPotentialTripletRxnidx) = model.ub(currentPotentialTripletRxnidx);
                    
                end
            
            %Reset bounds on idx reaction
            modelSL.lb(currentPotentialPair) = model.lb(currentPotentialPair);
            modelSL.ub(currentPotentialPair) = model.ub(currentPotentialPair);
    end
    
    modelSL.lb(currentSingleSL) = model.lb(currentSingleSL);
    modelSL.ub(currentSingleSL) = model.ub(currentSingleSL);
    waitbar(iRxn/length(targets),h,[num2str(round(iRxn*100/length(targets))) '% completed...']);

end
close(h);


%Eliminate double lethal reaction deletions in triple lethal reactions
%IS THIS REALLY NECESSARY IN THE COVID? THERE ARE FEW TARGETS, SO PROBABLY
%SOME OF THEM WILL APPEAR IN BOTH, TRIPLETS AND PAIRS.
temporary = [];
g = zeros(1,length(effectiveTargetPair));
for iRxn = 1:length(effectiveTargetTriplet)
    for jRxn = 1:length(effectiveTargetPair)
        g(jRxn) = sum(ismember(effectiveTargetTriplet(iRxn,:),effectiveTargetPair(jRxn,:)));
        if g(jRxn)>= 2
            break;
        end
    end
    if max(g)<2
        temporary=[temporary;effectiveTargetTriplet(iRxn,:)];
    end
end
effectiveTargetTriplet=temporary;

%Eliminate duplicates in triple reaction deletions
effectiveTargetTriplet=unique(sort(effectiveTargetTriplet,2),'rows');



efficientSingleTargets  = model.rxns(efficientSingleTargets);
effectiveTargetPair     = model.rxns(effectiveTargetPair);
effectiveTargetTriplet  = model.rxns(effectiveTargetTriplet);



    case 2
%% Applying ROOM
disp('Identifying Identifying effective target pairs & effective target triplets - Part 1 of 1...')

%Required for making parfor work
environment                       = getEnvironment();
SLTripletsIndexPart1              = cell(100000,1);
SLPairsIndexPart1                 = cell(100000,1);
numberOfPotentialTargets          = length(targetsIdx);

%Main Loop
for iRxn = 1:numberOfPotentialTargets
    targetsParfor = targetsIdx;
   %Restoring environment for making parfor work
    restoreEnvironment(environment);
    modelKO = modelPlacebo; 
    %CurrentSingleKO and get potential pairs WITHIN TARGETS
        currentSingleSLidx    = targetsParfor(iRxn,1);
        currentSingleSL       = modelKO.rxns(currentSingleSLidx);
        potentialPairRxns     = targetsParfor(~ismember(targetsParfor,targetsParfor(1:iRxn)));
    %For to loop through every potential target.
    %Perform a ROOM and compare the euclidean distance between
    %drugflux-ROOM to drugflux-fluxWithNoKOs
    %If the distance is less, then it is a essential target
       for jRxn = 1:length(potentialPairRxns)
            currentPotentialPairIdx     = potentialPairRxns(jRxn);
            currentPotentialPair        = modelPlacebo.rxns(currentPotentialPairIdx);
            tempPair = [currentSingleSL  currentPotentialPair];
            [~,~,currentDistMin]     = ROOMtolerance(modelKO, drugEuclMinSolution,tempPair,'delta', 0.1,'epsilon', 0.01,'printLevel',0);
            [~,~,currentDistMax]     = ROOMtolerance(modelKO, drugEuclMaxSolution,tempPair,'delta', 0.1,'epsilon', 0.01,'printLevel',0);
            
            if acceptedDistMax > currentDistMax && acceptedDistMin > currentDistMin
                    SLPairsIndexPart1{iRxn,1}{jRxn,1} = [currentSingleSLidx currentPotentialPairIdx];
            end
   
            %Finding triplets
            %Loop through potential triplets and repeat procedure as in
            %pairs.
                targetsTriplets = potentialPairRxns(~ismember(potentialPairRxns,targetsParfor(1:iRxn)));
                targetsTriplets = targetsTriplets(~ismember(targetsTriplets,potentialPairRxns(1:jRxn)));
                tempArraySLTriplets = zeros(length(targetsTriplets),3);
             
                for kRxn = 1:length(targetsTriplets)   
                    modelKO = modelPlacebo;
                    currentPotentialTripletRxnidx = targetsTriplets(kRxn);
                    currentPotentialTriplet       = modelKO.rxns(currentPotentialTripletRxnidx);
                    tempTriplet = [tempPair currentPotentialTriplet];
                    [~,~,currentDistMin] = ROOMtolerance(modelKO, drugEuclMinSolution,tempTriplet,'delta', 0.1,'epsilon', 0.01,'printLevel',0);
                    [~,~,currentDistMax] = ROOMtolerance(modelKO, drugEuclMaxSolution,tempTriplet,'delta', 0.1,'epsilon', 0.01,'printLevel',0);        
                    if acceptedDistMax > currentDistMax && acceptedDistMin > currentDistMin
                         tempIdxs = [currentSingleSLidx currentPotentialPairIdx currentPotentialTripletRxnidx];
                         tempArraySLTriplets(kRxn,:)  =  tempIdxs;
                    end
                end
                SLTripletsIndexPart1{iRxn,1}{jRxn,1} = tempArraySLTriplets;
       end
       
       disp(['Finished analyzing reaction number ' num2str(iRxn) ' out of ' num2str(numberOfPotentialTargets)])  
       %waitbar(iRxn/length(potentialSynLethals),h,[num2str(round(iRxn*100/length(potentialSynLethals))) '% completed...']);
end




%Delete empty cells
indexEmptyCellsPart1                    = cellfun(@isempty,SLPairsIndexPart1);
SLPairsIndexPart1(indexEmptyCellsPart1) = [];
SLPairsFinal                             = SLPairsIndexPart1;
%check again empty cells, just in case any of those was empty
indexEmptyCellsFinal               = cellfun(@isempty,SLPairsFinal);
SLPairsFinal(indexEmptyCellsFinal) = [];

%Pass cell array to a number format
for i = 1:size(SLPairsFinal,1)
    currentSLPairs = SLPairsFinal{i,1};
    currentSLPairs = cell2mat(currentSLPairs);
    if i == 1 
        effectiveTargetPair = currentSLPairs;
        continue
    end
    effectiveTargetPair = [effectiveTargetPair; currentSLPairs];
end
effectiveTargetPair = effectiveTargetPair(find(any(effectiveTargetPair,2)),:);

%Delete empty cells
    indexEmptyCellsPart1                       = cellfun(@isempty,SLTripletsIndexPart1);
    SLTripletsIndexPart1(indexEmptyCellsPart1) = [];
    SLTripletsCell                              = SLTripletsIndexPart1;
%check again empty cells, just in case any of those was empty
    indexEmptyCellsFinal                 = cellfun(@isempty,SLTripletsCell);
    SLTripletsCell(indexEmptyCellsFinal) = [];
    
%Pass cell array to a number format
for i = 1:size(SLTripletsCell,1)
    currentSetOfTripelts = SLTripletsCell{i,1};
    currentSetOfTripelts = cell2mat(currentSetOfTripelts);
    if i == 1
        effectiveTargetTriplet = currentSetOfTripelts;
        continue
    end
   effectiveTargetTriplet = [effectiveTargetTriplet; currentSetOfTripelts];
end
effectiveTargetTriplet = effectiveTargetTriplet(find(any(effectiveTargetTriplet,2)),:);

%Eliminate double lethal reaction deletions in triple lethal reactions
%IS THIS REALLY NECESSARY IN THE COVID? THERE ARE FEW TARGETS, SO PROBABLY
%SOME OF THEM WILL APPEAR IN BOTH, TRIPLETS AND PAIRS.
temporary = [];
g = zeros(1,length(effectiveTargetPair));
for iRxn = 1:size(effectiveTargetTriplet,1)
    for jRxn = 1:size(effectiveTargetPair,1)
        g(jRxn) = sum(ismember(effectiveTargetTriplet(iRxn,:),effectiveTargetPair(jRxn,:)));
        if g(jRxn)>= 2
            break;
        end
    end
    if max(g)<2
        temporary=[temporary;effectiveTargetTriplet(iRxn,:)]; %#ok<AGROW>
    end
end
effectiveTargetTriplet = temporary;

%Eliminate duplicates in triple reaction deletions
effectiveTargetTriplet = unique(sort(effectiveTargetTriplet,2),'rows');



efficientSingleTargets = modelDrug.rxns(efficientSingleTargets);
effectiveTargetPair    = modelDrug.rxns(effectiveTargetPair);
effectiveTargetTriplet = modelDrug.rxns(effectiveTargetTriplet);
end
end
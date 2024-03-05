function [efficientSingleTargets,effectiveTargetPair,effectiveTargetTriplet] = fastSLROOMCOVID(targets,modelDrug,modelPlacebo, placeboCentralSolution,placeboFurthestSolution,drugFurthestSolution,drugCentralSolution,threshold)
%This finds the most efficient drug targets, when comparing the metabolism of
%a person treated with drug and a person treated with a placebo. 
%The analysis is done in singles, pairs, and triplets.
%This function looks for those lethals that reduce the euclidean distance, between the
%placebo flux distribution and the drug flux distribution, the most.  For
%doing this it uses the central and furthest solution as inputs for the
%ROOM.
% targets: Cell array containing three cells. Each cell will include a
    % vertical cell array corresponding to the reactions that form
    % single lethals, pair lethals, and triplet lethals.
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
%Threshold: Value that will help to determine which distance is enough to consider a reaction as essential            
            
%OUTPUT
% EffectiveSingleTargets        Targets that manage to reduce the euclidean
            % distance between the drugWT and the original Placebo WT
% EffectivePairTargets        Indices of pair of reactions  that manage to reduce the euclidean
            % distance between the drugWT and the original Placebo WT
% EffectiveTripletsTargets        Indices of triplets of reactions  that manage to reduce the euclidean
            % distance between the drugWT and the original Placebo WT

%%
%The input could be reaction indexes but also just the reaction names
if iscell(targets)
    targetsIdx = find(ismember(modelPlacebo.rxns,targets));
end

if or(or(~exist('placeboCentralSolution','var'),~exist('placeboFurthestSolution','var')),or(~exist('drugFurthestSolution','var'),~exist('drugCentralSolution','var')))
    error('For ROOM fastSL algorithm, one of the following, variables is missing: placeboCentralSolution placeboFurthestSolution drugFurthestSolution drugCentralSolution')
end
if ~exist('threshold','var') ||isempty(threshold)
     disp('No threshold set for ROOM fastSL, default threshold will be used...threshold = .2')
     threshold = 0.2;
end

%% Step1 Identify Single Lethal Reactions...
 %Compute original distances
      origDistFurthestSolutions = norm(placeboFurthestSolution - drugFurthestSolution);
      acceptedDistFurthest      = origDistFurthestSolutions * threshold;
      origDistCentralSolutions  = norm(placeboCentralSolution - drugCentralSolution);
      acceptedDistCentral       = origDistCentralSolutions * threshold;
 %Recovers efficient single targets to not use them again
      [efficientSingleTargets, ~] = singleSLCovid(targets, modelDrug,modelPlacebo,placeboCentralSolution,placeboFurthestSolution,drugFurthestSolution,drugCentralSolution,threshold);
      efficientSingleTargets      = find(ismember(modelDrug.rxns,efficientSingleTargets));
      %Eliminate those rxns which already have been classified as essential rxns
      %in singleSL
      %CHECK THIS
        targets = targets(~ismember(targets,efficientSingleTargets));
       
                                                
effectiveTargetPair = [];
effectiveTargetTriplet = [];

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
            [~,~,currentDistMin]     = ROOMtolerance(modelKO, drugCentralSolution,tempPair,'delta', 0.1,'epsilon', 0.01,'printLevel',0);
            [~,~,currentDistMax]     = ROOMtolerance(modelKO, drugFurthestSolution,tempPair,'delta', 0.1,'epsilon', 0.01,'printLevel',0);
            
            if acceptedDistFurthest > currentDistMax && acceptedDistCentral > currentDistMin
                    SLPairsIndexPart1{iRxn,1}{jRxn,1} = [currentSingleSLidx currentPotentialPairIdx];
            end
   
            %Finding triplets
            %Loop through potential triplets and repeat procedure as in
            %pairs.
                targetsTriplets = potentialPairRxns(~ismember(potentialPairRxns,targetsParfor(1:iRxn)));
                targetsTriplets = targetsTriplets(~ismember(targetsTriplets,potentialPairRxns(1:jRxn)));
                tempArraySLTriplets = zeros(length(targetsTriplets),3);
             
                parfor kRxn = 1:length(targetsTriplets)   
                    modelKO = modelPlacebo;
                    currentPotentialTripletRxnidx = targetsTriplets(kRxn);
                    currentPotentialTriplet       = modelKO.rxns(currentPotentialTripletRxnidx);
                    tempTriplet = [tempPair currentPotentialTriplet];
                    [~,~,currentDistMin] = ROOMtolerance(modelKO, drugCentralSolution,tempTriplet,'delta', 0.1,'epsilon', 0.01,'printLevel',0);
                    [~,~,currentDistMax] = ROOMtolerance(modelKO, drugFurthestSolution,tempTriplet,'delta', 0.1,'epsilon', 0.01,'printLevel',0);        
                    if acceptedDistFurthest > currentDistMax && acceptedDistCentral > currentDistMin
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
SLPairsFinal                            = SLPairsIndexPart1;
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

%CHECK THIS
%Eliminate double lethal reaction deletions in triple lethal reactions
%IS THIS REALLY NECESSARY IN THE COVID? THERE ARE FEW TARGETS, SO PROBABLY
%SOME OF THEM WILL APPEAR IN BOTH, TRIPLETS AND PAIRS.
%temporary = [];
%g = zeros(1,length(effectiveTargetPair));
%for iRxn = 1:size(effectiveTargetTriplet,1)
%    for jRxn = 1:size(effectiveTargetPair,1)
%        g(jRxn) = sum(ismember(effectiveTargetTriplet(iRxn,:),effectiveTargetPair(jRxn,:)));
%        if g(jRxn)>= 2
%            break;
%        end
%    end
%    if max(g)<2
%        temporary=[temporary;effectiveTargetTriplet(iRxn,:)]; %#ok<AGROW>
%    end
%end
%effectiveTargetTriplet = temporary;

%Eliminate duplicates in triple reaction deletions
effectiveTargetTriplet = unique(sort(effectiveTargetTriplet,2),'rows');


efficientSingleTargets = modelDrug.rxns(efficientSingleTargets);
effectiveTargetPair    = modelDrug.rxns(effectiveTargetPair);
effectiveTargetTriplet = modelDrug.rxns(effectiveTargetTriplet);
end
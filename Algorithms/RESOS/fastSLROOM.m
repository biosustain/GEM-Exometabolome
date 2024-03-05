function [Jsl,Jdl,Jtl] = fastSLROOM(model,cutoff,eliList,furthestSolution,centralSolution)
%THIS FUNCTION WAS NOT USED IN THE FINAL PIPELINE, THE FUNCTION USED THERE
%WAS fastSLROOM2

%This function finds lethals based on the biomass reduction. It differs
%from the traditional fastSL as this one uses the furthest and central
%solution and the ROOM algorithm. 
    %The reduction in the biomass must be lower than the cutoff in both
    %fluxes produced by ROOM (furthest and central solutions)
%This function also eliminates the lethals that produce aditive effect
%Works but it is too slow

% INPUT
% model (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
%   rxns         Reaction Names
%OPTIONAL
% cutoff         cutoff percentage value for lethality.Default is 0.01.
% eliList        List of reactions to be ignored for lethality
    % analysis:Exchange Reactions, ATPM etc.
% furthestSolution = After a sampling and percentile filtration, this
% solution is the one that has the maximum euclidean distance value within the
% valid solutions.
% centralSolution = After a sampling and percentile filtration, this
% solution is the one that has minimum euclidean distance value within the
% valid solutions.
            
%OUTPUT
% Jsl        Indices of single lethal reactions identified
% Jdl        Indices of double lethal reactions identified
% Jtl        Indices of triple lethal reactions identified
% Fernando Silva-Lance 2021
%%
if ~exist('furthestSolution','var') || ~exist('centralSolution','var') || isempty(furthestSolution) || isempty(centralSolution)  
        error('Missing variables: furthestSolution or centralSolution')
end
if exist('cutoff', 'var')  || isempty(cutoff)
        cutoff = 0.01;
end

if exist('eliList', 'var') || isempty(eliList)
        eliList = model.rxns(find(model.c)); 
end


        nonZeroRxnsIdx = find(~eq(centralSolution, 0) & ~eq(furthestSolution, 0));
        BioMassIndex = find(model.c);
        BioMassMinSolution = centralSolution(find(model.c)); %#ok<*FNDSB>
        BioMassMaxSolution = furthestSolution(find(model.c));
        if (~isempty(eliList))
            eliIdx = find(ismember(model.rxns,eliList)); %Index of reactions not considered for lethality analysis
            potentialSynLethals = nonZeroRxnsIdx(~ismember(nonZeroRxnsIdx,eliIdx)); %nonZeroFluxes
        else
            potentialSynLethals = nonZeroRxnsIdx;
        end
        [Jsl, BMSinglePercentReduction] = singleSLModified(model,cutoff,eliList,furthestSolution,centralSolution);
         %Eliminate those rxns which already have been classified as essential rxns
         %in singleSL
        potentialSynLethals = potentialSynLethals(~ismember(potentialSynLethals,Jsl));
        BMSinglePercentReduction = BMSinglePercentReduction(ismember(BMSinglePercentReduction(:,1),potentialSynLethals),:);
        
        %Discard Pairs or Triplets that present aditive effect
        counterPairs    = 1;
        counterTriplets = 1;
        for i = 1:size(BMSinglePercentReduction,1)
            currentSingleKOidx = BMSinglePercentReduction(i,1);
            currentSingleKOmin = BMSinglePercentReduction(i,2);
            currentSingleKOmax = BMSinglePercentReduction(i,3);
            for j = 1:size(BMSinglePercentReduction,1)
                if i <= j
                    continue
                end
                currentPotentialPairidx = BMSinglePercentReduction(j,1);
                currentPotentialPairMin = BMSinglePercentReduction(j,2);
                currentPotentialPairMax = BMSinglePercentReduction(j,3);
        %CheckPair for aditive effect and put in discarded SLPairs
                if (currentSingleKOmin*currentPotentialPairMin <= cutoff)&& (currentSingleKOmax*currentPotentialPairMax <= cutoff)
                    discardedSLPairs(counterPairs,:) = [currentSingleKOidx currentPotentialPairidx];
                    counterPairs = counterPairs + 1;
                continue
                end
                for k = 1:size(BMSinglePercentReduction,1)
                    if k <= i || k <= j
                    continue
                    end
                    currentPotentialTripletIdx       = BMSinglePercentReduction(k,1);
                    currentPotentialTripletMin = BMSinglePercentReduction(k,2);
                    currentPotentialTripletMax = BMSinglePercentReduction(k,3);
            %CheckTriplet for aditive effect and put in discarded
            %SLTriplets
                    if (currentSingleKOmin*currentPotentialPairMin*currentPotentialTripletMin <= cutoff)&& (currentSingleKOmax * currentPotentialPairMax * currentPotentialTripletMax<= cutoff)
                        discardedSLTriplets(counterTriplets,:) = [currentSingleKOidx currentPotentialPairidx currentPotentialTripletIdx];
                        counterTriplets = counterTriplets + 1;
                    end
                end
            end
        end
        
if ~exist('discardedSLTriplets', 'var')
    discardedSLTriplets = [0 0 0];
end
if ~exist('discardedSLPairs', 'var')
    discardedSLPairs = [0 0 0];
end

                                                
Jdl=[];
Jtl=[];

%% Applying ROOM
disp('Identifying Pairs & Triplets - Part 1 of 1...')
environment = getEnvironment();
SLTripletsIndexPart1              = cell(100000,1);
SLPairsIndexPart1                 = cell(100000,3);
numberOfPotentialSynLethals       = length(potentialSynLethals);
usedRxns =  [];
for iRxn = 1:numberOfPotentialSynLethals
   %Restoring environment for making parfor work
    modelKO = model; 
    currentSingleSLidx    = potentialSynLethals(iRxn);
    currentSingleSL       = modelKO.rxns(currentSingleSLidx);
    [fluxMinSingleSL,~,~] = ROOMtolerance(modelKO, centralSolution,currentSingleSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
    [fluxMaxSingleSL,~,~] = ROOMtolerance(modelKO, furthestSolution,currentSingleSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
    %Fluxes different from 0 in the results could be pairs that could lower
    %the biomass value
    discardRxns           = find(fluxMinSingleSL == 0 | fluxMaxSingleSL == 0);
    potentialPairRxns     = potentialSynLethals(~ismember(potentialSynLethals,discardRxns));

    if (~isempty(eliList))
        potentialPairRxns = potentialPairRxns(~ismember(potentialPairRxns,eliIdx));
    end
    
    potentialPairRxns = potentialPairRxns(~ismember(potentialPairRxns,potentialSynLethals(1:iRxn)));
    
       for jRxn = 1:length(potentialPairRxns)
            currentPotentialPairIdx     = potentialPairRxns(jRxn);
            idxsPair                    = [currentSingleSLidx currentPotentialPairIdx];
            matchPairWithDiscardedPairs = sum(ismember(discardedSLPairs,idxsPair),2);
            if ~isempty(matchPairWithDiscardedPairs(matchPairWithDiscardedPairs >= 2))
                continue
            end
            currentPotentialPair    = modelKO.rxns(currentPotentialPairIdx);
            tempPairSL = [currentSingleSL currentPotentialPair];
            [fluxMinPairSL,~,~]     = ROOMtolerance(modelKO, centralSolution,tempPairSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
            [fluxMaxPairSL,~,~]     = ROOMtolerance(modelKO, furthestSolution,tempPairSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
            fluxMinBMreduction      = (fluxMinPairSL(BioMassIndex)/BioMassMinSolution)  ;
            fluxMaxBMreduction      = (fluxMaxPairSL(BioMassIndex)/BioMassMaxSolution)  ;
            if (fluxMinBMreduction + fluxMaxBMreduction) <= 2*cutoff || ((sum(isnan(fluxMinPairSL)) ~= 0) && (sum(isnan(fluxMaxPairSL)) ~= 0))
                SLPairsIndexPart1{iRxn,1}{jRxn,1} = [currentSingleSLidx currentPotentialPairIdx];
                usedRxns = [usedRxns;currentSingleSLidx currentPotentialPairIdx];
                continue
            elseif (sum(isnan(fluxMinBMreduction)) ~= 0 && fluxMaxBMreduction <= cutoff) ||  (sum(isnan(fluxMaxBMreduction)) ~= 0 && fluxMinBMreduction <= cutoff)
                SLPairsIndexPart1{iRxn,1}{jRxn,1} = [currentSingleSLidx currentPotentialPairIdx];
                usedRxns = [usedRxns;currentSingleSLidx currentPotentialPairIdx];
                continue
            end
            
            %Finding triplets
                discardRxns = find((fluxMaxPairSL == 0) | (fluxMinPairSL == 0));
                potentialTriplets   = potentialSynLethals(~ismember(potentialSynLethals,discardRxns));
                
                if (~isempty(eliList))
                    potentialTriplets = potentialTriplets(~ismember(potentialTriplets,eliIdx)); 
                end
                if (~isempty(usedRxns))
                    potentialTriplets = potentialTriplets(~ismember(potentialTriplets,usedRxns(:,1)));
                    potentialTriplets = potentialTriplets(~ismember(potentialTriplets,usedRxns(:,2)));
                end
                
               %Eliminate Aditive effect triplets
                if (~isempty(discardedSLTriplets))
                    discardAditiveEffect = zeros(length(discardedSLTriplets),1);
                    for i =1:length(potentialTriplets) 
                        currentPotentialTripletRxnidx = potentialTriplets(i);
                        idxsTriplet                       = [currentSingleSLidx currentPotentialPairIdx currentPotentialTripletRxnidx];
                        matchTripletWithDiscardedTriplets = sum(ismember(discardedSLTriplets,idxsTriplet),2);
                        if ~isempty(matchTripletWithDiscardedTriplets(matchTripletWithDiscardedTriplets >= 3))
                            discardAditiveEffect(i,1) = currentPotentialTripletRxnidx;
                         continue
                        end
                    end
                    if (~isempty(discardAditiveEffect))
                        potentialTriplets = potentialTriplets(~ismember(potentialTriplets,discardAditiveEffect));
                    end
                end
                
                potentialTriplets = potentialTriplets(~ismember(potentialTriplets,potentialSynLethals(1:iRxn)));
                potentialTriplets = potentialTriplets(~ismember(potentialTriplets,potentialPairRxns(1:jRxn)));
  
                tempArraySLTriplets = zeros(length(potentialTriplets),3);
                numberOfPotentialTriplets = length(potentialTriplets) ;
                parfor kRxn = 1:numberOfPotentialTriplets
                    restoreEnvironment(environment);
                    modelKO1 = model;
                    currentPotentialTripletRxnidx = potentialTriplets(kRxn);
                    currentPotentialTriplet = modelKO1.rxns(currentPotentialTripletRxnidx);
                    tempTripletSL           = [tempPairSL currentPotentialTriplet];
                    [fluxMinTripletSL,~,~]  = ROOMtolerance(modelKO1, centralSolution,tempTripletSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
                    [fluxMaxTripletSL,~,~]  = ROOMtolerance(modelKO1, furthestSolution,tempTripletSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
                    fluxMinBMreduction  = (fluxMinTripletSL(BioMassIndex)/BioMassMinSolution)  ;
                    fluxMaxBMreduction  = (fluxMaxTripletSL(BioMassIndex)/BioMassMaxSolution)  ;         
                    if or((fluxMinBMreduction + fluxMaxBMreduction) <= 2*cutoff,((sum(isnan(fluxMinTripletSL)) ~= 0) && (sum(isnan(fluxMaxTripletSL)) ~= 0)))
                         tempIdxs = [currentSingleSLidx currentPotentialPairIdx currentPotentialTripletRxnidx];
                         tempArraySLTriplets(kRxn,:)  =  tempIdxs;
                    elseif (sum(isnan(fluxMinBMreduction)) ~= 0 && fluxMaxBMreduction <= cutoff) ||  (sum(isnan(fluxMaxBMreduction)) ~= 0 && fluxMinBMreduction <= cutoff)
                         tempIdxs = [currentSingleSLidx currentPotentialPairIdx currentPotentialTripletRxnidx];
                         tempArraySLTriplets(kRxn,:)  =  tempIdxs;
                    end
                    %disp(kRxn)
                end
                SLTripletsIndexPart1{iRxn,1}{jRxn,1} = tempArraySLTriplets;
       end
       
       disp(['Finished analyzing reaction number ' num2str(iRxn) ' out of ' num2str(numberOfPotentialSynLethals)])  
       %waitbar(iRxn/length(potentialSynLethals),h,[num2str(round(iRxn*100/length(potentialSynLethals))) '% completed...']);
end


% Join pairs indexes
    indexEmptyCellsPart1 = cellfun(@isempty,SLPairsIndexPart1);
    SLPairsIndexPart1(indexEmptyCellsPart1) = [];
    SLPairsFinal = [SLPairsIndexPart1']; 
%check again empty cells, just in case any of those was empty
    indexEmptyCellsFinal = cellfun(@isempty,SLPairsFinal);
    SLPairsFinal(indexEmptyCellsFinal) = [];
for i = 1:size(SLPairsFinal,1)
    currentSLPairs = SLPairsFinal{i,1};
    currentSLPairs = cell2mat(currentSLPairs);
    if i == 1 
        Jdl = currentSLPairs;
        continue
    end
    Jdl = [Jdl; currentSLPairs];
end
Jdl = Jdl(find(any(Jdl,2)),:);

%Join SLTripletsIndexPart1 and SLTripletsIndexPart2
    indexEmptyCellsPart1                       = cellfun(@isempty,SLTripletsIndexPart1);
    SLTripletsIndexPart1(indexEmptyCellsPart1) = [];
    SLTripletsCell                             = [SLTripletsIndexPart1];
%check again empty cells, just in case any of those was empty
    indexEmptyCellsFinal                 = cellfun(@isempty,SLTripletsCell);
    SLTripletsCell(indexEmptyCellsFinal) = [];
for i = 1:size(SLTripletsCell,1)
    currentSetOfTripelts = SLTripletsCell{i,1};
    currentSetOfTripelts = cell2mat(currentSetOfTripelts);
    if i == 1
        Jtl = currentSetOfTripelts;
        continue
    end
   Jtl = [Jtl; currentSetOfTripelts];
end
Jtl = Jtl(find(any(Jtl,2)),:);

%Eliminate double lethal reaction deletions in triple lethal reactions
temporary = [];
g = zeros(1,length(Jdl));
for iRxn = 1:length(Jtl)
    for jRxn = 1:length(Jdl)
        g(jRxn) = sum(ismember(Jtl(iRxn,:),Jdl(jRxn,:)));
        if g(jRxn)>= 2
            break;
        end
    end
    if max(g)<2
        temporary=[temporary;Jtl(iRxn,:)]; %#ok<AGROW>
    end
end
Jtl = temporary;

%Eliminate duplicates in triple reaction deletions
Jtl=unique(sort(Jtl,2),'rows');



Jsl = model.rxns(Jsl);
Jdl = model.rxns(Jdl);
Jtl = model.rxns(Jtl);
save('BackupTripleSL')
end
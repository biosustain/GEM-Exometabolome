function [Jsl,Jdl,Jtl] = tripleSLModified(model,cutoff,eliList,atpm, method,euclMaxSolution,euclMinSolution)
%%  [slist_id,dlist_id,tlist_id]=tripleSL(model,cutoff,eliList,atpm)
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
% is true.
% atpm           ATPM Reaction Id in model.rxns if other than 'ATPM'
% method
    %1 = Original fastSL algorithm
    %2 =  ROOM (Requires two extra inputs: maxSolution minSolution)
% maxSolution = After a sampling and euclidean distance filtration, this
% solution is the one that has the maximum euclidean distance value within the
% valid solutions.
% minSolution = After a sampling and euclidean distance filtration, this
% solution is the one that has minimum euclidean distance value within the
% valid solutions.
            
%OUTPUT
% Jsl        Indices of single lethal reactions identified
% Jdl        Indices of double lethal reactions identified
% Jtl        Indices of triple lethal reactions identified
% Aditya Pratapa       7/1/14.
%%
if ~exist('method','var') ||  isempty(method)
    method = 1;
    disp('Variable "method" not set. Using traditional fastSL algorithm...')
end

if or(method == 2, method == 3)
    if ~exist('euclMaxSolution','var') || ~exist('euclMinSolution','var') || isempty(euclMaxSolution) || isempty(euclMinSolution)  
        error('Missing variables: euclMaxSolution or euclMinSolution')
    end
else
    euclMinSolution = 0;
    euclMaxSolution = 0;
end


if exist('cutoff', 'var')
    if isempty(cutoff)
        cutoff = 0.01;
    end
else
    cutoff = 0.01;
end


if exist('atpm', 'var')
    if isempty(atpm)
        atpm = 'ATPM'; %Reaction Id of ATP maintenance reaction- by default it takes 'ATPM'
    end
else
    atpm = 'ATPM';
end


if exist('eliList', 'var')
    if isempty(eliList)
        eliList = model.rxns(ismember(model.rxns,atpm)); %To eliminate ATPM.
    end
else
    eliList = model.rxns(ismember(model.rxns,atpm));
end

%Wildtype FBA solution
%Step1 Identify Single Lethal Reactions...
switch method 
    
    case 1
        %Identify minNorm flux distribution
        solWT             = optimizeCbModel(model,'max','one');
        ObjectiveValueWT  = solWT.f;
        %Reactions with fluxes different to 0
        nonZeroRxnsIdx    = find(~eq(solWT.x,0));
        if (~isempty(eliList))
            eliIdx = find(ismember(model.rxns,eliList)); %Index of reactions not considered for lethality analysis
            potentialSynLethals    = nonZeroRxnsIdx(~ismember(nonZeroRxnsIdx,eliIdx)); %nonZeroFluxes
        else
            potentialSynLethals = nonZeroRxnsIdx;
        end
        [Jsl, BMSinglePercentReduction] = singleSLModified(model,cutoff,eliList,[],method,euclMaxSolution,euclMinSolution);
        potentialSynLethals      = potentialSynLethals(~ismember(potentialSynLethals,Jsl));
        BMSinglePercentReduction = BMSinglePercentReduction(ismember(BMSinglePercentReduction(:,1),potentialSynLethals),:);
        %discardedSLPairs
        counterPairs = 1;
        counterTriplets = 1;
        for i = 1:size(BMSinglePercentReduction,1)
            currentSingleKO = BMSinglePercentReduction(i,2);
            singleKOIdx     = BMSinglePercentReduction(i,1);
            for j = 1:size(BMSinglePercentReduction,1)
                if i <= j
                    continue
                end
                currentPotentialPair = BMSinglePercentReduction(j,2);
                potentialPairIdx     = BMSinglePercentReduction(j,1);
        %CheckPair for aditive effect and put in discarded SLPairs
                if (currentSingleKO*currentPotentialPair <= cutoff)
                    discardedSLPairs(counterPairs,:) = [singleKOIdx potentialPairIdx];
                    counterPairs = counterPairs + 1;
                continue
                end
                for k = 1:size(BMSinglePercentReduction,1)
                    if k <= i || k <= j
                    continue
                    end
                    currentpotentialTripletIdx = BMSinglePercentReduction(k,1);
                    curretPotentialTriplet = BMSinglePercentReduction(k,2);
            %CheckTriplet for aditive effect and put in discarded
            %SLTriplets
                    if (currentSingleKO*currentPotentialPair*curretPotentialTriplet <= cutoff)
                        discardedSLTriplets(counterTriplets,:) = [singleKOIdx potentialPairIdx currentpotentialTripletIdx];
                        counterTriplets = counterTriplets + 1;
                    end
                end
            end
        end
        
    case 2 
        nonZeroRxnsIdx = find(~eq(euclMinSolution, 0) & ~eq(euclMaxSolution, 0));
        BioMassIndex = find(model.c);
        BioMassMinSolution = euclMinSolution(find(model.c)); %#ok<*FNDSB>
        BioMassMaxSolution = euclMaxSolution(find(model.c));
        if (~isempty(eliList))
            eliIdx = find(ismember(model.rxns,eliList)); %Index of reactions not considered for lethality analysis
            potentialSynLethals = nonZeroRxnsIdx(~ismember(nonZeroRxnsIdx,eliIdx)); %nonZeroFluxes
        else
            potentialSynLethals = nonZeroRxnsIdx;
        end
        [Jsl, BMSinglePercentReduction] = singleSLModified(model,cutoff,eliList,[],method,euclMaxSolution,euclMinSolution);
         %Eliminate those rxns which already have been classified as essential rxns
         %in singleSL
        potentialSynLethals = potentialSynLethals(~ismember(potentialSynLethals,Jsl));
        BMSinglePercentReduction = BMSinglePercentReduction(ismember(BMSinglePercentReduction(:,1),potentialSynLethals),:);
        
        %discardedSLPairs
        counterPairs = 1;
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
        
end
if ~exist('discardedSLTriplets', 'var')
    discardedSLTriplets = [0 0 0];
end
if ~exist('discardedSLPairs', 'var')
    discardedSLPairs = [0 0 0];
end

                                                
Jdl=[];
Jtl=[];
%%
switch method
    case 1 
disp('Original fastSL algortihm')
h = waitbar(0,'0.00','Name','Identifying Jdl & Jtl - Part 1 of 2...');
modelSL=model;

for iRxn = 1:length(potentialSynLethals)
    currentSingleSL = potentialSynLethals(iRxn);
    modelSL.lb(currentSingleSL) = 0; 
    modelSL.ub(currentSingleSL) = 0;
    solKO_i = optimizeCbModel(modelSL,'max','one'); %It can't be a single lethal so we can proceed further
    potentialPairRxns = find(~eq(solKO_i.x,0));
    potentialPairRxns = potentialPairRxns(~ismember(potentialPairRxns,nonZeroRxnsIdx));
    
    if (~isempty(eliList))
        potentialPairRxns = potentialPairRxns(~ismember(potentialPairRxns,eliIdx)); %Eliminate Exchange and ATP Maintenance reactions
    end
    
    for jRxn = 1:length(potentialPairRxns)
            currentPotentialPair  = potentialPairRxns(jRxn);
            idxsPair                    = [currentSingleSL currentPotentialPair];
            matchPairWithDiscardedPairs = sum(ismember(discardedSLPairs,idxsPair),2);
            if ~isempty(matchPairWithDiscardedPairs(matchPairWithDiscardedPairs >= 2))
                continue
            end
            modelSL.lb(currentPotentialPair) = 0;
            modelSL.ub(currentPotentialPair) = 0;
            testPairSL = optimizeCbModel(modelSL,'max','one');
            
            if (testPairSL.f < cutoff*ObjectiveValueWT && ~eq(testPairSL.stat,0)) 
                Jdl = [Jdl;currentSingleSL currentPotentialPair];
            else
                  if eq(testPairSL.stat,0)
                    testPairSL = optimizeCbModel(modelSL);
                    if (testPairSL.f < cutoff*ObjectiveValueWT || isnan(testPairSL.f)) 
                        Jdl = [Jdl;currentSingleSL currentPotentialPair];
                        %Reset bounds on idx reaction
                            modelSL = model;
                        continue;
                    end
                  end
            %Finding triplets
                %Find those rxns with a flux different to 0 and they must
                %not belong to the initial nonZeroRxns, as the algorithm
                %tries to avoid repeating rxns in singleKO, SLpairs, and
                %SLtriplets
                NonZeroFluxesPairSLidx = find(~eq(testPairSL.x,0));
                NonZeroFluxesPairSLidx = NonZeroFluxesPairSLidx(~ismember(NonZeroFluxesPairSLidx,nonZeroRxnsIdx));
                
                if (~isempty(eliList))
                    NonZeroFluxesPairSLidx = NonZeroFluxesPairSLidx(~ismember(NonZeroFluxesPairSLidx,eliIdx)); %Eliminate Exchange and ATPM reactions
                end
                
                for kRxn = 1:length(NonZeroFluxesPairSLidx)           
                    
                    currentPotentialTripletRxnidx = NonZeroFluxesPairSLidx(kRxn);
                    idxsTriplet = [currentSingleSL currentPotentialPair currentPotentialTripletRxnidx];
                    matchTripletWithDiscardedTriplets = sum(ismember(discardedSLTriplets,idxsTriplet),2);
                    if ~isempty(matchTripletWithDiscardedTriplets(matchTripletWithDiscardedTriplets >= 3))
                        continue
                    end
                    modelSL.lb(currentPotentialTripletRxnidx) = 0;
                    modelSL.ub(currentPotentialTripletRxnidx) = 0;
                    testTripletSL = optimizeCbModel(modelSL);
                    if (testTripletSL.f < cutoff*ObjectiveValueWT || isnan(testTripletSL.f))
                        Jtl = [Jtl;currentSingleSL currentPotentialPair currentPotentialTripletRxnidx];
                    end
                    %Reset bounds on idx reactION
                    modelSL.lb(currentPotentialTripletRxnidx) = model.lb(currentPotentialTripletRxnidx);
                    modelSL.ub(currentPotentialTripletRxnidx) = model.ub(currentPotentialTripletRxnidx);
                    
                end
            end
            %Reset bounds on idx reaction
            modelSL.lb(currentPotentialPair) = model.lb(currentPotentialPair);
            modelSL.ub(currentPotentialPair) = model.ub(currentPotentialPair);
    end
    
    modelSL.lb(currentSingleSL) = model.lb(currentSingleSL);
    modelSL.ub(currentSingleSL) = model.ub(currentSingleSL);
    waitbar(iRxn/length(potentialSynLethals),h,[num2str(round(iRxn*100/length(potentialSynLethals))) '% completed...']);

end
close(h);

%
h = waitbar(0,'0.00','Name','Identifying Jdl & Jtl - Part 2 of 2...');

for iRxn = 1:length(potentialSynLethals)
    for jRxn = 1:length(potentialSynLethals)
        if (jRxn<iRxn)
            modelSL = model;
            currentSingleSL      = potentialSynLethals(iRxn);
            currentPotentialPair = potentialSynLethals(jRxn);
            idxsPair                    = [currentSingleSL currentPotentialPair];
            matchPairWithDiscardedPairs = sum(ismember(discardedSLPairs,idxsPair),2);
            if ~isempty(matchPairWithDiscardedPairs(matchPairWithDiscardedPairs >= 2))
                continue
            end
            modelSL.lb(currentSingleSL) = 0;
            modelSL.ub(currentSingleSL) = 0;
            modelSL.lb(currentPotentialPair) = 0;
            modelSL.ub(currentPotentialPair) = 0;
            %Optimization minimizing taxicab distance (Similar to euclidean
            %but done according to the cartesian coordenades)
            testPairSL = optimizeCbModel(modelSL,'max','one');
            if (testPairSL.f < cutoff*ObjectiveValueWT && ~eq(testPairSL.stat,0))
                Jdl=[Jdl;currentSingleSL currentPotentialPair];
            else
                  if eq(testPairSL.stat,0)
                    testPairSL = optimizeCbModel(modelSL);
                    %Por qué al dar NAN acepta el synthetic lethal?
                    if (testPairSL.f<cutoff*ObjectiveValueWT || isnan(testPairSL.f)) 
                        Jdl=[Jdl;currentSingleSL currentPotentialPair];
                        modelSL.lb(currentPotentialPair)=model.lb(currentPotentialPair);
                        modelSL.ub(currentPotentialPair)=model.ub(currentPotentialPair);
                        continue;
                    end
                   end
                NonZeroFluxesPairSLidx = find(~eq(testPairSL.x,0));
                NonZeroFluxesPairSLidx = NonZeroFluxesPairSLidx(~ismember(NonZeroFluxesPairSLidx,nonZeroRxnsIdx));
                
                if (~isempty(eliList))
                    NonZeroFluxesPairSLidx = NonZeroFluxesPairSLidx(~ismember(NonZeroFluxesPairSLidx,eliIdx)); %Eliminate Exchange and ATPM reactions
                end
                
                for kRxn=1:length(NonZeroFluxesPairSLidx)
                    currentPotentialTripletRxnidx=NonZeroFluxesPairSLidx(kRxn);
                    idxsTriplet = [currentSingleSL currentPotentialPair currentPotentialTripletRxnidx];
                    matchTripletWithDiscardedTriplets = sum(ismember(discardedSLTriplets,idxsTriplet),2);
                    if ~isempty(matchTripletWithDiscardedTriplets(matchTripletWithDiscardedTriplets >= 3))
                        continue
                    end
                    modelSL.lb(currentPotentialTripletRxnidx)=0;
                    modelSL.ub(currentPotentialTripletRxnidx)=0;
                    testTripletSL = optimizeCbModel(modelSL);
                    if (testTripletSL.f<cutoff*ObjectiveValueWT || isnan(testTripletSL.f))
                        Jtl=[Jtl;currentSingleSL currentPotentialPair currentPotentialTripletRxnidx];
                    end
                    
                    modelSL.lb(currentPotentialTripletRxnidx)=model.lb(currentPotentialTripletRxnidx);
                    modelSL.ub(currentPotentialTripletRxnidx)=model.ub(currentPotentialTripletRxnidx);
                    
                end
                
                for kRxn=1:length(potentialSynLethals)
                    
                    if (kRxn<jRxn)
                        currentPotentialTripletRxnidx = potentialSynLethals(kRxn);
                        idxsTriplet                       = [currentSingleSL currentPotentialPair currentPotentialTripletRxnidx];
                        matchTripletWithDiscardedTriplets = sum(ismember(discardedSLTriplets,idxsTriplet),2);
                    if ~isempty(matchTripletWithDiscardedTriplets(matchTripletWithDiscardedTriplets >= 3))
                        continue
                    end
                        modelSL.lb(currentPotentialTripletRxnidx)=0;modelSL.ub(currentPotentialTripletRxnidx)=0;
                        testTripletSL=optimizeCbModel(modelSL);
                        if (testTripletSL.f<cutoff*ObjectiveValueWT ||isnan(testTripletSL.f))
                            Jtl=[Jtl;currentSingleSL currentPotentialPair currentPotentialTripletRxnidx ];
                        end
                        
                        modelSL.lb(currentPotentialTripletRxnidx)=model.lb(currentPotentialTripletRxnidx);
                        modelSL.ub(currentPotentialTripletRxnidx)=model.ub(currentPotentialTripletRxnidx);
                        
                    else
                        break;
                    end
                end
            end
            modelSL.lb(currentPotentialPair)=model.lb(currentPotentialPair);
            modelSL.ub(currentPotentialPair)=model.ub(currentPotentialPair);
        else
            break;
        end
        
    end
    modelSL.lb(currentSingleSL)=model.lb(currentSingleSL);
    modelSL.ub(currentSingleSL)=model.ub(currentSingleSL);
    waitbar(iRxn*(iRxn-1)*(iRxn-2)/(length(potentialSynLethals)*(length(potentialSynLethals)-1)*(length(potentialSynLethals)-2)),h,[num2str(round(iRxn*(iRxn-1)*(iRxn-2)*100/(length(potentialSynLethals)*(length(potentialSynLethals)-1)*(length(potentialSynLethals)-2)))) '% completed...']);

end
close(h);

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
        temporary=[temporary;Jtl(iRxn,:)];
    end
end
Jtl=temporary;

%Eliminate duplicates in triple reaction deletions
Jtl=unique(sort(Jtl,2),'rows');



Jsl=model.rxns(Jsl);
Jdl=model.rxns(Jdl);
Jtl=model.rxns(Jtl);



    case 2
%% Applying ROOM
disp('Identifying Jdl & Jtl - Part 1 of 1...')
environment = getEnvironment();
%Part 1 Starts with doing a ROOM with an individual KO, then, the
%NonZeroFluxes in the result determine the possible pairs for SL. Once
%determined the potential pairs, the process is repeated with pairs,
%a ROOM is computed with the PairKO and the nonZeroFluxes
%determine the potential reactions that could form the triplet.
SLTripletsIndexPart1              = cell(100000,1);
SLPairsIndexPart1                 = cell(100000,3);
numberOfPotentialSynLethals       = length(potentialSynLethals);
usedRxns =  [];
for iRxn = 1:numberOfPotentialSynLethals
   %Restoring environment for making parfor work
    modelKO = model; 
    currentSingleSLidx    = potentialSynLethals(iRxn);
    currentSingleSL       = modelKO.rxns(currentSingleSLidx);
    [fluxMinSingleSL,~,~] = ROOMtolerance(modelKO, euclMinSolution,currentSingleSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
    [fluxMaxSingleSL,~,~] = ROOMtolerance(modelKO, euclMaxSolution,currentSingleSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
    %Fluxes different from 0 in the results could be pairs that could lower
    %the biomass value
    discardRxns           = find(fluxMinSingleSL == 0 | fluxMaxSingleSL == 0);
    potentialPairRxns     = potentialSynLethals(~ismember(potentialSynLethals,discardRxns));

    if (~isempty(eliList))
        potentialPairRxns = potentialPairRxns(~ismember(potentialPairRxns,eliIdx)); %Eliminate Exchange and ATP Maintenance reactions
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
            [fluxMinPairSL,~,~]     = ROOMtolerance(modelKO, euclMinSolution,tempPairSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
            [fluxMaxPairSL,~,~]     = ROOMtolerance(modelKO, euclMaxSolution,tempPairSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
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
                    potentialTriplets = potentialTriplets(~ismember(potentialTriplets,eliIdx)); %Eliminate Exchange and ATPM reactions
                end
                if (~isempty(usedRxns))
                    potentialTriplets = potentialTriplets(~ismember(potentialTriplets,usedRxns(:,1)));
                    potentialTriplets = potentialTriplets(~ismember(potentialTriplets,usedRxns(:,2)));
                end
                
               
                if (~isempty(discardedSLTriplets))
                    discardAditiveEffect = zeros(length(discardedSLTriplets),1);
                    parfor i =1:length(potentialTriplets) 
                        restoreEnvironment(environment);
                       potTriplets = potentialTriplets;
                       currentPotentialTripletRxnidx = potTriplets(i);
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
                    [fluxMinTripletSL,~,~] = ROOMtolerance(modelKO1, euclMinSolution,tempTripletSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
                    [fluxMaxTripletSL,~,~] = ROOMtolerance(modelKO1, euclMaxSolution,tempTripletSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
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


%disp('Identifying Jdl & Jtl - Part 2 of 2...')
%Part 2 Starts with doing a ROOM with a Pair, then, the
%NonZeroFluxes in the result determine the possible triplets for SL. In the
%meantime, another triplet prediction is done using those reactions not
%located in the list of essential reactions generated by singleSL.
%Later Synthetic Lethals that involve reactions appearing in pairs and triplets will be deleted from the triplets.
%SLTripletsIndexPart2         = cell(100000,1);
%SLPairsIndexPart2            = cell(100000,1);

%for iRxn = 1:length(potentialSynLethals)
    %Restoring environment for making parfor work
%    restoreEnvironment(environment);
%    for jRxn = 1:length(potentialSynLethals)
%        if ~(jRxn<iRxn)
%           break; 
%        else
%         currentSingleSLidx      = potentialSynLethals(iRxn);
%         currentPotentialPairIdx = potentialSynLethals(jRxn);
%         idxsPair                    = [currentSingleSLidx currentPotentialPairIdx];
%         matchPairWithDiscardedPairs = sum(ismember(discardedSLPairs,idxsPair),2);
%            if ~isempty(matchPairWithDiscardedPairs(matchPairWithDiscardedPairs >= 2))
%                continue
%            end
%         currentSingleSL         = model.rxns(currentSingleSLidx);
%         currentPotentialPair    = model.rxns(currentPotentialPairIdx);
%         tempPairSL = [currentSingleSL currentPotentialPair];
%         [fluxMinPairSL,~,~] = ROOMtolerance(model, euclMinSolution,tempPairSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
%         [fluxMaxPairSL,~,~] = ROOMtolerance(model, euclMaxSolution,tempPairSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
%         fluxMinBMreduction  = (fluxMinPairSL(BioMassIndex)/BioMassMinSolution)  ;
%         fluxMaxBMreduction  = (fluxMaxPairSL(BioMassIndex)/BioMassMaxSolution)  ;
%            if (fluxMinBMreduction + fluxMaxBMreduction) <= 2*cutoff || ((sum(isnan(fluxMinPairSL)) ~= 0) && (sum(isnan(fluxMaxPairSL)) ~= 0))
%                SLPairsIndexPart2{iRxn,1}{jRxn,1} = [currentSingleSLidx currentPotentialPairIdx];
%                continue
%            elseif (sum(isnan(fluxMinBMreduction)) ~= 0 && fluxMaxBMreduction <= cutoff) ||  (sum(isnan(fluxMaxBMreduction)) ~= 0 && fluxMinBMreduction <= cutoff)
%                SLPairsIndexPart2{iRxn,1}{jRxn,1} = [currentSingleSLidx currentPotentialPairIdx];
%                continue
%            end
%                
%              NonZeroFluxesPairSLidx = find((fluxMaxPairSL ~= 0) & (fluxMinPairSL ~= 0));
%              NonZeroFluxesPairSLidx = NonZeroFluxesPairSLidx(~ismember(NonZeroFluxesPairSLidx,nonZeroRxnsIdx));
%                
%                if (~isempty(eliList))
%                    NonZeroFluxesPairSLidx = NonZeroFluxesPairSLidx(~ismember(NonZeroFluxesPairSLidx,eliIdx)); %Eliminate Exchange and ATPM reactions
%                end
%                
%                tempArraySLTriplets1 = zeros(length(NonZeroFluxesPairSLidx),3);
%                parfor kRxn=1:length(NonZeroFluxesPairSLidx) 
%                    restoreEnvironment(environment);
%                    modelKO2 = model;
%                    currentPotentialTripletRxnidx = NonZeroFluxesPairSLidx(kRxn);
%                    idxsTriplet                       = [currentSingleSLidx currentPotentialPairIdx currentPotentialTripletRxnidx];
%                    matchTripletWithDiscardedTriplets = sum(ismember(discardedSLTriplets,idxsTriplet),2);
%                    if ~isempty(matchTripletWithDiscardedTriplets(matchTripletWithDiscardedTriplets >= 3))
%                        continue
%                    end                   
%                    currentPotentialTripletRxn = modelKO2.rxns(currentPotentialTripletRxnidx);
%                    tempTripletSL = [tempPairSL currentPotentialTripletRxn];
%                     [fluxMinTripletSL,~,~] = ROOMtolerance(modelKO2, euclMinSolution,tempTripletSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
%                     [fluxMaxTripletSL,~,~] = ROOMtolerance(modelKO2, euclMaxSolution,tempTripletSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);       
%                     fluxMinBMreduction  = (fluxMinTripletSL(BioMassIndex)/BioMassMinSolution)  ;
%                     fluxMaxBMreduction  = (fluxMaxTripletSL(BioMassIndex)/BioMassMaxSolution)  ;         
%                     
%                     if or((fluxMinBMreduction + fluxMaxBMreduction) <= 2*cutoff,((sum(isnan(fluxMinTripletSL)) ~= 0) && (sum(isnan(fluxMaxTripletSL)) ~= 0)))
%                         tempIdxs = [currentSingleSLidx currentPotentialPairIdx currentPotentialTripletRxnidx];
%                         tempArraySLTriplets1(kRxn,:)  =  tempIdxs;
%                         
%                    elseif (sum(isnan(fluxMinTripletSL)) ~= 0 && fluxMaxBMreduction <= cutoff) ||  (sum(isnan(fluxMinTripletSL)) ~= 0 && fluxMaxBMreduction <= cutoff)
%                         tempIdxs = [currentSingleSLidx currentPotentialPairIdx currentPotentialTripletRxnidx];
%                         tempArraySLTriplets1(kRxn,:)  =  tempIdxs;
%                    end
%                  
%                end
%                %SLTripletsIndexPart2{iRxn,1}{jRxn,1} = tempArraySLTriplets;
%                %clear idx tempArraySLTriplets tempIdxs modelKO currentPotentialTripletRxn currentPotentialTripletRxnidx
%                
%                
%                tempArraySLTriplets2 = zeros(length(potentialSynLethals),3);
%                parfor kRxn=1:length(potentialSynLethals)
%                    restoreEnvironment(environment);
%                    modelKO3 = model;
%                    if ~(kRxn<jRxn)
%                        continue
%                    else
%                        currentPotentialTripletRxnidx = potentialSynLethals(kRxn);
%                        idxsTriplet                       = [currentSingleSLidx currentPotentialPairIdx currentPotentialTripletRxnidx];
%                        matchTripletWithDiscardedTriplets = sum(ismember(discardedSLTriplets,idxsTriplet),2);
%                       if ~isempty(matchTripletWithDiscardedTriplets(matchTripletWithDiscardedTriplets >= 3))
%                         continue
%                      end                       
%                        currentPotentialTripletRxn = modelKO3.rxns(currentPotentialTripletRxnidx);
%                        tempTripletSL = [tempPairSL currentPotentialTripletRxn];
%                        [fluxMinTripletSL,~,~] = ROOMtolerance(modelKO3, euclMinSolution,tempTripletSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
%                        [fluxMaxTripletSL,~,~] = ROOMtolerance(modelKO3, euclMaxSolution,tempTripletSL,'delta', 0.15,'epsilon', 0.015,'printLevel',0); 
%                        fluxMinBMreduction  = (fluxMinTripletSL(BioMassIndex)/BioMassMinSolution)  ;
%                        fluxMaxBMreduction  = (fluxMaxTripletSL(BioMassIndex)/BioMassMaxSolution)  ;       
%                        
%                        if or((fluxMinBMreduction + fluxMaxBMreduction) <= 2*cutoff,((sum(isnan(fluxMinTripletSL)) ~= 0) && (sum(isnan(fluxMaxTripletSL)) ~= 0)))
%                         tempIdxs = [currentSingleSLidx currentPotentialPairIdx currentPotentialTripletRxnidx];
%                         tempArraySLTriplets2(kRxn,:)  =  tempIdxs;
%                         
%                        elseif (sum(isnan(fluxMinTripletSL)) ~= 0 && fluxMaxBMreduction <= cutoff) ||  (sum(isnan(fluxMinTripletSL)) ~= 0 && fluxMaxBMreduction <= cutoff)
%                         tempIdxs = [currentSingleSLidx currentPotentialPairIdx currentPotentialTripletRxnidx];
%                         tempArraySLTriplets2(kRxn,:)  =  tempIdxs;
%                        end
%                        
%                    end
%                    
%                end
%                tempArraySLTripletsBoth = [tempArraySLTriplets1;tempArraySLTriplets2];
%                SLTripletsIndexPart2{iRxn,1}{jRxn,1} = tempArraySLTripletsBoth;
               %clear idx tempArraySLTriplets tempIdxs modelKO
%        end 
%    end
%    disp(['Finished analyzing reaction number ' num2str(iRxn) ' out of ' num2str(length(potentialSynLethals))])  
%end

%Join pairs indexes
indexEmptyCellsPart1 = cellfun(@isempty,SLPairsIndexPart1);
%indexEmptyCellsPart2 = cellfun(@isempty,SLPairsIndexPart2);
SLPairsIndexPart1(indexEmptyCellsPart1) = [];
%SLPairsIndexPart2(indexEmptyCellsPart2) = [];
SLPairsFinal = [SLPairsIndexPart1']; %SLPairsIndexPart2];
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
    indexEmptyCellsPart1 = cellfun(@isempty,SLTripletsIndexPart1);
    %indexEmptyCellsPart2 = cellfun(@isempty,SLTripletsIndexPart2);
    SLTripletsIndexPart1(indexEmptyCellsPart1) = [];
    %SLTripletsIndexPart2(indexEmptyCellsPart2) = [];
    SLTripletsCell = [SLTripletsIndexPart1];%SLTripletsIndexPart2];
%check again empty cells, just in case any of those was empty
    indexEmptyCellsFinal = cellfun(@isempty,SLTripletsCell);
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
end
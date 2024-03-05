function [Jsl, BMReductionValues]=singleSLModified(model,cutoff,eliList,furthestSolution,centralSolution)
%% [Jsl]=singleSL(model,cutoff,eliList,atpm)
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
% atpm           ATPM Reaction Id in model.rxns if other than 'ATPM'
% method
    %1 = Original fastSL algorithm
    %2 =  MOMA (Requires two extra inputs: maxSolution minSolution)
    %3 =  ROOM (Requires two extra inputs: maxSolution minSolution)
% maxSolution = After a sampling and euclidean distance filtration, this
% solution is the one that has the maximum euclidean distance value within the
% valid solutions.
% minSolution = After a sampling and euclidean distance filtration, this
% solution is the one that has minimum euclidean distance value within the
% valid solutions.
%OUTPUT
% Jsl            Single lethal reactions identified
% Aditya Pratapa       6/26/14. 

    if ~exist('furthestSolution','var') || ~exist('centralSolution','var') || isempty(furthestSolution) || isempty(centralSolution)  
        error('Missing variables: furthestSolution or centralSolution')
    end


if ~exist('cutoff', 'var') || isempty(cutoff)
        cutoff = 0.01;
end

if exist('eliList', 'var') || isempty(eliList)
        eliList = model.rxns(find(model.c)); 
end


Jsl=[];

%Step1 Identify Single Lethal Reactions...

        nonZeroFluxesIdxs = find(furthestSolution ~= 0 & centralSolution ~= 0);
        %solWT             = optimizeCbModel(model,'max','one');
        %nonZeroFluxesIdxs = find(~eq(solWT.x,0));
        BioMassIndex = find(model.c);
        BioMassMinSolution = centralSolution(find(model.c)); %#ok<*FNDSB>
        BioMassMaxSolution = furthestSolution(find(model.c));
        if (~isempty(eliList))
            eliIdx = find(ismember(model.rxns,eliList)); %Index of reactions not considered for lethality analysis
            nonZeroFluxesIdxs    = nonZeroFluxesIdxs(~ismember(nonZeroFluxesIdxs,eliIdx)); %nonZeroFluxes
        end
        
   
 %Identify Single Lethal Reaction Deletions...

        BMReductionValues       = zeros (length(nonZeroFluxesIdxs),3);
        model.lb(find(model.c)) = 0;
        essentialRxnsIdx        = zeros(size(nonZeroFluxesIdxs,1),1);
        environment             = getEnvironment();
        potentialEssentialRxns  = length(nonZeroFluxesIdxs);
        parfor iRxn = 1:potentialEssentialRxns
            restoreEnvironment(environment);
            modelKO  = model;
            currentKOidx  = nonZeroFluxesIdxs(iRxn);
            currentKO     = modelKO.rxns(currentKOidx);
            %Delta and epsilon values are recovered from the ROOM paper in
            %which, for lethality purposes they use the values used here. 
                [fluxMin,~,~] =  ROOMtolerance(modelKO, centralSolution,currentKO,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
                [fluxMax,~,~] =  ROOMtolerance(modelKO, furthestSolution,currentKO,'delta', 0.15,'epsilon', 0.015,'printLevel',0);
            %Compute biomass reduction porcentage(without multiplying by
            %100)
                fluxMinBMreduction        = (fluxMin(BioMassIndex)/BioMassMinSolution) ;
                fluxMaxBMreduction        = (fluxMax(BioMassIndex)/BioMassMaxSolution) ;
                BMReductionValues(iRxn,:) = [currentKOidx fluxMinBMreduction fluxMaxBMreduction];
            %If the sum of both porcentages is equal or less than 2 times
            %the cutoff value, then it is an essential reaction.
                if ((fluxMinBMreduction + fluxMaxBMreduction) <= 2*cutoff) || ((sum(isnan(fluxMin)) ~= 0) && (sum(isnan(fluxMax)) ~= 0))
                    essentialRxnsIdx(iRxn) = currentKOidx; 
                elseif sum(isnan(fluxMinBMreduction)) ~= 0 && fluxMaxBMreduction <= cutoff
                    essentialRxnsIdx(iRxn) = currentKOidx;
                elseif sum(isnan(fluxMaxBMreduction)) ~= 0 && fluxMinBMreduction <= cutoff
                    essentialRxnsIdx(iRxn) = currentKOidx;
                end
                    
           disp(['Finished analyzing reaction number ' num2str(iRxn) ' out of ' num2str(potentialEssentialRxns)])
        end
        essentialRxnsIdx = essentialRxnsIdx(essentialRxnsIdx~=0);
        Jsl = essentialRxnsIdx;
         
end
   
function [Jsl,Jdl,Jtl,BMReduction] = fastSLROOM2(model,cutoff,ArrayWithLethals,furthestSolution,centralSolution)
%This function tests the lethals returned by the fastSL algorithm. The
%lethal is validated as lethal if it reduces the biomass lower than the
%cutoff in both fluxes, furthest and central solution. It simulates the
%knockout with ROOM algorithm
%Inputs:
    %ArrayWithLethals: Cell array containing three cells. Each cell will include a
    % vertical cell array corresponding to the reactions that form
    % single lethals, pair lethals, and triplet lethals.

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
% Jsl        single lethal reactions identified
% Jdl        double lethal reactions identified
% Jtl        triple lethal reactions identified
% Fernando Silva-Lance 2021
%%
if ~exist('furthestSolution','var') || ~exist('centralSolution','var') || isempty(furthestSolution) || isempty(centralSolution)  
        error('Missing variables: furthestSolution or centralSolution')
end
        BioMassIndex = find(model.c);
        BioMassMinSolution = centralSolution(find(model.c)); %#ok<*FNDSB>
        BioMassMaxSolution = furthestSolution(find(model.c));

%% Applying ROOM
        model.lb(find(model.c)) = 0;
        
        environment             = getEnvironment();
        singleLethals = ArrayWithLethals{1,1};
            essentialSingles       = cell(size(singleLethals,1),1);
        doubleLethals = ArrayWithLethals{2,1};
            essentialPairs        = cell(size(doubleLethals,1),2);
        tripleLethals = ArrayWithLethals{3,1};
            essentialTriplets     = cell(size(tripleLethals,1),3);
            
        BMReductionValuesSingles   = cell (size(singleLethals,1),3);
        BMReductionValuesPairs     = cell (size(doubleLethals,1),4);
        BMReductionValuesTriplets  = cell (size(tripleLethals,1),5);
        parfor iRxn = 1:size(singleLethals,1)
            restoreEnvironment(environment);
            modelKO  = model;
            currentLethal     = singleLethals(iRxn,:);
            %Delta and epsilon values are recovered from the ROOM paper in
            %which, for lethality purposes they use the values used here. 
                [fluxMin,~,~] =  ROOMtolerance(modelKO, centralSolution,currentLethal,'delta', 0.1,'epsilon', 0.01,'printLevel',0);
                [fluxMax,~,~] =  ROOMtolerance(modelKO, furthestSolution,currentLethal,'delta', 0.1,'epsilon', 0.01,'printLevel',0);
            %Compute biomass reduction porcentage(without multiplying by
            %100)
                fluxMinBMreduction        = (fluxMin(BioMassIndex)/BioMassMinSolution) ;
                fluxMaxBMreduction        = (fluxMax(BioMassIndex)/BioMassMaxSolution) ;
                BMReductionValuesSingles(iRxn,:) = [currentLethal num2cell(fluxMinBMreduction) num2cell(fluxMaxBMreduction)];
            %If the sum of both porcentages is equal or less than 2 times
            %the cutoff value, then it is an essential reaction.
                if ((fluxMinBMreduction + fluxMaxBMreduction) == 0) || ((sum(isnan(fluxMin)) ~= 0) && (sum(isnan(fluxMax)) ~= 0))
                    essentialSingles(iRxn) = currentLethal; 
                elseif sum(isnan(fluxMinBMreduction)) ~= 0 && fluxMaxBMreduction <= cutoff
                    essentialSingles(iRxn) = currentLethal;
                elseif sum(isnan(fluxMaxBMreduction)) ~= 0 && fluxMinBMreduction <= cutoff
                    essentialSingles(iRxn) = currentLethal;
                end
                    
           disp(['Finished analyzing singleLethal number ' num2str(iRxn) ' out of ' num2str(size(singleLethals,1))])
        end
        

        parfor iRxn = 1:size(doubleLethals,1)
            restoreEnvironment(environment);
            modelKO  = model;
            currentLethal     = doubleLethals(iRxn,:);
            %Delta and epsilon values are recovered from the ROOM paper in
            %which, for lethality purposes they use the values used here. 
                [fluxMin,~,~] =  ROOMtolerance(modelKO, centralSolution,currentLethal,'delta', 0.1,'epsilon', 0.01,'printLevel',0);
                [fluxMax,~,~] =  ROOMtolerance(modelKO, furthestSolution,currentLethal,'delta', 0.1,'epsilon', 0.01,'printLevel',0);
            %Compute biomass reduction porcentage(without multiplying by
            %100)
                fluxMinBMreduction        = (fluxMin(BioMassIndex)/BioMassMinSolution) ;
                fluxMaxBMreduction        = (fluxMax(BioMassIndex)/BioMassMaxSolution) ;
                BMReductionValuesPairs(iRxn,:) = [currentLethal num2cell(fluxMinBMreduction) num2cell(fluxMaxBMreduction)];
            %If the sum of both porcentages is equal or less than 2 times
            %the cutoff value, then it is an essential reaction.
                if ((fluxMinBMreduction + fluxMaxBMreduction) == 0) || ((sum(isnan(fluxMin)) ~= 0) && (sum(isnan(fluxMax)) ~= 0))
                    essentialPairs(iRxn,:) = currentLethal; 
                elseif sum(isnan(fluxMinBMreduction)) ~= 0 && fluxMaxBMreduction <= cutoff
                    essentialPairs(iRxn,:) = currentLethal;
                elseif sum(isnan(fluxMaxBMreduction)) ~= 0 && fluxMinBMreduction <= cutoff
                    essentialPairs(iRxn,:) = currentLethal;
                end
                    
           disp(['Finished analyzing PairLethal number ' num2str(iRxn) ' out of ' num2str(size(doubleLethals,1))])
        end
       
        
        parfor iRxn = 1:size(tripleLethals,1)
            restoreEnvironment(environment);
            modelKO  = model;
            currentLethal     = tripleLethals(iRxn,:);
            %Delta and epsilon values are recovered from the ROOM paper in
            %which, for lethality purposes they use the values used here. 
                [fluxMin,~,~] =  ROOMtolerance(modelKO, centralSolution,currentLethal,'delta', 0.1,'epsilon', 0.01,'printLevel',0);
                [fluxMax,~,~] =  ROOMtolerance(modelKO, furthestSolution,currentLethal,'delta', 0.1,'epsilon', 0.01,'printLevel',0);
            %Compute biomass reduction porcentage(without multiplying by
            %100)
                fluxMinBMreduction        = (fluxMin(BioMassIndex)/BioMassMinSolution) ;
                fluxMaxBMreduction        = (fluxMax(BioMassIndex)/BioMassMaxSolution) ;
                BMReductionValuesTriplets(iRxn,:) = [currentLethal num2cell(fluxMinBMreduction) num2cell(fluxMaxBMreduction)];
            %If the sum of both porcentages is equal or less than 2 times
            %the cutoff value, then it is an essential reaction.
                if ((fluxMinBMreduction + fluxMaxBMreduction) == 0) || ((sum(isnan(fluxMin)) ~= 0) && (sum(isnan(fluxMax)) ~= 0))
                    essentialTriplets(iRxn,:) = currentLethal; 
                elseif sum(isnan(fluxMinBMreduction)) ~= 0 && fluxMaxBMreduction <= cutoff
                    essentialTriplets(iRxn,:) = currentLethal;
                elseif sum(isnan(fluxMaxBMreduction)) ~= 0 && fluxMinBMreduction <= cutoff
                    essentialTriplets(iRxn,:) = currentLethal;
                end
                    
           disp(['Finished analyzing TripleLethal number ' num2str(iRxn) ' out of ' num2str(size(tripleLethals,1))])
        end
 Jtl = essentialTriplets;  
    indexEmptyCells = cellfun(@isempty,Jtl(:,1));
    Jtl(indexEmptyCells,:) = [];
  Jdl = essentialPairs;
     indexEmptyCells = cellfun(@isempty,Jdl(:,1));
     Jdl(indexEmptyCells,:) = [];
  Jsl = essentialSingles;
      indexEmptyCells = cellfun(@isempty,Jsl(:,1));
     Jsl(indexEmptyCells,:) = [];
  BMReduction = [{BMReductionValuesSingles};{BMReductionValuesPairs};{BMReductionValuesTriplets}];
save('BackupTripleSL')
end
function [effectiveDrugTargets, dist2WT] = singleSLCovid(targets, modelDrug,modelPlacebo,placeboCentralSolution,placeboFurthestSolution,drugFurthestSolution,drugCentralSolution,threshold)
%% [Jsl]=singleSL(model,cutoff,eliList,atpm)
% INPUT
% model (the following fields are required - others can be supplied)       
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
%   rxns         Reaction Names
% euclMaxSolution = After a sampling and euclidean distance filtration, this
%   solution is the one that has the maximum euclidean distance value within the
%   valid solutions.
% euclMinSolution = After a sampling and euclidean distance filtration, this
%   solution is the one that has minimum euclidean distance value within the
%   valid solutions.
%OUTPUT
% Jsl            Single lethal reactions identified
% Aditya Pratapa       6/26/14. 
if iscell(targets)
    targetsIdx = find(ismember(modelPlacebo.rxns,targets));
end

    if or(or(~exist('placeboCentralSolution','var'),~exist('placeboFurthestSolution','var')),or(~exist('drugFurthestSolution','var'),~exist('drugCentralSolution','var')))
        error('For ROOM fastSL algorithm, one of the following, variables is missing: placeboFurthestSolution drugFurthestSolution drugCentralSolution placeboCentralSolution')
    end
    if ~exist('threshold','var') ||isempty(threshold)
        disp('No threshold set for ROOM fastSL, default threshold will be used...threshold = .5')
        threshold = 0.5;
    end
effectiveDrugTargets=[];

%Step1 Identify Single Lethal Reactions...
   
 %Identify Single Lethal Reaction Deletions...
%dist2WT = Vector that saves the distance to the placeboWT for each KO
        %Distance from WTplacebo to WTdrug
        origDistFurthest     = norm(placeboFurthestSolution - drugFurthestSolution);
        acceptedDistFurthest = origDistFurthest * threshold;
        origDistCentral      = norm(placeboCentralSolution - drugCentralSolution);
        acceptedDistCentral  = origDistCentral * threshold;
        dist2WT = zeros(size(targets,1),3);
        modelPlacebo.lb(find(modelPlacebo.c)) = 0;
        effectiveTargetsRxnsIdx               = zeros(size(targets,1),1);
        environment = getEnvironment();
        parfor iRxn = 1:length(targets)
            restoreEnvironment(environment);
            modelKO       = modelPlacebo;
            currentKOidx  = targetsIdx(iRxn);
            currentKO     = modelKO.rxns(currentKOidx);
            %Delta and epsilon values are recovered from the ROOM paper in
            %which, for lethality purposes they use the values used here. 
            %QUE PASA CON AQUELLOS KOs que devuelven infeasible?
                [~,~,currentDistMin] = ROOMtolerance(modelKO, drugCentralSolution,currentKO,'delta', 0.1,'epsilon', 0.01);
                [~,~,currentDistMax] = ROOMtolerance(modelKO, drugFurthestSolution,currentKO,'delta', 0.1,'epsilon', 0.01);
                dist2WT(iRxn,:) = [currentKOidx currentDistMin currentDistMax];
            %If the sum of both porcentages is equal or less than 2 times
            %the cutoff value, then it is an essential reaction.
                if acceptedDistFurthest >= currentDistMax && acceptedDistCentral >= currentDistMin
                    effectiveTargetsRxnsIdx(iRxn) = currentKOidx; 
                end
            %waitbar(iRxn/length(nonZeroFluxesIdxs),h,[num2str(round(iRxn*100/length(nonZeroFluxesIdxs))) '% completed...']);
        end
        effectiveTargetsRxnsIdx = effectiveTargetsRxnsIdx(effectiveTargetsRxnsIdx~=0);
        effectiveDrugTargets    = modelDrug.rxns(effectiveTargetsRxnsIdx);      
end
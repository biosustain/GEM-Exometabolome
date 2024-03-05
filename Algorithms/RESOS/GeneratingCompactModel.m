function modelCompacted = GeneratingCompactModel(modelLP,BioMassIndex,NotCompressedIndex)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Modified by Fernando Silva-Lance 2021 %%%%%%%%%%%%%%%%%
%inputs:
    %modelLP: Cobra Model in linear structure
    %BiomassIndex:Index of the biomass
    %NotCompressedIndex: Index of those reactions that the user wishes to
    %not compared
%Create reaction sets
[Sset,SetLb,SetUb,ReactionSet,MetaboliteSet] = findReactionSetsNotCompressingSeveralReactions(full(modelLP.S),modelLP.lb,modelLP.ub,BioMassIndex,NotCompressedIndex);

% Match reaction names and only include rxns not associated to biomass
numRxns = size(SetLb,1);
rxns{numRxns,1} = [];
counter = 1;

%Generate CompactModel
while counter <= numRxns
    rxnIdx = find(ReactionSet(:,1) == counter);
    
    % If the reaction has not been lumped
    if length(rxnIdx) == 1
        rxns{counter} = modelLP.rxns{rxnIdx};
    else
        % Find reactions with the same identifier
        for j = 1:length(rxnIdx)
            if j == 1
                tempName = modelLP.rxns{rxnIdx(j)};
            else
%                tempName = [tempName,'_',modelLP.rxns{rxnIdx(j)}];
                tempName = [tempName,'@',modelLP.rxns{rxnIdx(j)}];
            end
        end
        rxns{counter} = tempName;
        
    end
    counter = counter + 1;
end

% Build final model structure
modelCompacted.S        = sparse(Sset);
modelCompacted.lb       = SetLb;
modelCompacted.ub       = SetUb;
modelCompacted.c        = zeros(size(SetLb));
modelCompacted.b        = zeros(size(Sset,1),1);
modelCompacted.rxns     = rxns;
modelCompacted.rxnNames = rxns;
modelCompacted.mets     = modelLP.mets(MetaboliteSet~=-1);
modelCompacted.metNames = modelCompacted.mets;
modelCompacted.description = [modelLP.description,'_compacted'];
modelCompacted.numRxns     = numel(modelCompacted.c); 
end 
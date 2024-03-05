function modelFinal = compactCbModelFinal(model,biomass, NotCompressed)
% Compacts constrained-based model structure
% Inputs:  model sampling structure
    %biomass      :  Cell array with the biomass reaction
    %NotCompressed: Reaction Names which will be tried to not be
                    %compressed.
% Outputs: compacted model structure
%%%%%%%%%%%%%%%%%%%%%% based on Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Fernando Silva-Lance 2021 %%%%%%%%%%%%%%%%%%%%%%%%%%

%Linearize the model in an LP structure
    modelLP = linearStructure(model);

disp('Performing initial FVA reduction...')
%Perform an FVA using the modified version
[modelLP,EliminatedRxns1] = FVA(modelLP,NotCompressed,biomass);

BioMassIndex       = find(ismember(modelLP.rxns, biomass));
NotCompressedIndex = find(ismember(modelLP.rxns, NotCompressed));

%While loop that tries compacting the model,if the outcome is an infeasible model, 
%then tries compacting less reactions by randomly excluding them of the
%compaction.
    flag = 0;
    percentRxnsCompacted = 100;
    numberOfRxns = size(modelLP.rxns,1);
while flag == 0
    %%%MAIN FUNCTION
    modelCompacted  = GeneratingCompactModel(modelLP,BioMassIndex,NotCompressedIndex);
    %%%
    indexBM         = find(ismember(modelCompacted.rxns,biomass));
    modelCompacted.c(indexBM) = 1;
    test = optimizeCbModel(modelCompacted);
    if isnan(test.f)
       disp('Error in optimization, compacting less reactions');
        %Compacting 10% less reactions
            percentRxnsCompacted   = percentRxnsCompacted - 10;
            newNumberNotCompressed = round(numberOfRxns * ((100-percentRxnsCompacted)/100));
        %Choose random reactions
            randomRxnsIndex        = randsample(numberOfRxns,newNumberNotCompressed)';
            newNotCompressedRxns   = modelLP.rxns(randomRxnsIndex);
            NotCompressed          = [NotCompressed;newNotCompressedRxns];
            %Just in case some are repeated
                NotCompressed      = unique(NotCompressed);
            %Checking that biomass is not being compressed
                NotCompressed      = NotCompressed(~ismember(NotCompressed,biomass));
                NotCompressedIndex = find(ismember(modelLP.rxns, NotCompressed));
        continue
    end
    % Run final consistency check using FVA
    modelFinal = linearStructure(modelCompacted);
    disp('Performing final FVA reduction...')
    [modelFinal,~]               = FVA(modelFinal,NotCompressed,biomass);
    modelFinal.c                 = zeros(size(modelFinal.rxns,1),1);
    BioMassIndex2                = find(ismember(modelFinal.rxns, biomass));
    modelFinal.c(BioMassIndex2)  = 1;
    test = optimizeCbModel(modelFinal);
    if isnan(test.f)
       disp('Error in optimization, compacting less reactions');
        %Compacting 10% less reactions
            percentRxnsCompacted   = percentRxnsCompacted - 10;
            newNumberNotCompressed = round(numberOfRxns * ((100-percentRxnsCompacted)/100));
        %Choose random reactions
            randomRxnsIndex        = randsample(numberOfRxns,newNumberNotCompressed)';
            newNotCompressedRxns   = modelLP.rxns(randomRxnsIndex);
            NotCompressed          = [NotCompressed;newNotCompressedRxns];
            %Just in case some are repeated
                NotCompressed      = unique(NotCompressed);
            %Checking that biomass is not being compressed
                NotCompressed      = NotCompressed(~ismember(NotCompressed,biomass));
                NotCompressedIndex = find(ismember(modelLP.rxns, NotCompressed));
        continue
    end
flag = 1;
end

     
end
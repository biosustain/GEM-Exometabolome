function [modelSimplificated,modelLoopless] = modelSimplification(model,NotEliminate)
%%%%%%%%%%%%%%%%%%%Fernando Silva-Lance 2021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simplifies a model by performing a loopless and a compaction.
%It does not eliminate exchange reactions if no other reactions are specified. 
%Inputs:
    %model:COBRA Structure model
    %NotEliminate: Vertical cell array containing reactions that you wish
    %to preserve
%Outputs
    %modelLoopless
    %modelSimplificated: Loopless + compact
if ~exist('NotEliminate', 'var') || isempty(NotEliminate)
    disp('Not compressing exchange reactions')
           [ExcRxns, ~]   = findExcRxns(model);
           NotEliminate   = model.rxns(ExcRxns);
end

%% Loopless

           tic
            [modelLoopless,RemovedRxns] = LooplessModel(model, NotEliminate);
           toc
           

 %% Compact Model

    %In a vertical cell array, add every reaction that you wish not to be
    %compressed. Biomass should be at the bottom of the array
        biomass        = checkObjective(modelLoopless);
    %Compact CbModel and retrieve the reaction sets
            modelSimplificated    = compactCbModelFinal(modelLoopless,biomass, NotEliminate);
            modelSimplificated.modelsense = 'max';
            modelSimplificated.lb(find(modelSimplificated.c)) = 0;
end
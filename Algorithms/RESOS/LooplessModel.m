function [modelLoopless2,EliminatedRxns] = LooplessModel(model, NotEliminate, toleranceBiomass, toleranceLoops)
%%%%%%%%%%%%%%%%%%%%%%%%%%% Fernando Silva-Lance 2021 %%%%%%%%%%%%%%%%%%%%
%Creates model without loops
%Inputs
    %model: CobraModel 
    %NotEliminate: Vertical Cell arrays with the name of those reactions
    %that the user wishes to keep. 
    %toleranceBiomass: Tolerance referent to the maximum variation the objective could have
            %in a loopless optimization.
    %toleranceLoops: Tolerance that determines if , after loopless
        %optimization, that reaction is indeed part of a loop.
%Outputs
    %modelLoopless2: LooplessModel
    %EliminatedRxns: Those reactions that were removed from the model.
    
%Find if toleranceBiomass exist 
if ~exist('toleranceBiomass', 'var') || isempty(toleranceBiomass)
  toleranceBiomass = .01;
  disp('Setting default value of toleranceBiomass. This value ensures that biomass production is not lower than optimal - tolerance');
end

%Find if toleranceLoops exist 
if ~exist('toleranceLoops', 'var') || isempty(toleranceLoops)
  toleranceLoops = .001;
  disp('Setting default value of toleranceLoops.'); %Identify and remove loop reactions
end

 %Determine the optimal BM production and the fluxes allowing loops
    BioMassIndex  = find(model.c);
    BiomassRxn    = model.rxns( BioMassIndex);
    SolLoop       = optimizeCbModel(model);
    FluxSetLoop   = SolLoop.x;
    Biomass       = SolLoop.f;
    modelIni      = model;
    
%Ensure that biomass production is not lower than the optimal-tolerance
    modelIni.lb(BioMassIndex) = Biomass*(1-toleranceBiomass);  
    notEliminateIndex         = find(ismember(modelIni.rxns,NotEliminate));
    
%Optimize Model with Loopless option
    SolLoopless     = optimizeCbModel(modelIni,'max',0,0);
    FluxSetLoopless = SolLoopless.x;

%Identify and remove loop reactions

%( abs(FluxSetLoop)>toleranceLoops ) : Returns every reaction whose absolute value
% is over the tolerance in the normal optimization. 

% abs(FluxSetLoopless)<toleranceLoops) : Returns every reaction whose absolute value
% after loopless optimization is under de tolerannce. 
%If the reaction appears in both, then it is a loopless reaction.
%In other word a reaction is part of a loop when it changes its flux from a value higher
%than the tolerance to a value below the tolerance.
try
LoopRxn            = find( ( ( abs(FluxSetLoop)>toleranceLoops ) + ( abs(FluxSetLoopless)<toleranceLoops) )==2);
LoopRxn            = LoopRxn( ~ismember(LoopRxn,notEliminateIndex) );
EliminatedRxns     = modelIni.rxns(LoopRxn);

%Test, if removing the reactions causes infeasibility, try removing
%individually
modelLoopless      = removeRxns(modelIni,modelIni.rxns(LoopRxn));
test = optimizeCbModel(modelLoopless);
if isnan(test.f)
    disp('Removing loopless reactions caused infeasibility')
    modelLoopless = modelIni;
    counter       = 1;
    
    for i = 1:size(EliminatedRxns,1)
         currentRxn     = EliminatedRxns(i);
         modelLoopless = removeRxns(modelLoopless,currentRxn);
         test = optimizeCbModel(modelLoopless);
        if isnan(test.f)
            modelLoopless          = modelIni;
            errorMoving(counter,1) = model.rxns(LoopRxn(i));
            counter                = counter +1;
            continue
        end
        modelIni = modelLoopless;
    end
    
end

catch
    modelLoopless = modelIni;
    EliminatedRxns = {};
end

%FVA for not leaving isolated reactions.
    modelLP                            = linearStructure(modelLoopless);
    [modelLoopless2,FVAeliminatedRxns] = FVA(modelLP,NotEliminate,BiomassRxn);
    EliminatedRxns                     = [EliminatedRxns ;FVAeliminatedRxns]; 
end
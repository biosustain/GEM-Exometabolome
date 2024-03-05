function [modelRmvdRxns,EliminatedRxns] = FVA(model,NotEliminate,BiomassRxn)
% Performs FVA
% Inputs:  model structure
            %BiomassRxn in a cell array
            %NotEliminate: Indexes of those reactions which are wished to
            %not elimiante
% Outputs: model structure with new bounds and directionalities

%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Modified by Fernando Silva-Lance 2021 %%%%%%%%%%%%%%%%%%%%%%
%This modification deals with solvers tolerance issues when having
%experimental boundaries that could cause infeasibiity.
%If you want to check which modifications were done, look for the word:
    %MODIFICATION
tol = 1e-9;

% Define and run FVA problem
f   = model.obj;
FVA = zeros(model.numRxns,2);
progIndicator = fix(model.numRxns/10);

% Main loop
%Counters used to save those reactions that produce infeasibiity when being
%the objective. Change the lower or upper boundary to 0 when infeasibility
%is the outcome.
counterUpper = 1;
counterLower = 1;
for i = 1:model.numRxns   
    
    % Show progress
    if ~mod(i,progIndicator) && i>0
        disp(['Current progress... ',num2str(10*i/progIndicator),'%']);
    end
    model.obj    = f;
    model.obj(i) = 1;
    model.c = model.obj;
    % Mimization problem
    %MODIFICATION
   try
        sol = optimizeCbModel(model,'min');
        FVA(i,1)     = sol.x(i);
   catch
       %Change lower boundary to 0 
       FVA(i,1) = 0;
       checkLowerBoundaries(counterLower,1) = i;
       %disp(model.rxns(i));
       counterLower = counterLower +1;
   end
   
    % Maximization problem;
     %MODIFICATION
    try
    sol = optimizeCbModel(model,'max');
    FVA(i,2)     = sol.x(i);
    catch
       %Change upper boundary to 0 
       %disp(model.rxns(i));
       FVA(i,2) = 0;
       checkUpperBoundaries(counterUpper,1) = i;
       counterUpper = counterUpper +1;
    end
    
end

 %MODIFICATION
%Check reactions that produce infeasibility when minimizing
    %If the modified lower boundary (0) is equal to the upper Bound, then, not
    %eliminate reaction because that modified boundary could not necessary be a
    %0. If the modified boundary is higher than the upper Bound, then eliminate
    %the reaction by changing both bounds to 0.
counter = 1;
if exist('checkLowerBoundaries', 'var')
    
    for i = 1: size(checkLowerBoundaries,1)
        currentRxn = checkLowerBoundaries(i);
        modifiedLB = FVA(currentRxn,1);
        currentUB  = FVA(currentRxn,2);
        if modifiedLB == currentUB
            notEliminate(counter) = currentRxn;
            counter               = counter + 1;
            continue
        elseif modifiedLB > currentUB 
            modifiedLB        = currentUB;
            FVA(currentRxn,1) = modifiedLB;
        end
    end
    
end

%MODIFICATION
%Check reactions that produce infeasibility when maximizing
    %If the modified upper boundary (0) is equal to the lower Bound, then, not
    %eliminate reaction because that modified boundary could not necessary be a
    %0. If the lower boundary is higher than the modified upper Bound, then eliminate
    %the reaction by changing both bounds to 0.
if exist('checkUpperBoundaries', 'var')
    
    for i = 1: size(checkUpperBoundaries,1)
        currentRxn = checkUpperBoundaries(i);
        currentLB  = FVA(currentRxn,1);
        modifiedUB = FVA(currentRxn,2);
        if currentLB == modifiedUB
            notEliminate(counter) = currentRxn;
            counter               = counter + 1;
            continue
        elseif currentLB > modifiedUB 
            modifiedUB        = currentLB;
            FVA(currentRxn,2) = modifiedUB;
        end
    end
    
end

 %MODIFICATION    
% Find zero reactions and avoid eliminating those reactions determined by
% the user.
zeroRxns = find((abs(FVA(:,1))<tol)+(abs(FVA(:,2))<tol)==2);
if exist('notEliminate', 'var')
    zeroRxns = zeroRxns(~ismember(zeroRxns,notEliminate));
end
zeroRxns       = zeroRxns(~ismember(model.rxns(zeroRxns),NotEliminate));
EliminatedRxns = model.rxns(zeroRxns);


%MODIFICATION
if isempty(zeroRxns)
    % Reset objective function
        model.obj  = f;
        model.c    = model.obj;
        model.c(:) = 0;
        biomassIndex          = find(ismember(model.rxns,BiomassRxn));
        model.c(biomassIndex) = 1;
    % Redefine model boundaries
        modelOri = model;
        model.lb(1:model.numRxns) = FVA(:,1).*(abs(FVA(:,1))>tol);
        model.ub(1:model.numRxns) = FVA(:,2).*(abs(FVA(:,2))>tol);
        %If the change causes infeasibility, return to the last stable
        %version.
        test = optimizeCbModel(model);
        if isnan(test.f)
            model = modelOri;
        end
        modelRmvdRxns = model;
    % Re-format the problem
    %modelFVA = linearStructure(model);
else
    
%MODIFICATION
    model.c(:)            = 0;
    biomassIndex          = find(ismember(model.rxns,BiomassRxn));
    model.c(biomassIndex) = 1;
    modelIni   = model;
    modeltest = removeRxns(modelIni,EliminatedRxns);
    %Test the model for infeasibility, if it is infeasible, return to the
    %last feasible model and remove and test each reaction individually.
    %If, when removing the reaction there is an infeasible model, then do
    %not remove.
    test = optimizeCbModel(modeltest);
    if isnan(test.f)
        disp('Removing zero reactions caused infeasibility, removing individually...')
        modelRmvdRxns = modelIni;
        counter       = 1;
    for i = 1:size(EliminatedRxns,1)
        currentRxn     = EliminatedRxns(i);
         modelRmvdRxns = removeRxns(modelRmvdRxns,currentRxn);
         test = optimizeCbModel(modelRmvdRxns);
        if isnan(test.f)
            modelRmvdRxns          = modelIni;
            errorMoving(counter,1) = model.rxns(zeroRxns(i));
            counter                = counter +1;
            continue
        end
        modelIni = modelRmvdRxns;
    end
    else
    modelRmvdRxns = modeltest;
    end
end

% Assign reversibilities to thermodynamically allowable rxns
numberOfRxns      = size(modelRmvdRxns.rxns,1);
modelRmvdRxns.rev = (modelRmvdRxns.lb(1:numberOfRxns)<0).*(modelRmvdRxns.ub(1:numberOfRxns)>0);

end
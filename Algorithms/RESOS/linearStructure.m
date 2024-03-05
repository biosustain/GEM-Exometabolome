function modelLP = linearStructure(model)
% Translates MILP problem into LP problem
% Inputs:  model   (MILP formulation)
% Outputs: modelLP (LP formulation)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelLP.numRxns     = numel(model.rxns);
modelLP.obj         = zeros(modelLP.numRxns,1);
modelLP.A           = sparse(model.S);
modelLP.S           = sparse(model.S);
modelLP.rhs         = zeros(size(model.S,1),1);
modelLP.lb          = model.lb(1:modelLP.numRxns);
modelLP.ub          = model.ub(1:modelLP.numRxns);
modelLP.vtype       = 'C';
modelLP.modelsense  = 'min';
modelLP.sense       = '=';
modelLP.rxns        = model.rxns;
modelLP.rxnNames    = model.rxnNames;
modelLP.mets        = model.mets;
modelLP.metNames    = model.metNames;
modelLP.description = model.description;
modelLP.rev = (modelLP.lb<0).*(modelLP.ub>0);
end
function ADSB_parser(modelArg, looplessArg, compactArg)

    %%% STEP 0: Load CobraToolbox and change Solver to Gurobi
    tic
    ti = cputime;
    %cd /zhome/2e/2/164651/cobratoolbox
    %initCobraToolbox()
    %cd /zhome/2e/2/164651/ADSB_sampler_implementation-main/work/modified_code/final_scripts
    cd /home/igor/Documents/Exometabolomic_data_integration/Eric/ADSB_sampler_implementation-main/work/modified_code/final_scripts
    %changeCobraSolver('gurobi')
    addpath('/home/igor/Documents/Exometabolomic_data_integration/Eric/ADSB_sampler_implementation-main/work/modified_code/looplessFluxSampler-master')
    addpath('/home/igor/Documents/Exometabolomic_data_integration/Eric/ADSB_sampler_implementation-main/work/modified_code/looplessFluxSampler-master/fxns')
    addpath('/home/igor/Documents/Exometabolomic_data_integration/Eric/ADSB_sampler_implementation-main/work/modified_code/models')
    addpath('/home/igor/Documents/Exometabolomic_data_integration/Eric/ADSB_sampler_implementation-main/work/modified_code/FinalFunctions')
    addpath('/home/igor/Documents/Exometabolomic_data_integration/Eric/ADSB_sampler_implementation-main/work/modified_code/final_scripts/results')
    %cd /zhome/2e/2/164651/ADSB_sampler_implementation-main/work/modified_code/final_scripts
    cd /home/igor/Documents/Exometabolomic_data_integration/Eric/ADSB_sampler_implementation-main/work/modified_code/final_scripts  
    modelName = flip(split(modelArg,'/')); % New added by Igor
    modelName = split(modelName{1},'.'); % New added by Igor
    modelName = modelName{1}; % New added by Igor
    times = []; % Empty list to store times

    % Measure time for STEP 0
    times = [times append('Time for STEP 0: ', num2str(toc)), ';\n'];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Characterize the Baseline model   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic
    %%% STEP 1: Load model
    myModel = load(modelArg); % Original model load
    fieldNames = fieldnames(myModel); % Get field names to reach nested structure
    model = myModel.(fieldNames{1}); % Get nested model structure
    originalModel = model; % Save original model

    % Convert arguments ASCII --> double
    looplessArg = str2double(looplessArg);
    compactArg = str2double(compactArg);

    % Debug
    disp(modelArg)
    disp(looplessArg)
    disp(compactArg)

    % Measure time for STEP 1
    times = [times append('Time for STEP 1: ', num2str(toc)), ';\n'];

    tic
    %%% STEP 2: Model pre-processing + Sampling using the new algorithm (Adaptive Directions Sampling on a Box - ADSB)
    FBA = optimizeCbModel(originalModel); % This FBA is part of the random sampling protocol it is used to establish the set the boundary of the possible flux predictions to consider
    model.lb(find(model.c)) = 0.5*FBA.f; % Here we find half the values of the maximal predicted fluxes

    % Set up sampling parameters
    options.numSamples = length(originalModel.rxns)*2;
    options.stepsPerPoint = 2e1;
    options.loopless = looplessArg; % Turn on/off loopless option
    options.compact = compactArg; % Turn on/off compact option
    biomass = printObjective(originalModel);

    % Here we take the half maximal model and set up the random sampling
    [ADSB_sampleMetaOutC, newModel] = looplessFluxSampler_v3(model,biomass,options); % NEW FUNCTION
    pointsC = ADSB_sampleMetaOutC.points;

    % Save results
    filename1 = sprintf('results_IM/ADSB_sample_%s_Loopless_%d_Compact_%d.mat', modelName, looplessArg, compactArg); % Modified by Igor modelArg -> modelName
    filename2 = sprintf('results_IM/ADSB_model_%s_Loopless_%d_Compact_%d.mat', modelName, looplessArg, compactArg); % Modified by Igor modelArg -> modelName
    save(filename1, 'ADSB_sampleMetaOutC');
    save(filename2, 'newModel');

    % Measure time for STEP 2
    times = [times append('Time for STEP 2: ', num2str(toc)), ';\n'];

    tic
    %%% STEP 3: Get Fluxes Stats %%%

    % Set model after random sampling
    model = ADSB_sampleMetaOutC;

    % Get the stats about the reaction fluxes
    [FluxesPreclinical, fluxPoints] = getFluxStats(model);

    % Measure time for STEP 3
    times = [times append('Time for STEP 3: ', num2str(toc)), ';\n'];

    fprintf(['Times:\n', times])
    fprintf(['Total time: ',num2str(cputime-ti), ' cputime (s)']);

end
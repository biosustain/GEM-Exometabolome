%% Example ADSB

% load model
load('IMM_iCHO_R_33.mat');

% Prepare parameters for sampling
    options.loopless = 0; % Turn off loopless option
    options.compact = 0; % Turn off compact option
    options.numSamples = 10000; % Sampling Points
    biomass = printObjective(model3_3);
% Perform sampling and save sampling points in a variable
     [results_pool, temp] = looplessFluxSampler_v3(model3_3,biomass,options);
    pointsEcoli = results_pool.points;

% Analyze sampling points with Continuous Method
figure
    disp("For Continuous Method");
    tic
        [minSolutionCon,maxSolutionCon] = EuclideanDistancesMaxMin(pointsEcoli);
    toc
    % Finish Plot
        title("Dynamic Percentile for Continuous Method")

% Analyze sampling points with Discrete Method
figure
    disp("For Continuous Method");
    tic
        [minSolutionDis,maxSolutionDis] = FrequencyHistogramsMaxMin(pointsEcoli);
    toc
    % Finish Plot
        title("Dynamic Percentile for Discrete Method")
%% Example ADSB

% load model
load('RECON1.mat');

% Prepare parameters for sampling
    options.loopless = 0; % Turn off loopless option
    options.compact = 0; % Turn off compact option
    options.numSamples = 10000; % Sampling Points
    
% Perform sampling and save sampling points in a variable
     [results_pool] = looplessFluxSampler(RECON1,options);
    pointsRECON1 = results_pool.points;

% Analyze sampling points with Continuous Method
figure
    disp("For Continuous Method");
    tic
        [minSolutionCon,maxSolutionCon] = EuclideanDistancesMaxMin(pointsRECON1);
    toc
    % Finish Plot
        title("Dynamic Percentile for Continuous Method")

% Analyze sampling points with Discrete Method
figure
    disp("For Continuous Method");
    tic
        [minSolutionDis,maxSolutionDis] = FrequencyHistogramsMaxMin(pointsRECON1);
    toc
    % Finish Plot
        title("Dynamic Percentile for Discrete Method")

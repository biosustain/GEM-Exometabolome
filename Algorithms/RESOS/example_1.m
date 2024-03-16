%% Example gpSampler

% load model
load('RECON1.mat');

% Prepare parameters for sampling
    sampling_time   = 120;
    sampling_points = 10000;

% Perform sampling and save sampling points in a variable
    [results_pool, temp] = gpSampler(RECON1,sampling_points,[],sampling_time);
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
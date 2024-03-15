%% Example gpSampler

% load model
model = load('IMM_iCHO_R_33.mat');

% Prepare parameters for sampling
    sampling_time   = 120;
    sampling_points = 10000;

% Perform sampling and save sampling points in a variable
    [results_pool, temp] = gpSampler(model3_3,sampling_points,[],sampling_time);
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
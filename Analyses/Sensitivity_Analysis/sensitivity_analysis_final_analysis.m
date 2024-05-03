
sample_size = 95;
number_of_tasks = 16;
iterations = 100;

load('Ex2Task_Work_Space.mat');

if ~exist('Plots', 'dir')
   mkdir('Plots')
end

number_of_metabolites = length(ithIndex);

summary_all_patients = cell(2,number_of_metabolites);

for i = 1:sample_size
    filename = sprintf('WS_%d.mat', i); % create the filename
    data = load(filename); % load the file
    data = data.SummaryData;

    modelname =  sprintf('modelPatient_Sampled_%d_Mean.mat', i); % create the model name
    model = load(modelname); % load the file
    model = model.sampleMetaOutC;

    summary = zeros(number_of_metabolites,number_of_tasks);  
    boundaries = [model.lb(ithIndex) model.ub(ithIndex)];

    for j = 1:number_of_metabolites
	    ith_sum = data{j};
	    ith_sum = ith_sum(:,2:end);	        
	    ith_grad_conc = linspace(boundaries(j,1),boundaries(j,2),iterations);
	    [~,ith_gradient] = gradient(ith_sum);
	    ith_gradient(isnan(ith_gradient)) = 0;
	    ith_gradient(isinf(ith_gradient)) = 0;
	    ith_gradient = sum(ith_gradient);
	    ith_gradient = ith_gradient/(boundaries(j,2)-boundaries(j,1)); %slope
	    summary(j,:) = ith_gradient;
    end
    summary_all_patients{1,i} = summary;
    summary(isnan(summary)) = 0;
    summary = normc(summary);
    summary = rescale(summary,-1,1);
    summary_all_patients{2,i} = summary;

    % Create the clustergram  
    cgo=clustergram(summary,'Standardize','Row', 'Colormap', 'redbluecmap');
    %uniq_task_list=unique(join(TasksTable{:,2:4}," - "));
    uniq_task_list=unique(TasksTable{:,4});   
    set(cgo,'Linkage','complete','Dendrogram',3,'RowLabels', TXTUrxn1, 'ColumnLabels', uniq_task_list, 'ColumnLabelsRotate', 16);
    title = addTitle(cgo,sprintf('Patient %d', i),'Color','black');
    title.FontSize = 10;
    
    % Adjust the size of the figure window
    %set(cgo.figureHandle, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 1]);
    set(cgo.figureHandle, 'Position', [0 0 1000 800]); % Adjust the figure size to your preference
    
    % Adjust the paper size to match the figure size
    %set(cgo.figureHandle, 'PaperPositionMode', 'auto');

    % Define the filename base
    filename_base = sprintf('Plots/Sens_An_Patient_%d', i);
    
    % Save the figure in SVG, and PNG formats in the 'Plots' subfolder
    saveas(cgo.figureHandle, [filename_base, '.svg']);
    saveas(cgo.figureHandle, [filename_base, '.png']);
end

%% Prepare data table for further analyses
% Assume 'summary_all_patients' is your 2x95 cell array
numPatients = size(summary_all_patients, 2);
numMetabolites = size(summary_all_patients{1,1}, 1);
numTasks = size(summary_all_patients{1,1}, 2);

% Reshape each matrix into a vector and stack all vectors into a single matrix
dataMatrix = zeros(numPatients, numMetabolites*numTasks);
for i = 1:numPatients
    dataMatrix(i, :) = reshape(summary_all_patients{1,i}, 1, []);
end

dataMatrix(isnan(dataMatrix)) = 0;
dataMatrix = normc(dataMatrix);
dataMatrix = rescale(dataMatrix,-1,1);

patients_group =["C", "A", "C", "D", "D", "D", "A", "A", "B", "D", "C", "D", "D", "C", "D", "C", "D", "C", "C", "B", "B", "D", "A", "D", "D", "D", "D", "D", "B", "C", "C", "D", "D", "C", "D", "B", "B", "C", "C", "C", "C", "C", "B", "D", "C", "D", "D", "B", "C", "D", "B", "A", "A", "D", "D", "B", "D", "A", "C", "D", "C", "C", "C", "C", "D", "D", "D", "D", "D", "A", "B", "D", "A", "D", "D", "A", "C", "D", "C", "C", "D", "C", "D", "D", "B", "C", "D", "D", "D", "A", "D", "D", "C", "A", "C"];

%% Clustering analysis
% Determine the optimal number of clusters using the Elbow Method
maxNumClusters = 10; % Maximum number of clusters to consider
sumOfSquaredDistances = zeros(maxNumClusters, 1);
for k = 1:maxNumClusters
    [~, ~, sumd] = kmeans(dataMatrix, k);
    sumOfSquaredDistances(k) = sum(sumd);
end

optimalNumClusters = 4;

% Define the colors for each group
%groupColors = containers.Map({'A', 'B', 'C', 'D'}, {'orange', 'green', 'black', 'blue'});
%groupColors = containers.Map({'A', 'B', 'C', 'D'}, {'red', 'blue', 'purple', 'darkorange'});
groupColors = containers.Map({'A', 'B', 'C', 'D'}, {'red', [1, 0.5, 0], 'magenta','blue' });

% Perform k-means clustering
[idx, ~] = kmeans(dataMatrix, optimalNumClusters);

% Perform PCA for dimensionality reduction
coeff = pca(dataMatrix);
reducedData = dataMatrix * coeff(:, 1:2);

% Create a scatter plot
fig = figure;
hold on;
colors = {'red', [1, 0.5, 0], 'magenta','blue'}; % colors for the clusters
for i = 1:numPatients
    % Plot points and labels
    scatter(reducedData(i,1), reducedData(i,2), [], groupColors(patients_group(i)), 'filled');
    text(reducedData(i,1), reducedData(i,2), num2str(i), 'Color', groupColors(patients_group(i)));
end


% Draw a convex hull around each cluster
for i = 1:optimalNumClusters
    % Find the points in the cluster
    clusterPoints = reducedData(idx==i,:);
    
    % Calculate the convex hull of the cluster
    K = convhull(clusterPoints(:,1), clusterPoints(:,2));
    
    % Draw the convex hull
    plot(clusterPoints(K,1), clusterPoints(K,2), 'Color', colors{i}, 'LineWidth', 2);
end
hold off;
xlabel('Principal Component 1');
ylabel('Principal Component 2');
%legend('Location','best');
grid on;

% Save the figure in SVG, and PNG formats in the 'Plots' subfolder
filename_base = sprintf('Plots/Clustering_Analysis');
saveas(fig, [filename_base,'.svg']);
saveas(fig, [filename_base,'.png']);


%% PLS-DA
% Convert group labels to numeric
[~,~,groups_numeric] = unique(patients_group);

% Perform PLS-DA
[XL, YL, XS, YS, BETA, PCTVAR] = plsregress(dataMatrix, groups_numeric, 20);  % 10 is the number of components, adjust as needed

% In this code, plsregress function is used to perform PLS-DA. 
% The scores plot is a scatter plot of the first two latent variables and it 
% can be used to visualize the separation between the groups. 
% Each point in the plot corresponds to a patient, and the color of the 
% point corresponds to the group of the patient. If the PLS-DA model is 
% good, patients from the same group should cluster together in the plot.
% The loadings (Xloadings) can be used to understand which 
% “metabolite concentration”/“metabolic task response” are important for 
% the separation of the groups. A high absolute value of a loading for a 
% certain “metabolite concentration”/“metabolic task response” means that 
% this variable is important for the separation of the groups.

% Plot
% Define the color map
%colors = ['r', 'g', 'b', 'y'];  % adjust as needed
colors = [[0,138,216]; [255,109,106]; [68,215,168]; [255,130,0]]/255;  % adjust as needed

% Create a new figure 2D for Scores: C1 vs C2
fig = figure;

% Plot each group with a different color
hold on
for i = 1:4
    scatter(XS(groups_numeric == i, 1), XS(groups_numeric == i, 2), 50, colors(i,:), 'filled')
end
hold off

% Set the labels and title
xlabel(sprintf('1st component (%.2f%% variance)', 100*PCTVAR(1,1)))
ylabel(sprintf('2nd component (%.2f%% variance)', 100*PCTVAR(1,2)))
title('PLS-DA Scores Plot')
legend('Group A', 'Group B', 'Group C', 'Group D') % Create the legend

% Save the figure in SVG, and PNG formats in the 'Plots' subfolder
filename_base = sprintf('Plots/PLS-DA_Scores');
saveas(fig, [filename_base,'.svg']);
saveas(fig, [filename_base,'.png']);



% Create a new figure 2D for Scores: C1 vs C3
fig = figure;

% Plot each group with a different color
hold on
for i = 1:4
    scatter(XS(groups_numeric == i, 1), XS(groups_numeric == i, 5), 50, colors(i,:), 'filled')
end
hold off

% Set the labels and title
xlabel(sprintf('1st component (%.2f%% variance)', 100*PCTVAR(1,1)))
ylabel(sprintf('5th component (%.2f%% variance)', 100*PCTVAR(1,5)))
title('PLS-DA Scores Plot')
legend('Group A', 'Group B', 'Group C', 'Group D') % Create the legend

% Save the figure in SVG, and PNG formats in the 'Plots' subfolder
filename_base = sprintf('Plots/PLS-DA_Scores_C1-C5');
saveas(fig, [filename_base,'.svg']);
saveas(fig, [filename_base,'.png']);






% Determine what are the "metabolite concentration"/"metabolic task
% response" that contributes the most to separate groups
% Find the index of the most important "metabolite concentration"/"metabolic task response"

% Number of most important "metabolite concentration"/"metabolic task response" to find
N = 5;  % adjust as needed
PC = 20;
XL2 = XL(:,1:PC);
% Initialize arrays to store the most important indices and loading values for each component
most_important_indices = zeros(N, size(XL2, 2));
loading_values = zeros(N, size(XL2, 2));

% For each component...
for j = 1:size(XL2, 2)
    % Find the indices of the N most important "metabolite concentration"/"metabolic task response"
    [loading_values(:, j), most_important_indices(:, j)] = maxk(abs(XL2(:, j)), N);
end

% Convert the indices into metabolite and task indices
most_important_metabolites = mod(most_important_indices-1, numMetabolites) + 1;
most_important_tasks = ceil(most_important_indices / numMetabolites);

% Print the most important "metabolite concentration"/"metabolic task response" for each component
for j = 1:size(XL2, 2)
    fprintf('For component %d:\n', j);
    for i = 1:N
        fprintf('The %dth most important metabolite is %d and the task is %d\n', i, most_important_metabolites(i, j), most_important_tasks(i, j));
    end
    fprintf('\n');
end

% Calculate the sum of the absolute values of the loadings across all components
sum_abs_loadings = sum(abs(XL2), 2);

% Find the indices of the N most important "metabolite concentration"/"metabolic task response"
[loading_values, most_important_indices] = maxk(sum_abs_loadings, N);

% Convert the indices into metabolite and task indices
most_important_metabolites = mod(most_important_indices-1, numMetabolites) + 1;
most_important_tasks = ceil(most_important_indices / numMetabolites);

% Print the most important "metabolite concentration"/"metabolic task response"
for i = 1:N
    fprintf('The %dth most important overall metabolite is %d and the task is %d\n', i, most_important_metabolites(i), most_important_tasks(i));
end

% Create labels for the x-axis
labels = cell(N, 1);
for i = 1:N
    labels{i} = sprintf('Metabolite %d, Task %d', most_important_metabolites(i), most_important_tasks(i));
end

% Create a bar plot
fig = figure;
bar(loading_values)
xlabel('Metabolite, Task')
ylabel('Loading Value')
title('Importance of Each Metabolite Concentration/Metabolic Task Response')
xticklabels(labels)
xtickangle(45)  % rotate x-axis labels for better visibility

% Save the figure in SVG, and PNG formats in the 'Plots' subfolder
filename_base = sprintf('Plots/PLS-DA_Relevant_Metabolite-Task');
saveas(fig, [filename_base,'.svg']);
saveas(fig, [filename_base,'.png']);

% Create a new figure 2D for Loadings
% Define a color map with a unique color for each task
colorMap = lines(numTasks);  % adjust as needed

% Create a new figure
fig = figure;

% Plot each "metabolite concentration"/"metabolic task response" with a different color
hold on
for i = 1:numMetabolites*numTasks
    % Calculate the task index
    taskIndex = mod(i-1, numTasks) + 1;
    
    % Scatter plot
    scatter(XL(i, 1), XL(i, 2), 50, colorMap(taskIndex,:), 'filled');
end

% Add labels for the N most important "metabolite concentration"/"metabolic task response"
for i = 1:N
    text(XL(most_important_indices(i), 1), XL(most_important_indices(i), 2), ...
         sprintf('Metabolite %d', most_important_metabolites(i)), ...
         'Color', 'k', 'FontSize', 10);
end

% % Add labels for the N most important "metabolite concentration"/"metabolic task response"
% for i = 1:N
%     text(XL(most_important_indices(i), 1), XL(most_important_indices(i), 2), ...
%          sprintf('Metabolite %d, Task %d', most_important_metabolites(i), most_important_tasks(i)), ...
%          'Color', 'k', 'FontSize', 10);
% end

% Set the labels and title
xlabel(sprintf('1st loading (%.2f%% variance)', 100*PCTVAR(2,1)))
ylabel(sprintf('2nd loading (%.2f%% variance)', 100*PCTVAR(1,1)))
title('PLS-DA Loadings Plot')

% Create a legend
legend(arrayfun(@(x) sprintf('Task %d', x), 1:numTasks, 'UniformOutput', false))

% Save the figure in SVG, and PNG formats in the 'Plots' subfolder
filename_base = sprintf('Plots/PLS-DA_Loadings');
saveas(fig, [filename_base,'.svg']);
saveas(fig, [filename_base,'.png']);

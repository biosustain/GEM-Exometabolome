case_models = readtable('Patient_Models.xls','sheet','Sheet1');

for i=1:height(case_models)
    % Load case-specific data
    model_max = load(case_models{i,1}{1});
    model_mean = load(case_models{i,2}{1});
    model_min = load(case_models{i,3}{1});
    overall_points = [model_min.sampleMetaOutC.points,model_mean.sampleMetaOutC.points,model_max.sampleMetaOutC.points];
    
    % Perform PCA on the concatenated matrix
    [coeff, score, ~, ~, explained] = pca(overall_points');

    i_str = num2str(i);
    
    groups = {1:6012, 6013:12024, 12025:18036};
    group_names = {'min group', 'mean group', 'max group'};
    %colors = {'r', 'g', 'b'};
    colors = {[1.00 0.57 0.48], [0.46 1.00 0.57], [0.38 0.38 1.00]};  % softer colors. Look here to choose colors: https://rgbcolorcode.com/

    % Plot the first two principal components in a 2D plot
    figure_2D = figure('visible','off');
    hold on;
    for j = 1:length(groups)
        scatter(score(groups{j}, 1), score(groups{j}, 2), 10, colors{j}, 'filled');
    end
    xlabel(['PC1 (' num2str(explained(1), '%.2f') '%)']);
    ylabel(['PC2 (' num2str(explained(2), '%.2f') '%)']);
    title(['PCA of Patient ' i_str]);  % concatenate string variable to the title
    legend(group_names, 'Location', 'best');
    grid on;
    hold off;

    % Save the plot as SVG and PNG files 2D
    filename = ['Patient_' i_str '_PCA'];
    saveas(figure_2D, filename, 'svg');
    saveas(figure_2D, filename, 'png');
    
%     % Create a new figure without displaying it for 3D plot
%     figure_3D = figure('visible','off');
%     hold on;
%     for j = 1:length(groups)
%         scatter3(score(groups{j}, 1), score(groups{j}, 2), score(groups{j}, 3), 10, colors{j}, 'filled');
%     end
%     xlabel(['PC1 (' num2str(explained(1), '%.2f') '%)']);
%     ylabel(['PC2 (' num2str(explained(2), '%.2f') '%)']);
%     zlabel(['PC3 (' num2str(explained(3), '%.2f') '%)']);
%     title(['3D PCA of Patient ' i_str]);
%     legend(group_names, 'Location', 'best');
%     grid on;
%     hold off;
%     
%     % Save the 3D plot as SVG and PNG files
%     filename = ['Patient_' i_str '_3D_PCA'];
%     saveas(figure_3D, filename, 'svg');
%     saveas(figure_3D, filename, 'png');

end
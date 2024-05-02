num_files = 95;

% Loop over each file
for i = 1:num_files

    % Create the file name
    file_name = sprintf('sensitivity_analysis_%d', i);
    
    % Execute the file
    run(file_name);    
end

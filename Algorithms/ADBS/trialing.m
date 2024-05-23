% function trialing(arg1, arg2)

%     cd /zhome/2e/2/164651/ADSB_sampler_implementation-main/work/modified_code/final_scripts
%     addpath('/zhome/2e/2/164651/ADSB_sampler_implementation-main/work/modified_code/final_scripts')

%     % Print the current working directory
%     disp(pwd);

%     % Convert input arguments from string to whatever type you need
%     arg1 = str2double(arg1);
%     arg2 = str2double(arg2);
    
%     % Display
%     disp(['Argument 1: ', num2str(arg1)]);
%     disp(['Argument 2: ', num2str(arg2)]);
% end

function output = trialing(modelName, looplessOption, compactOption)
    
    % Your function code
    looplessOption = str2double(looplessOption);
    compactOption = str2double(compactOption);

    % Generate output filename based on input arguments
    filename = sprintf('%s_Loopless_%d_Compact_%d.mat', modelName, looplessOption, compactOption);
    output = 1234
    
    % Save the output to a .mat file
    save(filename, 'output');
end
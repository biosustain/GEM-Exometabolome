function [centralSolution,furthestSolution] = SamplingWithNoKOsMaxMin(model,samplingTime,ContinuousOrDiscrete,MILPproblem)
%This function performs a normal sampling, returning the central and
%furthest solution. If there is a MILP problem, it is possible to use it.
%ContinuosOrDiscrete is a variable used to determine which analysis
%will be executed
        % 'Continuous' for Continuos distance
        % 'Discrete' for Discrete distance
cellArrayStrings2Compare = [{'Continuous'};{'Discrete'}];
[ContinuousOrDiscrete,~] = bestMatchString(cellArrayStrings2Compare,ContinuousOrDiscrete);
ContinuousOrDiscrete = ContinuousOrDiscrete{1};
%Find if sampling is done with or without MILPproblem
if ~exist('MILPproblem', 'var') || isempty(MILPproblem)
  MILPproblem = model;
  disp('MILP Problem not detected, sampling will be done with the original model');
else
    MILPproblem.S = full(MILPproblem.A);    
end
        
 %1. Generate a pool of solutions and save the solutions given in a cell
    %array
        %Minimum points required is always the amount of reactions * 2
            minimumPoints = length(MILPproblem.lb)*2;
            try
                test = zeros(length(MILPproblem.lb),minimumPoints);
                clear test
            catch
               disp('Array using MILP problem will lead to out of memory error...Using less warming points')
               try
                   test = zeros(length(MILPproblem.lb),minimumPoints/2);
                   minimumPoints = minimumPoints/2;
                   clear test
               catch
                   test = zeros(length(MILPproblem.lb),round(minimumPoints/4));
                   minimumPoints = round(minimumPoints/4);
                   clear test
               end
            end
            
        %Create pool, the sampling time is defined in the input of the
            %function
            %[~,coreNumbers] = evalc('feature(''numcores'')');

        results_pool   = gpSamplerTolerance(MILPproblem,minimumPoints,[],samplingTime*10,50);%,coreNumbers);
        Relevant_Array = results_pool.points(1:length(model.rxns),:);
        templogic      = Relevant_Array(find(model.c),:) >= 0;
        Relevant_Array = Relevant_Array(:,templogic);
        clear results_pool
   %Function that calculates distance according to the input
           switch ContinuousOrDiscrete
               case 'Continuous'
                   
                disp('Continuous Distance...')
                    [centralSolution,furthestSolution] = EuclideanDistancesMaxMin(Relevant_Array); 
                disp('...')
                clear Relevant_ArrayNormal
               case 'Discrete'  
                    %Relevant_ArrayTransposed = tall(Relevant_Array');
                disp('Discrete Distance...')
                    [centralSolution,furthestSolution] = FrequencyHistogramsMaxMin(Relevant_Array); 
                disp('...')
                %columns_relevant_results = columns_relevant_results';
                clear Relevant_ArrayTransposed
           end
        
end
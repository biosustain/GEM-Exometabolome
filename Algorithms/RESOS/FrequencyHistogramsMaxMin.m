function [centricSolution,furthestSolution] = FrequencyHistogramsMaxMin(ArrayOfSolutions,userDefinedPercentile)   
%%%%%%%%%%%%%%%%%%% Fernando Silva-Lance 2021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate total frequency of the value of each reaction in 
          %every  solution. Summing the frequencies per reaction value will
          %give the "proximity" of one reaction conpared to the rest.
          %Then, a percentile is used to delimit the space of solutions and
          %eliminate outliers.
   %Input:
        %UserDefinedPercentile: preselected value that will be used to
            %delimit the space of solutions. If this value does not exist,a dynamic
            %percentile will be computed.
        %ArrayOfSolutions: Results of a sampling, just the points.
   %Output:
        %centricSolution: Most centric solution within the space of solutions
        %after filtrating with the percentile
        %furthestSolution: Furthest solution within the space of solutions
        %after filtrating with the percentile   


%UserDefinedPercentile is a preselected value that will be used to
   %delimit the space of solutions. If this value does not exist,a dynamic
   %percentile will be computed.
   if ~exist('userDefinedPercentile','var') || isempty(userDefinedPercentile)
       userDefinedPercentile = 0;
       disp('Calculating a dynamic percentile');
   end
   
   
%Histogram with sum of the bins = 1.
     validSolutions    = size(ArrayOfSolutions,2);
     reactionsInModel  = size(ArrayOfSolutions,1);
     %Columns are now rows and rows are now columns
        frequencyMatrix   = zeros(reactionsInModel,validSolutions,'single');

      %For every row/reaction, do an histogram and assign the corresponding
      %frequency value according to the bin each value belongs to. 
      for i = 1:reactionsInModel
          hist_Data                   = ArrayOfSolutions(i,:);
          hist_Data                   = rescale(hist_Data);
         [hist_frequency,hist_bins,~] = histcounts(hist_Data,'Normalization','pdf');
         currentReaction              = zeros(1,validSolutions,'single');
         
         for j = 2:size(hist_bins,2)
            lowEdge  = hist_bins(:,j-1);
            highEdge = hist_bins(:,j);
            indexValuesCurrentBin = hist_Data >= lowEdge & hist_Data <= highEdge;
            currentReaction(:,indexValuesCurrentBin) = hist_frequency(:,j-1);
        end
        
        frequencyMatrix(i,:) = currentReaction;
        
      end

    %4 Sum frequencies per solution and Filter with the desired percentile, could be 90 
        FrequenciesSumVector = sum(frequencyMatrix,1) * -1;
       
       if userDefinedPercentile == 0
            PercentileDenseSpace             = DynamicPercentileSampling(FrequenciesSumVector);
            valueForFiltering                = prctile(FrequenciesSumVector,PercentileDenseSpace);
            columns_relevant_results         = find(FrequenciesSumVector <= valueForFiltering); 
            
       else
            valueForFiltering                = prctile(FrequenciesSumVector,userDefinedPercentile);
            columns_relevant_results         = find(FrequenciesSumVector <= valueForFiltering);
       end
         
            
        %8. Create a new array only considering those solutions within the
        % most dense space of solutions
            ArrayOfSolutions     = ArrayOfSolutions(:,columns_relevant_results);
            FrequenciesSumVector = FrequenciesSumVector(columns_relevant_results);
            
        %9. Within the most dense space of solutions, find the solution
        %with the maximum euclidean distance and the one with the minimum
        %euclidean distance
            [~,maxDistanceSolution] = max(FrequenciesSumVector);
            [~,minDistanceSolution] = min(FrequenciesSumVector);
            
       %10. Recover the relevant solutions
            furthestSolution = ArrayOfSolutions(:,maxDistanceSolution);
            centricSolution  = ArrayOfSolutions(:,minDistanceSolution);
            distanceFurthestCentral = norm(furthestSolution - centricSolution);
          %This is Saved for future samplings, this structure is used as an
          %input of the function named: gpSamplerPercentile
            structureForNextSamplings.numberOfSolutions       = length(columns_relevant_results);
            structureForNextSamplings.distanceFurthestCentral = distanceFurthestCentral;
            structureForNextSamplings.modeOfReactions         = mode(ArrayOfSolutions,2);
            save('DataForNextSamplings','structureForNextSamplings')
            
end
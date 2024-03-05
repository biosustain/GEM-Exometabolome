function [centricSolution,furthestSolution] = EuclideanDistancesMaxMin(ArrayOfSolutions,userDefinedPercentile)
%%%%%%%%%%%%%%%%%%% Fernando Silva-Lance 2021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate euclidean distance for every solution against ALL of
          %the other solutions. The outcome should be a symmetrical matrix
          %with a diagonal of 0s.
          %Every solution´s euclidean distances will then be summed
          %up to know the distance with respect to every other solution.
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
        
   if ~exist('userDefinedPercentile','var') || isempty(userDefinedPercentile)
       userDefinedPercentile = 0;
       disp('Calculating a dynamic percentile');
   end
      numberofRelevantSolutions = size(ArrayOfSolutions,2);
      TotalEuclideanDistance    = zeros(numberofRelevantSolutions,1,'single');
      nonZeroElements           = (sum(0:numberofRelevantSolutions-1));
      
     %Main Loop, euclidean distances are calculated. 
                    %If you want to save memory, uncomment the 'single',
                    %however this is not accepted in some versions of
                    %matlab when calling the sparse matrix.
           columnIndexes = zeros(nonZeroElements,1);%,'single');
           rowIndexes    = zeros(nonZeroElements,1);%,'single');
        %calculedED can´t be single because sparse only admits double
           calculedED    = zeros(nonZeroElements,1);
           counter = 1;
           for column = 1:numberofRelevantSolutions
                currentSolution = ArrayOfSolutions(:,column);
                %Creating a factor for jumping results already calculated
                    factor = 1 + column;
                    row =  factor;
                %While to calcule the euclidean distance of the current
                %solution against every other solution
                while row <= numberofRelevantSolutions 
                    
                    calculedED(counter,1) = norm(currentSolution - ArrayOfSolutions(:,row));
                    columnIndexes(counter,1) = column;
                    rowIndexes(counter,1)    = row;
                    counter = counter + 1;
                    row     = row + 1;
        
                end
           end
           
           %Create SparseMatrix
           euclideanPerSolution = sparse(rowIndexes,columnIndexes,calculedED);
           
           %Making the sum for each solution. Sum column N, then Sum row N
           %to retrieve those values skipped in the loop (skipped to not make
           %the calculations twice)
           Index1 = 1;
           Index2 = 2;
           border = numberofRelevantSolutions;
           for i = 1: border
               columnIndex = 1:Index1;
               rowIndex = Index2:border;
               if i == 1
                   TotalEuclideanDistance(i) = full(sum(euclideanPerSolution(rowIndex,i)));
                   Index2 = Index2 + 1;
                   continue
               end
               
               if i == border
                   TotalEuclideanDistance(i) = full(sum(euclideanPerSolution(i,columnIndex)));
                   continue
               end
               TotalEuclideanDistance(i) = full(sum(euclideanPerSolution(rowIndex,i)))+ full(sum(euclideanPerSolution(i,columnIndex)));
               Index2 = Index2 + 1;
               Index1 = Index1 + 1;
           end
            
        %6. Normalizing the euclidean distance
            TotalSum                    = sum(TotalEuclideanDistance);
            NormalizedEuclideanDistance = TotalEuclideanDistance/TotalSum;
            
        %7. Filtering the distances with a percentile, those above the
        %percentile are not considered to be within the most dense space of
        %solutions.

       if userDefinedPercentile == 0
            PercentileDenseSpace             = DynamicPercentileSampling(NormalizedEuclideanDistance);
            valueForFiltering                = prctile(NormalizedEuclideanDistance,PercentileDenseSpace);
            columns_relevant_results         = find(NormalizedEuclideanDistance <= valueForFiltering); 
            
       else
            valueForFiltering                = prctile(NormalizedEuclideanDistance,userDefinedPercentile);
            columns_relevant_results         = find(NormalizedEuclideanDistance <= valueForFiltering);
       end
            
        %8. Create a new array only considering those solutions within the
        % most dense space of solutions
            ArrayOfSolutions            = ArrayOfSolutions(:,columns_relevant_results);
            NormalizedEuclideanDistance = NormalizedEuclideanDistance(columns_relevant_results);
            
        %9. Within the most dense space of solutions, find the solution
        %with the maximum euclidean distance and the one with the minimum
        %euclidean distance
            [~,maxDistanceSolution] = max(NormalizedEuclideanDistance);
            [~,minDistanceSolution] = min(NormalizedEuclideanDistance);
            
       %10. Recover the relevant solutions
            furthestSolution        = ArrayOfSolutions(:,maxDistanceSolution);
            centricSolution         = ArrayOfSolutions(:,minDistanceSolution);
            distanceFurthestCentral = norm(furthestSolution - centricSolution);
            
          %This is Saved for future samplings, this structure is used as an
          %input of the function named: gpSamplerPercentile
            structureForNextSamplings.numberOfSolutions       = length(columns_relevant_results);
            structureForNextSamplings.distanceFurthestCentral = distanceFurthestCentral;
            structureForNextSamplings.modeOfReactions         = mode(ArrayOfSolutions,2);
            
            save('DataForNextSamplings','structureForNextSamplings')
            
end
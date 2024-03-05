function [FiltratedResults,columns_relevant_results] = EuclideanDistancesSmallModels(ArrayOfSolutions,userDefinedPercentile)
%%%%%%%%%%%%%%%%%%% Fernando Silva-Lance 2021 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate euclidean distance for every solution against ALL of
          %the other solutions. The outcome should be a simetrical matrix
          %with a diagonal of 0s.
          %Every solution´s euclidean distances will then be summed
          %up to know the distance with respect to every other solution.
          %Then, a percentile is used to delimit the space of solutions and
          %eliminate outliers.
   %Input:
        %UserDefinedPercentile: preselected value that will be used to
            %delimit the space of solutions. If this value does not exist,a dynamic
            %percentile will be computed.
        %ArrayOfSolutions: Results of a sampling, just the points.
   %Output:
        %FiltratedResults : Array containing every result within the
            %percentile.
        %columns_relevant_results :  Column index of the original array that are within the
            %percentile.
        
   if ~exist('userDefinedPercentile','var') || isempty(userDefinedPercentile)
       userDefinedPercentile = 0;
       disp('Calculating a dynamic percentile');
   end
          
          %As array of solutions is a tall array, gather function needs to be used
            numberofRelevantSolutions = size(ArrayOfSolutions,2);
            TotalEuclideanDistance    = zeros(numberofRelevantSolutions,1,'single');
            nonZeroElements           = (sum(0:numberofRelevantSolutions-1));
            %tempArrayOfSolutions      = gather(ArrayOfSolutions);
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
                    rowIndexes(counter,1) = row;
                    counter = counter + 1;
                    row = row + 1;
        
                end
           end
           %Create SparseMatrix
           euclideanPerSolution = sparse(rowIndexes,columnIndexes,calculedED);
           
           %Making the sum for each solution Sum column N, then Sum row N
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
            
       if userDefinedPercentile == 0
            PercentileDenseSpace             = DynamicPercentileSampling(NormalizedEuclideanDistance);
            valueForFiltering                = prctile(NormalizedEuclideanDistance,PercentileDenseSpace);
            columns_relevant_results         = find(NormalizedEuclideanDistance <= valueForFiltering); 
            
       else
            valueForFiltering                = prctile(NormalizedEuclideanDistance,userDefinedPercentile);
            columns_relevant_results         = find(NormalizedEuclideanDistance <= valueForFiltering);
       end
            
        %7. Filter solutions through percentile and create a new table
            FiltratedResults = ArrayOfSolutions(:,columns_relevant_results);

end
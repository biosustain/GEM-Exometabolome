function [modelNewBoundaries] = IntegrateRelativeMetabolomicData(model, percent25, percent75, RxnNames,RelExValues)
%Incorporates boundaries based on exametabolomic data
% Relative ExametabolomicFilename: Filename of the excel spreadsheet which
% includes the rxn names in the first column and exametabolomic values in
% the second column.
%%NOTE: Sometimes errors can be returned depending on the format of the
%%reactions located in RxnNames.
    %Make sure they are all in a vertical cell array, if there is n error
    %of the type:
        %"Input A of class cell and input B of class cell must be cell
        %arrays of character vectors, unless one is a character vector.%
     %Try changing The parenthesis for curly braces or the other way around
     %in those lines where RxnNames is indexed, for example:
                    %RxnNames{i}  change to RxnNames(i)
                    %RxnNames(i)  change to RxnNames{i}
%%%%%%%%%%%%%%%%%%%         Igor Marin de Más         %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Function by: Fernando Silva-Lance 2021 %%%%%%%%
    upperB = zeros(size(RxnNames,1),2);
    lowerB = zeros(size(RxnNames,1),2);
%1. For loop to change and test that boundaries difference is not 0
counter = 1;
    for i = 1:size(RxnNames,1)
        index       = find(ismember(model.rxns,RxnNames(i,1)));
        if isempty(index)
            notInModel(counter) = i;
            counter = counter + 1;
            disp(['Reaction not in model: ' RxnNames{i,1}])
            continue
        else
            notInModel(counter) = 0;
        end
        upperB(i,1) = percent75(index);
        lowerB(i,1) = percent25(index);
        if (abs(lowerB(i,1)) - abs(upperB(i,1)) == 0) 
            upperB(i,1) = 10;
            lowerB(i,1) = -10;
        end
        if  isinf(RelExValues(i)) | isnan(RelExValues(i))
            RelExValues(i) = 1;
        end
        upperB(i,2) = upperB(i,1) + (upperB(i,1)/abs(upperB(i,1)))*(upperB(i,1)*(RelExValues(i)-1));
        lowerB(i,2) = lowerB(i,1) + (lowerB(i,1)/abs(lowerB(i,1)))*(lowerB(i,1)*(RelExValues(i)-1));
     
        if lowerB(i,1)<0 && lowerB(i,2)>0
            lowerB(i,3) = lowerB(i,1)/RelExValues(i);
        end
        if upperB(i,1)>0 && upperB(i,2)<0
            upperB(i,2) = upperB(i,1)*RelExValues(i);
        end  
    end
    
%Create cell array with an uniform series of numbers which will be used to
%relax the new boundaries in case the first ones boundaries dont work
    Int = (cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);

    modelNewBoundaries = model;
%Forloop to change and test every new boundary
    for i = 1:size(RxnNames,1)
        if sum(i == notInModel) ~= 0
            continue
        end
        %Save an initial model which wont suffer any changes until the
        %boundaries are tested and approved
        modelInitial = modelNewBoundaries;
        
        %Change boundaries and test
        modelNewBoundaries = changeRxnBounds(modelNewBoundaries,RxnNames(i,1),lowerB(i,2),'l');
        modelNewBoundaries = changeRxnBounds(modelNewBoundaries,RxnNames(i,1),upperB(i,2),'u');
        testBoundaries     = optimizeCbModel(modelNewBoundaries);
    
        %If there is not an Nan in the objective function, save boundaries in
        %both models
        if ~isnan(testBoundaries.f)
            modelInitial = modelNewBoundaries;
        else
        %Relax the boundaries using a factor located in Int, then test
        %again
            modelNewBoundaries = modelInitial;
            flag = 0;
            j = 1;
        
        %While loop to keep relaxing boundaries until they give a solution
        %or until the factors in Int are all used
            while flag == 0
                
                modelNewBoundaries = changeRxnBounds(modelNewBoundaries,RxnNames(i,1),(lowerB(i,2)-(lowerB(i,2)/abs(lowerB(i,2)))*(lowerB(i,2)*Int(j))),'l');
                modelNewBoundaries = changeRxnBounds(modelNewBoundaries,RxnNames(i,1),(upperB(i,2)+(upperB(i,2)/abs(upperB(i,2)))*(upperB(i,2)*Int(j))),'u');                   
                testBoundaries = optimizeCbModel(modelNewBoundaries);
%               for testing purposes, otherwise comment 
%               i
%               j
%               TestTmpModel.f   

            %If there is not an Nan in the objective function, save boundaries in
            %both models and leave loop
            if ~isnan(testBoundaries.f)==1
                modelInitial = modelNewBoundaries;
                flag=1;
            
            %Elseif the j is equal to the length of Int, then stop relaxing the
            %boundaries and leave loop
            elseif j==length(Int)
                disp(['Boundaries not relaxed succesfully, reaction: ' RxnNames{i,1}])
                modelNewBoundaries = modelInitial;
                flag = 1;  
            else
                modelNewBoundaries=modelInitial;
                j = j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
            end 
        end
    end
finalTest = optimizeCbModel(modelNewBoundaries);
end
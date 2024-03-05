function [mostSimilarString,idxBestOption] = bestMatchString(cellArrayStrings2Compare,baseString)
%Compares one base string to a set of possible similar strings. 
%The best fit is chosen based on a points system.
%The base string is compared to the possible matches by testing substrings
%of each, which start from one letter until n-1 letters. 
%Inputs 
    %cellArrayStrings2Compare: Vertical cell array including the possible
        %matches in a char format.
    %Base string: Main word, it must be in a char format.
%Outputs
    %mostSimilarString: Cell array containing the string that is more
        %similar to the base string
    %idxBestOption: Index of the best match in the variable
    %"cellArrayStrings2Compare".

%%%%%%%%%%%%%%%%%%%%%Fernando Silva-Lance 2021 %%%%%%%%%%%%%%%%%%%%%%%%%%%


%Test variables
testCellArray = iscell(cellArrayStrings2Compare);
if testCellArray == 1
    if size(cellArrayStrings2Compare,2) > 1
       error('Words should be arranged in a vertical way in the cell array');
    end
    if ~ischar(cellArrayStrings2Compare{1,1})
        error('The words are not in char format')
    end      
end
if ~ischar(baseString)
    error('Base string is not in char format')
end

pointsStrings       = zeros(size(cellArrayStrings2Compare,1),1);
baseStringOriginal  = baseString;

%Loop that goes through every word to compare the base string with 
for numberOfStrings = 1:size(cellArrayStrings2Compare,1)
    currentString2Compare = cellArrayStrings2Compare{numberOfStrings,1};
    %Lengths must be the same
    try
        currentString2Compare = currentString2Compare(1:length(baseString));
    catch
        baseString = baseString(1:length(currentString2Compare));
    end
    
    %Loop for changing the number of letters compared in each string
    for i = 1:length(baseString) - 1
        
        counter = 1;
        %Loop for changing the substring of the base string
        for j = i:length(baseString)
            currentBaseStringLetters = baseString(counter:j);
            counter = counter + 1;
            
            %Loop for making the comparison of each possible substring using 
            %the current number of letters compared
            %Stores the points
            counter2 = 1;
            for k = i:length(baseString)
                subString = currentString2Compare(counter2:k);
                counter2 = counter2 + 1;
                pointsStrings(numberOfStrings) = pointsStrings(numberOfStrings) + strncmpi(subString,currentBaseStringLetters,i);  
            end
            
        end
        
    end
    baseString = baseStringOriginal;
    
end
    %Which one has the most points
    [~,idxBestOption]    = max(pointsStrings);
    mostSimilarString    = cellArrayStrings2Compare(idxBestOption);
end
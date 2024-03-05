function [modelNewExchangeReactions,FinalRelativeTableData] = IntegrateCovidData(model)
%%%%%%%%%%%%%%%%%%%%%% COVID DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%
%IMPORTANT: THIS SCRIPT WORKS ONLY WITH THE FILES NAMED 'PlaceboDrug.csv' and 'ControlPatients.csv','keggIDS.xlsx'
%CHANGE THEM THE LEAST POSSIBLE
%NOTE IT IS POSSIBLE THAT READCELL FUNCTION IS NOT AVAILABLE IN MATLAB 2017
%VERSION
%% Recover Drug and patients Data and Means
    PlaceboDrugData     = readcell('PlaceboDrug.csv');
    
    for i = 5:size(PlaceboDrugData,2)
    for j = 2:size(PlaceboDrugData,1)
        currentNumber = PlaceboDrugData{j,i};
        if ischar(currentNumber)
            currentNumber = strrep(currentNumber,',','.');
            PlaceboDrugData{j,i} = str2double(currentNumber);
        else
            PlaceboDrugData{j,i} = currentNumber;
        end
    end
    end
%Check number of NaN per column, if the number is at least one third of the
%total rows, then leave the reaction unconstraint (here I will put every
%value as NaN, when ratios are being introduced in the model, that is when
%the reaction will be left unconstraint)
    for i = 5:size(PlaceboDrugData,2)
         tempVar = cell2table(PlaceboDrugData(2:end,i));
         tempVar = tempVar{:,1};
         if sum(isnan(tempVar)) > round(size(tempVar,1)/3)
             PlaceboDrugData(2:end,i) = num2cell(NaN(size(tempVar,1),1));
         end
    end
 %Logical vectors
 % 0 = ALIVE
 % 1 = DEATH
 PlaceboDrugData = cell2table(PlaceboDrugData(2:end,:), 'VariableNames',PlaceboDrugData(1,:));
 logicalVectorTreated         = ismember(PlaceboDrugData{:,4},{'Ilomedin'}) & ismember(PlaceboDrugData{:,3},{'24H'});
 logicalVectorPlacebo         = ismember(PlaceboDrugData{:,4},{'Placebo'})  & ismember(PlaceboDrugData{:,3},{'24H'});
 logicalVectorPlaceboSurvived = ismember(PlaceboDrugData{:,4},{'Placebo'}) & PlaceboDrugData{:,2} == 0 & ismember(PlaceboDrugData{:,3},{'24H'});
 logicalVectorDeath           = PlaceboDrugData{:,2} == 1 & ismember(PlaceboDrugData{:,3},{'BASELINE'});
 logicalVectorTreatedDeath    = ismember(PlaceboDrugData{:,4},{'Ilomedin'}) & PlaceboDrugData{:,2} == 1 & ismember(PlaceboDrugData{:,3},{'BASELINE'});
 logicalVectorPlaceboDeath    = ismember(PlaceboDrugData{:,4},{'Placebo'}) & PlaceboDrugData{:,2} == 1 & ismember(PlaceboDrugData{:,3},{'BASELINE'});
 %Means
 meanTreated         =  mean(PlaceboDrugData{logicalVectorTreated,5:end},1,'omitnan');
 meanPlacebo         =  mean(PlaceboDrugData{logicalVectorPlacebo,5:end},1,'omitnan');
 meanPlaceboSurvived =  mean(PlaceboDrugData{logicalVectorPlaceboSurvived,5:end},1,'omitnan');
 meanDeath           =  mean(PlaceboDrugData{logicalVectorDeath,5:end},1,'omitnan');
 meanTreatedDeath    =  mean(PlaceboDrugData{logicalVectorTreatedDeath,5:end},1,'omitnan');
 meanPlaceboDeath    =  mean(PlaceboDrugData{logicalVectorPlaceboDeath,5:end},1,'omitnan');





%% Recover Control Patients Data and Mean
controlPatients = readcell('ControlPatients.csv');
controlPatients = controlPatients(:,4:end);
for i = 1:size(controlPatients,2)
    for j = 2:size(controlPatients(:,i),1)
        currentNumber = controlPatients{j,i};
        currentNumber = strrep(currentNumber,',','.');
        controlPatients{j,i} = str2double(currentNumber);
    end
end
controlPatients = cell2table(controlPatients(2:end,:), 'VariableNames',controlPatients(1,:));
meanControl     = mean(controlPatients{:,:},1,'omitnan');


 
 %% Compute ratio
    RatioTreated         =  (meanTreated         ./ meanControl)';
    RatioPlacebo         =  (meanPlacebo         ./ meanControl)';                
    RatioPlaceboSurvived =  (meanPlaceboSurvived ./ meanControl)';
    RatioDeath           =  (meanDeath           ./ meanControl)';
    RatioTreatedDeath    =  (meanTreatedDeath    ./ meanControl)';
    RatioPlaceboDeath    =  (meanPlaceboDeath    ./ meanControl)';
 

 
 %% Recover keggIDstrings
 %KEGGIds of the metabolites measured in the data
keggIDsData             = readcell('keggIDS.xlsx')';
    excRxns             = model.rxns(findExcRxns(model));
    metListExcRxns      = findMetsFromRxns(model, excRxns);
    metNamesListExcRxns = model.metNames(ismember(model.mets,metListExcRxns));
    

%The following file is located in human1 repository, it includes relevant information
%about the metabolites of the model, for example, kegg ids.
    filename = 'metsHuman1.txt';
    fid      = fopen(filename);
    C        = textscan(fid,'%s','delimiter','\n');
    fclose(fid);
    SeparatedColumns = cell(size(C{1,1},1),1);
    for i=2:size(C{1},1)
        SeparatedColumns{i,1} = split(C{1,1}(i))';
        for j = [1,4]
            if j == 4
                keggIDsModel{i-1,2} = SeparatedColumns{i,1}{j}(2:end-1); 
            else
                keggIDsModel{i-1,1} = SeparatedColumns{i,1}{j}(2:end-1);
            end
        end  % find empty spaces
    end
    
%Update Metabolite names so they all are in the same format
    for i=1:size(keggIDsModel,1)
        TempName = keggIDsModel{i,1}; 
        if contains(TempName, 'p')
            keggIDsModel{i,1} = strrep(TempName,'p','p[p]');
        end
        if contains(TempName, 'c')
            keggIDsModel{i,1} = strrep(TempName,'c','c[c]');
        end
        if contains(TempName, 's')
            keggIDsModel{i,1} = strrep(TempName,'s','s[s]');
        end
        if contains(TempName, 'l')
            keggIDsModel{i,1} = strrep(TempName,'l','l[l]');
        end
        if contains(TempName, 'm')
            keggIDsModel{i,1} = strrep(TempName,'m','m[m]');
        end
        if contains(TempName, 'r')
            keggIDsModel{i,1} = strrep(TempName,'r','r[r]');
        end
        if contains(TempName, 'g')
            keggIDsModel{i,1} = strrep(TempName,'g','g[g]');
        end   
        if contains(TempName, 'n')
            keggIDsModel{i,1} = strrep(TempName,'n','n[n]');
        end           
    end
    
%Create a table only including exchange reactions
    keggIDsModelExcRxns = keggIDsModel(ismember(keggIDsModel(:,1),metListExcRxns),:);

    
    Data(:,1)           = keggIDsData';
    for i=1:size(Data,1)
        TempName = Data{i,1}; 
        if contains(TempName, '_')
            Data{i,1} = strrep(TempName,'_','');
        end
    end
    

% Map KeggIDs To metabolites
    for i = 1:size(Data,1)
        currentID = Data(i,1);
        %First Map with ExcRxns mets
            currentMets = keggIDsModelExcRxns(ismember(keggIDsModelExcRxns(:,2),currentID),1);
            if size(currentMets,1) == 1
                Data(i,2)   = currentMets;
                continue
            elseif size(currentMets,1) > 1
                Data{i,2}   = currentMets;
                continue
            end
          %Now rest of reactions mets
           currentMets = keggIDsModel(ismember(keggIDsModel(:,2),currentID),1);
           Data{i,2}   = currentMets;
    end

% Map Metabolites to reactions
    for i = 1:size(Data,1)
        if isempty(Data{i,2})
            continue
        end
        if size(Data{i,2},1)>1
            for j = 1:size(Data{i,2},1)
                temptMets = Data{i,2};
                [rxnList, ~] = findRxnsFromMets(model, temptMets(j,1));
                rxnList = rxnList(ismember(rxnList,excRxns));
                Data{i,3}{j,1} = rxnList;
            end
            continue
        end
       [rxnList, ~] = findRxnsFromMets(model, Data{i,2});
       rxnList = rxnList(ismember(rxnList,excRxns));
       Data{i,3} = rxnList;
    end


%% As KeggIDs are not available for every reaction
% Match metabolites using Common metabolite name 
commonNames = PlaceboDrugData.Properties.VariableNames;
commonNames = commonNames(5:end);
    for i = 1:size(Data,1)
      if ~isempty(Data{i,2})
          continue
      end
      currentString = commonNames{i};
      [bestMatch,idx]   = bestMatchString(metNamesListExcRxns,currentString);
      possibleMets{i,1} = currentString;
      possibleMets(i,3) = metListExcRxns(idx,1);
      possibleMets(i,2) = bestMatch;
      possibleMets{i,4} = idx;
    end

 %The following reactions did not return satisfying reactions, a more restrictive mapping is done here. 
    checkAgain   = [2;4;15;34;38];
    tempMetNames = model.metNames(~ismember(model.metNames,metNamesListExcRxns));
    tempMets     =  model.mets(~ismember(model.metNames,metNamesListExcRxns));
    for i = 1:size(checkAgain,1)
      currentString  = commonNames{checkAgain(i)};
      [bestMatch,idx] = bestMatchString(tempMetNames,currentString);
      possibleMets{checkAgain(i),1} = currentString;
      possibleMets(checkAgain(i),3) = tempMets(idx,1);
      possibleMets(checkAgain(i),2) = bestMatch;
      possibleMets{checkAgain(i),4} = idx;
    end
    
 %% Find reactions
    for i = 1:size(Data,1)
       if ~isempty(Data{i,2})
          continue
       end
       currentMet   = possibleMets(i,3);
       [rxnList, ~] = findRxnsFromMets(model, currentMet);
       rxnList      = rxnList(ismember(rxnList,excRxns));
       if ~isempty(rxnList)
            Data{i,3} = rxnList;
       end
      Data{i,2} = currentMet;
    end

% Mets with two reaction matches
    mets       = Data{16,2};
    Data{16,2} = Data{16,2}{1};
    Data{16,3} = Data{16,3}{1};
    
% Change formatting of Human1 Reactions
for i = 1:size(Data,1)
    if ~isempty(Data{i,3})
        currentReactions = Data{i,3}{1};
        Data{i,3} = currentReactions;
    end
end
    
% Using alternative name for linolenic acid
    [bestMatch,idx] = bestMatchString(model.metNames,'linolenate [Extracellular]');
    tempMet         =  model.mets(ismember(model.metNames,bestMatch));
    [rxnList, ~]    = findRxnsFromMets(model, tempMet);
    rxnList         = rxnList(ismember(rxnList,excRxns));
    Data(34,3)      = rxnList;
% Reactions with no exchange reactions will be created an exchange reaction in the model as
% they currently dont have one
  %Decanoyl Carnitine
     [bestMatch,idx] = bestMatchString(model.metNames,'Decanoyl-coa[Cytosol]');
     previousMetID   = model.mets(idx);
     ExchngeMetName      = strrep(bestMatch{1,1},'Cytosol','Extracellular');
     newMetID        = strrep(previousMetID{1,1},'c[c]','s[s]');
     model2          = addMetabolite(model,newMetID, ExchngeMetName);
     model2          = addReaction(model2,'MAR99991','reversible',1,'reactionFormula', [previousMetID{1} ' <=> ' newMetID],'lowerBound',-10,'upperBound',10);
     model2          = addExchangeRxn(model2, newMetID, -10, 10);
     [bestMatch,~]   = bestMatchString(model2.rxns,['EX_' newMetID]);
     Data{15,2}      = newMetID;
     Data{15,3}      = bestMatch{1,1};
  %Trimethylamine N-Oxide
     [bestMatch,idx] = bestMatchString(model.metNames,'Trimethylamine N-Oxide[Cytosol]');
     previousMetID   = model.mets(idx);
     ExchngeMetName      = strrep(bestMatch{1,1},'Cytosol','Extracellular');
     newMetID        = strrep(previousMetID{1,1},'c[c]','s[s]');
     model3          = addMetabolite(model2,newMetID, ExchngeMetName);
     model3          = addReaction(model3,'MAR99992','reversible',1,'reactionFormula', [previousMetID{1} ' <=> ' newMetID],'lowerBound',-10,'upperBound',10);
     model3          = addExchangeRxn(model3, newMetID, -10, 10);
      [bestMatch,~] = bestMatchString(model3.rxns,['EX_' newMetID]);
     Data{53,2}      = newMetID;
     Data{53,3}      = bestMatch{1,1};
  %Myriostiocarnitine
     [bestMatch,idx] = bestMatchString(model.metNames,'O-tetradecanoylcarnitine[Cytosol]');
     previousMetID   = model.mets(idx);
     ExchngeMetName      = strrep(bestMatch{1,1},'Cytosol','Extracellular');
     newMetID        = strrep(previousMetID{1,1},'c[c]','s[s]');
     model4          = addMetabolite(model3,newMetID, ExchngeMetName);
     model4          = addReaction(model4,'MAR99993','reversible',1,'reactionFormula', [previousMetID{1} ' <=> ' newMetID],'lowerBound',-10,'upperBound',10);
     model4          = addExchangeRxn(model4, newMetID, -10, 10);
     [bestMatch,~] = bestMatchString(model4.rxns,['EX_' newMetID]);
     Data{38,2}      = newMetID;
     Data{38,3}      = bestMatch{1,1};     
%% Final Table
sz = [size(Data,1) 10];
varTypes = ["cell","cell","cell","cell","double","double","double","double","double","double"];
FinalRelativeTableData = table('Size',sz,'VariableTypes',varTypes,'VariableNames',[{'MetNames'},{'KeggIDs'},{'MetsHuman1'},{'ReactionsHuman1'},{'Placebo'},{'Treated'},{'Death'},{'PlaceboSurvived'},{'PlaceboDeath'},{'TreatedDeath'}]);

FinalRelativeTableData.MetNames          = PlaceboDrugData.Properties.VariableNames(5:end)';
FinalRelativeTableData.KeggIDs           = Data(:,1);
FinalRelativeTableData.MetsHuman1        = Data(:,2);
FinalRelativeTableData.ReactionsHuman1   = Data(:,3);
FinalRelativeTableData.Placebo           = RatioPlacebo;
FinalRelativeTableData.Treated           = RatioTreated;
FinalRelativeTableData.Death             = RatioDeath;
FinalRelativeTableData.PlaceboSurvived   = RatioPlaceboSurvived;
FinalRelativeTableData.PlaceboDeath      = RatioPlaceboDeath;
FinalRelativeTableData.TreatedDeath      = RatioTreatedDeath;
modelNewExchangeReactions = model4;

end
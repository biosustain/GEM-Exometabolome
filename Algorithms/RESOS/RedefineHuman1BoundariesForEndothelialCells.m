%Redefines human1 Boundaries so that they represent the capacity of an
%endothelial cell
%clear
load('Human1.mat')
load('Recon1Endothelial.mat')
load('Recon3D.mat')
% Create tables of exchange rxns(mets and names) for Recon1 Recon 3D and
% Human 1
    excRxnsHuman1 = findExcRxns(model);
    excRxnsRecon1 = findExcRxns(Recon1Endothelial);
    excRxnsRecon3 = findExcRxns(Recon3D);
%Table for Human 1 exchange reactions  
    Human1(:,1) = model.rxns(excRxnsHuman1);
        [metListHuman, ~] = findMetsFromRxns(model, model.rxns(excRxnsHuman1));
        %Pass to cell array format
            for i = 1:size(metListHuman)
                newMetListHuman(i,1) =  metListHuman{i}(1);
            end
            for i = 1:size(newMetListHuman,1)
                Human1(i,2) =  model.metFormulas(ismember(model.mets,newMetListHuman{i}));
            end    

 %Table for RECON 1 exchange reactions
    Recon1IDs = UpdateRxnNamesRecon1(Recon1Endothelial);
    Recon1(:,1) = Recon1IDs(excRxnsRecon1);
    [metListRecon, ~] = findMetsFromRxns(Recon1Endothelial, Recon1Endothelial.rxns(excRxnsRecon1));
    %Pass to cell array format
        for i = 1:size(metListRecon)
            newMetListRecon(i,1) =  metListRecon{i}(1);
        end
        for i = 1:size(newMetListRecon,1)
           Recon1(i,2) =  Recon1Endothelial.metFormulas(ismember(Recon1Endothelial.mets,newMetListRecon{i}));
        end
   
%Table for RECON3D exchange reactions
    Recon3IDs         =  Recon3D.rxns;
    excRxnsRecon3     = findExcRxns(Recon3D);
    Recon3(:,1)       = Recon3IDs(excRxnsRecon3);
     Recon3(:,1)      = UpdateRxnNamesRecon3(Recon3(:,1));
    [metListRecon, ~] = findMetsFromRxns(Recon3D, Recon3D.rxns(excRxnsRecon3));
        %Pass to cell array format
            for i = 1:size(metListRecon)
                newMetListRecon(i,1) =  metListRecon{i}(1);
            end
            for i = 1:size(newMetListRecon,1)
                Recon3(i,2) =  Recon3D.metFormulas(ismember(Recon3D.mets,newMetListRecon{i}));
            end
  
   %NOTE: THERE ARE REPEATED EXCHANGE REACTIONS in Recon1
     
%FIRST mapping: Map Human1 metabolites to Recon1 metabolites. Add
%potential matches to Recon1 Table
for i = 1: size(Recon1,1)
        if  isempty(Recon1{i,2})
            continue
        end
        metMappedIdx = find(ismember(Human1(:,2),Recon1(i,2)));
        if isempty(metMappedIdx)
            continue
        end
        Recon1{i,3} = Human1(metMappedIdx,1);
end


%In Human 1 repository, there is a table that gives Recon3D reaction names
%equivalents in HUMAN 1. Recover that table.
        filename = 'rxnsHuman1.txt';
        fid      = fopen(filename);
        C        = textscan(fid,'%s','delimiter','\n');
        fclose(fid);
        SeparatedColumns = cell(size(C{1,1},1),1);
        for i=2:size(C{1},1)
            SeparatedColumns{i,1} = split(C{1,1}(i))';
            for j = [1,7]
                if j == 7
                  Human1Recon3DID{i-1,1} = SeparatedColumns{i,1}{j}(2:end-1); 
                else
                  Human1IDS{i-1,1} = SeparatedColumns{i,1}{j}(2:end-1);
                end
            end  % find empty spaces
        end
        
        tableSimilarsR3H1(:,1) = Human1Recon3DID;
        tableSimilarsR3H1(:,2) = Human1IDS;


%SECOND Mapping: Mapping Recon1 metabolites with Recon3 exchange reactions metabolites and
%recover the Human1 equivalent by comparing the Recon3 name to Human1 Name
    for i = 1: size(Recon1,1)
        numberOfHuman1Matches = size(Recon1{i,3},1);
        if numberOfHuman1Matches >=1
            continue
        end
        currentMet = Recon1(i,2);
        %MapMetabolite with recon 3
        mapMetRecon3 = find(ismember(Recon3(:,2),currentMet));
        for j = 1: length(mapMetRecon3)
            human1Equivalent = tableSimilarsR3H1(find(ismember(tableSimilarsR3H1(:,1),Recon3(mapMetRecon3(j,1)))),2);
            if ~isempty(human1Equivalent)
                idxInHuman1Table = find(ismember(Human1(:,1),human1Equivalent));
                equivalentExcRxn = Human1(idxInHuman1Table,1);
                if ~isempty(equivalentExcRxn)
                    Recon1{i,3} = equivalentExcRxn;
                end
            end
            
        end
    end

%THIRD mapping: When no possible match was returned in the first mapping, Map Recon1 Rxns
%to Recon3 rxn names, using exact name match (There are name differences between recon 3 and recon 1)
% to then recover human 1 equivalent.
    for i = 1: size(Recon1,1)
        numberOfHuman1Matches = size(Recon1{i,3},1);
        if numberOfHuman1Matches >=1
            continue
        end
        currentRxn    = Recon1(i,1);
        mapNameRecon3 = find(ismember(Recon3(:,1),currentRxn));
        for j = 1: length(mapNameRecon3)
            human1Equivalent = tableSimilarsR3H1(find(ismember(tableSimilarsR3H1(:,1),Recon3(mapNameRecon3(j,1)))),2);
            if ~isempty(human1Equivalent)
                idxInHuman1Table = find(ismember(Human1(:,1),human1Equivalent));
                equivalentExcRxn = Human1(idxInHuman1Table,1);
                if ~isempty(equivalentExcRxn)
                    Recon1{i,3} = equivalentExcRxn;
                end
            end
        end
    end

% FOURTH mapping: When the first mapping with human 1 returned more than one possible Human1 match, 
% map Recon1 names with Recon3 (using the table of equivalencies between Recon3 and Human1)
    % names to retrieve human1 reaction equivalent
        for k=1:size(Recon1,1)
            numberOfHuman1Matches     = size(Recon1{k,3},1);
            if numberOfHuman1Matches <=1
                continue
            end
            potentialMatches = Recon1{k,3};
            EqtoRecon3 = tableSimilarsR3H1(find(ismember(tableSimilarsR3H1(:,2),potentialMatches)),1);
            [~,matchRxnIdx] = bestMatchString(EqtoRecon3, Recon1{k,1});
            Recon1(k,3) = Recon1{k,3}(matchRxnIdx);
        end

%Eliminate apostrophes in Human1 reaction names
    tempTable   = tableSimilarsR3H1(ismember(tableSimilarsR3H1(:,1),Recon3(:,1)),:);
    for i =1:size(Recon1,1)
        value = Recon1{i,3};
        value = char(value);
        tempCellArray{i,1} = value;
    end

%Manual Mapping of reactions with names that change a lot between Recon1
%and Recon 3
    %DM_dctp_m -> DCTPtm
     Recon1(1,3) = tableSimilarsR3H1(9702,2);
    %DM_4hrpo -> '4HPRO_LTte'
     Recon1(210,3) = tableSimilarsR3H1(9778,2);
    %EX5MTP   -> EX_3mtp_e (MOST SIMILAR MATCH)
     Recon1(246,3) = tableSimilarsR3H1(9951,2);
     Recon1(299,3) = tableSimilarsR3H1(8591,2);
    %sink-ser/Thr[g] -> Asn_X_Ser_Thrtr'
    Recon1(301,3) = tableSimilarsR3H1(7360,2);
    equivalencyNotFound = [6,213,214,215,219,273,274,275,308];
%FINAL MAPPING: Mapping once again names of Recon1 with Recon3 without exact match.
   %As the name format is different between Recon1 and Recon 3, some
   %reaction names cannot be mapped exactly. Therefore, here we use a
   %function, bestMatchString,that compares strings based on a points system
                rxnsNotUsed = find(~ismember(tempTable(:,2),tempCellArray(:,1)));
                tempTable   =  tempTable(rxnsNotUsed,:);
            for k=1:size(Recon1,1)
                
                if sum(ismember(k,equivalencyNotFound))==1
                    continue
                end
                numberOfHuman1Matches     = size(Recon1{k,3},1);
                if numberOfHuman1Matches == 1 
                    continue
                end
                currentName = Recon1{k,1};
                [a,matchRxnIdx] = bestMatchString(tempTable(:,1), currentName);
                Recon1(k,3) = tempTable(matchRxnIdx,2);
            end

%Pass potential Boundaries
    excRxnsRecon1 = findExcRxns(Recon1Endothelial);
    Recon1(:,4) = num2cell(Recon1Endothelial.lb(excRxnsRecon1));
%Recon1(1,8) = {'Upper Boundaries'};
    Recon1(:,5) = num2cell(Recon1Endothelial.ub(excRxnsRecon1));

testmodel = model;
for i=1:size(Recon1(:,1),1)
    currentUB = cell2mat(Recon1(i,5));
    currentLB = cell2mat(Recon1(i,4));
    currentHuman1Rxn = Recon1{i,3};
    if isempty(currentHuman1Rxn)
        continue
    end
        modelIni = testmodel;
        testmodel = changeRxnBounds(testmodel,currentHuman1Rxn,currentUB,"u");
        testmodel = changeRxnBounds(testmodel,currentHuman1Rxn,currentLB,"l");
        if isnan(optimizeCbModel(testmodel).f)
            testmodel = modelIni;
            continue
        end
        modelIni = testmodel;
       
end
model = testmodel;
clearvars -except model;



  function Recon1IDs =  UpdateRxnNamesRecon1(Recon1Endothelial)
    Recon1IDs = Recon1Endothelial.rxns;
    for i = 1:size(Recon1IDs,1)

        if ~isempty(strfind(Recon1IDs{i,1}, '(e)'))
             Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'(e)','[e]');
        elseif ~isempty(strfind(Recon1IDs{i,1}, '24,25'))
             Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'24,25','24_25');
        elseif ~isempty(strfind(Recon1IDs{i,1}, '25HVITD2tin_m'))
            Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'25HVITD2tin_m','25HVITD2tin_m;25HVITD2tm');
        elseif ~isempty(strfind(Recon1IDs{i,1}, '(NADP)'))
            Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'(NADP)','_NADP_');        
        elseif ~isempty(strfind(Recon1IDs{i,1}, '(2)r'))
            Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'(2)r','_2_r');        
        elseif ~isempty(strfind(Recon1IDs{i,1}, '(211)r'))
            Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'(211)r','_211_r');          
        elseif ~isempty(strfind(Recon1IDs{i,1}, '(311)r'))
            Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'(211)r','_211_r');      
        elseif ~isempty(strfind(Recon1IDs{i,1}, 'u10m'))
            Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'u10m','u10mi');               
        elseif ~isempty(strfind(Recon1IDs{i,1}, 'EX_3h5chola'))
            Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'EX_3h5chola','EX_cholate[e]');                
        elseif ~isempty(strfind(Recon1IDs{i,1}, 'Ex_glc-D[e]'))
            Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'Ex_glc-D[e]','EX_glc_D[e]');             
        elseif ~isempty(strfind(Recon1IDs{i,1}, 'Ex_lac-L[e]'))
            Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'Ex_lac-L[e]','EX_lac_L[e]');              
        elseif ~isempty(strfind(Recon1IDs{i,1}, 'Ex_cys_L[e]'))
            Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'Ex_cys_L[e]','EX_cys_L[e]');              
        elseif ~isempty(strfind(Recon1IDs{i,1}, 'EX_adp'))
           Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'EX_adp','EX_adp[e]');   
        elseif ~isempty(strfind(Recon1IDs{i,1}, 'Excspg_a'))
             Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'Excspg_a','EX_cspg_a[e]');
        elseif ~isempty(strfind(Recon1IDs{i,1}, 'Excspg_b'))
             Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'Excspg_b','EX_cspg_b[e]');
        elseif ~isempty(strfind(Recon1IDs{i,1}, 'Excspg_c'))
             Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'Excspg_c','EX_cspg_c[e]');
        elseif ~isempty(strfind(Recon1IDs{i,1}, 'Excspg_d'))
             Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'Excspg_d','EX_cspg_d[e]');
        elseif ~isempty(strfind(Recon1IDs{i,1}, 'Excspg_e'))
             Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'Excspg_e','EX_cspg_e[e]');             
        elseif ~isempty(strfind(Recon1IDs{i,1}, 'Ex_'))
             Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'Ex_','EX_'); 
        elseif ~isempty(strfind(Recon1IDs{i,1}, 'tp(m)'))
             Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'(m)','_m');
        elseif ~isempty(strfind(Recon1IDs{i,1}, 'sig(er)'))
             Recon1IDs{i,1} = strrep(Recon1IDs{i,1},'sig(er)','sig_r');
        end         
    end
  end
  
  function Recon3IDs =  UpdateRxnNamesRecon3(Recon3IDs)
    %excRxns = findExcRxns(Recon3);
    %Recon3IDs = Recon3.rxns(excRxns);
    patDM = 'DM_' + alphanumericsPattern + '_' + alphanumericsPattern;
    for i=1:size(Recon3IDs)
        
        if ~isempty(strfind(Recon3IDs{i,1}, '_e'))
         Recon3IDs{i,1} = strrep(Recon3IDs{i,1},'_e','[e]');
        end
        if ~isempty(strfind(Recon3IDs{i,1}, '__'))
        Recon3IDs{i,1} = strrep(Recon3IDs{i,1},'__','_');
        end
        %if i == 374
        %    i
        %end
        %if contains(Recon3IDs{i,1}, patDM)
        %     str = Recon3IDs{i,1};
        %     pat = "_[abcdefghijklmnopqrstuvwxyz]_[abcdefghijklmnopqrstuvwxyz]";
        %     characters = extract(str,pat);
        %     Recon3IDs{i,1} = strrep(Recon3IDs{i,1},'DM_\w_\w','DM_[a-z]+([a-z]+)');
        %end
    end
  end
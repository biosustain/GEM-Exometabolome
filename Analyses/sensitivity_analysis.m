%% Initial parameters
sample_n = 1;
Iterations=10;
N = 8;  % Number of workers

% Get the 'Processes' cluster
c = parcluster('Processes');

% Change the NumWorkers property
c.NumWorkers = N;  % Or any number up to 512

% Save the changes
saveProfile(c);

%% Start a parallel pool with N workers
if isempty(gcp('nocreate'))
    parpool(N);
end

%% Loading
load('Ex2Task_Work_Space');

%% Routine
% Create a temporary structure to store the Summary data
SummaryTemp = cell(1, sample_n);

for p = sample_n:sample_n
    input_model = sprintf('modelPatient_Sampled_%d_Mean.mat', p);
    output_WS = sprintf('WS_%d', p);    
    pth_model = load(input_model);
    sampleMetaOutC = pth_model.sampleMetaOutC;
    % Table, Headers, and colums:
    Tasks=unique(TasksTable(:,1));
    headers=split(unique(join(TasksTable{:,2:4},'_'),'stable'),'_');
    DataCorr=cell(height(Tasks)+1,length(ithIndex)+2);
    for i=1:length(headers)
        DataCorr(i+1,1)=headers(i,1);      
        DataCorr(i+1,2)=headers(i,2);    
    end
    for i=1:length(ithIndex)
        DataCorr(1,i+2)=Model_Annotation.rxns(ithIndex(i)); 
    end
    DataP = DataCorr;
    %Warm up
    model=sampleMetaOutC;
    model.rxns=Model_Annotation.rxns;
    model.mets=Model_Annotation.mets;
    Summary={};
    %Calc Correlations
    spmd
        lockFilePath = sprintf('/home/igor/Documents/Trauma_Ind_sensitivity/cobratoolbox_2/.git/modules/papers/index.lock', labindex);
        if isfile(lockFilePath)
            system(['rm ', lockFilePath]);
        end
        addpath('/home/igor/Documents/Trauma_Ind_sensitivity/cobratoolbox_2/');
        initCobraToolbox_2
        lockFilePath = sprintf('/home/igor/Documents/Trauma_Ind_sensitivity/cobratoolbox_2/.git/modules/papers/index.lock', labindex);
        if isfile(lockFilePath)
            system(['rm ', lockFilePath]);
        end
        metaboliteIndices = labindex:numlabs:length(ithIndex);  % Indices of metabolites for this worker
        for j = 1:length(metaboliteIndices) %labindex:numlabs:length(ithIndex) % Loop for extracellular measured metabolites
            metaboliteIndex = ithIndex(metaboliteIndices(j));  % Actual index of the metabolite in ithIndex           
            JthIndexLb=min(model.points(metaboliteIndex,:));
            JthIndexUb=max(model.points(metaboliteIndex,:));
            JthIndexFractions=abs(JthIndexLb-JthIndexUb)/Iterations;
            JthCorrTable=zeros(Iterations,height(Tasks)+1);	
            for l=1:Iterations % loop for range of extrallular flux values
                modelb=model;
                %Change the jth nutrient uptake-secretion boundarie
                LthBoundary=JthIndexLb+JthIndexFractions*(l-1);
                modelb.lb(metaboliteIndex)=LthBoundary;
                modelb.ub(metaboliteIndex)=LthBoundary;
                JthCorrTable(l,1)=LthBoundary; 
                for k=1:height(Tasks) % Loop for Tasks
                    IthIndex=find(ismember(TasksTable(:,1),Tasks(k,1))); 
                    IthSubstrate=table2cell(TasksTable(IthIndex,5));
                    IthProduct=table2cell(TasksTable(IthIndex,8));
                    [Flux, ~,~] = testPathway(modelb,IthSubstrate,IthProduct);
                    disp(['j: ', num2str(j), ', l: ', num2str(l), ', k: ', num2str(k)]);                	
                    JthCorrTable(l,k+1)=Flux;
                end
            end                
            Summary{metaboliteIndices(j)} = JthCorrTable;  % Store the result at the position corresponding to the worker's %Summary{ithIndex(j)}=JthCorrTable;
        end    
    end
end

% Extract data from the Composite object
composite_copy = cell(1, length(Summary));
for h = 1:length(Summary)
    composite_copy{h} = Summary{h};
end

% Extract data from the Composite object
SummaryData = cell(1, length(ithIndex));  % Initialize SummaryData with the correct size
for i = 1:length(ithIndex)
    for j = 1:length(composite_copy)
        if i <= length(composite_copy{j})
            tmp_summary = composite_copy{j};
            if ~isempty(tmp_summary{i}) 
                SummaryData{i} = tmp_summary{i};  % Get the value of Summary on the ith metabolite
            end
        end
    end
end

% Save the extracted data
save(output_WS, 'SummaryData');

delete(gcp('nocreate'));

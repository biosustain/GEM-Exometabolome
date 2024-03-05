function [pointsCellArrayFurthest, pointsCellArrayCentral]= SamplingROOMSynLethals(centralSolution,furthestSolution,model,SLarray,samplingTime)  
%Samples the lethals returned by the analysis. It does 2 ROOMs, one for
%central and the other for furthest solution

try
    load('DataForNextSamplings')
catch
    structureForNextSamplings = [];
end
            numberOfReactions       = length(model.rxns);
            numberOfSL = size(SLarray,1);
            pointsCellArrayFurthest = cell(numberOfSL,1);
            pointsCellArrayCentral = cell(numberOfSL,1);
            %%Required for making parfor work.
            environment = getEnvironment();
            mr = mapreducer(0);
            parfor i = 1:numberOfSL
                %Restoring environment for making parfor work
                 restoreEnvironment(environment);
                 currentSL = SLarray(i,:);
                %ROOM and Recover MILP problem
                [fluxesMin, ~, ~,MILPproblemMin] = ROOMtolerance(model, centralSolution, currentSL,'delta', 0.1,'epsilon', 0.01);
                [fluxesMax, ~, ~,MILPproblemMax] = ROOMtolerance(model, furthestSolution, currentSL,'delta', 0.1,'epsilon', 0.01);
                points_samplingMin = length(MILPproblemMax.lb)*2;
                points_samplingMax = length(MILPproblemMin.lb)*2;
                mapreducer(mr)
                if and(sum(isnan(fluxesMin)) ~= 0,sum(isnan(fluxesMax)) ~= 0)
                    pointsCellArrayFurthest{i,1} = 0;
                    pointsCellArrayCentral{i,1} = 0;
                    continue
                end
                if and(sum(isnan(fluxesMin)) ~= 0,sum(isnan(fluxesMax)) == 0)
                     MILPproblemMax.S = MILPproblemMax.A;
                     pool_ROOM_Max = gpSamplerPercentile(MILPproblemMax,points_samplingMin,[],samplingTime,[],[],[],structureForNextSamplings);
                     pointsCellArrayFurthest{i,1} = tall(pool_ROOM_Max.points(1:numberOfReactions,:));
                     pointsCellArrayCentral{i,1} = 0;                    
                    continue
                elseif and(sum(isnan(fluxesMin)) == 0,sum(isnan(fluxesMax)) ~= 0)
                     MILPproblemMin.S = MILPproblemMin.A;
                     pool_ROOM_Min = gpSamplerPercentile(MILPproblemMin,points_samplingMin,[],samplingTime,[],[],[],structureForNextSamplings);
                     pointsCellArrayCentral{i,1} = tall(pool_ROOM_Min.points(1:numberOfReactions,:));
                     pointsCellArrayFurthest{i,1} = 0;
                    continue
                end
                MILPproblemMin.S = MILPproblemMin.A;
                MILPproblemMax.S = MILPproblemMax.A;
                %Pool of results 
                pool_ROOM_Min = gpSamplerPercentile(MILPproblemMin,points_samplingMin,[],samplingTime,[],[],[],structureForNextSamplings);
                pool_ROOM_Max = gpSamplerPercentile(MILPproblemMax,points_samplingMax,[],samplingTime,[],[],[],structureForNextSamplings);
                
                pointsCellArrayFurthest{i,1} = tall(pool_ROOM_Max.points(1:numberOfReactions,:));
                pointsCellArrayCentral{i,1} = tall(pool_ROOM_Min.points(1:numberOfReactions,:));

            end
end
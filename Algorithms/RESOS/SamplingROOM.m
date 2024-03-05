function pointsCellArray = SamplingROOM(WTRelevant_Array,model,KO,points_sampling,samplingTime)  
%Cambiar por el del server
            numberOfReactions       = length(model.rxns);
            numberOfRelevantResults = gather(size(WTRelevant_Array,2));
            WTRelevant_Array = gather(WTRelevant_Array);
            %ROOMfluxes = zeros(numberOfReactions,numberOfRelevantResults);
            pointsCellArray = cell(numberOfRelevantResults,1);
            mr = mapreducer(0);
            %%Required for making parfor work.
            environment = getEnvironment();
            parfor i = 1:numberOfRelevantResults
                %Restoring environment for making parfor work
                    restoreEnvironment(environment);                

                %ROOM and Recover MILP problem
                %Con delta y epsilon 0
                %delta = 0;
                %epsilon = 0;
                [fluxes,~,~,MILPproblem] = ROOM_2(model,WTRelevant_Array(:,i),KO,'delta', 0.15,'epsilon', 0.035);
                if sum(isnan(fluxes)) ~= 0
                    pointsCellArray{i,1} = zeros(numberOfReactions,numberOfRelevantResults,'single');
                    continue
                end
                MILPproblem.S = MILPproblem.A;
                %ROOMfluxes(:,i) = fluxes;
                %Pool of results 
                pool_ROOM_Results = gpSampler(MILPproblem,points_sampling,[],samplingTime);
                mapreducer(mr)
                pointsCellArray{i,1} = tall(pool_ROOM_Results.points(1:numberOfReactions,:));
            end
end
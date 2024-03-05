function [Summary_table,percentile25_vector,percentile75_vector] = SummaryTable(RelevantResults,model)
       %Statistic Summary of a sampling
           tempArray = gather(RelevantResults);
           IDs             = string(model.rxns);
           mean_vector     = zeros(length(IDs),1);
           mode_vector     = zeros(length(IDs),1);
           quantile_vector = zeros(length(IDs),1);
           sd_vector       = zeros(length(IDs),1);
           min_vector      = zeros(length(IDs),1);
           max_vector      = zeros(length(IDs),1);
           UPstd           = zeros(length(IDs),1);
           LOWstd          = zeros(length(IDs),1);
           percentile25_vector = zeros(length(IDs),1);
           percentile75_vector = zeros(length(IDs),1);
           quantileValue = 0.5;
           
           for i = 1:length(IDs)
              dataArray            = tempArray(i,:);
              mean_vector(i,1)     = mean(dataArray);
              quantile_vector(i,1) = quantile(dataArray,quantileValue);
              sd_vector(i,1)       = std(dataArray);
              min_vector(i,1)      = min(dataArray);
              UPstd(i,1)                = mean(dataArray) + (std(dataArray)*2);
              LOWstd(i,1)               = mean(dataArray) - (std(dataArray)*2);
              mode_vector(i,1)     = mode(dataArray);
              percentile25_vector(i,1) = prctile(dataArray,25);
              percentile75_vector(i,1) = prctile(dataArray,75);
              max_vector(i,1) = max(dataArray);
           end
           
           Summary_table = table(IDs,mean_vector,mode_vector,sd_vector,UPstd,LOWstd,min_vector,max_vector,quantile_vector,percentile75_vector,percentile25_vector);
end
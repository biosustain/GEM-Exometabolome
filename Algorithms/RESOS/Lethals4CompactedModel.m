function [JslCompactedFinal,JdlCompactedFinal,JtlCompactedFinal] = Lethals4CompactedModel(modelLoopless,model_compacted,Jsl,Jdl,Jtl)
%Using the lethals returned by fastSL when it was done with the original
%model, this function finds the corresponding targets in the compacted
%model.
%Eliminate possible Lethals that were in loops
JslLoopless      = Jsl(ismember(Jsl,modelLoopless.rxns),:);
tempVar(:,1) = ismember(Jdl(:,1),modelLoopless.rxns);
tempVar(:,2) = ismember(Jdl(:,2),modelLoopless.rxns);
JdlLoopless  = Jdl(sum(tempVar,2) == 2,:);
tempVar2(:,1) = ismember(Jtl(:,1),modelLoopless.rxns);
tempVar2(:,2) = ismember(Jtl(:,2),modelLoopless.rxns);
tempVar2(:,3) = ismember(Jtl(:,3),modelLoopless.rxns);
JtlLoopless   = Jtl(sum(tempVar2,2) == 3,:);
%Find Lethals in compacted model     
    %For loop to save index and name of the reactions in each set
    
            for i = 1:size(model_compacted.rxns,1)
                currentSet = split(model_compacted.rxns(i),'@');
                if sum(ismember(currentSet,JslLoopless(:,1))) > 0
                    JslCompacted(ismember(JslLoopless(:,1),currentSet),1) = model_compacted.rxns(i);
                end
                if sum(ismember(currentSet,JdlLoopless(:,1))) > 0
                    JdlCompacted(ismember(JdlLoopless(:,1),currentSet),1) = model_compacted.rxns(i);
                end
                if sum(ismember(currentSet,JdlLoopless(:,2))) > 0
                    JdlCompacted(ismember(JdlLoopless(:,2),currentSet),2) = model_compacted.rxns(i);
                end   
                if sum(ismember(currentSet,JtlLoopless(:,1))) > 0
                    JtlCompacted(ismember(JtlLoopless(:,1),currentSet),1) = model_compacted.rxns(i);
                end  
                if sum(ismember(currentSet,JtlLoopless(:,2))) > 0
                    JtlCompacted(ismember(JtlLoopless(:,2),currentSet),2) = model_compacted.rxns(i);
                end
                if sum(ismember(currentSet,JtlLoopless(:,3))) > 0
                    JtlCompacted(ismember(JtlLoopless(:,3),currentSet),3) = model_compacted.rxns(i);
                end                 
            end  
        JtlCompacted = JtlCompacted(~cellfun(@isempty,JtlCompacted(:,1)),:);
        JtlCompacted = JtlCompacted(~cellfun(@isempty,JtlCompacted(:,2)),:);
        JtlCompacted = JtlCompacted(~cellfun(@isempty,JtlCompacted(:,3)),:);
        JdlCompacted = JdlCompacted(~cellfun(@isempty,JdlCompacted(:,1)),:);
        JdlCompacted = JdlCompacted(~cellfun(@isempty,JdlCompacted(:,2)),:);
        JslCompacted = JslCompacted(~cellfun(@isempty,JslCompacted),:);
        
        %Eliminate possible Duplicates
        tempVar1 = cell2table(JslCompacted);
        tempVar2 = cell2table(JdlCompacted);
        tempVar3 = cell2table(JtlCompacted);
        
        JslCompactedFinal = table2cell(unique(tempVar1,'rows'));
        JdlCompactedFinal = table2cell(unique(tempVar2,'rows'));     
        JtlCompactedFinal = table2cell(unique(tempVar3,'rows')); 

        end
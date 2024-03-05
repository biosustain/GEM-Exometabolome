function [Sset,SetLb,SetUb,ReactionSet,MetaboliteSet] = findReactionSetsNotCompressingSeveralReactions(S,lb,ub,BioMassIndex,NotCompressedIndex)
% Find minimum reaction set
% Inputs:  stoichiometry, lower and upper bounds
% Outputs: minimal reaction and bounds sets
%%%%%%%%%%%%%%%%%%%%%% Lars Nielsen UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Modified by Fernando Silva-Lance 2021 %%%%%%%%%%%%%%%%%%%%%%%
            %Modification was done by adding conditions that do not
            %compress any desired reactions.
            
[NoMetabolites,NoReactions] = size(S);
ReactionSet   = zeros(NoReactions,4);
MetaboliteSet = zeros(NoMetabolites,1);
K = null(S);
for i=1:NoReactions
    if isempty(find(K(i,:)))
        ReactionSet(i,1)=-1;
        ReactionSet(i,2)= 0;
        ReactionSet(i,3)= 0;
        ReactionSet(i,4)= 0;            
    else
        [C,idx]=max(abs(K(i,:))); %C=max val, idx=col where the max is
        ReactionSet(i,2)=K(i,idx);
        K(i,:)=K(i,:)./K(i,idx); %Normalize the row by the max val => 1=max
    end
end
ReactionSetNo = 1;
for i=1:NoReactions
    if ~ReactionSet(i,1) && ~sum(ismember(BioMassIndex,i))
        ReactionSet(i,1) = ReactionSetNo;
        CurrentSetIdx = i;
        for j = i+1:NoReactions
             if sum(ismember(NotCompressedIndex,i)) 
                 continue
             end
            if ~sum(ismember(NotCompressedIndex,j)) && ~sum(ismember(BioMassIndex,j))
                if min(norm(K(i,:)-K(j,:)),norm(K(i,:)+K(j,:)))<1e-10 %Frobenius norm
                    ReactionSet(j,1) = ReactionSetNo;
                    CurrentSetIdx = [CurrentSetIdx,j];
                end
            end
        end
        CurrentSet = ReactionSet(CurrentSetIdx,:);
        [MinAbsWeight,idx] = min(abs(CurrentSet(:,2)));
        ReactionSet(CurrentSetIdx,2) = CurrentSet(:,2)/CurrentSet(idx,2);        
        CurrentLb = -1e9;
        CurrentUb = 1e9;
        for j = CurrentSetIdx
            if ReactionSet(j,2)>0
                CurrentLb = max(CurrentLb,lb(j)/ReactionSet(j,2));
                CurrentUb = min(CurrentUb,ub(j)/ReactionSet(j,2));
            else
                CurrentLb = max(CurrentLb,ub(j)/ReactionSet(j,2));
                CurrentUb = min(CurrentUb,lb(j)/ReactionSet(j,2));
            end
        end
        if abs(CurrentUb-CurrentLb)<1e-9
            %UNCOMMENT IF YOU WANT TO MAKE SURE EVERY EXC RXN IS KEPT EVEN
            %THOUGH THEIR UPPER AND LOWER BOUNDARY IS 0
             if sum(ismember(NotCompressedIndex,i)) 
                    ReactionSet(CurrentSetIdx,3) = CurrentLb;
                    ReactionSet(CurrentSetIdx,4) = CurrentUb; 
                    ReactionSetNo = ReactionSetNo + 1;
                 continue
             end 
            ReactionSet(CurrentSetIdx,1) = -2;
            ReactionSet(CurrentSetIdx,3) = CurrentLb;
            ReactionSet(CurrentSetIdx,4) = CurrentUb;            
        elseif abs(CurrentUb) <1e-9 % Change direction
            ReactionSet(CurrentSetIdx,2) = -ReactionSet(CurrentSetIdx,2);
            ReactionSet(CurrentSetIdx,3) = 0;
            ReactionSet(CurrentSetIdx,4) = -CurrentLb;            
            ReactionSetNo = ReactionSetNo + 1;
        else
            ReactionSet(CurrentSetIdx,3) = CurrentLb;
            ReactionSet(CurrentSetIdx,4) = CurrentUb;            
            ReactionSetNo = ReactionSetNo + 1;
        end
    %Entra aquí cuando es la biomasa
    elseif ~ReactionSet(i,1) && sum(ismember(BioMassIndex,i)) 
        ReactionSet(i,1) = ReactionSetNo;
        if CurrentUb <= ub(i)
            TempCurrentUb = max(CurrentUb,ub(i));
        else 
            TempCurrentUb = CurrentUb;
        end
        CurrentSetIdx = i;
        CurrentSet = ReactionSet(CurrentSetIdx,:);
        [MinAbsWeight,idx] = min(abs(CurrentSet(:,2)));
        ReactionSet(CurrentSetIdx,1) = ReactionSetNo;
        ReactionSet(CurrentSetIdx,2) = CurrentSet(:,2)/CurrentSet(idx,2);
        ReactionSet(CurrentSetIdx,3) = 0;
        ReactionSet(CurrentSetIdx,4) = TempCurrentUb;    
        ReactionSetNo = ReactionSetNo + 1;
    end
end
ReactionSetNo = max(ReactionSet(:,1));
Sset = zeros(NoMetabolites,ReactionSetNo);
SetLb = zeros(ReactionSetNo,1);
SetUb = zeros(ReactionSetNo,1);
for col = 1:ReactionSetNo
    CurrentSetIdx = find(ReactionSet(:,1)==col);
    for j= CurrentSetIdx
        Sset(:,col) = Sset(:,col) + S(:,j)*ReactionSet(j,2);
    end
    SetLb(col) = ReactionSet(CurrentSetIdx(1),3);
    SetUb(col) = ReactionSet(CurrentSetIdx(1),4);
end
Sset(find(abs(Sset)<1e-9))=0;
CurrentSetMetaboliteNo = 1;
for row = 1:NoMetabolites
    if isempty(find(Sset(row,:)))
        MetaboliteSet(row) = -1;
    else
        MetaboliteSet(row) = CurrentSetMetaboliteNo;
        CurrentSetMetaboliteNo = CurrentSetMetaboliteNo+1;
    end
end
Sset = Sset(find(MetaboliteSet>0),:);
end
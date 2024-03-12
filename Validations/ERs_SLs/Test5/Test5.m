%%TEST 5 -> LOOPLESS+COMPACTED+ROOM(All rxns)|JUNE 28 2022

load('EcoLooplessCompacted_28Jun.mat') %Workspace from tests 1 and 4 - here we use the compacted model to test all rxns with ROOM

%%Arrange model
model=EcoCompacted
JslCompactedFinal2=model.rxns(~ismember(model.rxns,IgnoreForLethality)) %only for singles
JdlCompactedFinal={}
JtlCompactedFinal={}


ArrayWithLethals = [{JslCompactedFinal2};{JdlCompactedFinal};{JtlCompactedFinal}];


%%Note on this ROOM: Uses a modified version to get only the ERs. 
%The mod version, was used because the calculation of the SLs takes too
%long
disp ('SLROOM') 
[JslROOM_All,BMReduction_All]= fastSLROOM2ERs(model,cutoff,ArrayWithLethals,maxSolution,minSolution);



function [NewNameData] = MatchRxnNamesPC3(filename)
%This function mathces the reaction names in PC3 Data with the format in
%the model 'Recon1 tunning'
fid =fopen(filename);
C=textscan(fid,'%s','delimiter','\n');
fclose(fid);
SeparatedColumns = cell(size(C{1,1},1),1);
for i=1:size(C{1},1)
    SeparatedColumns{i,1} = split(C{1,1}(i))'; % find empty spaces
end
TempName = cell(size(SeparatedColumns,1),1);
NewRxnName = cell(size(SeparatedColumns,1),1);
for i=1:size(C{1},1)
    TempName{i,1}   = split(SeparatedColumns{i,1}(1),'R_')'; % find empty spaces
    NewRxnName{i,1} = strrep(TempName{i,1}{2},'_LPAREN_e_RPAREN_','(e)');
    if ~isempty(strfind(TempName{i,1}{2}, '_LPAREN_m_RPAREN_'))
    NewRxnName{i,1} = strrep(TempName{i,1}{2},'_LPAREN_m_RPAREN_','(m)');
    end
    if ~isempty(strfind(TempName{i,1}{2}, '_LPAREN_2_RPAREN_r'))
    NewRxnName{i,1} = strrep(TempName{i,1}{2},'_LPAREN_2_RPAREN_r','(2)r');
    end
    if ~isempty(strfind(TempName{i,1}{2}, '_LPAREN_311_RPAREN_r'))
    NewRxnName{i,1} = strrep(TempName{i,1}{2},'_LPAREN_311_RPAREN_r','(311)r');
    end
    if ~isempty(strfind(TempName{i,1}{2}, '_LPAREN_211_RPAREN_r'))
    NewRxnName{i,1} = strrep(TempName{i,1}{2},'_LPAREN_211_RPAREN_r','(211)r');
    end
    if ~isempty(strfind(TempName{i,1}{2}, '_LPAREN_NADP_RPAREN_'))
    NewRxnName{i,1} = strrep(TempName{i,1}{2},'_LPAREN_NADP_RPAREN_','(NADP)');
    end
    if ~isempty(strfind(TempName{i,1}{2}, '_pnto_'))
    NewRxnName{i,1} = strrep(TempName{i,1}{2},'_pnto_','_pnto_R(e)');
    end  
end
NewNameData = cell(size(SeparatedColumns,1),3);

for i = 1:size(SeparatedColumns,1)
    NewNameData{i,1} = char(NewRxnName{i,1});
    NewNameData{i,2} = str2double(SeparatedColumns{i,1}{2});
    NewNameData{i,3} = str2double(SeparatedColumns{i,1}{3});
end
end

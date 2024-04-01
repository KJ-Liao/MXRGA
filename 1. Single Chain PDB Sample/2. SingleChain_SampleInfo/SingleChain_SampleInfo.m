% Load Single Chain protein ID/Sequence
SCH_Sample=table2struct(readtable('SCH_Seq.xlsx'));

SCH_Sample(1).Res=[];
SCH_Sample(1).SG=[];
SCH_Sample(1).Unit_Cell=[];
for n=1:size(SCH_Sample,1)  
    try
        File=fopen(strcat(lower(SCH_Sample(n).ID),'.pdb'));
        % Resolution
        while(1)
            line=fgetl(File);
                if isequal(sscanf(line(1:21),'%c'),'REMARK   2 RESOLUTION'),
                    SCH_Sample(n).Res=line(27:30);
                    break,
                end
        end

        % Space group
        while(1)
            line=fgetl(File);
            if isequal(sscanf(line(1:3),'%c'),'CRY'), 
                SCH_Sample(n).Unit_Cell=sscanf(line(8:55),'%f %f %f');
                SCH_Sample(n).SG=strcat(line(56:66),'');
                break, 
            end
        end
        
        fclose(File);
    end
end

writetable(struct2table(SCH_Sample), 'SCH_SampleInfo.xlsx')

%%%% Load C2.DiPER TopTank Result (25)
load('C2_TopRank_Result.mat')

%%%% DockQ Calculation
Sample(1).ID=[];
Sample(1).Output_Idx=[];

Sample(1).Fnat=[];
Sample(1).iRMS =[];
Sample(1).LRMS=[];
Sample(1).DockQ=[];

for r=1:size(C2_TopRank_Result,2)
    % Initillization
    Fnat_Result=[];
    iRMS_Result=[];
    LRMS_Result=[];
    DockQ_Result=[];

    parfor n=1:size(C2_TopRank_Result(r).Output_Idx,1)
        Xtal_Filename=['C2.Top25_Pose/', lower(C2_TopRank_Result(r).ID), '/', lower(C2_TopRank_Result(r).ID), '.002.', num2str(C2_TopRank_Result(r).Output_Idx(n)), '_Ans.pdb'];
        AF_Filename=['C2.Top25_Pose/', lower(C2_TopRank_Result(r).ID), '/', lower(C2_TopRank_Result(r).ID), '.002.', num2str(C2_TopRank_Result(r).Output_Idx(n)), '.pdb'];
         
        cmd_line=['./DockQ.py ', AF_Filename,' ', Xtal_Filename, ' -native_chain1 A -short'];
        [~, output] = system(cmd_line);

        Anchor_1=strfind(output, 'DockQ');
        Anchor_2=strfind(output, 'Fnat');
        Anchor_3=strfind(output, 'iRMS');
        Anchor_4=strfind(output, 'LRMS');
        Anchor_5=strfind(output, 'Fnonnat');

        if length(Anchor_1)==1 && length(Anchor_2)==1 && length(Anchor_3)==1 && length(Anchor_4)==1
            DockQ=str2double(output(Anchor_1+6:Anchor_2-2));
            Fnat=str2double(output(Anchor_2+5:Anchor_3-2));
            iRMS=str2double(output(Anchor_3+5:Anchor_4-2));
            LRMS=str2double(output(Anchor_4+5:Anchor_5-2));
        else
            Fnat=NaN;
            iRMS=NaN;
            LRMS=NaN;
            DockQ=NaN;
        end

        % Result
        Fnat_Result=[Fnat_Result; Fnat];
        iRMS_Result=[iRMS_Result; iRMS];
        LRMS_Result=[LRMS_Result; LRMS];
        DockQ_Result=[DockQ_Result; DockQ];
    end
    Sample(r).ID         = C2_TopRank_Result(r).ID;
    Sample(r).Output_Idx = C2_TopRank_Result(r).Output_Idx;
    Sample(r).Fnat       = Fnat_Result;
    Sample(r).iRMS       = iRMS_Result;
    Sample(r).LRMS       = LRMS_Result;
    Sample(r).DockQ      = DockQ_Result;

    r
end

save('Dimer_DockQ_Statistic.mat','Sample');
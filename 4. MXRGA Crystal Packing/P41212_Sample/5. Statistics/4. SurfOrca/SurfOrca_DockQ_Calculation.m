%%%% Load C4.DiPER TopTank Result (Top 300, Approximated RMSD_95<12A)
load('Sum_Table.mat')

%%%% DockQ Calculation
Sample(1).ID=[];
Sample(1).Output_Idx=[];

Sample(1).Fnat=[];
Sample(1).iRMS =[];
Sample(1).LRMS=[];
Sample(1).DockQ=[];

for r=1:size(Sum_Table,2)

    % FileName
    Xtal_Filename=['TopRank_300_Pose/', lower(Sum_Table(r).ID(1:4)), '/', Sum_Table(r).ID, '_Ans.pdb'];
    AF_Filename=['Refined_TopRank_300_Pose/', lower(Sum_Table(r).ID(1:4)), '/', Sum_Table(r).ID, '_Refined.pdb'];

    % AlphaFold2
    try
        AF=pdbread(AF_Filename);
        AF_status=1;
    catch
        AF_status=0;
    end

    if AF_status==1
        % Extract Coordinates of Alpha Carbon (CA)
        AF_Model=[AF.Model(1).Atom.X; AF.Model(1).Atom.Y; AF.Model(1).Atom.Z]';

        AF_resSeq_Idx=[AF.Model(1).Atom.resSeq]';
        pseudo_AF_CA_Idx=strcmp({AF.Model(1).Atom.AtomName}, 'CA')';
        AF_resSeq_Unique_Idx=unique(AF_resSeq_Idx(pseudo_AF_CA_Idx),'stable');

        AF_CA_Idx=zeros(size(AF_resSeq_Unique_Idx,1),1);
        for RSUidx=1:size(AF_resSeq_Unique_Idx,1)
            AF_CA_Idx(RSUidx,1)=find((AF_resSeq_Idx==AF_resSeq_Unique_Idx(RSUidx))&pseudo_AF_CA_Idx, 1, 'first');
        end

        Copy_Num=length(unique({AF.Model(1).Atom.chainID}));
        AF2_CA_Idx=zeros(length(AF_resSeq_Idx)/Copy_Num,1);
        AF2_CA_Idx(AF_CA_Idx)=1;
        AF_Unit=AF_Model(logical(repmat(AF2_CA_Idx,Copy_Num,1)),:);

        % [toXYZ; toXYZ_pair; fromXYZ; fromXYZ_pair]
        len_AFU=length(AF_Unit)/4;
        AF_Unit_to=AF_Unit(len_AFU*0+1:len_AFU*1,:);
        AF_Unit_tp=AF_Unit(len_AFU*1+1:len_AFU*2,:);
        AF_Unit_fr=AF_Unit(len_AFU*2+1:len_AFU*3,:);
        AF_Unit_fp=AF_Unit(len_AFU*3+1:len_AFU*4,:);

        DISM=squareform(pdist([mean(AF_Unit_to); mean(AF_Unit_tp); mean(AF_Unit_fr); mean(AF_Unit_fp)]));
        [~,Ind]=min(sum(DISM));

        if Ind==1
            ChainID='A';
        elseif Ind==2
            ChainID='B';
        elseif Ind==3
            ChainID='C';
        else
            ChainID='D';
        end

        % Command Line
        cmd_line=['./DockQ.py ', AF_Filename,' ', Xtal_Filename, ' -native_chain1 ', ChainID, ' -short'];
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
        Sample(r).Fnat       = Fnat;
        Sample(r).iRMS       = iRMS;
        Sample(r).LRMS       = LRMS;
        Sample(r).DockQ      = DockQ;
    else
        % Result
        Sample(r).Fnat       = NaN;
        Sample(r).iRMS       = NaN;
        Sample(r).LRMS       = NaN;
        Sample(r).DockQ      = NaN;
    end

    % Result
    Sample(r).ID         = Sum_Table(r).ID(1:4);
    Sample(r).Output_Idx = Sum_Table(r).ID;

    r
end

save('SurfOrca_DockQ_Statistic.mat','Sample');
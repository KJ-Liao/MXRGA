%%%% Load C4.DiPER TopTank Result (Top 300, Approximated RMSD_95<12A)
load('Sum_Table.mat')

%%%% RMSD95 Calculation
% Path Setting
Path1='/home/juju/Desktop/MXRGA_Sample Log/4. MXRGA Crystal Packing/P41212_Sample/6. SurfOrca_P41212/TopRank_300_Pose/';
Path2='/home/juju/Desktop/MXRGA_Sample Log/4. MXRGA Crystal Packing/P41212_Sample/6. SurfOrca_P41212/Refined_TopRank_300_Pose/';

% RMSD95 Calculation
Sample(1).ID=[];
Sample(1).Output_Idx=[];

Sample(1).Nat_Mod_OVLP=[];
Sample(1).Nat_Int_Res=[];
Sample(1).Fnat=[];
Sample(1).iRMS =[];
Sample(1).LRMS=[];
Sample(1).RRMS=[];
Sample(1).DockQ=[];

Sample(1).C4_RMSD95=[];
Sample(1).eRMSD=[];
Sample(1).RMSD_95=[];

for r=1:size(Sum_Table,2)

    % Xtal
    Xtal_Filename=[Path1, lower(Sum_Table(r).ID(1:4)), '/', Sum_Table(r).ID, '_Ans.pdb'];
    PDB=pdbread(Xtal_Filename);

    % Extract Coordinates of Alpha Carbon (CA)
    PDB_Model=[PDB.Model(1).Atom.X; PDB.Model(1).Atom.Y; PDB.Model(1).Atom.Z]';

    PDB_resSeq_Idx=[PDB.Model(1).Atom.resSeq]';
    pseudo_PDB_CA_Idx=strcmp({PDB.Model(1).Atom.AtomName}, 'CA')';
    PDB_resSeq_Unique_Idx=unique(PDB_resSeq_Idx(pseudo_PDB_CA_Idx),'stable');

    PDB_CA_Idx=zeros(size(PDB_resSeq_Unique_Idx,1),1);
    for RSUidx=1:size(PDB_resSeq_Unique_Idx,1)
        PDB_CA_Idx(RSUidx,1)=find((PDB_resSeq_Idx==PDB_resSeq_Unique_Idx(RSUidx))&pseudo_PDB_CA_Idx, 1, 'first');
    end

    % Sequence
    PDB_Seq=aminolookup(strrep(strrep([PDB.Model(1).Atom(PDB_CA_Idx).resName], 'SEC', 'CYS'), 'UNK', 'GLY'));

    %

    % AlphaFold2
    AF_Filename=[Path2, lower(Sum_Table(r).ID(1:4)), '/', Sum_Table(r).ID, '_Refined.pdb'];
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

        % Sequence
        AF_Seq=aminolookup(strrep(strrep([AF.Model(1).Atom(AF_CA_Idx).resName], 'SEC', 'CYS'), 'UNK', 'GLY'));

        %

        % Global Sequence Alignment to Extract Superimposed CA for RMSD Calculation
        Copy_Num=length(unique({PDB.Model(1).Atom.chainID}));
        [~, Alignment] = nwalign(repmat(PDB_Seq,1,Copy_Num), repmat(AF_Seq,1,Copy_Num));
        alignidex=find((Alignment(1,:)~='-')&(Alignment(3,:)~='-'));
        aseq1=find(Alignment(1,:)=='-');
        aseq2=find(Alignment(3,:)=='-');

        aseq1_idx=[]; aseq2_idx=[];
        for a=1:length(alignidex)
            aseq1_idx=[aseq1_idx, alignidex(a)-sum(aseq1<alignidex(a))];
            aseq2_idx=[aseq2_idx, alignidex(a)-sum(aseq2<alignidex(a))];
        end

        % RMSD 95 Calculation
        PDB2_CA_Idx=zeros(length(PDB_resSeq_Idx)/Copy_Num,1);
        PDB2_CA_Idx(PDB_CA_Idx)=1;
        Asy_Unit=PDB_Model(logical(repmat(PDB2_CA_Idx,Copy_Num,1)),:);

        AF2_CA_Idx=zeros(length(AF_resSeq_Idx)/Copy_Num,1);
        AF2_CA_Idx(AF_CA_Idx)=1;
        AF_Unit=AF_Model(logical(repmat(AF2_CA_Idx,Copy_Num,1)),:);

        [~,~,eRMSD,~,Sum_Square]=CoordiExam(Asy_Unit(aseq1_idx,:), AF_Unit(aseq2_idx,:));
        RMSD_95=sqrt(sum(Sum_Square(1:ceil(length(alignidex)*0.95)))/ceil(length(alignidex)*0.95));

        %

        %%%% Docking Quality Assessment
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

        % Docking Quality Assessment
        [Nat_Mod_OVLP, Nat_Int_Res, Fnat, iRMS, LRMS, RRMS, DockQ] = QualiAssess (Xtal_Filename,AF_Filename,ChainID);

        % Result
        Sample(r).eRMSD      = eRMSD;
        Sample(r).RMSD_95    = RMSD_95;
        Sample(r).Nat_Mod_OVLP = Nat_Mod_OVLP;
        Sample(r).Nat_Int_Res  = Nat_Int_Res;
        Sample(r).Fnat       = Fnat;
        Sample(r).iRMS       = iRMS;
        Sample(r).LRMS       = LRMS;
        Sample(r).RRMS       = RRMS;
        Sample(r).DockQ      = DockQ;
    else
        % Result
        Sample(r).eRMSD      = NaN;
        Sample(r).RMSD_95    = NaN;
        Sample(r).Nat_Mod_OVLP = NaN;
        Sample(r).Nat_Int_Res  = NaN;
        Sample(r).Fnat       = NaN;
        Sample(r).iRMS       = NaN;
        Sample(r).LRMS       = NaN;
        Sample(r).RRMS       = NaN;
        Sample(r).DockQ      = NaN;
    end

    % Result
    Sample(r).ID         = Sum_Table(r).ID(1:4);
    Sample(r).Output_Idx = Sum_Table(r).ID;
    Sample(r).C4_RMSD95  = Sum_Table(r).C4_RMSD95;

    r
end

save('SurfOrca_RMSD_Statistic.mat','Sample');
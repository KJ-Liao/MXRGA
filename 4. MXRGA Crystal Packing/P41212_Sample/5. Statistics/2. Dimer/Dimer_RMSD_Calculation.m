%%%% Load C2.DiPER TopTank Result (25)
load('C2_TopRank_Result.mat')

%%%% RMSD95 Calculation
% Path Setting
Path='/home/juju/Desktop/MXRGA_Sample Log/4. MXRGA Crystal Packing/P41212_Sample/4. C2.DIPER_P41212/4. C2_TopRank/';

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

Sample(1).eRMSD=[];
Sample(1).RMSD_95=[];

for r=1:size(C2_TopRank_Result,2)
    % Initillization
    Nat_Mod_OVLP_Result=[];
    Nat_Int_Res_Result=[];
    Fnat_Result=[];
    iRMS_Result=[];
    LRMS_Result=[];
    RRMS_Result=[];
    DockQ_Result=[];
    eRMSD_Result=[];
    RMSD_95_Result=[];

    parfor n=1:size(C2_TopRank_Result(r).Output_Idx,1)
        % Xtal
        Xtal_Filename=[Path, 'C2.Top25_Pose/', lower(C2_TopRank_Result(r).ID), '/', lower(C2_TopRank_Result(r).ID), '.002.', num2str(C2_TopRank_Result(r).Output_Idx(n)), '_Ans.pdb'];
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
        AF_Filename=[Path, 'C2.Top25_Pose/', lower(C2_TopRank_Result(r).ID), '/', lower(C2_TopRank_Result(r).ID), '.002.', num2str(C2_TopRank_Result(r).Output_Idx(n)), '.pdb'];
        AF=pdbread(AF_Filename);

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
        [~, Alignment] = nwalign([PDB_Seq,PDB_Seq], [AF_Seq,AF_Seq]);
        alignidex=find((Alignment(1,:)~='-')&(Alignment(3,:)~='-'));
        aseq1=find(Alignment(1,:)=='-');
        aseq2=find(Alignment(3,:)=='-');

        aseq1_idx=[]; aseq2_idx=[];
        for a=1:length(alignidex)
            aseq1_idx=[aseq1_idx, alignidex(a)-sum(aseq1<alignidex(a))];
            aseq2_idx=[aseq2_idx, alignidex(a)-sum(aseq2<alignidex(a))];
        end

        % RMSD 95 Calculation
        PDB2_CA_Idx=zeros(length(PDB_resSeq_Idx)/2,1);
        PDB2_CA_Idx(PDB_CA_Idx)=1;
        Asy_Unit=PDB_Model(logical([PDB2_CA_Idx; PDB2_CA_Idx]),:);

        AF2_CA_Idx=zeros(length(AF_resSeq_Idx)/2,1);
        AF2_CA_Idx(AF_CA_Idx)=1;
        AF_Unit=AF_Model(logical([AF2_CA_Idx; AF2_CA_Idx]),:);

        [~,~,eRMSD,~,Sum_Square]=CoordiExam(Asy_Unit(aseq1_idx,:), AF_Unit(aseq2_idx,:));
        RMSD_95=sqrt(sum(Sum_Square(1:ceil(length(alignidex)*0.95)))/ceil(length(alignidex)*0.95));

        %

        %%%% Docking Quality Assessment
        [Nat_Mod_OVLP, Nat_Int_Res, Fnat, iRMS, LRMS, RRMS, DockQ] = QualiAssess (Xtal_Filename,AF_Filename,'A');

        % Result
        Nat_Mod_OVLP_Result=[Nat_Mod_OVLP_Result; Nat_Mod_OVLP];
        Nat_Int_Res_Result=[Nat_Int_Res_Result; Nat_Int_Res];
        Fnat_Result=[Fnat_Result; Fnat];
        iRMS_Result=[iRMS_Result; iRMS];
        LRMS_Result=[LRMS_Result; LRMS];
        RRMS_Result=[RRMS_Result; RRMS];
        DockQ_Result=[DockQ_Result; DockQ];
        eRMSD_Result=[eRMSD_Result; eRMSD];
        RMSD_95_Result=[RMSD_95_Result; RMSD_95];

    end
    Sample(r).ID         = C2_TopRank_Result(r).ID;
    Sample(r).Output_Idx = C2_TopRank_Result(r).Output_Idx;
    Sample(r).Nat_Mod_OVLP = Nat_Mod_OVLP_Result;
    Sample(r).Nat_Int_Res  = Nat_Int_Res_Result;
    Sample(r).Fnat       = Fnat_Result;
    Sample(r).iRMS       = iRMS_Result;
    Sample(r).LRMS       = LRMS_Result;
    Sample(r).RRMS       = RRMS_Result;
    Sample(r).DockQ      = DockQ_Result;
    Sample(r).eRMSD      = eRMSD_Result;
    Sample(r).RMSD_95    = RMSD_95_Result;

    r
end

save('Dimer_RMSD_Statistic.mat','Sample');
%%%% Load MXRGA Result (Top 300, Approximated RMSD_95<12A)
load('MXRGA_Result.mat')

%%%% RMSD95 Calculation
% Path Setting
Path='/home/juju/Desktop/MXRGA_Sample Log/4. MXRGA Crystal Packing/P41212_Sample/6. SurfOrca_P41212/Packing_Poses/';

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

for r=1:size(MXRGA,2)

    % Read Files
    Xtal_Filename=[Path, MXRGA(r).PDB_ID, '/', MXRGA(r).SiXt_ID, '_Xtal_Packing.pdb'];
    AF_Filename=[Path, MXRGA(r).PDB_ID, '/', MXRGA(r).SiXt_ID, '_Sim_Packing.pdb'];

    try
        PDB=pdbread(Xtal_Filename);
        AF=pdbread(AF_Filename);
        File_status=1;
    catch
        File_status=0;
    end

    if File_status==1
        % Xtal
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
        % Docking Quality Assessment
        [Nat_Mod_OVLP, Nat_Int_Res, Fnat, iRMS, LRMS, RRMS, DockQ] = QualiAssess (Xtal_Filename,AF_Filename,'A');

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
    Sample(r).ID         = MXRGA(r).PDB_ID;
    Sample(r).Output_Idx = MXRGA(r).SiXt_ID;

    r
end

save('MXRGA_RMSD_Statistic.mat','Sample');
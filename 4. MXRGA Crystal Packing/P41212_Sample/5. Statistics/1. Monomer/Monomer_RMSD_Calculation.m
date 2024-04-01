%%%% Load PDB ID
List=fopen('Input_P41212_Sample_190.txt');
% 5A2C is excluded
Sample(1).ID=[];s=1;
while (1)
    line=fgetl(List);
    if line==-1, break, end
    Sample(s).ID=sscanf(line, '%c');
    s=s+1;
end
fclose(List);

%%%% RMSD95 Calculation
% Path Setting
Path='/home/juju/Desktop/MXRGA_Sample Log/4. MXRGA Crystal Packing/P41212_Sample/4. C2.DIPER_P41212/2. P41212_reduced AF_sample/';

% RMSD95 Calculation
Sample(1).eRMSD=[];
Sample(1).RMSD_95=[];

for r=1:size(Sample,2)
    % Xtal
    Xtal_Filename=[Path, 'PDB/', lower(Sample(r).ID), '.pdb'];
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
    Asy_Unit=PDB_Model(PDB_CA_Idx,:);

    % Sequence
    PDB_Seq=aminolookup(strrep(strrep([PDB.Model(1).Atom(PDB_CA_Idx).resName], 'SEC', 'CYS'), 'UNK', 'GLY'));

    %

    % AlphaFold2
    AF_Filename=[Path, lower(Sample(r).ID), '_trun_AF.pdb'];
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
    AF_Unit=AF_Model(AF_CA_Idx,:);

    % Sequence
    AF_Seq=aminolookup(strrep(strrep([AF.Model(1).Atom(AF_CA_Idx).resName], 'SEC', 'CYS'), 'UNK', 'GLY'));
    
    %

    % Global Sequence Alignment to Extract Superimposed CA for RMSD Calculation
    [~, Alignment] = nwalign(PDB_Seq,AF_Seq);
    alignidex=find((Alignment(1,:)~='-')&(Alignment(3,:)~='-'));
    aseq1=find(Alignment(1,:)=='-');
    aseq2=find(Alignment(3,:)=='-');

    aseq1_idx=[]; aseq2_idx=[];
    for a=1:length(alignidex)
        aseq1_idx=[aseq1_idx, alignidex(a)-sum(aseq1<alignidex(a))];
        aseq2_idx=[aseq2_idx, alignidex(a)-sum(aseq2<alignidex(a))];
    end

    % RMSD 95 Calculation
    [~,~,eRMSD,~,Sum_Square]=CoordiExam(Asy_Unit(aseq1_idx,:), AF_Unit(aseq2_idx,:));
    RMSD_95=sqrt(sum(Sum_Square(1:ceil(length(alignidex)*0.95)))/ceil(length(alignidex)*0.95));

    %

    % Result
    Sample(r).eRMSD=eRMSD;
    Sample(r).RMSD_95=RMSD_95;

    r
end

save('Monomer_RMSD_Statistic.mat','Sample');

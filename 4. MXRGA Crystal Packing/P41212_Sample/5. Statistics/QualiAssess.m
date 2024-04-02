function [Nat_Mod_OVLP, Nat_Int_Res, Fnat, iRMS, LRMS, RRMS, DockQ] = QualiAssess (Native_Filename,Model_Filename,ChainID)

    %%%% Load Corresponding PDB File
    % Parse Info of Naive PDB (Xtal)
    Native_PDB=PDBParse(Native_Filename);
    Native_Coord=[Native_PDB.X, Native_PDB.Y, Native_PDB.Z];

    Native_ChainA=strcmp(Native_PDB.chainID,ChainID);
    Native_Receptor_Coord=Native_Coord(Native_ChainA,:);
    Native_Receptor_resNum=Native_PDB.resNum(Native_ChainA);
    Native_Receptor_atomName=Native_PDB.atomName(Native_ChainA);

    Refined_Nat_resNum=min(Native_Receptor_resNum)-1;
    Native_Receptor_resNum=Native_Receptor_resNum-Refined_Nat_resNum;
    
    Native_CA_Idx=strcmp(Native_Receptor_atomName,'CA');
    Native_Receptor_resName=Native_PDB.resName(Native_CA_Idx);
    Native_Receptor_Seq=aminolookup(strrep(strrep(Native_Receptor_resName,'SEC','CYS'),'UNK','GLY'));

    % Parse Info of Model PDB (Simulated)
    Model_PDB=PDBParse(Model_Filename);
    Model_Coord=[Model_PDB.X, Model_PDB.Y, Model_PDB.Z];

    Model_ChainA=strcmp(Model_PDB.chainID,ChainID);
    Model_Receptor_Coord=Model_Coord(Model_ChainA,:);
    Model_Receptor_resNum=Model_PDB.resNum(Model_ChainA);
    Model_Receptor_atomName=Model_PDB.atomName(Model_ChainA);

    Refined_Mod_resNum=min(Model_Receptor_resNum)-1;
    Model_Receptor_resNum=Model_Receptor_resNum-Refined_Mod_resNum;

    Model_CA_Idx=strcmp(Model_Receptor_atomName,'CA');
    Model_Receptor_resName=Model_PDB.resName(Model_CA_Idx);
    Model_Receptor_Seq=aminolookup(strrep(strrep(Model_Receptor_resName,'SEC','CYS'),'UNK','GLY'));

    %

    %%%% Paired Ligand Coordinate (Other than Chain A)
    % Native_Ligand
    len_native_mono=size(Native_Receptor_resNum,1);
    Native_Copy=length(Native_ChainA)/len_native_mono-1;
    if ChainID=='A'
        Native_Ligand_Idx=[false(len_native_mono,1); repmat(true(len_native_mono,1),Native_Copy,1)];
    elseif ChainID=='B'
        Native_Ligand_Idx=[true(len_native_mono,1); false(len_native_mono,1); repmat(true(len_native_mono,1),Native_Copy-1,1)];
    elseif ChainID=='C'
        Native_Ligand_Idx=[repmat(true(len_native_mono,1),Native_Copy-1,1); false(len_native_mono,1); true(len_native_mono,1)];
    else
        Native_Ligand_Idx=[repmat(true(len_native_mono,1),Native_Copy,1); false(len_native_mono,1)];
    end

    Native_Ligand_Coord=Native_Coord(Native_Ligand_Idx,:);
    Native_Ligand_atomName=Native_PDB.atomName(Native_Ligand_Idx);

    Native_Ligand_resNum=[];
    for i=1:Native_Copy
        Native_Ligand_resNum=[Native_Ligand_resNum; Native_Receptor_resNum+(i-1)*max(Native_Receptor_resNum)]; 
    end

    % Model_Ligand
    len_model_mono=size(Model_Receptor_resNum,1);
    Model_Copy=length(Model_ChainA)/len_model_mono-1;
    if ChainID=='A'
        Model_Ligand_Idx=[false(len_model_mono,1); repmat(true(len_model_mono,1),Model_Copy,1);];
    elseif ChainID=='B'
        Model_Ligand_Idx=[true(len_model_mono,1); false(len_model_mono,1); repmat(true(len_model_mono,1),Model_Copy-1,1)];
    elseif ChainID=='C'
        Model_Ligand_Idx=[repmat(true(len_model_mono,1),Model_Copy-1,1); false(len_model_mono,1); true(len_model_mono,1)];
    else
        Model_Ligand_Idx=[repmat(true(len_model_mono,1),Model_Copy,1); false(len_model_mono,1)];
    end

    Model_Ligand_Coord=Model_Coord(Model_Ligand_Idx,:);
    Model_Ligand_atomName=Model_PDB.atomName(Model_Ligand_Idx);

    Model_Ligand_resNum=[];
    for j=1:Model_Copy
        Model_Ligand_resNum=[Model_Ligand_resNum; Model_Receptor_resNum+(j-1)*max(Model_Receptor_resNum)]; 
    end

    %

    %%%% Sequence Alignment to Identify Paired Residue No
    % Global Sequence Alignment to Extract Paired Residue Numbers
    [~, Alignment] = nwalign(Native_Receptor_Seq{1},Model_Receptor_Seq{1});
    aseq1=find(Alignment(1,:)=='-');
    aseq2=find(Alignment(3,:)=='-');

    alignidex=find((Alignment(1,:)~='-'));
    Realignment=Alignment(:,alignidex);

    len_alignidex=length(alignidex);
    aseq1_idx=NaN(1,len_alignidex);
    aseq2_idx=NaN(1,len_alignidex);
    for o=1:len_alignidex
        % Native Seq
        if Realignment(1,o)~='-'
            aseq1_idx(1,o)=alignidex(o)-sum(aseq1<alignidex(o));
        end
        % Model Seq
        if Realignment(3,o)~='-'
            aseq2_idx(1,o)=alignidex(o)-sum(aseq2<alignidex(o));
        end
    end

    %
    %
    %

    %%%% Model Quality Assessment
    %%% Fnat
    % Distance Map of Rcpt and Lgnd (Native)
    Native_Rcpt_CA_Idx = strcmp(Native_Receptor_atomName,'CA');
    Native_Rcpt_CA=Native_Receptor_Coord(Native_Rcpt_CA_Idx,:);
    Native_Lgnd_CA_Idx = strcmp(Native_Ligand_atomName,'CA');
    Native_Lgnd_CA=Native_Ligand_Coord(Native_Lgnd_CA_Idx,:);

    DISM_Nat=pdist2(Native_Rcpt_CA,Native_Lgnd_CA);

    % Distance Map of Rcpt and Lgnd (Model)
    Model_Rcpt_CA_Idx = strcmp(Model_Receptor_atomName,'CA');
    Model_Rcpt_CA=Model_Receptor_Coord(Model_Rcpt_CA_Idx,:);
    Model_Lgnd_CA_Idx = strcmp(Model_Ligand_atomName,'CA');
    Model_Lgnd_CA=Model_Ligand_Coord(Model_Lgnd_CA_Idx,:);

    DISM_Mod=pdist2(Model_Rcpt_CA,Model_Lgnd_CA);

    % Build up the Distance Map
    % Native_Rcpt_list=unique(Native_Receptor_resNum,'stable');
    Native_Rcpt_list=unique(Native_Receptor_resNum(Native_Rcpt_CA_Idx,:),'stable');
    % Model_Rcpt_list=unique(Model_Receptor_resNum,'stable');
    Model_Rcpt_list=unique(Model_Receptor_resNum(Model_Rcpt_CA_Idx,:),'stable');

    Native_Rcpt_No=length(Native_Rcpt_CA);
    Model_Rcpt_No=length(Model_Rcpt_CA);
    len_native_rcpt=length(Native_Receptor_resNum);
    len_model_rcpt=length(Model_Receptor_resNum);

    for c=1:Native_Copy
        % Atom Paired-wise Map
        DISM_Nat_Sub_Atom=pdist2(Native_Receptor_Coord,Native_Ligand_Coord(len_native_rcpt*(c-1)+1:len_native_rcpt*c,:));
        DISM_Mod_Sub_Atom=pdist2(Model_Receptor_Coord,Model_Ligand_Coord(len_model_rcpt*(c-1)+1:len_model_rcpt*c,:));

        % Native: Identify Atom Interactions (<5) within CA < 25A
        DISM_Nat_Sub=DISM_Nat(:,Native_Rcpt_No*(c-1)+1:Native_Rcpt_No*c);
        [Nat_row,Nat_col]=find(DISM_Nat_Sub<25);
        for n=1:length(Nat_row)
            nr=Nat_row(n);
            nc=Nat_col(n);
            DISM_Nat_Sub(nr,nc)=min(DISM_Nat_Sub_Atom(Native_Receptor_resNum==Native_Rcpt_list(nr),Native_Receptor_resNum==Native_Rcpt_list(nc)),[],'all');
        end

        % Model: Identify Atom Interactions (<5) within CA < 25A
        DISM_Mod_Sub=DISM_Mod(:,Model_Rcpt_No*(c-1)+1:Model_Rcpt_No*c);
        [Mod_row,Mod_col]=find(DISM_Mod_Sub<25);
        for m=1:length(Mod_row)
            mr=Mod_row(m);
            mc=Mod_col(m);
            DISM_Mod_Sub(mr,mc)=min(DISM_Mod_Sub_Atom(Model_Receptor_resNum==Model_Rcpt_list(mr),Model_Receptor_resNum==Model_Rcpt_list(mc)),[],'all');
        end

        DISM_Nat(:,Native_Rcpt_No*(c-1)+1:Native_Rcpt_No*c)=DISM_Nat_Sub;
        DISM_Mod(:,Model_Rcpt_No*(c-1)+1:Model_Rcpt_No*c)=DISM_Mod_Sub;
    end

    DISM_Mod_Aligned=NaN(len_alignidex,len_alignidex*Model_Copy);
    Mod_Idx=ismember(1:length(Model_Rcpt_CA), aseq2_idx(~isnan(aseq2_idx)));
    DISM_Mod_Aligned(~isnan(aseq2_idx),repmat(~isnan(aseq2_idx),1,Native_Copy))=DISM_Mod(Mod_Idx,repmat(Mod_Idx',1,Native_Copy));
    Nat_Int_Res=sum(DISM_Nat<5,'all');
    Nat_Mod_OVLP=sum(DISM_Mod_Aligned(DISM_Nat<5)<5,'all');
    if Nat_Int_Res==0
        Fnat=NaN;
    else
        Fnat=Nat_Mod_OVLP/Nat_Int_Res;
    end  

    %%% iRMS
    % Aligned reNum
    alignidex2=find((Alignment(1,:)~='-')&(Alignment(3,:)~='-'));
    aligned_aseq1_idx=[];
    aligned_aseq2_idx=[];
    for o2=1:length(alignidex2)
        aligned_aseq1_idx=[aligned_aseq1_idx, alignidex2(o2)-sum(aseq1<alignidex2(o2))]; 
        aligned_aseq2_idx=[aligned_aseq2_idx, alignidex2(o2)-sum(aseq2<alignidex2(o2))]; 
    end

    aligned_nat_rcpt_resnum_list=Native_Rcpt_list(aligned_aseq1_idx);
    aligned_mod_rcpt_resnum_list=Model_Rcpt_list(aligned_aseq2_idx);

    Native_Rcpt_Idx=zeros(len_native_rcpt,1);
    Model_Rcpt_Idx=zeros(len_model_rcpt,1);

    for aidx=1:length(alignidex2)
        aligned_nat_atom_no=sum(Native_Receptor_resNum==aligned_nat_rcpt_resnum_list(aidx));
        aligned_mod_atom_no=sum(Model_Receptor_resNum==aligned_mod_rcpt_resnum_list(aidx));
        if aligned_nat_atom_no == aligned_mod_atom_no
            Native_Rcpt_Idx=Native_Rcpt_Idx+(Native_Receptor_resNum==aligned_nat_rcpt_resnum_list(aidx));
            Model_Rcpt_Idx=Model_Rcpt_Idx+(Model_Receptor_resNum==aligned_mod_rcpt_resnum_list(aidx));
        end
    end

    % Paired Native_Receptor (Cahin A)
    Paired_Native_Receptor_Idx=logical(Native_Rcpt_Idx);
    Paired_Native_Receptor_Coord=Native_Receptor_Coord(Paired_Native_Receptor_Idx,:);
    Paired_Native_Receptor_resNum=Native_Receptor_resNum(Paired_Native_Receptor_Idx);
    Paired_Native_Receptor_atomName=Native_Receptor_atomName(Paired_Native_Receptor_Idx);

    % Paired Model_Receptor (Cahin A)
    Paired_Model_Receptor_Idx=logical(Model_Rcpt_Idx);
    Paired_Model_Receptor_Coord=Model_Receptor_Coord(Paired_Model_Receptor_Idx,:);
    Paired_Model_Receptor_resNum=Model_Receptor_resNum(Paired_Model_Receptor_Idx);
    Paired_Model_Receptor_atomName=Model_Receptor_atomName(Paired_Model_Receptor_Idx);

    % Paired Native_Ligand
    Native_Lgnd_Idx=repmat(Native_Rcpt_Idx,Native_Copy,1);
    Paired_Native_Ligand_Idx=logical(Native_Lgnd_Idx);
    Paired_Native_Ligand_Coord=Native_Ligand_Coord(Paired_Native_Ligand_Idx,:);
    Paired_Native_Ligand_resNum=Native_Ligand_resNum(Paired_Native_Ligand_Idx);
    Paired_Native_Ligand_atomName=Native_Ligand_atomName(Paired_Native_Ligand_Idx);

    % Paired Model_Ligand
    Model_Lgnd_Idx=repmat(Model_Rcpt_Idx,Model_Copy,1);
    Paired_Model_Ligand_Idx=logical(Model_Lgnd_Idx);
    Paired_Model_Ligand_Coord=Model_Ligand_Coord(Paired_Model_Ligand_Idx,:);
    Paired_Model_Ligand_resNum=Model_Ligand_resNum(Paired_Model_Ligand_Idx);
    Paired_Model_Ligand_atomName=Model_Ligand_atomName(Paired_Model_Ligand_Idx);

    % Backbone Extraction
    Native_Rcpt_BB_Idx = [find(strcmp(Paired_Native_Receptor_atomName,'N')); find(strcmp(Paired_Native_Receptor_atomName,'CA'));
                         find(strcmp(Paired_Native_Receptor_atomName,'C')); find(strcmp(Paired_Native_Receptor_atomName,'O'))];
    Native_Rcpt_BB=Paired_Native_Receptor_Coord(Native_Rcpt_BB_Idx,:);
    Nat_RBB_resNum=Paired_Native_Receptor_resNum(Native_Rcpt_BB_Idx,:);

    Native_Lgnd_BB_Idx = [find(strcmp(Paired_Native_Ligand_atomName,'N')); find(strcmp(Paired_Native_Ligand_atomName,'CA'));
                         find(strcmp(Paired_Native_Ligand_atomName,'C')); find(strcmp(Paired_Native_Ligand_atomName,'O'))];
    Native_Lgnd_BB=Paired_Native_Ligand_Coord(Native_Lgnd_BB_Idx,:);
    Nat_LBB_resNum=Paired_Native_Ligand_resNum(Native_Lgnd_BB_Idx,:);

    Model_Rcpt_BB_Idx = [find(strcmp(Paired_Model_Receptor_atomName,'N')); find(strcmp(Paired_Model_Receptor_atomName,'CA'));
                        find(strcmp(Paired_Model_Receptor_atomName,'C')); find(strcmp(Paired_Model_Receptor_atomName,'O'))];
    Model_Rcpt_BB=Paired_Model_Receptor_Coord(Model_Rcpt_BB_Idx,:);
    Mod_RBB_resNum=Paired_Model_Receptor_resNum(Model_Rcpt_BB_Idx,:);

    Model_Lgnd_BB_Idx = [find(strcmp(Paired_Model_Ligand_atomName,'N')); find(strcmp(Paired_Model_Ligand_atomName,'CA'));
                        find(strcmp(Paired_Model_Ligand_atomName,'C')); find(strcmp(Paired_Model_Ligand_atomName,'O'))];
    Model_Lgnd_BB=Paired_Model_Ligand_Coord(Model_Lgnd_BB_Idx,:);
    Mod_LBB_resNum=Paired_Model_Ligand_resNum(Model_Lgnd_BB_Idx,:);

    % Interacting Residues (Atom-Paired Distance <10)
    DISM_Idx=ismember(Native_Rcpt_list, unique(Nat_RBB_resNum,'stable'));
    DISM_Int = DISM_Nat(DISM_Idx,repmat(DISM_Idx',1,Native_Copy))<10;


    if sum(DISM_Int,'all')==0
        iRMS=NaN;
    else
        Rcpt_Idx=sum(DISM_Int,2)~=0;
        uniq_Nat_RBB_resNum=unique(Nat_RBB_resNum,'stable');
        uniq_Mod_RBB_resNum=unique(Mod_RBB_resNum,'stable');

        Nat_RBB_resNum_list=uniq_Nat_RBB_resNum(Rcpt_Idx);
        Mod_RBB_resNum_list=uniq_Mod_RBB_resNum(Rcpt_Idx);
        Coord_Rcpt=[];
        for p=1:sum(Rcpt_Idx)
            Coord_Rcpt=[Coord_Rcpt; [Native_Rcpt_BB(Nat_RBB_resNum==Nat_RBB_resNum_list(p),:), Model_Rcpt_BB(Mod_RBB_resNum==Mod_RBB_resNum_list(p),:)]];
        end

        Lgnd_Idx=sum(DISM_Int,1)~=0;
        uniq_Nat_LBB_resNum=unique(Nat_LBB_resNum,'stable');
        uniq_Mod_LBB_resNum=unique(Mod_LBB_resNum,'stable');

        Nat_LBB_resNum_list=uniq_Nat_LBB_resNum(Lgnd_Idx);
        Mod_LBB_resNum_list=uniq_Mod_LBB_resNum(Lgnd_Idx);
        Coord_Lgnd=[];
        for q=1:sum(Lgnd_Idx)
            Coord_Lgnd=[Coord_Lgnd; [Native_Lgnd_BB(Nat_LBB_resNum==Nat_LBB_resNum_list(q),:), Model_Lgnd_BB(Mod_LBB_resNum==Mod_LBB_resNum_list(q),:)]];
        end

        % iRMS
        Combined_Coord=[Coord_Rcpt; Coord_Lgnd];
        [~,~,iRMS]=CoordiExam(Combined_Coord(:,1:3), Combined_Coord(:,4:6));
    end

    %%% LRMS
    % Larger Proteins as Recptors
    if length(Native_Rcpt_BB)>=length(Native_Lgnd_BB)
        [R,T]=CoordiExam(Model_Rcpt_BB,Native_Rcpt_BB);
        Trans=Model_Lgnd_BB*R+T-Native_Lgnd_BB;
    else
        [R,T]=CoordiExam(Model_Lgnd_BB,Native_Lgnd_BB);
        Trans=Model_Rcpt_BB*R+T-Native_Rcpt_BB;
    end

    % LRMS
    LRMS=sqrt(sum((Trans).^2,'all')/length(Trans));

    %%% RRMS
    [R,T]=CoordiExam(Model_Rcpt_BB,Native_Rcpt_BB);
    Trans=Model_Lgnd_BB*R+T-Native_Lgnd_BB;
    RRMS=sqrt(sum((Trans).^2,'all')/length(Trans));

    %%% DockQ
    try
        DockQ=(Fnat+1/(1+(iRMS/1.5)*(iRMS/1.5))+1/(1+(LRMS/8.5)*(LRMS/8.5)))/3;
    catch
        DockQ=NaN;
    end
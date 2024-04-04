%%%% Extract C2_TopRank_Result
% Load C2.DIPER Result
load('Sample_600_1800.mat');

C2_TopRank_Result(1).ID=[];
C2_TopRank_Result(1).Lowest_RMSD_95=[];
C2_TopRank_Result(1).Output_Idx=[];
C2_TopRank_Result(1).Output_RMSD_95=[];
C2_TopRank_Result(1).Best_Idx=[];
C2_TopRank_Result(1).Best_RMSD_95=[];

for r=1:size(Sample,2)
    % Extract Top 25 Results
    Pose_No=size(Sample(r).Output_Idx,1);
    if Pose_No>25, Pose_No=25; end

    % Record ID, Idx, RMSD95 & Best Result
    C2_TopRank_Result(r).ID=Sample(r).ID;
    C2_TopRank_Result(r).Lowest_RMSD_95=Sample(r).Lowest_RMSD_95;
    C2_TopRank_Result(r).Output_Idx=Sample(r).Output_Idx(1:Pose_No);
    C2_TopRank_Result(r).Output_RMSD_95=Sample(r).Output_RMSD_95(1:Pose_No);

    % Best Result
    [val, row]=min(C2_TopRank_Result(r).Output_RMSD_95);
    C2_TopRank_Result(r).Best_Idx=C2_TopRank_Result(r).Output_Idx(row);
    C2_TopRank_Result(r).Best_RMSD_95=val;
end

save('C2_TopRank_Result.mat', 'C2_TopRank_Result');

%
%
%

%%%% Construct C2_TopRank_Result Homodimers
% load('C2_TopRank_Result.mat')

% Load PIPER C2-Symmetry/Rotation Matrix
Rots_list=table2array(readtable('C2_rots.txt'));

% Default Parameters
prm=0; set=['ft.00', num2str(prm), '.00']; Limit=1800;

% Combination List for Translation Table
itr=0; l=20; List=zeros((l*2+1)^3,3);
for a=-l:l
    for b=-l:l
        for c=-l:l
            itr=itr+1;
            List(itr,1:3)=[a,b,c];
        end
    end
end

% Establish Deposited File of C2.DIPER Top 25 Outputs
mkdir('C2.Top25_Pose');

for r=1:size(C2_TopRank_Result,2)

    %%%% Cull Crystal Homodimers for RMSD Examination
    % Load PDB File
    PDB=pdbread(['PDB/',lower(C2_TopRank_Result(r).ID),'.pdb']);
    File=fopen(['PDB/',lower(C2_TopRank_Result(r).ID),'.pdb']);
    Ans.Model.Atom=repmat(PDB.Model.Atom,1,2);

    % Symmetry Operator
    Operator=[];
    while(1)
        line=fgetl(File);
        if isequal(sscanf(line(1:10),'%c'),'REMARK 300'), break, end

        if isequal(sscanf(line(1:10),'%c'),'REMARK 290')
            Operator=[Operator; sscanf(line(25:79),'%f %f %f')'];
        end
    end

    % Unit Cell Parameter
    Unit_Cell=[];
    while(1)
        line=fgetl(File);
        if isequal(sscanf(line(1:3),'%c'),'CRY')
            Unit_Cell=sscanf(line(8:55),'%f %f %f');
            break,
        end
    end
    fclose(File);

    % Extract Coordinates of Alpha Carbon (CA)
    PDB_Model=[PDB.Model(1).Atom.X; PDB.Model(1).Atom.Y; PDB.Model(1).Atom.Z]';

    PDB_resSeq_Idx=[PDB.Model(1).Atom.resSeq]';
    pseudo_PDB_CA_Idx=strcmp({PDB.Model(1).Atom.AtomName}, 'CA')';
    PDB_resSeq_Unique_Idx=unique(PDB_resSeq_Idx(pseudo_PDB_CA_Idx));

    PDB_CA_Idx=zeros(size(PDB_resSeq_Unique_Idx,1),1);
    for RSUidx=1:size(PDB_resSeq_Unique_Idx,1)
        PDB_CA_Idx(RSUidx,1)=find((PDB_resSeq_Idx==PDB_resSeq_Unique_Idx(RSUidx))&pseudo_PDB_CA_Idx, 1, 'first');
    end
    Asy_Unit=PDB_Model(PDB_CA_Idx,:);

    %%% Construct Reduced Unit Box
    % Basis: (a, b, c) to (x, y, z)
    basis=eye(3);
    basis(1,2)=cosd(Unit_Cell(6));
    basis(2,2)=sind(Unit_Cell(6));
    basis(1,3)=cosd(Unit_Cell(5));
    basis(2,3)=(cosd(Unit_Cell(4))-basis(1,2)*basis(1,3))/basis(2,2);
    basis(3,3)=sqrt(1-basis(1,3)^2 -basis(2,3)^2);
    Unit_Cell_extent=basis.*repmat(Unit_Cell(1:3)',3,1);

    % Set of Symmertic Mates
    Unit=[];
    Sym_No=size(Operator,1)/3;
    for j=1:Sym_No
        SMTRY_Unit=Asy_Unit*Operator(3*j-2:3*j,1:3)'+ Operator(3*j-2:3*j,4)';

        % Shift
        Shift=List*Unit_Cell_extent';
        Translation_Table=mean(Asy_Unit)-(mean(SMTRY_Unit)+Shift);

        % Store Shifted SMTRY_Unit
        [~,idx]=min(sqrt(sum(Translation_Table.^2,2)));
        SMTRY_Unit=SMTRY_Unit+Shift(idx,:);
        Unit=[Unit; SMTRY_Unit];
    end

    % Unit Box
    Unit_Box=Unit;
    for x=-2:2
        for y=-2:2
            for z=-2:2
                if ~(x==0 && y==0 && z==0)
                    NB_Unit=Unit+[x,y,z]*Unit_Cell_extent';
                    Unit_Box=[Unit_Box; NB_Unit];
                end
            end
        end
    end

    % Identify Closest Packing Neighbors (PDB_NBs)
    NB_Cutoff=12;
    while(1)
        i=0; PDB_NB=struct('Coord', []);
        for j=2:size(Unit_Box,1)/size(Asy_Unit,1)
            DISM=pdist2(Asy_Unit, Unit_Box(size(Asy_Unit,1)*(j-1)+1:size(Asy_Unit,1)*j,:));
            if min(min(DISM))<NB_Cutoff
                i=i+1;
                PDB_NB(i).Coord=Unit_Box(size(Asy_Unit,1)*(j-1)+1:size(Asy_Unit,1)*j,:);
            end
        end

        NeiB_Idx=[];
        RMSD_Table=zeros(size(PDB_NB,2));
        for k=1:size(PDB_NB,2)
            for l=1:size(PDB_NB,2)
                [~,~, eRMSD_1]=CoordiExam([Asy_Unit; PDB_NB(k).Coord], [Asy_Unit; PDB_NB(l).Coord]);
                [~,~, eRMSD_2]=CoordiExam([Asy_Unit; PDB_NB(k).Coord], [PDB_NB(l).Coord; Asy_Unit]);
                RMSD_Table(k,l)=min(eRMSD_1, eRMSD_2);
            end
            NeiB_Idx=[NeiB_Idx, find(RMSD_Table(k,:)<1, 1)];
        end

        % Generate Closest/"Unique" Packing Neighbors (PDB_UNBs)
        UNB=unique(NeiB_Idx);
        m=0; PDB_UNB=struct('Coord',[]);
        for n=1:length(UNB)
            R=CoordiExam(Asy_Unit, PDB_NB(UNB(n)).Coord);
            Eigen=sortrows([real(eig(R)), imag(eig(R))])*[1;sqrt(-1)];

            % Examine C2 Symmetry (Angle)
            if norm(Eigen-[-1; -1; 1])<6*pi/180
                [V,D]=eig(R);
                [eig_v,eig_idx]=max(real(diag(D)'));

                % Examine C2 Symmetry (Axis: loose criteria (<5) for tolerance of CA mass center shift)
                if abs(acosd((mean(Asy_Unit)-mean(PDB_NB(UNB(n)).Coord))/norm(mean(Asy_Unit)-mean(PDB_NB(UNB(n)).Coord))*V(:,eig_idx))-90)<5
                    m=m+1;
                    PDB_UNB(m).Coord=[Asy_Unit; PDB_NB(UNB(n)).Coord];
                end
            end
        end

        % Termination
        if m>0
            break,
        else
            % Expand Searching Region
            NB_Cutoff=NB_Cutoff+6;
        end
    end

    %

    %%%% Construct C2-DIPER Homodimers
    % Establish Deposited File
    mkdir(['C2.Top25_Pose/', lower(C2_TopRank_Result(r).ID)]);

    % Load PIPER Docking Result
    DIPER=fopen([lower(C2_TopRank_Result(r).ID),'/',set]);
    PIPER_result=zeros(Limit, 10);
    for p=1:Limit
        line=fgetl(DIPER);
        PIPER_result(p,1:10)=sscanf(line,'%f')';
    end
    fclose(DIPER);

    % Extract CA Coordinate of Truncated Alphafold2 Monomer
    AF=pdbread([lower(C2_TopRank_Result(r).ID),'_trun_AF_pnon.pdb']);
    Model=[AF.Model.Atom.X; AF.Model.Atom.Y; AF.Model.Atom.Z]';
    AF_CA_Idx=strcmp({AF.Model.Atom.AtomName}, 'CA')';
    Rcpt=Model(AF_CA_Idx,:);

    % Construct Representative Pose (RePo)
    AF.Model.Atom=repmat(AF.Model.Atom,1,2);
    Pose_No=size(C2_TopRank_Result(r).Output_Idx,1);

    % Extract 1-letter Sequences
    PDB_Seq=aminolookup(strrep(strrep([PDB.Model(1).Atom(PDB_CA_Idx).resName], 'SEC', 'CYS'), 'UNK', 'GLY'));
    AF_Seq=aminolookup([AF.Model.Atom(AF_CA_Idx).resName]);

    seq1=[PDB_Seq, PDB_Seq]; seq2=[AF_Seq, AF_Seq];

    % Global Sequence Alignment to Extract Superimposed CA for RMSD Calculation
    [~, Alignment] = nwalign(seq1,seq2);
    alignidex=find((Alignment(1,:)~='-')&(Alignment(3,:)~='-'));
    aseq1=find(Alignment(1,:)=='-');
    aseq2=find(Alignment(3,:)=='-');

    aseq1_idx=[]; aseq2_idx=[];
    for o=1:length(alignidex)
        aseq1_idx=[aseq1_idx, alignidex(o)-sum(aseq1<alignidex(o))];
        aseq2_idx=[aseq2_idx, alignidex(o)-sum(aseq2<alignidex(o))];
    end

    %

    % Construct C2-DIPER Homodimers
    for q=1:Pose_No
        Result_Idx=C2_TopRank_Result(r).Output_Idx(q);
        Trans_M=PIPER_result(Result_Idx,2:4);
        Rot_Info=Rots_list(PIPER_result(Result_Idx,1)+1,:);
        Rot_M=[Rot_Info(2:4); Rot_Info(5:7); Rot_Info(8:10)];
        Lgnd=(Rcpt-mean(Rcpt))*Rot_M'+repmat(Trans_M+mean(Rcpt),size(Rcpt,1),1);

        % Output C2-Symmetric PDB File
        [R_AF,T_AF]=CoordiExam(Rcpt, Lgnd);
        Model_Lgnd=Model*R_AF+repmat(T_AF,size(Model,1),1);

        len_AF=size(Model_Lgnd,1);
        for s=1:len_AF
            AF.Model.Atom(s).chainID='A';
            AF.Model.Atom(len_AF+s).chainID='B';
            AF.Model.Atom(s).AtomSerNo=s;
            AF.Model.Atom(len_AF+s).AtomSerNo=len_AF+s;
            AF.Model.Atom(len_AF+s).X=Model_Lgnd(s,1);
            AF.Model.Atom(len_AF+s).Y=Model_Lgnd(s,2);
            AF.Model.Atom(len_AF+s).Z=Model_Lgnd(s,3);
        end

        warning('off','all');
        pdbwrite(['C2.Top25_Pose/', lower(C2_TopRank_Result(r).ID), '/', lower(C2_TopRank_Result(r).ID), set(3:7), num2str(Result_Idx), '.pdb'], AF);

        %

        % RMSD 95 calculation (ATTTENTION: ORDER OF MONOMERS)
        PDB_1=[Rcpt; Lgnd]; PDB_CA_1=PDB_1(aseq2_idx,:);
        PDB_2=[Lgnd; Rcpt]; PDB_CA_2=PDB_2(aseq2_idx,:);
        len_CA=size(PDB_CA_1,1);

        RMSD_95_Table=[];
        for w=1:size(PDB_UNB,2)
            [~,~,~,~, Sort_Sum_Square_1]=CoordiExam(PDB_UNB(w).Coord(aseq1_idx,:), PDB_CA_1);
            [~,~,~,~, Sort_Sum_Square_2]=CoordiExam(PDB_UNB(w).Coord(aseq1_idx,:), PDB_CA_2);
            RMSD_95_1=sqrt(sum(Sort_Sum_Square_1(1:ceil(len_CA*0.95)))/ceil(len_CA*0.95));
            RMSD_95_2=sqrt(sum(Sort_Sum_Square_2(1:ceil(len_CA*0.95)))/ceil(len_CA*0.95));
            RMSD_95_Table(w)=min(RMSD_95_1, RMSD_95_2);
        end
        [~, row]=min(RMSD_95_Table,[],'all');

        % Output Ans PDB File
        [R_PDB,T_PDB]=CoordiExam(Asy_Unit, PDB_UNB(row).Coord(length(Asy_Unit)+1:end,:));
        PDB_Model_Pair=PDB_Model*R_PDB+repmat(T_PDB,size(PDB_Model,1),1);

        len_PDB=size(PDB_Model_Pair,1);
        for t=1:len_PDB
            Ans.Model.Atom(t).chainID='A';
            Ans.Model.Atom(len_PDB+t).chainID='B';
            Ans.Model.Atom(t).AtomSerNo=t;
            Ans.Model.Atom(len_PDB+t).AtomSerNo=len_PDB+t;
            Ans.Model.Atom(len_PDB+t).X=PDB_Model_Pair(t,1);
            Ans.Model.Atom(len_PDB+t).Y=PDB_Model_Pair(t,2);
            Ans.Model.Atom(len_PDB+t).Z=PDB_Model_Pair(t,3);
        end

        warning('off','all');
        pdbwrite(['C2.Top25_Pose/', lower(C2_TopRank_Result(r).ID), '/', lower(C2_TopRank_Result(r).ID), set(3:7), num2str(Result_Idx), '_Ans.pdb'], Ans);
    end
    r
end

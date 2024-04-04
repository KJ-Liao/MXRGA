%%%% Extract C4_TopRank_Result
% Load C2.TopRank Result
load('C2_TopRank_Result.mat')

% Load C4.DIPER Result
load('C4_DIPER_Result_P41212.mat');

% Default Setting
Top_Pose_No=25;

C4_TopRank_Result(1).ID=[];
C4_TopRank_Result(1).Lowest_C4_RMSD_95=[];
C4_TopRank_Result(1).Lowest_Overall_RMSD_95=[];
C4_TopRank_Result(1).Output_Idx=[];
C4_TopRank_Result(1).Output_C4_RMSD_95=[];
C4_TopRank_Result(1).Output_Overall_RMSD_95=[];
C4_TopRank_Result(1).Best_Idx=[];
C4_TopRank_Result(1).Best_C4_RMSD_95=[];
C4_TopRank_Result(1).Best_Overall_RMSD_95=[];

for r=1:size(Sample,2)
    % Extract Top Results
    Pose_No=size(Sample(r).Output_Idx,1);
    if Pose_No>Top_Pose_No, Pose_No=Top_Pose_No; end

    % Record ID, Idx, RMSD95 & Best Result
    C4_TopRank_Result(r).ID=Sample(r).ID;
    C4_TopRank_Result(r).Lowest_C4_RMSD_95=Sample(r).Lowest_C4_RMSD_95;
    C4_TopRank_Result(r).Lowest_Overall_RMSD_95=Sample(r).Lowest_Overall_RMSD_95;
    C4_TopRank_Result(r).Output_Idx=Sample(r).Output_Idx(1:Pose_No);
    C4_TopRank_Result(r).Output_C4_RMSD_95=Sample(r).Output_C4_RMSD_95(1:Pose_No);
    C4_TopRank_Result(r).Output_Overall_RMSD_95=Sample(r).Output_Overall_RMSD_95(1:Pose_No);

    % Best Result
    [~, row_20]=min(C4_TopRank_Result(r).Output_Overall_RMSD_95(1:Pose_No));
    [~, row_15]=min(C4_TopRank_Result(r).Output_Overall_RMSD_95(1:min(Pose_No,15)));
    [~, row_12]=min(C4_TopRank_Result(r).Output_Overall_RMSD_95(1:min(Pose_No,12)));

    C4_TopRank_Result(r).Best_Idx=[C4_TopRank_Result(r).Output_Idx(row_20), C4_TopRank_Result(r).Output_Idx(row_15), C4_TopRank_Result(r).Output_Idx(row_12)];
    C4_TopRank_Result(r).Best_C4_RMSD_95=[C4_TopRank_Result(r).Output_C4_RMSD_95(row_20), C4_TopRank_Result(r).Output_C4_RMSD_95(row_15), C4_TopRank_Result(r).Output_C4_RMSD_95(row_12)];
    C4_TopRank_Result(r).Best_Overall_RMSD_95=[C4_TopRank_Result(r).Output_Overall_RMSD_95(row_20), C4_TopRank_Result(r).Output_Overall_RMSD_95(row_15), C4_TopRank_Result(r).Output_Overall_RMSD_95(row_12)];
end

save('C4_TopRank_Result.mat', 'C4_TopRank_Result');

%
%
%

%%%% Rank Top 300 from 25*25 Table
% Overall Result
clear Result;
Result(1).ID=[];
Result(1).Idx_Table=[];
Result(1).C4_RMSD_Table=[];
Result(1).All_RMSD_Table=[];
Result(1).Linear_Table=[];

C2_No=25; C4_No=25;

for m=1:size(C2_TopRank_Result,2)

    % Index of Corresponding ID
    Idx_list=find(contains({C4_TopRank_Result.ID}', lower(C2_TopRank_Result(m).ID)));
    if size(Idx_list,1)>C2_No
        Idx_list=Idx_list(1:C2_No,:);
    end

    % Pre-allocate Memorry for Restoration
    Idx_M=zeros(C2_No,C4_No);
    C4_RMSD_M=ones(C2_No,C4_No)*20;
    All_RMSD_M=ones(C2_No,C4_No)*20;

    for n=1:size(Idx_list,1)
        % Extract Top Results
        r=Idx_list(n);
        Pose_No=size(Sample(r).Output_Idx,1);
        if Pose_No>C4_No, Pose_No=C4_No; end

        % Record Idx & RMSD95
        Idx_M(n,1:Pose_No)=Sample(r).Output_Idx(1:Pose_No);
        C4_RMSD_M(n,1:Pose_No)=Sample(r).Output_C4_RMSD_95(1:Pose_No);
        All_RMSD_M(n,1:Pose_No)=Sample(r).Output_Overall_RMSD_95(1:Pose_No);
    end

    Result(m).ID=lower(C2_TopRank_Result(m).ID);
    Result(m).Idx_Table=Idx_M;
    Result(m).C4_RMSD_Table=C4_RMSD_M;
    Result(m).All_RMSD_Table=All_RMSD_M;
end

% Anti-Diagonal Rank
Input=zeros(C2_No,C4_No);

[r,c]=size(Input);
Mat=repmat((1:r)',1,c)+repmat(1:c,r,1);

t=0;
Idx_list=unique(Mat);
Linear_Idx=zeros(r*c,3);
for m=1:size(Idx_list,1)
    [row,col]=find(Mat==Idx_list(m));
    % rc_Ind=sortrows([row,col],1,'ascend');
    rc_Ind=sortrows([row,col],2,'ascend');
    for rc=1:size(rc_Ind,1)
        t=t+1;
        Linear_Idx(t,:)=[t, rc_Ind(rc,:)];
        Input(rc_Ind(rc,1),rc_Ind(rc,2))=t;
    end
end

Linear_Table=[Linear_Idx, zeros(r*c,3)];
for i=1:size(Result,2)
    Result_Idx_Table=Result(i).Idx_Table;
    Result_C4_RMSD_Table=Result(i).C4_RMSD_Table;
    Result_All_RMSD_Table=Result(i).All_RMSD_Table;
    for j=1:size(Linear_Table,1)
        Linear_Table(j,4)=Result_Idx_Table(Linear_Table(j,2), Linear_Idx(j,3));
        Linear_Table(j,5)=Result_C4_RMSD_Table(Linear_Table(j,2), Linear_Idx(j,3));
        Linear_Table(j,6)=Result_All_RMSD_Table(Linear_Table(j,2), Linear_Idx(j,3));
    end
    Result(i).Linear_Table=Linear_Table(Linear_Table(:,4)~=0,:);
end

save('TopRank_Result.mat', 'Result');

% Statistics
Boundary=300;
minResult=ones(size(Result,2),1)*20;
for g=1:size(Result,2)
    Table=Result(g).Linear_Table;
    Limit=size(Table,1);
    if Limit>Boundary
        Limit=Boundary;
    end
    minResult(g,1)=min(Table(1:Limit,6));
end

%

% format shortG
Count=size(Result,2);

Brief_Sum=[[5.0, sum(minResult<5.0)/Count*100];
    [5.5, sum(minResult<5.5)/Count*100];
    [6.0, sum(minResult<6.0)/Count*100];
    [6.5, sum(minResult<6.5)/Count*100];
    [7.0, sum(minResult<7.0)/Count*100]]

%%%% P41212
%      300; 250; 200; 150
% 5       28  26 25 18
% 5.5     34  31 29 22
% 6       43  40 39 33
% 6.5     45  43 42 36
% 7       53  50 48 42

%
%
%

%%%% Construct C4_TopRank_Result Homodimers
% Parallel Computatoin
% Activate the Parallel Cores (default: 6 cores)
parpool(6);
clear all

% Establish Deposited File
mkdir('3. TopRank_300_Pose');

% Load C4.DIPER Result
load('2. C4_TopRank/TopRank_Result.mat');
load('2. C4_TopRank/C2_TopRank_Result.mat');

% Load PIPER C4-Symmetry/Rotation Matrix
Destn='0. C4.DIPER_P41212/';
Rots_list=table2array(readtable([Destn,'C4_rots.txt']));

% Default Parameters
prm=0; set=['ft.00', num2str(prm), '.00'];

% Combination List for Translation Table
itr=0; l=20; List=zeros((l*2+1)^3,3);
for p=-l:l
    for q=-l:l
        for r=-l:l
            itr=itr+1;
            List(itr,1:3)=[p,q,r];
        end
    end
end

%%%% Construct C4-DIPER Homodimers
parfor r=1:size(Result,2)

    % Establish Deposited File
    mkdir(['3. TopRank_300_Pose/',Result(r).ID]);

    %%%% Cull Crystal Homodimers for RMSD Examination
    % Load PDB File
    PDB=pdbread([Destn,'PDB/',Result(r).ID,'.pdb']);
    File=fopen([Destn,'PDB/',Result(r).ID,'.pdb']);

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

    % Sequence
    PDB_Seq=aminolookup(strrep(strrep([PDB.Model(1).Atom(PDB_CA_Idx).resName], 'SEC', 'CYS'), 'UNK', 'GLY'));

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

    % Identify Closest Packing Neighbors (PDB_NBs, Clostest CA Dist.<12A)
    NB_Cutoff=12;
    while(1)
        i=0; PDB_NB=struct('Coord',[]);
        for j=2:size(Unit_Box,1)/size(Asy_Unit,1)
            DISM=pdist2(Asy_Unit, Unit_Box(size(Asy_Unit,1)*(j-1)+1:size(Asy_Unit,1)*j,:));
            if min(DISM,[],'all')<NB_Cutoff
                i=i+1;
                PDB_NB(i).Coord=Unit_Box(size(Asy_Unit,1)*(j-1)+1:size(Asy_Unit,1)*j,:);
            end
        end

        NeiB_Idx=zeros(size(PDB_NB));
        for k=1:size(PDB_NB,2)
            l=0;
            while(1)
                l=l+1;
                [~,~,eRMSD_1]=CoordiExam([Asy_Unit; PDB_NB(k).Coord], [Asy_Unit; PDB_NB(l).Coord]);
                [~,~,eRMSD_2]=CoordiExam([Asy_Unit; PDB_NB(k).Coord], [PDB_NB(l).Coord; Asy_Unit]);
                if min(eRMSD_1,eRMSD_2)<1
                    NeiB_Idx(k)=l; break,
                end
            end
        end

        % Generate Closest/Unique "NON-C2" Packing Neighbors (PDB_nC2_UNB)
        UNB=unique(NeiB_Idx);
        p_C2=0; PDB_C2_UNB=struct('Coord',[]);
        p_nC2=0; PDB_nC2_UNB=struct('Coord',[]);
        for q=1:length(UNB)
            R=CoordiExam(Asy_Unit, PDB_NB(UNB(q)).Coord);
            [V,D]=eig(R); [eig_v,eig_idx]=max(real(diag(D)'));
            Eigen=sortrows([real(eig(R)), imag(eig(R))])*[1;sqrt(-1)];

            % Examine C2 Symmetry (Angle & Axis: loose criteria (<5) for tolerance of CA mass center shift)
            if ~(norm(Eigen-[-1; -1; 1])<6*pi/180 && abs(acosd((mean(Asy_Unit)-mean(PDB_NB(UNB(q)).Coord))/norm(mean(Asy_Unit)-mean(PDB_NB(UNB(q)).Coord))*V(:,eig_idx))-90)<5)
                p_nC2=p_nC2+1;
                PDB_nC2_UNB(p_nC2).Coord=[Asy_Unit; PDB_NB(UNB(q)).Coord];
            else
                p_C2=p_C2+1;
                PDB_C2_UNB(p_C2).Coord=[Asy_Unit; PDB_NB(UNB(q)).Coord];
            end
        end

        % Termination
        if p_C2>0
            break,
        else
            % Expand Searching Region
            NB_Cutoff=NB_Cutoff+6;
        end
    end

    %
    %
    %

    %%%% Tetramer Generation
    Total_Pose_No=size(Result(r).Linear_Table,1);
    if Total_Pose_No>300, Total_Pose_No=300; end

    Res_Cutoff=12;
    Qualified_Idx=find(Result(r).Linear_Table(1:Total_Pose_No,end)<Res_Cutoff);

    if ~isempty(Qualified_Idx)
        for QI=1:size(Qualified_Idx,1)
            rp=Qualified_Idx(QI,1);

            %%% Refine the Homodimer (HMD) into C2-Symmetry
            % Extract Coordinates of Alpha Carbon (CA)
            HMD_ID=[Result(r).ID, '.000.', num2str(C2_TopRank_Result(r).Output_Idx(Result(r).Linear_Table(rp,2)))];
            HMD=pdbread([Destn,HMD_ID,'_pnon.pdb']);
            HMD_CA_Idx=strcmp({HMD.Model.Atom.AtomName}, 'CA')';

            % Sequence
            AF_Seq=aminolookup(strrep(strrep([HMD.Model.Atom(HMD_CA_Idx).resName], 'SEC', 'CYS'), 'UNK', 'GLY'));

            % Fit the Model_Rcpt/Lgnd into C2-Symmetry
            Model=[HMD.Model.Atom.X; HMD.Model.Atom.Y; HMD.Model.Atom.Z]';
            Rcpt=Model(HMD_CA_Idx,:);
            Rcpt_A=Rcpt(1:size(Rcpt,1)/2,:);
            Rcpt_B=Rcpt(size(Rcpt,1)/2+1:end,:);

            % Calculate the Rot Axis & Spin Angle of Model_Rcpt/Lgnd (for Angle Correction)
            Model_Rcpt=Model(1:size(Model,1)/2,:);
            Model_Lgnd=Model(size(Model,1)/2+1:end,:);
            Model_R=CoordiExam(Model_Rcpt, Model_Lgnd);
            Model_Axis_Rot=quat2axang(rotm2quat(Model_R));

            % Compute C2 Correction Matrix for Model_Rcpt (C2_Rot_M_R) & Model_Lgnd (C2_Rot_M_L)
            C2_Rot_M_R=axang2rotm([Model_Axis_Rot(1:3), -(pi-Model_Axis_Rot(4))/2]);
            C2_Rot_M_L=axang2rotm([Model_Axis_Rot(1:3), pi-(pi-Model_Axis_Rot(4))/2]);

            % Calibirate the Angle Error of C2 Symmetry
            pseudo_C2_Rcpt=Model_Rcpt*C2_Rot_M_R+repmat(mean(Model_Rcpt)-mean(Model_Rcpt*C2_Rot_M_R),size(Model_Rcpt,1),1);
            pseudo_C2_Lgnd=Model_Rcpt*C2_Rot_M_L+repmat(mean(Model_Lgnd)-mean(Model_Rcpt*C2_Rot_M_L),size(Model_Lgnd,1),1);

            % Calculate the Rot Axis & Spin Angle of pseudo_C2_Rcpt/Lgnd (for Translation Correction)
            pseudo_C2_R=CoordiExam(pseudo_C2_Rcpt, pseudo_C2_Lgnd);
            pseudo_C2_Axis_Rot=quat2axang(rotm2quat(pseudo_C2_R));

            % Calibirate the Translation Error of C2 Symmetry
            t=dot(mean(pseudo_C2_Rcpt)-mean(pseudo_C2_Lgnd), pseudo_C2_Axis_Rot(1:3));
            MassC_lgnd=mean(pseudo_C2_Lgnd)+pseudo_C2_Axis_Rot(1:3)*t;
            fitted_pseudo_C2_Lgnd=(Model_Rcpt-mean(Model_Rcpt))*C2_Rot_M_L+repmat(MassC_lgnd,size(pseudo_C2_Lgnd,1),1);

            % Divide Fitted Rcpt Homodimer (f_Rcpt) into two Monomers: f_Rcpt_A, f_Rcpt_B
            f_Rcpt=[pseudo_C2_Rcpt; fitted_pseudo_C2_Lgnd];
            f_Rcpt_A=pseudo_C2_Rcpt;
            f_Rcpt_B=fitted_pseudo_C2_Lgnd;
            Direct_Point=mean([f_Rcpt_A(1:ceil(size(f_Rcpt_A,1)/2),:); f_Rcpt_B(1:ceil(size(f_Rcpt_B,1)/2),:)]);

            % Adjust Fitted Rcpt_A/B Mass Center (MC) to Z=0 Plane
            XY_dist=norm(mean(f_Rcpt_A)-mean(f_Rcpt_B))/2/sqrt(2);
            New_Direct_Point=sqrt(norm(mean(f_Rcpt_B)-Direct_Point)^2-(norm(mean(f_Rcpt_A)-mean(f_Rcpt_B))/2)^2)/sqrt(2);

            Origin_MC_Pair=[mean(f_Rcpt_A); mean(f_Rcpt_B); Direct_Point];
            New_MC_Pair=[XY_dist, -XY_dist, 0; -XY_dist, XY_dist, 0; New_Direct_Point, New_Direct_Point, 0];

            [Rot,Trn]=CoordiExam(Origin_MC_Pair, New_MC_Pair);
            r_Rcpt=f_Rcpt*Rot+repmat(Trn, size(f_Rcpt,1), 1);

            % Refined/Normalized r_Rcpt_A, r_Rcpt_B
            r_CA_Rcpt=r_Rcpt(HMD_CA_Idx,:);
            r_CA_Rcpt_A=r_CA_Rcpt(1:size(r_CA_Rcpt,1)/2,:);
            r_CA_Rcpt_B=r_CA_Rcpt_A*[0, 1, 0; 1, 0, 0; 0, 0, -1];

            %

            %%%% Determine C2-Homodimers for Further RMSD 95 Calculation
            toXYZ=r_CA_Rcpt_A; toXYZ_pair=r_CA_Rcpt_B;

            % Global Sequence Alignment to Extract Superimposed CA for RMSD Calculation
            [~, Alignment] = nwalign([PDB_Seq,PDB_Seq], AF_Seq);
            alignidex=find((Alignment(1,:)~='-')&(Alignment(3,:)~='-'));
            aseq1=find(Alignment(1,:)=='-');
            aseq2=find(Alignment(3,:)=='-');

            aseq1_idx=[]; aseq2_idx=[];
            for alex=1:length(alignidex)
                aseq1_idx=[aseq1_idx, alignidex(alex)-sum(aseq1<alignidex(alex))];
                aseq2_idx=[aseq2_idx, alignidex(alex)-sum(aseq2<alignidex(alex))];
            end

            % RMSD 95 calculation
            HMD_RMSD_Table=zeros(1,size(PDB_C2_UNB,2));
            AF_PDB=[toXYZ; toXYZ_pair];

            for HRT=1:size(PDB_C2_UNB,2)
                [~,~,~,~,Sort_Sum_Square]=CoordiExam(PDB_C2_UNB(HRT).Coord(aseq1_idx,:), AF_PDB(aseq2_idx,:));
                HMD_RMSD_Table(HRT)=sqrt(sum(Sort_Sum_Square(1:ceil(length(aseq2_idx)*0.95)))/ceil(length(aseq2_idx)*0.95));
            end

            [~,C2_col]=min(HMD_RMSD_Table);
            C2_Homo=PDB_C2_UNB(C2_col).Coord;
            C2_Homo_A=C2_Homo(1:size(C2_Homo,1)/2,:);
            C2_Homo_B=C2_Homo(size(C2_Homo,1)/2+1:end,:);

            %

            %%%% Construct Representative Pose (RePo)
            Result_C4=fopen([Destn,HMD_ID,'/',set]);
            PIPER_result=zeros(size(Rots_list,1), 10);
            for p=1:size(Rots_list,1)
                line=fgetl(Result_C4);
                PIPER_result(p,1:10)=sscanf(line,'%f')';
            end
            fclose(Result_C4);

            % Construct Top Outputs
            Output_C2C4=struct('Model',[]);
            Output_SiXt=struct('Model',[]);
            Output_Xtal=struct('Model',[]);
            Output_C2C4.Model.Atom=repmat(HMD.Model.Atom,1,2);
            Output_SiXt.Model.Atom=repmat(HMD.Model.Atom,1,2);
            Output_Xtal.Model.Atom=repmat(PDB.Model(1).Atom,1,4);

            Result_Idx=Result(r).Linear_Table(rp,4);
            Trans_M=PIPER_result(Result_Idx,2:4);
            Rot_Info=Rots_list(PIPER_result(Result_Idx,1)+1,:);
            Rot_M=[Rot_Info(2:4); Rot_Info(5:7); Rot_Info(8:10)];
            Lgnd_Model=(Model-mean(Model))*Rot_M'+repmat(Trans_M+mean(Model),size(Model,1),1);

            % Divide Lgnd Homodimer into two Monomers: Lgnd_A, Lgnd_B
            Lgnd=Lgnd_Model(HMD_CA_Idx,:);
            Lgnd_A=Lgnd(1:size(Lgnd,1)/2,:);
            Lgnd_B=Lgnd(size(Lgnd,1)/2+1:end,:);

            % Determine the Interaction Pair (Force into 1-5/7-4 Interaction Pair)
            if norm(mean(Rcpt_A)-mean(Lgnd))<norm(mean(Rcpt_B)-mean(Lgnd))
                % 1-5 Interaction Pair
                [r_R,r_T]=CoordiExam(Rcpt_A, r_CA_Rcpt_A);
                fromXYZ=Lgnd_B*r_R+repmat(r_T, size(r_CA_Rcpt_A,1), 1);
                fromXYZ_pair=Lgnd_A*r_R+repmat(r_T, size(r_CA_Rcpt_A,1), 1);
            else
                % 7-3 Interaction Pair
                [r_R,r_T]=CoordiExam(Rcpt_B, r_CA_Rcpt_A);
                fromXYZ=Lgnd_A*r_R+repmat(r_T, size(r_CA_Rcpt_A,1), 1);
                fromXYZ_pair=Lgnd_B*r_R+repmat(r_T, size(r_CA_Rcpt_A,1), 1);
            end

            % Determine the Closest C4 Interacting-Dimer
            C4_ASU=struct('Coord',[]);
            if norm(mean(toXYZ)-mean(fromXYZ))>norm(mean(toXYZ)-mean(fromXYZ_pair))
                C4_ASU(1).Coord=[toXYZ; fromXYZ_pair];
                C4_ASU(2).Coord=[fromXYZ_pair; toXYZ];
            else
                C4_ASU(1).Coord=[toXYZ; fromXYZ];
                C4_ASU(2).Coord=[fromXYZ; toXYZ];
            end

            len=size(aseq2_idx,2);
            RMSD_Table=zeros(size(PDB_nC2_UNB,2),2);
            for v=1:2
                for u=1:size(PDB_nC2_UNB,2)
                    [~,~,~,~,Sort_Sum_Square]=CoordiExam(PDB_nC2_UNB(u).Coord(aseq1_idx,:), C4_ASU(v).Coord(aseq2_idx,:));
                    RMSD_Table(u,v)=sqrt(sum(Sort_Sum_Square(1:ceil(len*0.95)))/ceil(len*0.95));
                end
            end

            % Pose Re-Generation
            [row, col]=find(RMSD_Table==min(RMSD_Table,[],'all'));

            [Rot_1,Trn_1]=CoordiExam(C2_Homo_A, PDB_nC2_UNB(row).Coord(1:size(C2_Homo_A,1),:));
            [Rot_2,Trn_2]=CoordiExam(C2_Homo_A, PDB_nC2_UNB(row).Coord(size(C2_Homo_A,1)+1:end,:));
            Remained_Part_Xtal=[C2_Homo_B*Rot_1+Trn_1; C2_Homo_B*Rot_2+Trn_2];
            Xtal=[PDB_nC2_UNB(row).Coord; Remained_Part_Xtal];

            [Rt_1,Tr_1]=CoordiExam(r_CA_Rcpt_A, C4_ASU(col).Coord(1:size(r_CA_Rcpt_A,1),:));
            [Rt_2,Tr_2]=CoordiExam(r_CA_Rcpt_A, C4_ASU(col).Coord(size(r_CA_Rcpt_A,1)+1:end,:));
            Remained_Part_SiXt=[r_CA_Rcpt_B*Rt_1+Tr_1; r_CA_Rcpt_B*Rt_2+Tr_2];
            SiXt=[C4_ASU(col).Coord; Remained_Part_SiXt];

            % Rearrange SiXt/Xtal in Order of [toXYZ; toXYZ_pair; fromXYZ; fromXYZ_pair]
            len_SiXt=size(SiXt,1); len_Xtal=size(Xtal,1);
            if norm(mean(toXYZ)-mean(fromXYZ))>norm(mean(toXYZ)-mean(fromXYZ_pair))
                if col==1
                    % SiXt=[toXYZ; fromXYZ_pair; toXYZ_pair; fromXYZ];
                    f_SiXt=[SiXt(1:len_SiXt/4,:); SiXt(len_SiXt/4*2+1:len_SiXt/4*3,:); SiXt(len_SiXt/4*3+1:len_SiXt/4*4,:); SiXt(len_SiXt/4+1:len_SiXt/4*2,:)];
                    f_Xtal=[Xtal(1:len_Xtal/4,:); Xtal(len_Xtal/4*2+1:len_Xtal/4*3,:); Xtal(len_Xtal/4*3+1:len_Xtal/4*4,:); Xtal(len_Xtal/4+1:len_Xtal/4*2,:)];
                else
                    % SiXt=[fromXYZ_pair; toXYZ; fromXYZ; toXYZ_pair];
                    f_SiXt=[SiXt(len_SiXt/4+1:len_SiXt/4*2,:); SiXt(len_SiXt/4*3+1:len_SiXt/4*4,:); SiXt(len_SiXt/4*2+1:len_SiXt/4*3,:); SiXt(1:len_SiXt/4,:)];
                    f_Xtal=[Xtal(len_Xtal/4+1:len_Xtal/4*2,:); Xtal(len_Xtal/4*3+1:len_Xtal/4*4,:); Xtal(len_Xtal/4*2+1:len_Xtal/4*3,:); Xtal(1:len_Xtal/4,:)];
                end
            else
                if col==1
                    % SiXt=[toXYZ; fromXYZ; toXYZ_pair; fromXYZ_pair];
                    f_SiXt=[SiXt(1:len_SiXt/4,:); SiXt(len_SiXt/4*2+1:len_SiXt/4*3,:); SiXt(len_SiXt/4+1:len_SiXt/4*2,:); SiXt(len_SiXt/4*3+1:len_SiXt/4*4,:)];
                    f_Xtal=[Xtal(1:len_Xtal/4,:); Xtal(len_Xtal/4*2+1:len_Xtal/4*3,:); Xtal(len_Xtal/4+1:len_Xtal/4*2,:); Xtal(len_Xtal/4*3+1:len_Xtal/4*4,:)];
                else
                    % SiXt=[fromXYZ; toXYZ; fromXYZ_pair; toXYZ_pair];
                    f_SiXt=[SiXt(len_SiXt/4+1:len_SiXt/4*2,:); SiXt(len_SiXt/4*3+1:len_SiXt/4*4,:); SiXt(1:len_SiXt/4,:); SiXt(len_SiXt/4*2+1:len_SiXt/4*3,:)];
                    f_Xtal=[Xtal(len_Xtal/4+1:len_Xtal/4*2,:); Xtal(len_Xtal/4*3+1:len_Xtal/4*4,:); Xtal(1:len_Xtal/4,:); Xtal(len_Xtal/4*2+1:len_Xtal/4*3,:)];
                end
            end

            %

            %%%% Output SiXt/Xtal PDB File
            % C2C4/SiXt PDB File
            [R1,T1]=CoordiExam(Rcpt_A, f_SiXt(1:len_SiXt/4,:));
            [R2,T2]=CoordiExam(Rcpt_A, f_SiXt(len_SiXt/4+1:len_SiXt/4*2,:));
            [R3,T3]=CoordiExam(Rcpt_A, f_SiXt(len_SiXt/4*2+1:len_SiXt/4*3,:));
            [R4,T4]=CoordiExam(Rcpt_A, f_SiXt(len_SiXt/4*3+1:len_SiXt/4*4,:));
            SiXt_Model=[Model_Rcpt*R1+T1; Model_Rcpt*R2+T2; Model_Rcpt*R3+T3; Model_Rcpt*R4+T4];

            % SiXt PDB File
            for n=1:size(SiXt_Model,1)
                Output_SiXt.Model.Atom(n).AtomSerNo=n;
                Output_SiXt.Model.Atom(n).X=SiXt_Model(n,1);
                Output_SiXt.Model.Atom(n).Y=SiXt_Model(n,2);
                Output_SiXt.Model.Atom(n).Z=SiXt_Model(n,3);
                if n<size(SiXt_Model,1)/4+1
                    Output_SiXt.Model.Atom(n).chainID='A';
                elseif n>size(SiXt_Model,1)/4 && n<size(SiXt_Model,1)/4*2+1
                    Output_SiXt.Model.Atom(n).chainID='B';
                elseif n>size(SiXt_Model,1)/4*2 && n<size(SiXt_Model,1)/4*3+1
                    Output_SiXt.Model.Atom(n).chainID='C';
                else
                    Output_SiXt.Model.Atom(n).chainID='D';
                end
            end

            warning('off','all');
            pdbwrite(['3. TopRank_300_Pose/', Result(r).ID, '/', HMD_ID, set(3:7), num2str(Result_Idx), '.pdb'], Output_SiXt);
            
            %

            % C2C4 PDB File
            SiXt_CA_Model=[Rcpt_A*R1+T1; Rcpt_A*R2+T2; Rcpt_A*R3+T3; Rcpt_A*R4+T4];
            C2C4=struct('Coord',[]);
            C2C4(1).Coord=[Rcpt_A; Rcpt_B; Lgnd_A; Lgnd_B]; C2C4(2).Coord=[Rcpt_A; Rcpt_B; Lgnd_B; Lgnd_A];
            C2C4(3).Coord=[Rcpt_B; Rcpt_A; Lgnd_A; Lgnd_B]; C2C4(4).Coord=[Rcpt_B; Rcpt_A; Lgnd_B; Lgnd_A];
            C2C4(5).Coord=[Lgnd_A; Lgnd_B; Rcpt_A; Rcpt_B]; C2C4(6).Coord=[Lgnd_B; Lgnd_A; Rcpt_A; Rcpt_B];
            C2C4(7).Coord=[Lgnd_A; Lgnd_B; Rcpt_B; Rcpt_A]; C2C4(8).Coord=[Lgnd_B; Lgnd_A; Rcpt_B; Rcpt_A];

            C2C4_RMSD_Table=zeros(8,1);
            for C2C4_Idx=1:8
                [~,~,C2C4_RMSD]=CoordiExam(SiXt_CA_Model, C2C4(C2C4_Idx).Coord);
                C2C4_RMSD_Table(C2C4_Idx,1)=C2C4_RMSD;
            end
            [~,C2C4_ridx]=min(C2C4_RMSD_Table);

            C2C4_CA_Model=C2C4(C2C4_ridx).Coord;
            len_C2C4=length(C2C4_CA_Model);
            [Rc1,Tc1]=CoordiExam(Rcpt_A, C2C4_CA_Model(1:len_C2C4/4,:));
            [Rc2,Tc2]=CoordiExam(Rcpt_A, C2C4_CA_Model(len_C2C4/4+1:len_C2C4/4*2,:));
            [Rc3,Tc3]=CoordiExam(Rcpt_A, C2C4_CA_Model(len_C2C4/4*2+1:len_C2C4/4*3,:));
            [Rc4,Tc4]=CoordiExam(Rcpt_A, C2C4_CA_Model(len_C2C4/4*3+1:len_C2C4/4*4,:));
            C2C4_Model=[Model_Rcpt*Rc1+Tc1; Model_Rcpt*Rc2+Tc2; Model_Rcpt*Rc3+Tc3; Model_Rcpt*Rc4+Tc4];

            for w=1:size(C2C4_Model,1)
                Output_C2C4.Model.Atom(w).AtomSerNo=w;
                Output_C2C4.Model.Atom(w).X=C2C4_Model(w,1);
                Output_C2C4.Model.Atom(w).Y=C2C4_Model(w,2);
                Output_C2C4.Model.Atom(w).Z=C2C4_Model(w,3);
                if w<size(C2C4_Model,1)/4+1
                    Output_C2C4.Model.Atom(w).chainID='A';
                elseif w>size(C2C4_Model,1)/4 && w<size(C2C4_Model,1)/4*2+1
                    Output_C2C4.Model.Atom(w).chainID='B';
                elseif w>size(C2C4_Model,1)/4*2 && w<size(C2C4_Model,1)/4*3+1
                    Output_C2C4.Model.Atom(w).chainID='C';
                else
                    Output_C2C4.Model.Atom(w).chainID='D';
                end
            end

            warning('off','all');
            pdbwrite(['3. TopRank_300_Pose/', Result(r).ID, '/', HMD_ID, set(3:7), num2str(Result_Idx), '_C2C4.pdb'], Output_C2C4);
            
            %

            % Xtal PDB File
            [Rt1,Tn1]=CoordiExam(Asy_Unit, f_Xtal(1:len_Xtal/4,:));
            [Rt2,Tn2]=CoordiExam(Asy_Unit, f_Xtal(len_Xtal/4+1:len_Xtal/4*2,:));
            [Rt3,Tn3]=CoordiExam(Asy_Unit, f_Xtal(len_Xtal/4*2+1:len_Xtal/4*3,:));
            [Rt4,Tn4]=CoordiExam(Asy_Unit, f_Xtal(len_Xtal/4*3+1:len_Xtal/4*4,:));
            SiXt_Ans=[PDB_Model*Rt1+Tn1; PDB_Model*Rt2+Tn2; PDB_Model*Rt3+Tn3; PDB_Model*Rt4+Tn4];

            for m=1:size(SiXt_Ans,1)
                Output_Xtal.Model.Atom(m).AtomSerNo=m;
                Output_Xtal.Model.Atom(m).X=SiXt_Ans(m,1);
                Output_Xtal.Model.Atom(m).Y=SiXt_Ans(m,2);
                Output_Xtal.Model.Atom(m).Z=SiXt_Ans(m,3);
                if m<size(SiXt_Ans,1)/4+1
                    Output_Xtal.Model.Atom(m).chainID='A';
                elseif m>size(SiXt_Ans,1)/4 && m<size(SiXt_Ans,1)/4*2+1
                    Output_Xtal.Model.Atom(m).chainID='B';
                elseif m>size(SiXt_Ans,1)/4*2 && m<size(SiXt_Ans,1)/4*3+1
                    Output_Xtal.Model.Atom(m).chainID='C';
                else
                    Output_Xtal.Model.Atom(m).chainID='D';
                end
            end

            warning('off','all');
            pdbwrite(['3. TopRank_300_Pose/', Result(r).ID, '/', HMD_ID, set(3:7), num2str(Result_Idx), '_Ans.pdb'], Output_Xtal);
        end
    end
    r
end

% Shut Down Idle Cores
delete(gcp('nocreate'));

%
%
%

% Load TopRank Result
load('TopRank_Result.mat');
load('C2_TopRank_Result.mat')

% Analysis
t=0;
Res_Cutoff=12;
Sum_Table(1).ID=[];
Sum_Table(1).C4_RMSD95=[];
Sum_Table(1).Overall_RMSD95=[];

for r=1:size(Result,2)
    Total_Pose_No=size(Result(r).Linear_Table,1);
    if Total_Pose_No>300, Total_Pose_No=300; end
    Qualified_Idx=find(Result(r).Linear_Table(1:Total_Pose_No,end)<Res_Cutoff);

    if ~isempty(Qualified_Idx)
        for QI=1:size(Qualified_Idx,1)
            rp=Qualified_Idx(QI,1);
            t=t+1;
            Sum_Table(t).ID=[Result(r).ID, '.000.', num2str(C2_TopRank_Result(r).Output_Idx(Result(r).Linear_Table(rp,2))),'.000.',num2str(Result(r).Linear_Table(rp,4))];
            Sum_Table(t).C4_RMSD95=Result(r).Linear_Table(rp,5);
            Sum_Table(t).Overall_RMSD95=Result(r).Linear_Table(rp,6);
        end
    end
    r
end

save('Sum_Table.mat', 'Sum_Table');

% The Only Exception
% I 41/a Space Group

%%%% Default Settings
% Residue Table
Res={'GLY';'ALA';'VAL';'ILE';'LEU';
    'SER';'THR';'ASP';'ASN';'GLU';
    'GLN';'LYS';'ARG';'CYS';'MET';
    'PHE';'TYR';'TRP';'HIS';'PRO'};

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

% Distance Cutoff (A) of Interacting Residues
Threshold=10;

% Establish Deposited Folder
mkdir('Reduced_Sample_Int_Res_Results');

%

%%%% Cull Coordinates of Interacting Residues
% List of Sample ID
% Sample=readtable('Reduced_Sample_Training_Test_Dataset.xlsx');
Sample=readtable('Qualified_Reduced_Sample_3737(Exclude 4YKN).txt');
Scale_Table=table2struct(Sample);
Scale_Table(1).Scale_Factor=[];
Error_ID=[];

parfor t=1:size(Sample,1) % Example: t=1

    try
        %%% Load Corresponding PDB File
        PDB_ID=lower(Sample.ID{t});
        File=fopen(['PDB/', PDB_ID,'.pdb']);

        % Symmetry Operator
        Operator=[];
        while(1)
            line=fgetl(File);
            if isequal(sscanf(line(1:10),'%c'),'REMARK 300'), break, end

            if isequal(sscanf(line(1:10),'%c'),'REMARK 290')
                SMTRY=sscanf(line(25:79),'%f %f %f');
                Operator=[Operator; SMTRY'];
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
        PDB=pdbread(['PDB/', PDB_ID,'.pdb']);
        PDB_Model=[PDB.Model(1).Atom.X; PDB.Model(1).Atom.Y; PDB.Model(1).Atom.Z]';
        
        PDB_resSeq_Idx=[PDB.Model(1).Atom.resSeq]';
        pseudo_PDB_CA_Idx=strcmp({PDB.Model(1).Atom.AtomName}, 'CA')';
        PDB_resSeq_Unique_Idx=unique(PDB_resSeq_Idx(pseudo_PDB_CA_Idx));

        PDB_CA_Idx=zeros(size(PDB_resSeq_Unique_Idx,1),1);
        for RSUidx=1:size(PDB_resSeq_Unique_Idx,1)
            PDB_CA_Idx(RSUidx,1)=find((PDB_resSeq_Idx==PDB_resSeq_Unique_Idx(RSUidx))&pseudo_PDB_CA_Idx, 1, 'first');
        end
        Asy_Unit=PDB_Model(PDB_CA_Idx,:);

        %

        if PDB.Sequence.NumOfResidues<20
            Error_ID=[Error_ID; t];
        elseif size(Asy_Unit,1)<15
            Error_ID=[Error_ID; t];
        else

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

            % Reduced Unit Box
            Up_B=max(Asy_Unit,[],1)+Threshold;
            Low_B=min(Asy_Unit,[],1)-Threshold;
            X=Unit_Box(:,1)>Low_B(1)&Unit_Box(:,1)<Up_B(1);
            Y=Unit_Box(:,2)>Low_B(2)&Unit_Box(:,2)<Up_B(2);
            Z=Unit_Box(:,3)>Low_B(3)&Unit_Box(:,3)<Up_B(3);

            InBox=find(X+Y+Z==3);
            Reduced_Box=Unit_Box(InBox,:);
            Reduced_Neighbor_Box=Reduced_Box(size(Asy_Unit,1)+1:end,:);

            %

            %%% Identify Interacting Residues
            % Index of Interacting Residues
            DISM=pdist2(Asy_Unit, Reduced_Neighbor_Box);
            Asy_Unit_Idx=find(sum(DISM<Threshold,2));       % Row-wise
            Neighbor_Idx=find(sum(DISM<Threshold,1));       % Column-wise

            % Residue Type List of Interacting Residues
            PDB_Res={PDB.Model(1).Atom.resName}';
            Res_type=PDB_Res(PDB_CA_Idx,:);
            Unit_Box_Res_type=repmat(Res_type,5*5*5*size(Operator,1)/3,1);

            Reduced_Box_Res_type=Unit_Box_Res_type(InBox);
            Reduced_Neighbor_Box_Res_type=Reduced_Box_Res_type(size(Asy_Unit,1)+1:end);

            %

            %%% Output
            % Coordinate of Interacting Residues
            Asy_Coord=Asy_Unit(Asy_Unit_Idx,:);
            NB_Coord=Reduced_Neighbor_Box(Neighbor_Idx,:);

            % Residue Type of Interacting Residues
            Asy_Rtype=Res_type(Asy_Unit_Idx);
            NB_Rtype=Reduced_Neighbor_Box_Res_type(Neighbor_Idx);

            % Res_Type_Seperated Structure    try
        % Extract Coordinates of Alpha Carbon (CA)
        PDB_ID=lower(Sample.ID{t});
        PDB=pdbread(['PDB/', PDB_ID,'.pdb']);
        PDB_Model=[PDB.Model(1).Atom.X; PDB.Model(1).Atom.Y; PDB.Model(1).Atom.Z]';

        PDB_CA_Idx=strcmp({PDB.Model(1).Atom.AtomName}, 'CA')';
        PDB_resSeq_Idx=[PDB.Model(1).Atom.resSeq]';

        Res_ID_Idx=zeros(size(PDB_resSeq_Idx,1),1);
        Sub_Idx=sort(unique(PDB_resSeq_Idx));
        for Sidx=1:size(Sub_Idx,1)
            Res_ID_Idx(PDB_resSeq_Idx==Sub_Idx(Sidx))=Sidx;
        end

        % Load Protein Molecular Surface of Corresponding PDB ID
        Ply_file = fopen(['Surf/',PDB_ID,'/',PDB_ID, '.ply']);
        while(1)
            line=fgetl(Ply_file);
            if length(line)>14 && strcmp(line(1:14),'element vertex')
                Vertex_No=str2double(line(16:end));
            elseif length(line)>12 && strcmp(line(1:12),'element face')
                Face_No=str2double(line(14:end));
            elseif length(line)==10 && strcmp(line(1:10),'end_header')
                RawText = fread(Ply_file,inf,'*char');
                splitLines = cell2mat(textscan(RawText, '%f %f %f %f %*[^\n]'));
                break;
            else
                continue;
            end
        end
        fclose(Ply_file);
        Vertex_Table=splitLines(1:Vertex_No,1:3);

        % Extract Surface Residues (Dedault Threshold: 1.8A)
        DISM=pdist2(PDB_Model, Vertex_Table);
        Surf_Atom=find(min(DISM, [], 2)<1.8);
        Surf_Res=unique(Res_ID_Idx(Surf_Atom));
        ref_Surf_Res=intersect(Surf_Res, unique(Res_ID_Idx(PDB_CA_Idx)));

        Surf_Res_Sym_Idx=zeros(size(ref_Surf_Res,1),1);
        for RSRidx=1:size(ref_Surf_Res,1)
            Surf_Res_Sym_Idx(RSRidx,1)=find((Res_ID_Idx==ref_Surf_Res(RSRidx))&PDB_CA_Idx, 1, 'first');
        end
        PDB_Res={PDB.Model(1).Atom.resName}';
        Surf_Res_Sym=PDB_Res(Surf_Res_Sym_Idx,:);

        Surf_Res_Idx=zeros(size(ref_Surf_Res,1),1);
        for r=1:20
            Surf_Res_Idx=Surf_Res_Idx+strcmp(Surf_Res_Sym,Res{r})*r;
        end
        Surf_Res_Info=[PDB_Model(Surf_Res_Sym_Idx,:), Surf_Res_Idx];

        % non_Int_Surf_Res Coordiinates
        Int_File=table2array(readtable(['Reduced_Sample_Int_Res_Results/', PDB_ID,'.txt']));
        Int_Res_Info=Int_File(Int_File(:,5)==1,1:4);

        % Output All Surf Res
        % Union_Surf_Res_Info=union(Surf_Res_Info, Int_Res_Info, 'rows');     
        % Surf_Output=[];
        % Surf_Output.X=Union_Surf_Res_Info(:,1);
        % Surf_Output.Y=Union_Surf_Res_Info(:,2);
        % Surf_Output.Z=Union_Surf_Res_Info(:,3);
        % Surf_Output.Res=Union_Surf_Res_Info(:,4);
        % writetable(struct2table(Surf_Output),['Reduced_Sample_Surf_Res_Results/', PDB_ID,'.txt']);

        % Output n_Int Surf Res
        n_Int_Res_Info=setdiff([Surf_Res_Info; Int_Res_Info], Int_Res_Info, 'rows');

        n_Int_Output=[];
        n_Int_Output.X=n_Int_Res_Info(:,1);
        n_Int_Output.Y=n_Int_Res_Info(:,2);
        n_Int_Output.Z=n_Int_Res_Info(:,3);
        n_Int_Output.Res=n_Int_Res_Info(:,4);
        writetable(struct2table(n_Int_Output),['Reduced_Sample_nInt_Res_Results/', PDB_ID,'.txt']);

    catch
        Error_ID=[Error_ID; t];
    end
    t
            Coord=[Asy_Coord; NB_Coord];
            Mo_Idx=[ones(size(Asy_Co    try
        % Extract Coordinates of Alpha Carbon (CA)
        PDB_ID=lower(Sample.ID{t});
        PDB=pdbread(['PDB/', PDB_ID,'.pdb']);
        PDB_Model=[PDB.Model(1).Atom.X; PDB.Model(1).Atom.Y; PDB.Model(1).Atom.Z]';

        PDB_CA_Idx=strcmp({PDB.Model(1).Atom.AtomName}, 'CA')';
        PDB_resSeq_Idx=[PDB.Model(1).Atom.resSeq]';

        Res_ID_Idx=zeros(size(PDB_resSeq_Idx,1),1);
        Sub_Idx=sort(unique(PDB_resSeq_Idx));
        for Sidx=1:size(Sub_Idx,1)
            Res_ID_Idx(PDB_resSeq_Idx==Sub_Idx(Sidx))=Sidx;
        end

        % Load Protein Molecular Surface of Corresponding PDB ID
        Ply_file = fopen(['Surf/',PDB_ID,'/',PDB_ID, '.ply']);
        while(1)
            line=fgetl(Ply_file);
            if length(line)>14 && strcmp(line(1:14),'element vertex')
                Vertex_No=str2double(line(16:end));
            elseif length(line)>12 && strcmp(line(1:12),'element face')
                Face_No=str2double(line(14:end));
            elseif length(line)==10 && strcmp(line(1:10),'end_header')
                RawText = fread(Ply_file,inf,'*char');
                splitLines = cell2mat(textscan(RawText, '%f %f %f %f %*[^\n]'));
                break;
            else
                continue;
            end
        end
        fclose(Ply_file);
        Vertex_Table=splitLines(1:Vertex_No,1:3);

        % Extract Surface Residues (Dedault Threshold: 1.8A)
        DISM=pdist2(PDB_Model, Vertex_Table);
        Surf_Atom=find(min(DISM, [], 2)<1.8);
        Surf_Res=unique(Res_ID_Idx(Surf_Atom));
        ref_Surf_Res=intersect(Surf_Res, unique(Res_ID_Idx(PDB_CA_Idx)));

        Surf_Res_Sym_Idx=zeros(size(ref_Surf_Res,1),1);
        for RSRidx=1:size(ref_Surf_Res,1)
            Surf_Res_Sym_Idx(RSRidx,1)=find((Res_ID_Idx==ref_Surf_Res(RSRidx))&PDB_CA_Idx, 1, 'first');
        end
        PDB_Res={PDB.Model(1).Atom.resName}';
        Surf_Res_Sym=PDB_Res(Surf_Res_Sym_Idx,:);

        Surf_Res_Idx=zeros(size(ref_Surf_Res,1),1);
        for r=1:20
            Surf_Res_Idx=Surf_Res_Idx+strcmp(Surf_Res_Sym,Res{r})*r;
        end
        Surf_Res_Info=[PDB_Model(Surf_Res_Sym_Idx,:), Surf_Res_Idx];

        % non_Int_Surf_Res Coordiinates
        Int_File=table2array(readtable(['Reduced_Sample_Int_Res_Results/', PDB_ID,'.txt']));
        Int_Res_Info=Int_File(Int_File(:,5)==1,1:4);

        % Output All Surf Res
        % Union_Surf_Res_Info=union(Surf_Res_Info, Int_Res_Info, 'rows');     
        % Surf_Output=[];
        % Surf_Output.X=Union_Surf_Res_Info(:,1);
        % Surf_Output.Y=Union_Surf_Res_Info(:,2);
        % Surf_Output.Z=Union_Surf_Res_Info(:,3);
        % Surf_Output.Res=Union_Surf_Res_Info(:,4);
        % writetable(struct2table(Surf_Output),['Reduced_Sample_Surf_Res_Results/', PDB_ID,'.txt']);

        % Output n_Int Surf Res
        n_Int_Res_Info=setdiff([Surf_Res_Info; Int_Res_Info], Int_Res_Info, 'rows');

        n_Int_Output=[];
        n_Int_Output.X=n_Int_Res_Info(:,1);
        n_Int_Output.Y=n_Int_Res_Info(:,2);
        n_Int_Output.Z=n_Int_Res_Info(:,3);
        n_Int_Output.Res=n_Int_Res_Info(:,4);
        writetable(struct2table(n_Int_Output),['Reduced_Sample_nInt_Res_Results/', PDB_ID,'.txt']);

    catch
        Error_ID=[Error_ID; t];
    end
    tord,1),1); zeros(size(NB_Coord,1),1)];
            Rtype=[Asy_Rtype; NB_Rtype];
            Rtype_Idx=zeros(size(Rtype));

            for r=1:20
                Rtype_Idx=Rtype_Idx+strcmp(Rtype,Res{r})*r;
            end

            % Output
            Output=[];
            Output.X=Coord(:,1);
            Output.Y=Coord(:,2);
            Output.Z=Coord(:,3);
            Output.Res=Rtype_Idx;
            Output.Mo=Mo_Idx;
            writetable(struct2table(Output),['Reduced_Sample_Int_Res_Results/',lower(Sample.ID{t}),'.txt']);

            % Scale Record
            Asy_Unit_Scale=pdist2(mean(Asy_Unit,1),Asy_Unit);
            Asy_Coord_Scale=pdist2(mean(Asy_Coord,1),Asy_Coord);
            NB_Coord_Scale=pdist2(mean(NB_Coord,1),NB_Coord);
            Scale_Table(t).Scale_Factor=[max(Asy_Unit_Scale), max(Asy_Coord_Scale), max(NB_Coord_Scale)];
        end
    catch
        Error_ID=[Error_ID; t];
    end
    t
end

writetable(struct2table(Scale_Table),'Scale_Table.xlsx');
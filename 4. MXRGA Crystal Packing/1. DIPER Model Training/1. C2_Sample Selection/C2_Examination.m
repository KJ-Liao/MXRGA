%%%% C2 Homodimer Examination
% Load PDB ID
List=fopen('Reduced_Sample_Merged.txt');
Sample(1).ID=[];s=1;
while (1)
    line=fgetl(List);
    if line==-1, break, end
    Sample(s).ID=sscanf(line, '%c');
    s=s+1;
end
fclose(List);

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

% Parellel Computation (default core:10)
parpool(10)

Sample(1).C2=[];
for r=1:size(Sample,2)

    %%%% Cull Crystal Homodimers for RMSD Examination
    % Load PDB File
    try
        PDB=pdbread(['/home/juju/Desktop/SingleChain_PDB_2022-05-27/',lower(Sample(r).ID),'.pdb']);
        File=fopen(['/home/juju/Desktop/SingleChain_PDB_2022-05-27/',lower(Sample(r).ID),'.pdb']);

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

        %

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
        j=0; PDB_NB=struct('Coord', 1);
        for k=2:size(Unit_Box,1)/size(Asy_Unit,1)
            DISM=pdist2(Asy_Unit, Unit_Box(size(Asy_Unit,1)*(k-1)+1:size(Asy_Unit,1)*k,:));
            if min(min(DISM))<NB_Cutoff
                j=j+1;
                PDB_NB(j).Coord=Unit_Box(size(Asy_Unit,1)*(k-1)+1:size(Asy_Unit,1)*k,:);
            end
        end

        if j~=0
            NeiB_Idx=[];
            RMSD_Table=zeros(size(PDB_NB,2));
            for m=1:size(PDB_NB,2)
                for n=1:size(PDB_NB,2)
                    [~,~,eRMSD_1]=CoordiExam([Asy_Unit; PDB_NB(m).Coord], [Asy_Unit; PDB_NB(n).Coord]);
                    [~,~,eRMSD_2]=CoordiExam([Asy_Unit; PDB_NB(m).Coord], [PDB_NB(n).Coord; Asy_Unit]);
                    RMSD_Table(m,n)=min(eRMSD_1, eRMSD_2);
                end
                NeiB_Idx=[NeiB_Idx, min(find(RMSD_Table(m,:)<1))];
            end

            % Generate Closest/Unique Packing Neighbors (PDB_UNBs)
            p=0; UNB=unique(NeiB_Idx);
            for q=1:length(UNB)
                R=CoordiExam(Asy_Unit, PDB_NB(UNB(q)).Coord);
                Eigen=sortrows([real(eig(R)), imag(eig(R))])*[1;sqrt(-1)];
                % Examine C2 Symmetry
                if norm(Eigen-[-1; -1; 1])<6*pi/180
                    [V,D]=eig(R);
                    [eig_v,eig_idx]=max(real(diag(D)'));
                    if abs(acosd((mean(Asy_Unit)-mean(PDB_NB(UNB(q)).Coord))/norm(mean(Asy_Unit)-mean(PDB_NB(UNB(q)).Coord))*V(:,eig_idx))-90)<5
                        p=p+1;
                    end
                end
            end

            % Label
            if p>0
                Sample(r).C2=1;
            else
                Sample(r).C2=0;
            end
        end
    catch
        Sample(r).C2=2;
    end
    r
end

save('C2_Examination_Result.mat', 'Sample');
writetable(struct2table(Sample),'C2_Examination_Result.xlsx');

%

%%%% Output Set of Resulted PDB
% Load C2 Homodimer PDB ID
List=fopen('C2_Sample.txt');
Sample(1).ID=[];s=1;
while (1)
    line=fgetl(List);
    if line==-1, break, end
    Sample(s).ID=sscanf(line, '%c');
    s=s+1;
end
fclose(List);

% Copyfile to Destiniation (Result\C2_Sample)
for r=1:size(Sample,2)
    File=[lower(Sample(r).ID), '.pdb'];
    status = copyfile(File, 'Result\C2_Sample');
end

%

%%%% PIPER Script
% Create the Deposited File
File=fopen('C2_Sample.txt');
s=1; Sample_Info(1).ID=[];
while (1)
    line=fgetl(File);
    if line==-1, break, end 
    Sample_Info(s).ID=sscanf(line,'%c');
    mkdir(lower(Sample_Info(s).ID(1:4)));
    s=s+1;         
end
fclose(File);

%

% Sample Preparation: Sample_prep.txt
Script_Prep = fopen([pwd,'/Sample_prep.txt'],'w');
fprintf(Script_Prep,'#!/bin/bash');
fprintf(Script_Prep,'\n');
for r=1:size(Sample_Info,2)
    fprintf(Script_Prep,'\n');
    fprintf(Script_Prep,['../protein_prep/prepare.py',' ',lower(Sample_Info(r).ID(1:4)),'.pdb']);   
end  
fclose(Script_Prep);

%

% Script fot PIPER: Sample_PIPER.txt
% Directory Path of "piper_package"
Path='/home/sunmd33452/Desktop/piper_package/C2_Sample/';

% Write Sample_PIPER.txt
Script_PIPER = fopen([pwd,'/Sample_PIPER.txt'],'w');
fprintf(Script_PIPER,'#!/bin/bash');
fprintf(Script_PIPER,'\n');
for r=1:size(Sample_Info,2)
    fprintf(Script_PIPER,'\n');
    fprintf(Script_PIPER,['../piper -vv -c1.0 -k4 --msur_k=1.0 --maskr=1.0 -T FFTW_EXHAUSTIVE -R 3281 -t 1 -p ../prms/atoms.prm -f ../prms/coeffs_DIPERX.prm -r ../prms/C2_rots.prm',' ',lower(Sample_Info(r).ID(1:4)),'_pnon.pdb',' ',lower(Sample_Info(r).ID(1:4)),'_pnon.pdb',' ','--o ', Path, lower(Sample_Info(r).ID(1:4)),'/']);
end  
fclose(Script_PIPER);

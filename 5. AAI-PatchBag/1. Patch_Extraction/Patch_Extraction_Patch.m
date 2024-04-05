%%%% Default Setting
% Parallel Computation
parpool(6)


% Establish Deposited Folder
mkdir('Reduced_Sample_Int_Res_Results/Patch_Results_6plus1');

% Patch Size (No of Alpha Carbons)
Patch_Size=6;

% Pre-allocate RAM for Error ID & Patch Struct
Error_ID=[]; fields = {'Coord','Norm_V','Res_Type'};

% Load PDB ID of Qualified Reduced Samples
Filename='Reduced_Sample_Training_Test_Dataset.xlsx';
Sample_Info=table2struct(readtable(Filename));

% Patch Extraction
parfor i=1:size(Sample_Info,1) % Example: i=1
    try
        % Load Info of Interacting Residues
        PDB_ID=lower(Sample_Info(i).ID);
        Int_File=table2array(readtable(['Reduced_Sample_Int_Res_Results/', PDB_ID,'.txt']));
        Int_Asy_Unit=Int_File(Int_File(:,5)==1,1:4);

        if size(Int_Asy_Unit,1)<6
            Patch = cell2struct(cell(length(fields), []), fields);
            writetable(struct2table(Patch),['Reduced_Sample_Int_Res_Results/Patch_Results_6plus1/',PDB_ID, '_Patch.txt']);
        else
            % Read Molecule Surface
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

            % Normal Vector Calculation
            Vertex_Table=splitLines(1:Vertex_No,1:3);
            Residue_Surface_DISM=pdist2(Int_Asy_Unit(:,1:3), Vertex_Table);
            for k=1:size(Residue_Surface_DISM,1)
                [~, Nidx]=mink(Residue_Surface_DISM(k,:),10);
                V=Int_Asy_Unit(k,1:3)-mean(Vertex_Table(Nidx,:),1);
                Int_Asy_Unit(k,5:7)=V/norm(V);
            end

            %
            %
            %

            % Patch Extraction
            DISM=squareform(pdist(Int_Asy_Unit));
            Patch_Ridx=zeros(size(Int_Asy_Unit,1),Patch_Size);
            for p=1:size(Int_Asy_Unit,1)
                [~, Pidx]=mink(DISM(p,:), Patch_Size);
                Patch_Ridx(p,:)=sort(Pidx);
            end

            %%%% NEED TO BE IMPROVED %%%%
            Patch = cell2struct(cell(length(fields), []), fields);
            n_Int_File=Int_File(Int_File(:,5)==0,1:3);
            for j=1:size(Patch_Ridx,1)
                Nei_Dist=pdist2(Int_Asy_Unit(j,1:3), n_Int_File);
                [~, Nei_ID]=min(Nei_Dist);
                Patch(j).Coord=[Int_Asy_Unit(Patch_Ridx(j,:),1:3); n_Int_File(Nei_ID, :)];

                B=sum(Int_Asy_Unit(Patch_Ridx(j,:),5:7));
                Patch(j).Norm_V=B/norm(B);
                Patch(j).Res_Type=[Int_Asy_Unit(Patch_Ridx(j,:),4); Int_File(Nei_ID, 4)];
            end

            % Output
            writetable(struct2table(Patch),['Reduced_Sample_Int_Res_Results/Patch_Results_6plus1/',PDB_ID, '_Patch.txt']);
        end
    catch
        Error_ID=[Error_ID; i];
    end
    i
end

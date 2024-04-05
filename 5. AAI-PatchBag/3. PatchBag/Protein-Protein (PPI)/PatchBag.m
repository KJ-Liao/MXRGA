%%%% Patch to Vector
% Load Qualified Reduced Sample ID
Filename='Reduced_Sample_Training_Test_Dataset.xlsx';
Sample_Info=table2struct(readtable(Filename));

% Load Merged Error ID
load('Merged_Error_ID.mat')
Sample_Info=Sample_Info(setdiff([1:2472], Merged_Error_ID));

% Load Clust Info (Representative Patch Library)
Max_Cluster_No=300;
load(['Clust_medoids_', num2str(Max_Cluster_No), '_PPI.mat'])

% Patch Classification
Selected_Sample=Sample_Info;
Selected_Sample(1).Patch_code=[];
Selected_Sample(1).Patch_Res=[];
Selected_Sample(1).Bar_code=[];

Blank_ID=[];
P=perms(1:6);
for i=1:size(Selected_Sample,1)
    try
        % Patch Files
        PDB_ID=lower(Selected_Sample(i).ID);
        Patch_File=table2array(readtable([PDB_ID,'_Patch.txt']));
        [Unique_Patch, ~, Uidx]=unique(Patch_File, 'rows', 'stable');

        % Patch Codes
        Patch_code=zeros(size(Unique_Patch,1),1);

        parfor j=1:size(Unique_Patch,1)
            P1=reshape(Unique_Patch(j,1:21),7,3);
            RMSD_DISM=zeros(Max_Cluster_No, 1);

            for m=1:Max_Cluster_No
                RMSD_Table=zeros(size(P,1), 1);
                Angle_Table=zeros(size(P,1), 1);

                for r=1:size(P,1)
                    P2=Clust_medoids(m).Coord([P(r,:), 7],:);

                    [R,~,eRMSD] = CoordiExam_AC(P1, P2);
                    RMSD_Table(r,1)=eRMSD;
                    Angle_Table(r,1)=acosd(dot(Unique_Patch(j,22:24)*R, Clust_medoids(m).Norm_V));
                end

                if isempty(RMSD_Table(Angle_Table<90))
                    RMSD_DISM(m,1)=min(RMSD_Table,[],'all');
                else
                    RMSD_DISM(m,1)=min(RMSD_Table(Angle_Table<90),[],'all');
                end
            end
            [~,c]=min(RMSD_DISM);
            Patch_code(j,1)=c;
        end

        Unique_Count=hist(Uidx, (1:max(Uidx))');
        Selected_Sample(i).Patch_code=[Patch_code, Unique_Count'];
        Selected_Sample(i).Patch_Res=Unique_Patch(:,25:31);

        % Patch Vectorlization
        Bar_code=zeros(Max_Cluster_No,1);
        for u=1:Max_Cluster_No
            Bar_code(u,1)= sum(Unique_Count(Patch_code==u));
        end
        Selected_Sample(i).Bar_code=Bar_code;
    catch
        Blank_ID=[Blank_ID; i];
    end
    i
end

save('Blank_ID.mat', 'Blank_ID');
save('Selected_Sample.mat', 'Selected_Sample');

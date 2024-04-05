%%%% Patch Preprocessing
% Load Patch Merge
load('Unique_Patch_Merge_PSI.mat');

% Unique Res_Type Patches
Res_Type=[Unique_Patch_Merge.Res_Type];
Reshpaed_Res_Type=reshape(Res_Type,6,[])';
[~, Idx]=unique(Reshpaed_Res_Type, 'rows');
Unique_Res_Type_Patch=Unique_Patch_Merge(Idx);

Reshpaed_Unique_Res_Type=reshape([Unique_Res_Type_Patch.Res_Type],6,[])';
Unique_NonX_Res_Type_Patch=Unique_Res_Type_Patch(find(~sum(Reshpaed_Unique_Res_Type==0,2)));

%%%% % Patch Library
% Library Size
Count=8000;

% Patch Library
Seed=RandStream('mrg32k3a');
Idx=randsample(Seed, size(Unique_NonX_Res_Type_Patch,2), Count);

% Reshpaed_Norm_V=reshape(Norm_V,3,[])';
% plot3(Reshpaed_Norm_V(Idx,1), Reshpaed_Norm_V(Idx,2), Reshpaed_Norm_V(Idx,3), '.')
Patch_Library_8000=Unique_NonX_Res_Type_Patch(Idx);
clearvars -except Count Patch_Library_8000

% Distance Matrix (DISM) of Patch Library
% parpool (8)
P=perms(1:size(Patch_Library_8000(1).Coord,1)-1);
nm_RMSD_DISM=[];
Error_ID=[];

parfor n=1:(Count-1)
    P1=Patch_Library_8000(n).Coord;
    m_RMSD_DISM=zeros(Count-n, 1);

    for m=n+1:Count
        RMSD_Table=zeros(size(P,1), 1);
        Angle_Table=zeros(size(P,1), 1);

        for r=1:size(P,1)
            P2=Patch_Library_8000(m).Coord([P(r,:)],:);

            [R,~,eRMSD] = CoordiExam_AC(P1, P2);
            RMSD_Table(r,1)=eRMSD;
            Angle_Table(r,1)=acosd(dot(Patch_Library_8000(n).Norm_V*R, Patch_Library_8000(m).Norm_V));
        end

        if isempty(RMSD_Table(Angle_Table<90))
            Error_ID=[Error_ID; [n,m]];
            m_RMSD_DISM(m-n,1)=20;
        else
            m_RMSD_DISM(m-n,1)=min(RMSD_Table(Angle_Table<90),[],'all');
        end
    end
    nm_RMSD_DISM=[nm_RMSD_DISM; m_RMSD_DISM];
    n
end

Patch_Library_DISM=squareform(nm_RMSD_DISM);
save('Patch_Library_DISM.mat', 'Patch_Library_DISM');

%

%%%% Clustering of Patch Library
% Load PDB ID of Qualified Reduced  Samples
Filename='Qualified_Reduced_Sample_3738.txt';
Sample_Info=table2struct(readtable(Filename));

% Average Count No.
% mean(Count_Table)
Count_Table=zeros(size(Sample_Info));
for i=1:size(Sample_Info,1)
    PDB_ID=lower(Sample_Info(i).ID);
    File=table2array(readtable([PDB_ID, '_Patch.txt']));
    Count_Table(i)=size(File,1);
end


%%%% K-medoid Clustering
% Load Patch_Library_DISM.mat
load('Patch_Library_DISM.mat')

% Distance matrix
DISM=Patch_Library_DISM;
N = size(DISM,1);

Max_Cluster_No=300;
[idx, C] = kmedoids((1:N)', Max_Cluster_No, 'Distance', @(row, col) DISM(row, col), 'Replicates', 5);

% SH=zeros(Max_Cluster_No,2);
% for K=2:Max_Cluster_No
% idx = kmedoids((1:N)', K, 'Distance', @(row, col) DISM(row, col));
% s=silhouette((1:N)',idx, @(row, col) DISM(row, col));
% SH(K,1)=mean(s);
% SH(K,2)=Dunns(K,DISM,idx);
% K
% end
% [val, psn]=max(SH(:,2));
% [idx, C] = kmedoids((1:N)', psn, 'Distance', @(row, col) DISM(row, col));

% Medoid
load('Patch_Library_8000.mat')
Clust_medoids=Patch_Library_8000(C);

%

% Intra-distance
Clust_Table=zeros(Max_Cluster_No,3);
for t=1:Max_Cluster_No
    Tidx=find(idx==t);
    if size(Tidx,1)==1
        Clust_Table(t,1)=0;
        Clust_Table(t,2)=Tidx;
    else
        Cluster=Patch_Library_DISM(Tidx, Tidx);
        Clust_Table(t,1)=mean(squareform(Cluster));
        [~,c]=min(sum(Cluster,1));
        Clust_Table(t,2)=Tidx(c);
    end
end

% Inter_distance
P = perms(1:6);
Medoids_DISM=[];
parfor n=1:(Max_Cluster_No-1)
    P1=Clust_medoids(n).Coord;

    RMSD_Table=zeros(size(P,1), 1);
    Angle_Table=zeros(size(P,1), 1);
    for m=n+1:Max_Cluster_No
        for r=1:size(P,1)
            P2=Clust_medoids(m).Coord([P(r,:)],:);

            [R,~,eRMSD] = CoordiExam_AC(P1, P2);
            RMSD_Table(r,1)=eRMSD;
            Angle_Table(r,1)=acosd(dot(Clust_medoids(n).Norm_V*R, Clust_medoids(m).Norm_V));
        end
        
        % Default: 90
        Medoids_DISM=[Medoids_DISM, min(RMSD_Table(Angle_Table<120),[],'all')];
    end
    n
end
Medoids_DISM=squareform(Medoids_DISM);

% Clust_Table
Clust_Table(:,3)=sum(Medoids_DISM,1)/(Max_Cluster_No-1);
boxplot([Clust_Table(:,1), Clust_Table(:,3)], 'Labels', {'Intra-dist', 'Inter-dist'})

save('Clust_medoids_300.mat', 'Clust_medoids');
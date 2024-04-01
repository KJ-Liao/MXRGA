%%%% Load Space Group List
Data=readtable('SCH_SampleInfo.xlsx');
SpaceGroup_List=unique(Data.SG);
% SpaceGroup_List(23)={'I 41_a'};
% SpaceGroup_List(24)=[];

%%%% Pairwise sequence alignment matix (DISM) for Rplot
% Unadjusted dissimilarity as distance
for r=1:size(SpaceGroup_List,1)

    % Compute DISM
    Data=table2struct(readtable([SpaceGroup_List{r}, '.xlsx']));
    D=seqpdist({Data(1:size(Data,1)).Seq},'Method', 'p-distance');
    DISM=squareform(D);

    % Linkage Clustering
    Tree=seqlinkage(DISM);
    LeafClusters=cluster(Tree, 0.3);

    % Arrange Data in Tree Order
    LeafOrder=regexp(get(Tree,'LeafNames'), '\d*', 'match');
    LeafOrder=str2double(cat(2,LeafOrder{:})');

    % Select the Highest Resolution Structure as Representative PDB
    Reduced_sample_idx=[];
    Ordered_Data=Data(LeafOrder);
    for n=1:size(Ordered_Data,1)
        Ind=find(LeafClusters==n);
        [~,idx]=min(str2double({Ordered_Data(Ind).Res}'));
        Reduced_sample_idx=[Reduced_sample_idx; Ind(idx)];
    end

    % Output
    Reduced_Sample=Ordered_Data(Reduced_sample_idx);
    writetable(struct2table(Reduced_Sample),['Reduced Sample/', SpaceGroup_List{r}, '_Reduced_Sample.xlsx']);

end

% 56, 58, 61, 66, 67, 68, 69, 70, 71, 72

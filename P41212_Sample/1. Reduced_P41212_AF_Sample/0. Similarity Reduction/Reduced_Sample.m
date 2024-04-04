%%%% Sample reduction via dissimlarity clustering (<0.3)
% Dissimlarity clustering (<0.3)
Table=readtable('DISM.xlsx', 'VariableNamingRule', 'preserve');
DISM=table2array(Table(1:end,2:end));

Tree=seqlinkage(DISM);
LeafClusters=cluster(Tree, 0.3);

% Arrange data in tree order
LeafOrder=regexp(get(Tree,'LeafNames'), '\d*', 'match');
LeafOrder=str2double(cat(2,LeafOrder{:})');

Data=table2struct(readtable('P41212_Sample.xlsx'));
Ordered_Data=Data(LeafOrder);

% Select the lowest resolution structure as representative PDB
Reduced_sample_idx=[];
for n=1:size(Ordered_Data,1)
    Ind=find(LeafClusters==n);
    [~,idx]=min(str2double({Ordered_Data(Ind).Res}'));
    Reduced_sample_idx=[Reduced_sample_idx; Ind(idx)];
end

Reduced_Sample=Ordered_Data(Reduced_sample_idx);

writetable(struct2table(Reduced_Sample),'Reduced_Sample.xlsx');

%%%% Space Group statistics
Data=readtable('Accessible_SCH_SampleInfo.xlsx');
SpaceGroup_List=unique(Data.SG);

SpaceGroup(1).ID=[];
SpaceGroup(1).No=[];
for n=1:size(SpaceGroup_List,1)
    Index=ismember(Data.SG,SpaceGroup_List{n});
    SpaceGroup(n).ID=SpaceGroup_List{n};
    SpaceGroup(n).No=sum(Index);

    Sample=Data(ismember(Data.SG, SpaceGroup_List{n}),:);
    writetable(Sample, [SpaceGroup_List{n}, '.xlsx'])
end

writetable(struct2table(SpaceGroup), 'Space Group_statistics.xlsx')

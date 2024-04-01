%%%% Load Space Group List
Data=readtable('Accessible_SCH_SampleInfo.xlsx');
SpaceGroup_List=unique(Data.SG);
SpaceGroup_List(23)={'I 41_a'};
SpaceGroup_List(24)=[];

%%%% Record Reduced Unique Sample Info
Reduced_Sample(1).ID=[];
Reduced_Sample(1).Seq=[];
Reduced_Sample(1).SG=[];
m=1;

for r=1:size(SpaceGroup_List,1)

    % Record Reduced Unique Sample Info of ID & Seq
    Data=table2struct(readtable([SpaceGroup_List{r}, '_Reduced_Sample.xlsx']));
    for n=1:size(Data,1)
        Reduced_Sample(m).ID=Data(n).ID;
        Reduced_Sample(m).Seq=Data(n).Seq;
        Reduced_Sample(m).SG=Data(n).SG;
        m=m+1;
    end
    
end

% Output
writetable(struct2table(Reduced_Sample),'Reduced_Sample_Merged.xlsx');
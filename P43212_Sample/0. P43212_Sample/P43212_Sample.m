Data=readtable('Accessible_SCH_SampleInfo.xlsx');
P43212_Sample=Data(ismember(Data.SG,'P 43 21 2'),:);
writetable(P43212_Sample, 'P43212_Sample.xlsx')
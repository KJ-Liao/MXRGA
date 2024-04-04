Data=readtable('Accessible_SCH_SampleInfo.xlsx');
P41212_Sample=Data(ismember(Data.SG,'P 41 21 2'),:);
writetable(P41212_Sample, 'P41212_Sample.xlsx')
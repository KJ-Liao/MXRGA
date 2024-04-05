load('AA_Index.mat')
load('Selected_Sample.mat')

% Selected_Sample.Bar_code -> Spatial Distance
Count=[];
Merged_Table=[];
Selected_Sample=Selected_Sample';

parfor i=1:size(Selected_Sample,2)
    Patch_Res=Selected_Sample(i).Patch_Res;
    if ~isempty(Patch_Res)
        for j=1:size(Patch_Res,1)
            Idx=Patch_Res(j,:);
            Idx(Idx==0)=8;
            Merged_Table=[Merged_Table, mean(AA_Index(:,Idx),2)];
        end
        Count=[Count; Selected_Sample(i).Patch_code(:,2)];
    end
    i
end

Completed_Table=[];
parfor k=1:size(Count,1)
    Completed_Table=[Completed_Table, repmat(Merged_Table(:,k),1,Count(k))];
    k
end
% histogram(Completed_Table(1,:))

Ruler=[quantile(Completed_Table,0.2,2), quantile(Completed_Table,0.4,2), quantile(Completed_Table,0.6,2), quantile(Completed_Table,0.8,2)];
% MEAN=mean(Completed_Table,2);
% STD=std(Completed_Table,[],2);
% Ruler=[MEAN-3*STD, MEAN-2.5*STD, MEAN-2*STD, MEAN-1.5*STD, MEAN-1*STD, MEAN-0.5*STD, MEAN, MEAN+0.5*STD, MEAN+1*STD, MEAN+1.5*STD, MEAN+2*STD, MEAN+2.5*STD, MEAN+3*STD];
% Ruler=[MEAN-3*STD, MEAN-2*STD, MEAN-1*STD, MEAN, MEAN+1*STD, MEAN+2*STD, MEAN+3*STD];

Hist=zeros(size(Merged_Table));
for m=1:size(Merged_Table,2)
    Det=Merged_Table(:,m)<Ruler;
    for n=1:size(Merged_Table,1)
        if isempty(find(Det(n,:),1))
            Hist(n,m)=size(Ruler,2)+1;
        else
            Hist(n,m)=find(Det(n,:),1);
        end
    end
    m
end

Selected_Sample(1).Merged_AAIdx=[]; Acc=0;
for p=1:size(Selected_Sample,2)
    if ~isempty(Selected_Sample(p).Patch_code)
        Acc_Start=Acc+1;
        Acc=Acc+size(Selected_Sample(p).Patch_code,1);
        Selected_Sample(p).Merged_AAIdx=[Selected_Sample(p).Patch_code(:,2), Selected_Sample(p).Patch_code(:,1), Hist(:,Acc_Start:Acc)'];
    end
    p
end

Merged_AAIdx=Selected_Sample;
% save('Merged_AAIdx.mat','Merged_AAIdx')

% Patch Vectorlization
Max_Cluster_No=300; Bin=size(Ruler,2)+1;

clear AA_Idx_Patch_code;
AA_Idx_Patch_code(1).ID=[];
AA_Idx_Patch_code(1).Patch_code=[];

for u=1:158
    for v=1:size(Merged_AAIdx,2)
        if ~isempty(Merged_AAIdx(v).Patch_code)
            Data=Merged_AAIdx(v).Merged_AAIdx(:,[1,2,u+2]);
            Patch_code=zeros(Max_Cluster_No, Bin);
            for x=1:size(Data,1)
                Patch_code(Data(x,2), Data(x,3))=Patch_code(Data(x,2), Data(x,3))+Data(x,1);
            end
            AA_Idx_Patch_code(v).ID=Merged_AAIdx(v).ID;
            AA_Idx_Patch_code(v).Patch_code=sum(Patch_code,2);
            AA_Idx_Patch_code(v).(['Patch_code_', num2str(u)])=Patch_code;
        end
    end
    u
end
save('AA_Idx_Patch_code.mat','AA_Idx_Patch_code')

% Features of 2230 Reduced Samples (Including Outliers)
% MW_Patch Num_Int Matrix
load('MPI.mat')

% DISM of 2472 Reduced Samples
load('Seq_DISM_493.mat')
load('RMSD_DISM_493.mat')
load('Jerhoud_DISM_493.mat')
load('Kihara_DISM_493.mat')
load('Condition_DISM_050310_Temp_pH_493.mat')

Seq_Linear_DISM=reshape(Seq_DISM_493',1,493*2393);
RMSD_Linear_DISM=reshape(RMSD_DISM_493',1,493*2393);
Jehoud_Linear_DISM=reshape(Jerhoud_DISM_493',1,493*2393);
Kihara_Linear_DISM=reshape(Kihara_DISM_493',1,493*2393);
Condition_Linear_DISM=reshape(Condition_DISM_493',1,493*2393);

Protein_Sim=[Condition_Linear_DISM;
             Seq_Linear_DISM;
             RMSD_Linear_DISM;
             Jehoud_Linear_DISM;
             Kihara_Linear_DISM];

clearvars -except Protein_Sim MPI

% Combine Features of Reduced Samples (1+4+4 +69*2: 5 bin)
Result_Int=load('AA_Idx_Patch_Code/Int_Result_5bin.mat');
Result_nInt=load('AA_Idx_Patch_Code/nInt_Result_5bin.mat');
Feature=[Protein_Sim; MPI; Result_Int.Result; Result_nInt.Result];

%
%
%

%%%% Remove Outliers
% Load Unique_Problematic_Samples (Outliers)
load(['Outliers/', 'Unique_Problematic_Samples_3IQR.mat'])

%

% Load ID of Qualified_Sample_Condition_Table (RAW)
Data=readtable(['Package/', 'Qualified_Sample_Condition_Table(RAW).xlsx'],'VariableNamingRule','preserve');
Qualified_Sample_ID=table2struct(Data(:,1));

% Load PDB ID of Training/Test 493 Reduced Samples
Sample_Info=table2struct(readtable(['Package/', 'Reduced_Sample_493.txt']));

% Merge the Index of Problematic Samples
UPS_Idx=zeros(size(Unique_Problematic_Samples));
for r=1:size(Unique_Problematic_Samples,1)
    if any(strcmp({Sample_Info.ID}, Qualified_Sample_ID(Unique_Problematic_Samples(r)).ID))
        UPS_Idx(r,1)=find(strcmp({Sample_Info.ID}, Qualified_Sample_ID(Unique_Problematic_Samples(r)).ID));
    end
end
UPS_Idx=UPS_Idx(UPS_Idx~=0);

% Remove Outliers
Merged_UPS_Idx=unique(UPS_Idx);
OL_Idx1=setdiff([1:size(Sample_Info,1)], Merged_UPS_Idx,'stable');

%

% Load PDB ID of Database 2393 Reduced Samples
Sample_Info=table2struct(readtable(['Package/', 'Reduced_Sample_2393.txt']));

% Merge the Index of Problematic Samples
UPS_Idx=zeros(size(Unique_Problematic_Samples));
for r=1:size(Unique_Problematic_Samples,1)
    if any(strcmp({Sample_Info.ID}, Qualified_Sample_ID(Unique_Problematic_Samples(r)).ID))
        UPS_Idx(r,1)=find(strcmp({Sample_Info.ID}, Qualified_Sample_ID(Unique_Problematic_Samples(r)).ID));
    end
end
UPS_Idx=UPS_Idx(UPS_Idx~=0);

% Remove Outliers
Merged_UPS_Idx=unique(UPS_Idx);
OL_Idx2=setdiff([1:size(Sample_Info,1)], Merged_UPS_Idx,'stable');

%

%%%% Reduced_Feature
Idx=zeros(493,2393);
for r=1:size(OL_Idx1,2)
    for c=1:size(OL_Idx2,2)
        Idx(OL_Idx1(r),OL_Idx2(c))=1;
    end
end
Linear_Idx=reshape(Idx',1,2393*493);
Reduced_Feature=Feature(:,logical(Linear_Idx));

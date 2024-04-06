% Features of Reduced Samples
% Molecular Weight (MW) / Number of Used Patches / Interaction Matrix
load('MPI.mat')

% Crystallization Condition DISM and Conventional DISM of Reduced Samples
load('Seq_DISM_2393.mat')       % Pairwise Sequence Dissimilarity
load('RMSD_DISM_2393.mat')      % Pairwise RMSD
load('Jerhoud_DISM_2393.mat')   % Pairwise Jerhoud 3DZD
load('Kihara_DISM_2393.mat')    % Pairwise Kihara 3DZD
load('Crystallization_Condition_DISM/Condition_DISM_050310_Temp_pH.mat')

Seq_Linear_DISM=squareform(Seq_DISM_2393);
RMSD_Linear_DISM=squareform(RMSD_DISM_2393);
Jehoud_Linear_DISM=squareform(Jerhoud_DISM_2393);
Kihara_Linear_DISM=squareform(Kihara_DISM_2393);

Protein_Sim=[Condition_DISM;
             Seq_Linear_DISM;
             RMSD_Linear_DISM;
             Jehoud_Linear_DISM;
             Kihara_Linear_DISM];

clearvars -except Protein_Sim MPI

% Combine Features of Reduced Samples (1+4+4 +69*2: 5 bin)
Result_Int=load('AAI-PatchBag_Protein_Similarity_DISM/Int_Result_5bin.mat');
Result_nInt=load('AAI-PatchBag_Protein_Similarity_DISM/nInt_Result_5bin.mat');
Feature=[Protein_Sim; MPI; Result_Int.Result; Result_nInt.Result];

%
%
%

%%%% Remove Outliers
% Load ID of Qualified_Sample_Condition_Table (RAW)
Data=readtable(['Package/', 'Qualified_Sample_Condition_Table(RAW).xlsx'],'VariableNamingRule','preserve');
Qualified_Sample_ID=table2struct(Data(:,1));

% Load PDB ID of Training/Test 2393 Reduced Samples
Sample_Info=table2struct(readtable(['Package/', 'Reduced_Sample_2393.txt']));

% Load Unique_Problematic_Samples (Outliers)
load(['Outliers/', 'Unique_Problematic_Samples_3IQR.mat'])

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
OL_Idx=setdiff([1:size(Sample_Info,1)], Merged_UPS_Idx,'stable');

Reduced_Feature=zeros(size(Feature,1), size(OL_Idx,2)*(size(OL_Idx,2)-1)/2);
for f=1:size(Feature,1)
    SQF=squareform(Feature(f,:));
    Reduced_Feature(f,:)=squareform(SQF(OL_Idx,OL_Idx));
end

save('Reduced_Feature_Condition_DISM_050310_Temp_pH.mat','Reduced_Feature')

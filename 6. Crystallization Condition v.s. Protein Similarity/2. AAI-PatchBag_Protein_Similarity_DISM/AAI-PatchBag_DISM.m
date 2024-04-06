% Load PDB ID of Training Samples
% Load AA_Idx_Patch_code.mat
load('AA_Idx_Patch_code_Int_5bin_RA.mat')
load('AA_Idx_Patch_code_nInt_5bin_RA.mat')

Dataset=AA_Idx_Patch_code_RA;
clear AA_Idx_Patch_code_RA

Count=size(Dataset,2);
Result=zeros(length(fieldnames(Dataset))-2,Count*(Count-1)/2);

%
% Result
%
parpool(2)

parfor r=1:(length(fieldnames(Dataset))-2)
    t=0;
    localResult=zeros(1, Count*(Count-1)/2);
    for m=1:(Count-1)
        P1=Dataset(m).(['Patch_code_', num2str(r)]);
        if isempty(P1)
            P1=zeros(300,5);
        end
        for n=(m+1):Count
            P2=Dataset(n).(['Patch_code_', num2str(r)]);
            if isempty(P2)
                P2=zeros(300,5);
            end

            Diff=P2-P1; 

            t=t+1;
            localResult(t) = norm(Diff);
        end
    end
    Result(r, :) = localResult;
    r
end

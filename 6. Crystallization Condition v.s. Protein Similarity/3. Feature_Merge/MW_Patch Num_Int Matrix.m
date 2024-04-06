% Molecular Weight
load('Reduced_Sample_Sequence_Merged.mat')

MW=zeros(size(Sample_Info'));
for m=1:size(Sample_Info,2)
    S=Sample_Info(m).Sequence;
    Modified_S=strrep(strrep(S,'U','C'),'X','G');
    MW(m,1)=molweight(Modified_S);
    m
end

% No. of Int Patch
load('AA_Idx_Patch_code_Int_5bin_RA.mat')

Int_Patch_No=zeros(size(Sample_Info'));
for m=1:size(Sample_Info,2)
    Int_Patch_No(m,1)=sum(AA_Idx_Patch_code_RA(m).Patch_code);
    m
end

% No. of nInt Patch
load('AA_Idx_Patch_code_nInt_5bin_RA.mat')

nInt_Patch_No=zeros(size(Sample_Info'));
for m=1:size(Sample_Info,2)
    nInt_Patch_No(m,1)=sum(AA_Idx_Patch_code_RA(m).Patch_code);
    m
end

% MW/No_Int/No_nInt
Merge=[MW, Int_Patch_No, nInt_Patch_No];

Count=size(Merge,1);
Result=zeros(size(Merge,2),Count*(Count-1)/2);

for r=1:size(Merge,2)
    t=0;
    for m=1:(Count-1)
        P1=Merge(m,r);
        for n=(m+1):Count
            P2=Merge(n,r);

            Distance=abs(P1-P2);
            % Distance=P1-P2;

            t=t+1;
            Result(r,t)=Distance;
        end
    end
end

%
%
%

% Load PDB ID of Qualified Reduced Samples
% Filename='Reduced_Sample_2393.txt';
% Sample_Info=table2struct(readtable(Filename));

% Pre-allocate RAM for Error ID & Patch Struct
Sample_Info(1).Res_Matrix=[];

% Patch Extraction
Blank_ID=[];
for i=1:size(Sample_Info,2)

    % Load Info of Interacting Residues
    PDB_ID=lower(Sample_Info(i).Header);

    Int_File=table2array(readtable(['/home/juju/Desktop/PatchBag/Reduced_Sample_Int_Res_Results/', PDB_ID,'.txt']));
    Int_File_Res=Int_File(:,4); Int_File_Res(Int_File_Res==0)=1;

    if ~isempty(Int_File)
        Patch_File=table2array(readtable(['/home/juju/Desktop/PatchBag/Reduced_Sample_Int_Res_Results/Patch_Results_6plus1_2472/',PDB_ID,'_Patch.txt']));

        if ~isempty(Patch_File)
            Res_Info=Patch_File(:,25:31); Res_Info(Res_Info==0)=1;

            Res_Matrix=zeros(20);
            AA_Index_Matrix=zeros(1,158);
            for j=1:size(Patch_File,1)
                Res_Matrix(Int_File_Res(j,1),Res_Info(j,7))=Res_Matrix(Int_File_Res(j,1),Res_Info(j,7))+1;
            end
            Sample_Info(i).Res_Matrix=Res_Matrix;
        else
            Blank_ID=[Blank_ID; i];
        end
    else
        Blank_ID=[Blank_ID; i];
    end
    i
end

% Output
t=0;
for m=1:(Count-1)
    P1=Sample_Info(m).Res_Matrix;
    if isempty(P1)
        P1=zeros(20,20);
    end
    for n=(m+1):Count
        P2=Sample_Info(n).Res_Matrix;
        if isempty(P2)
            P2=zeros(20,20);
        end

        Distance=norm(P1-P2);

        t=t+1;
        Result(4,t)=Distance;
    end
    m
end

MPI=Result;

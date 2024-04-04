%%%% Load Acceptable Pose
load C2_TopRank_Result.mat;
Set_Idx='.000.';

mkdir('C4.DIPER_P43212');
for r=1:size(C2_TopRank_Result,2)
    for s=1:size(C2_TopRank_Result(r).Output_Idx,1)
        mkdir(['C4.DIPER_P43212/',lower(C2_TopRank_Result(r).ID),Set_Idx,num2str(C2_TopRank_Result(r).Output_Idx(s))]);
        File=['C2.Top25_Pose/',lower(C2_TopRank_Result(r).ID),'/',lower(C2_TopRank_Result(r).ID),Set_Idx,num2str(C2_TopRank_Result(r).Output_Idx(s)),'.pdb'];
        copyfile (File,'C4.DIPER_P43212');
    end
end          
                       
%

%%%% Sample Preparation: Sample_prep.txt
Script_Prep = fopen([pwd,'/Sample_prep.txt'],'w');
fprintf(Script_Prep,'#!/bin/bash');
fprintf(Script_Prep,'\n');
for r=1:size(C2_TopRank_Result,2)
    for s=1:size(C2_TopRank_Result(r).Output_Idx,1)
        fprintf(Script_Prep,'\n');
        fprintf(Script_Prep,['../protein_prep/prepare.py ',lower(C2_TopRank_Result(r).ID),Set_Idx,num2str(C2_TopRank_Result(r).Output_Idx(s)),'.pdb']);   
    end  
end
fclose(Script_Prep);

%

%%%% Script fot PIPER: Sample_PIPER.txt
% Directory Path of "piper_package"
Path='/home/sunmd33452/Desktop/piper_package/C4.DIPER_P43212/';

% Write Sample_PIPER.txt
Script_PIPER = fopen([pwd,'/Sample_PIPER.txt'],'w');
fprintf(Script_PIPER,'#!/bin/bash');
fprintf(Script_PIPER,'\n');
for r=1:size(C2_TopRank_Result,2)
    for s=1:size(C2_TopRank_Result(r).Output_Idx,1)
        fprintf(Script_PIPER,'\n');
        fprintf(Script_PIPER,['../piper -vv -c1.0 -k4 --msur_k=1.0 --maskr=1.0 -T FFTW_EXHAUSTIVE -R 3395 -t 1 -p ../prms/atoms.prm -f ../prms/Coeff_DIPER.prm -r ../prms/C4_rots.prm',' ',lower(C2_TopRank_Result(r).ID),Set_Idx,num2str(C2_TopRank_Result(r).Output_Idx(s)),'_pnon.pdb',' ',lower(C2_TopRank_Result(r).ID),Set_Idx,num2str(C2_TopRank_Result(r).Output_Idx(s)),'_pnon.pdb',' ','--o ',Path,lower(C2_TopRank_Result(r).ID),Set_Idx,num2str(C2_TopRank_Result(r).Output_Idx(s)),'/']);
    end  
end
fclose(Script_PIPER);

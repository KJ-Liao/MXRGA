%%%% Load ID of PDB Sample
% PDB ID line by line (if there were multiple samples)
% Example:
% 1B42
% 1C46
% ...
% Filename='P43212_Reduced_Sample.txt';
Filename='Example.txt';

%%%% Build up List of PDB Sample ID
File=fopen(Filename);
Sample(1).ID=[];m=1;
while (1)
    line=fgetl(File);
    if line==-1, break, end 
    Sample(m).ID=sscanf(line,'%c');
    m=m+1;         
end
fclose(File);

%%%% Download Files of PDB Sample from PDB Web
% Wrong or Changed(Renamed) PDB ID lead to 0 KB html logs
% If So, Please Correct the PDB ID or Remove them from Input_Sample.txt
for n=1:length(Sample)   
    URL=strcat('https://files.rcsb.org/download/',lower(Sample(n).ID),'.pdb');
    websave(lower(Sample(n).ID),URL);
    Reamined_No=length(Sample)-n
end

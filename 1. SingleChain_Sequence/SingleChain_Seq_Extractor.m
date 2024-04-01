%%%% Load ID of Available Single Chain PDB Sample (Before 05/27/2022)
File=fopen('Accessible_SingleChain_PDB_Sample_2022-05-27.txt');
m=1; SingleChain_Sample(1).ID=[];
while (1)
    line=fgetl(File);
    if line==-1, break, end 
    SingleChain_Sample(m).ID=sscanf(line,'%c');
    m=m+1;         
end
fclose(File);

%%%% Extract "Single Chain Human" proteins (SCHs)
SCH_Sample(1).ID=[];
SCH_Sample(1).Seq=[];
Error_Sample(1).ID=[];
i=1; j=1;
for n=1:size(SingleChain_Sample,2)  
    URL=strcat('https://www.rcsb.org/fasta/entry/',SingleChain_Sample(n).ID(1:4),'/display');
    try 
        Info=webread(URL);
        Seq_No=strfind(Info,'>');
        if size(Seq_No,2)==1
            SCH_Sample(i).ID=SingleChain_Sample(n).ID;
            Sta=sort([strfind(Info,'(9606)')+4, strfind(Info,'null')+2]);
            Sto=length(Info)+1;         
            SCH_Sample(i).Seq=Info(Sta(1)+3:Sto(1)-2);
            i=i+1;
        else
            SCH_Sample(i).ID=SingleChain_Sample(n).ID;
        end
    catch
        Error_Sample(j).ID=SingleChain_Sample(n).ID;
        j=j+1;
    end
    n
end
writetable(struct2table(SCH_Sample),'SCH_Seq.xlsx');

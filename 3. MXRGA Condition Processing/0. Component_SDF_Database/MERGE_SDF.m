MERGE=fopen('MERGE_SDF.txt','w');
File=fopen('CID_SDF.txt');
while(1)
    line=fgetl(File);
    if line==-1, break, end

    SDF=fopen(['Structure2D_CID_', line, '.sdf']);
    if SDF==-1, break, end

    while(1)
        SDF_line=fgetl(SDF);
        if SDF_line==-1, break, end
        
        fprintf(MERGE,SDF_line);
        fprintf(MERGE,'\n');
    end
    fclose(SDF);
end
fclose(File);
fclose(MERGE);

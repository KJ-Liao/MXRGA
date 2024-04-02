%%%% Merge Docking Results (Dimer/Tetramer/SurfOrca/Sim_Xtal)
% Load Statistics of Dimer/Tetramer
load('/home/juju/Desktop/MXRGA_Sample Log/4. MXRGA Crystal Packing/P41212_Sample/7. Statistics/3. Tetramer/TopRank_Result.mat')

load('/home/juju/Desktop/MXRGA_Sample Log/4. MXRGA Crystal Packing/P41212_Sample/7. Statistics/2. Dimer/Dimer_RMSD_Statistic.mat')
Dimer_Sample=Sample;

clearvars -except Result Dimer_Sample SurfOrca_Sample

% Load Statistics of MXRGA Simulated Xtal
load('SurfOrca_RMSD_Statistic.mat')

Sample(1).C2_RMSD95=[];
Sample(1).Rank=[];
for r=1:size(Sample,2)
    PDB_ID=upper(Sample(r).ID);
    C2_ID_Idx=strcmp({Dimer_Sample.ID}',PDB_ID);

    Idx=strfind(Sample(r).Output_Idx,'.');  
    C2_Idx=Sample(r).Output_Idx(Idx(2)+1:Idx(3)-1);
    C2_Rank=find(Dimer_Sample(C2_ID_Idx).Output_Idx==str2double(C2_Idx));
    C2_RMSD_95=Dimer_Sample(C2_ID_Idx).RMSD_95(C2_Rank);

    C4_Idx=Sample(r).Output_Idx(Idx(4)+1:end);
    C4_ID_Idx=strcmp({Result.ID}',Sample(r).ID);
    Rank=find(Result(C4_ID_Idx).Linear_Table(:,4)==str2double(C4_Idx) & Result(C4_ID_Idx).Linear_Table(:,2)==C2_Rank);
   
    Sample(r).Rank=Rank;
    Sample(r).C2_RMSD95=C2_RMSD_95;
end

%%%% Pack Quality Claasification
% Class Refined by KJ
Sample(1).Class=[];
for r=1:size(Sample,2)
    if Sample(r).Fnat < 0.3 || (Sample(r).RMSD_95 > 8.0 && Sample(r).iRMS > 4.0)
        % Incorrect
        Sample(r).Class=0;
    elseif (Sample(r).Fnat >= 0.3 && Sample(r).Fnat < 0.5) && (Sample(r).RMSD_95 <= 8.0 || Sample(r).iRMS <= 4.0) || (Sample(r).Fnat >= 0.5 && Sample(r).RMSD_95 > 5.0 && Sample(r).iRMS > 2.0)
        % Acceptable
        Sample(r).Class=1;
    elseif (Sample(r).Fnat >= 0.5 && Sample(r).Fnat < 0.7) && (Sample(r).RMSD_95 <= 5.0 || Sample(r).iRMS <= 2.0) || (Sample(r).Fnat >= 0.7 && Sample(r).RMSD_95 > 3.0 && Sample(r).iRMS > 1.0)
        % Medium
        Sample(r).Class=2;
    elseif Sample(r).Fnat >= 0.7 && (Sample(r).RMSD_95 <= 3.0 || Sample(r).iRMS <= 1.0)
        % High
        Sample(r).Class=3;
    else
        % Undef
        Sample(r).Class=-1;
    end
end

% Class Defined by CAPRI
Sample(1).CAPRI_Class=[];
for r=1:size(Sample,2)
    if Sample(r).Fnat < 0.1 || (Sample(r).LRMS > 10.0 && Sample(r).iRMS > 4.0)
        % Incorrect
        Sample(r).CAPRI_Class=0;
    elseif (Sample(r).Fnat >= 0.1 && Sample(r).Fnat < 0.3) && (Sample(r).LRMS <= 10.0 || Sample(r).iRMS <= 4.0) || (Sample(r).Fnat >= 0.3 && Sample(r).LRMS > 5.0 && Sample(r).iRMS > 2.0)
        % Acceptable
        Sample(r).CAPRI_Class=1;
    elseif (Sample(r).Fnat >= 0.3 && Sample(r).Fnat < 0.5) && (Sample(r).LRMS <= 5.0 || Sample(r).iRMS <= 2.0) || (Sample(r).Fnat >= 0.5 && Sample(r).LRMS > 1.0 && Sample(r).iRMS > 1.0)
        % Medium
        Sample(r).CAPRI_Class=2;
    elseif Sample(r).Fnat >= 0.5 && (Sample(r).LRMS <= 1.0 || Sample(r).iRMS <= 1.0)
        % High
        Sample(r).CAPRI_Class=3;
    else
        % Undef
        Sample(r).CAPRI_Class=-1;
    end
end


%%%% Merged Result
Merged_Result = struct('ID',{Sample.ID},...
                       'Output_Idx',{Sample.Output_Idx},...
                       'Fnat',{Sample.Fnat},...
                       'iRMS',{Sample.iRMS},...
                       'LRMS',{Sample.LRMS},...
                       'DockQ',{Sample.DockQ},...
                       'C2_RMSD95',{Sample.C2_RMSD95},...
                       'C4_RMSD95',{Sample.C4_RMSD95},...
                       'RMSD_95',{Sample.RMSD_95},...
                       'Class',{Sample.Class},...
                       'CAPRI_Class',{Sample.CAPRI_Class},...
                       'Rank',{Sample.Rank});

save('Merged_Result.mat','Merged_Result');

%
%
%

%%%% Statistics
% Load Files
load('Merged_Result.mat')
load('/home/juju/Desktop/MXRGA_Sample Log/4. MXRGA Crystal Packing/P41212_Sample/7. Statistics/3. Tetramer/TopRank_Result.mat')

% CAPRI Class
Unique_ID=unique({Merged_Result.ID}');
Table=nan(size(Unique_ID,1),5);
for r=1:size(Unique_ID,1)
    Corr_Ind=strcmp({Merged_Result.ID}',Unique_ID{r});
    CAPRI_Class=[Merged_Result(Corr_Ind).CAPRI_Class];
    Rank=[Merged_Result(Corr_Ind).Rank];
    for c=1:5
        try
            Table(r,c)=max(CAPRI_Class(Rank<=((c+1)*50)));
        catch
            Table(r,c)=NaN;
        end
    end
end

[size(Result,2)-(sum(Table(:,c)==1)+sum(Table(:,c)==2)+sum(Table(:,c)==3)),...
 (size(Result,2)-(sum(Table(:,c)==1)+sum(Table(:,c)==2)+sum(Table(:,c)==3)))/size(Result,2)*100;
 sum(Table(:,c)==1), sum(Table(:,c)==1)/size(Result,2)*100;
 sum(Table(:,c)==2), sum(Table(:,c)==2)/size(Result,2)*100;
 sum(Table(:,c)==3), sum(Table(:,c)==3)/size(Result,2)*100;
 sum(Table(:,c)==1)+sum(Table(:,c)==2)+sum(Table(:,c)==3),...
 (sum(Table(:,c)==1)+sum(Table(:,c)==2)+sum(Table(:,c)==3))/size(Result,2)*100]

% DockQ
Unique_ID=unique({Merged_Result.ID}');
Table=nan(size(Unique_ID,1),5);
for r=1:size(Unique_ID,1)
    Corr_Ind=strcmp({Merged_Result.ID}',Unique_ID{r});
    DockQ=[Merged_Result(Corr_Ind).DockQ];
    Rank=[Merged_Result(Corr_Ind).Rank];
    for c=1:5
        try
            Table(r,c)=max(DockQ(Rank<=((c+1)*50)));
        catch
            Table(r,c)=NaN;
        end
    end
end

[size(Result,2)-sum(Table(:,c)>=0.23&Table(:,c)<0.49)-sum(Table(:,c)>=0.49&Table(:,c)<0.8)-sum(Table(:,c)>0.8),...
 (size(Result,2)-sum(Table(:,c)>=0.23&Table(:,c)<0.49)-sum(Table(:,c)>=0.49&Table(:,c)<0.8)-sum(Table(:,c)>0.8))/size(Result,2)*100;
 sum(Table(:,c)>=0.23&Table(:,c)<0.49),sum(Table(:,c)>=0.23&Table(:,c)<0.49)/size(Result,2)*100;
 sum(Table(:,c)>=0.49&Table(:,c)<0.8),sum(Table(:,c)>=0.49&Table(:,c)<0.8)/size(Result,2)*100;
 sum(Table(:,c)>0.8),sum(Table(:,c)>0.8)/size(Result,2)*100;
 sum(Table(:,c)>=0.23&Table(:,c)<0.49)+sum(Table(:,c)>=0.49&Table(:,c)<0.8)+sum(Table(:,c)>0.8),...
 (sum(Table(:,c)>=0.23&Table(:,c)<0.49)+sum(Table(:,c)>=0.49&Table(:,c)<0.8)+sum(Table(:,c)>0.8))/size(Result,2)*100]

%
%
%

% Class
Unique_ID=unique({Merged_Result.ID}');
Table=nan(size(Unique_ID,1),5);
for r=1:size(Unique_ID,1)
    Corr_Ind=strcmp({Merged_Result.ID}',Unique_ID{r});
    Class=[Merged_Result(Corr_Ind).Class];
    Rank=[Merged_Result(Corr_Ind).Rank];
    for c=1:5
        try
            Table(r,c)=max(Class(Rank<=((c+1)*50)));
        catch
            Table(r,c)=NaN;
        end
    end
end

[size(Result,2)-(sum(Table(:,c)==1)+sum(Table(:,c)==2)+sum(Table(:,c)==3)),...
 (size(Result,2)-(sum(Table(:,c)==1)+sum(Table(:,c)==2)+sum(Table(:,c)==3)))/size(Result,2)*100;
 sum(Table(:,c)==1), sum(Table(:,c)==1)/size(Result,2)*100;
 sum(Table(:,c)==2), sum(Table(:,c)==2)/size(Result,2)*100;
 sum(Table(:,c)==3), sum(Table(:,c)==3)/size(Result,2)*100;
 sum(Table(:,c)==1)+sum(Table(:,c)==2)+sum(Table(:,c)==3),...
 (sum(Table(:,c)==1)+sum(Table(:,c)==2)+sum(Table(:,c)==3))/size(Result,2)*100]

% RMSD
C2_Ind_55=[Merged_Result.C2_RMSD95]<=5.5;
C4_Ind_55=[Merged_Result.C4_RMSD95]<=5.5;
SurfOrca_Ind=[Merged_Result.RMSD_95]<=8;
C2C4_Ind=([Merged_Result.C2_RMSD95]+[Merged_Result.C4_RMSD95])<=10;

Qualified_Merged_Result=Merged_Result(C2_Ind_55 & C4_Ind_55 & C2C4_Ind & SurfOrca_Ind);

Unique_ID=unique({Qualified_Merged_Result.ID}');
Table=nan(size(Unique_ID,1),5);
for r=1:size(Unique_ID,1)
    Corr_Ind=strcmp({Qualified_Merged_Result.ID}',Unique_ID{r});
    RMSD95=[Qualified_Merged_Result(Corr_Ind).RMSD_95];
    Rank=[Qualified_Merged_Result(Corr_Ind).Rank];
    for c=1:5
        try
            Table(r,c)=min(RMSD95(Rank<=((c+1)*50)));
        catch
            Table(r,c)=NaN;
        end
    end
end

m=3; n=5; o=8;
[size(Result,2)-(sum(Table(:,c)<=m)+sum(Table(:,c)>m&Table(:,c)<=n)+sum(Table(:,c)>n&Table(:,c)<=o)),...
 (size(Result,2)-(sum(Table(:,c)<=m)+sum(Table(:,c)>m&Table(:,c)<=n)+sum(Table(:,c)>n&Table(:,c)<=o)))/size(Result,2)*100;
 sum(Table(:,c)>n&Table(:,c)<=o),sum(Table(:,c)>n&Table(:,c)<=o)/size(Result,2)*100;
 sum(Table(:,c)>m&Table(:,c)<=n),sum(Table(:,c)>m&Table(:,c)<=n)/size(Result,2)*100;
 sum(Table(:,c)<=m),sum(Table(:,c)<=m)/size(Result,2)*100;
 sum(Table(:,c)<=m)+sum(Table(:,c)>m&Table(:,c)<=n)+sum(Table(:,c)>n&Table(:,c)<=o),...
 (sum(Table(:,c)<=m)+sum(Table(:,c)>m&Table(:,c)<=n)+sum(Table(:,c)>n&Table(:,c)<=o))/size(Result,2)*100]

%
%
%

% PackQ
best_d1=1.5;
best_d2=10.0; 
C=[0.36, 0.54, 0.67];

Merged_Result(1).PackQ=[];
for m=1:size(Merged_Result,2)
    Fnat=Merged_Result(m).Fnat;
    iRMS=Merged_Result(m).iRMS;
    RMSD_95=Merged_Result(m).RMSD_95;
    Merged_Result(m).PackQ=(Fnat + 1/(1+(iRMS/best_d1)^2) + 1/(1+(RMSD_95/best_d2)^2))/3;
end

Unique_ID=unique({Merged_Result.ID}');
Table=nan(size(Unique_ID,1),5);
for r=1:size(Unique_ID,1)
    Corr_Ind=strcmp({Merged_Result.ID}',Unique_ID{r});
    PackQ=[Merged_Result(Corr_Ind).PackQ];
    Rank=[Merged_Result(Corr_Ind).Rank];
    for c=1:5
        try
            Table(r,c)=max(PackQ(Rank<=((c+1)*50)));
        catch
            Table(r,c)=NaN;
        end
    end
end

[size(Result,2)-sum(Table(:,c)>=C(1)),...
 (size(Result,2)-sum(Table(:,c)>=C(1)))/size(Result,2)*100;
 sum(Table(:,c)>=C(1)&Table(:,c)<C(2)), sum(Table(:,c)>=C(1)&Table(:,c)<C(2))/size(Result,2)*100;
 sum(Table(:,c)>=C(2)&Table(:,c)<C(3)), sum(Table(:,c)>=C(2)&Table(:,c)<C(3))/size(Result,2)*100;
 sum(Table(:,c)>=C(3)), sum(Table(:,c)>=C(3))/size(Result,2)*100;
 sum(Table(:,c)>=C(1)), sum(Table(:,c)>=C(1))/size(Result,2)*100]
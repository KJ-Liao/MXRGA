Table=nan(size(Sample,2),5);
for r=1:size(Sample,2)
    for c=1:5
        try
            Table(r,c)=max(Sample(r).DockQ(1:c*5));
        catch
            Table(r,c)=max(Sample(r).DockQ(1:length(Sample(r).DockQ)));
        end
    end
end

sum(Table(:,c)<0.23)
sum(Table(:,c)>=0.23&Table(:,c)<0.49)
sum(Table(:,c)>=0.49&Table(:,c)<0.8)
sum(Table(:,c)>=0.8)

sum(Table(:,c)<0.23)/size(Sample,2)
sum(Table(:,c)>=0.23&Table(:,c)<0.49)/size(Sample,2)
sum(Table(:,c)>=0.49&Table(:,c)<0.8)/size(Sample,2)
sum(Table(:,c)>=0.8)/size(Sample,2)

%

Table=nan(size(Sample,2),5);
for r=1:size(Sample,2)
    for c=1:5
        try
            Table(r,c)=min(Sample(r).RMSD_95(1:c*5));
        catch
            Table(r,c)=min(Sample(r).RMSD_95(1:length(Sample(r).RMSD_95)));
        end
    end
end

sum(Table(:,c)<=1)
sum(Table(:,c)>1&Table(:,c)<=3)
sum(Table(:,c)>3&Table(:,c)<=5)
sum(Table(:,c)>5)

sum(Table(:,c)<=1)/size(Sample,2)
sum(Table(:,c)>1&Table(:,c)<=3)/size(Sample,2)
sum(Table(:,c)>3&Table(:,c)<=5)/size(Sample,2)
sum(Table(:,c)>5)/size(Sample,2)
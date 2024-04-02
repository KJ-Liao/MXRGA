ID=unique({Sample.ID}');

Result=zeros(size(ID,1),2);
for r=1:size(ID,1)
    Idx=strcmp({Sample.ID}',ID{r});
    Result(r,1)=max([Sample(Idx).DockQ]);
    Result(r,2)=min([Sample(Idx).RMSD_95]);
end

sum(Result(:,1)>0.3)
sum(Result(:,2)<8)
%%%% AA_Idx Reduction
% 1-R2 as Similarity Distance 
load('AA_Index.mat')
R2=corrcoef(AA_Index').^2;

% Hierachical Clustering (R2<0.3 as Cutoff)
Z=linkage(squareform(1-R2),'average');
T=cluster(Z,'cutoff',0.3,'criterion','distance');

% Representative AA_Idx Feature
Representative_Clust=zeros(max(T),1);
for t=1:max(T)
    Ind=find(T==t);
    if size(Ind,1)==1
        Representative_Clust(t,1)=Ind;
    elseif size(Ind,1)==2
        Ind1=sort(R2(:,Ind(1)),'descend');
        Ind2=sort(R2(:,Ind(2)),'descend');
        if Ind1(3)>Ind2(3)
            Representative_Clust(t,1)=Ind(2);
        else
            Representative_Clust(t,1)=Ind(1);
        end
    else
         [~,c]=min(sum(R2(Ind, Ind)));
         Representative_Clust(t,1)=Ind(c);
    end
end

%%%% AA_Idx_Patch_code Table Rearrangement
clear AA_Idx_Patch_code_RA;
AA_Idx_Patch_code_RA(1).ID=[];
AA_Idx_Patch_code_RA(1).Patch_code=[];

for u=1:size(AA_Idx_Patch_code,2)
    AA_Idx_Patch_code_RA(u).ID=AA_Idx_Patch_code(u).ID;
    AA_Idx_Patch_code_RA(u).Patch_code=AA_Idx_Patch_code(u).Patch_code;
    for v=1:max(T)
        AA_Idx_Patch_code_RA(u).(['Patch_code_', num2str(v)])=AA_Idx_Patch_code(u).(['Patch_code_', num2str(Representative_Clust(v))]);
    end
    u
end
save('AA_Idx_Patch_code_RA.mat','AA_Idx_Patch_code_RA')
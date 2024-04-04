%%%% Pairwise sequence alignment matix (DISM) for Rplot 
% Unadjusted dissimilarity as distance
Data=table2struct(readtable('P43212_Sample.xlsx'));
D=seqpdist({Data(1:size(Data,1)).Seq},'Method', 'p-distance');

% Output
Output=[{'PDB_ID'},{Data.ID}; {Data.ID}',num2cell(squareform(D))];
writecell(Output,'DISM.xlsx');
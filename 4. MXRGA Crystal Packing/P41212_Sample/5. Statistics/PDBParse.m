function PDB_Info = PDBParse (FileName)

%%% Load Corresponding PDB File
% Load Corresponding PDB File
FileID = fopen(FileName);
RawText= fread(FileID,inf,'*char');

%%% Parse Raw Text
% Parse lines by end-of-lines
SplitLines = textscan(RawText,'%s','delimiter',newline);

% Pre-allocate Memorie for Output Info
NumLines   = length(SplitLines{1});
recordName = cell(1,NumLines);
atomNum    = cell(1,NumLines);
atomName   = cell(1,NumLines);
altLoc     = cell(1,NumLines);
resName    = cell(1,NumLines);
chainID    = cell(1,NumLines);
resNum     = cell(1,NumLines);
X          = cell(1,NumLines);
Y          = cell(1,NumLines);
Z          = cell(1,NumLines);

% Read Each Line
for n=1:NumLines
    % Parse Line by Line
    line = cell2mat(SplitLines{1}(n));

    % Preserve Non-H Atoms with Firt altLoc
    if length(line)>4 && strcmp(line(1:4),'ATOM')...                                        % ATOM
            && strcmp(line(27),' ')...                                                      % No insertions
            && ~strcmp(line(13),'H') && ~strcmp(line(14),'H')...                            % Not H ATOM
            && ~strcmp(line(13),'X') && ~strcmp(line(14),'X')...                            % Not X ATOM
            && (strcmp(line(17),' ') || strcmp(line(17),'A') || strcmp(line(17),'1'))       % First altLoc

        recordName(n) = {line(1:6)};
        atomNum(n)    = {line(7:11)};
        atomName(n)   = {line(13:16)};
        altLoc(n)     = {line(17)};
        resName(n)    = {line(18:20)};
        chainID(n)    = {line(22)};
        resNum(n)     = {line(23:26)};
        X(n)          = {line(31:38)};
        Y(n)          = {line(39:46)};
        Z(n)          = {line(47:54)};      
    end
end

% Extract Info of Non-H Atoms with Firt altLoc
Extract_Idx = strcmp(recordName,'ATOM  ');
recordName  = recordName(Extract_Idx);
atomNum     = atomNum(Extract_Idx);
atomName    = atomName(Extract_Idx);
altLoc      = altLoc(Extract_Idx);
resName     = resName(Extract_Idx);
chainID     = chainID(Extract_Idx);
resNum      = resNum(Extract_Idx);
X           = X(Extract_Idx);
Y           = Y(Extract_Idx);
Z           = Z(Extract_Idx);

%%% Rearrange PDB Info in Struct
PDB_Info.recordName = (strtrim(recordName))';
PDB_Info.chainID    = (chainID)';
PDB_Info.altLoc     = (altLoc)';
PDB_Info.atomNum    = (str2double(atomNum))';
PDB_Info.atomName   = (strtrim(atomName))';
PDB_Info.resName    = (strtrim(resName))';
PDB_Info.resNum     = (str2double(resNum))';
PDB_Info.X          = (str2double(X))';
PDB_Info.Y          = (str2double(Y))';
PDB_Info.Z          = (str2double(Z))';

% ILE_Modification
PDB_Info.atomName(strcmp(PDB_Info.resName,'ILE') & strcmp(PDB_Info.atomName,'CD'))={'CD1'};

fclose(FileID);
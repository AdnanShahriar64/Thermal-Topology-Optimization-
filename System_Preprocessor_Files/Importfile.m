function Ab2Mat = Importfile(filename)
%IMPORTFILE Read an Abaqus .inp file into a cell array of text fields.
%   AB2MAT = IMPORTFILE(FILENAME) returns an N x 16 cell array of char:
%   one row per non-blank line, the line split at commas, each field with
%   surrounding spaces removed, and short lines padded with ''. Fields past
%   the 16th are dropped (as in the original MATLAB Import Tool version).
%   Works in both MATLAB and Octave (no string arrays).

nCol = 16;

fid = fopen(filename, 'r');
if fid < 0
    error('Importfile:notFound', ...
          'Cannot open "%s". Current folder is "%s".\nRun from the folder that contains the .inp file, or pass its full path.', ...
          filename, pwd);
end
txt = fread(fid, inf, '*char')';
fclose(fid);

lines = regexp(txt, '\r\n|\n|\r', 'split');          % any line ending
lines = lines(~cellfun(@(s) all(isspace(s)), lines)); % drop blank lines

Ab2Mat = repmat({''}, numel(lines), nCol);
for r = 1:numel(lines)
    f = strtrim(strsplit(lines{r}, ','));
    n = min(numel(f), nCol);
    Ab2Mat(r, 1:n) = f(1:n);
end

end

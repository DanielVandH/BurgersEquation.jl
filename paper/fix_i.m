% fix_i(filename)
%
% Edits the text in the numeric data file at `filename', replacing
% instances of 'i' with 'im'. This allows the data to be used in Julia.
function fix_i(filename)
fid = fopen(filename);
f = fread(fid, '*char');
fclose(fid);
f = strrep(f', 'i', 'im')';
fid = fopen(filename, 'w');
fprintf(fid, '%s', f);
fclose(fid);
end
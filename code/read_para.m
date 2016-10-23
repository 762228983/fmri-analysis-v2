function P = read_para(para_file, varargin)


%%
try
fid = fopen(para_file, 'r');
file_contents = textscan(fid,'%f%d%f%f%s');
fclose(fid);
catch
    fprintf('Problem in read_para.m with file:\n%s\n', para_file);
    keyboard
end

[P.onsets, P.condition_indices, P.durs, ~, P.conds] = file_contents{:}; %#ok<*NASGU>

if optInputs(varargin, 'remove-NULL')
    xi = ~ismember(P.conds, 'NULL');
    P.onsets = P.onsets(xi);
    P.condition_indices = P.condition_indices(xi);
    P.durs = P.durs(xi);
    P.conds = P.conds(xi);
end


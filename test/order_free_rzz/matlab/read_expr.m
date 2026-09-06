function [X, genes] = read_expr(fn)
%READ_EXPR  Read the repo's tab-delimited expression matrix (genes x cells).
%   Works in MATLAB and Octave (no readtable dependency).
    fid = fopen(fn, 'r'); assert(fid > 0, 'cannot open %s', fn);
    fgetl(fid);                                   % header line
    genes = {}; rows = {};
    while true
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        if isempty(strtrim(ln)), continue; end
        t = regexp(strtrim(ln), '\t', 'split');
        genes{end+1,1} = t{1};                    %#ok<AGROW>
        rows{end+1,1} = str2double(t(2:end));     %#ok<AGROW>
    end
    fclose(fid);
    X = cell2mat(rows);
end

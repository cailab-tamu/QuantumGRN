function [genes, pObs, act, ncells, B] = qgrn_read_data(fn, threshold)
%QGRN_READ_DATA  Preprocess the expression matrix for the QuantumGRN model.
%
%   [genes, pObs, act, ncells] = QGRN_READ_DATA(FN, THRESHOLD) reads the
%   tab-delimited genes-by-cells matrix FN and reproduces three steps of the
%   Python package: qsc_order_gene, qsc_distribution and qsc_activation_ratios.
%
%   Outputs
%     genes  : 1xn cellstr, ordered by activation frequency (descending)
%     pObs   : 2^n x 1 observed distribution over binary activation states
%     act    : n x 1 activation ratio per gene, in the ordered gene order
%     ncells : number of cells
%     B      : n x ncells binarized matrix, rows in the ordered gene order
%
%   BIT ORDERING (this is the detail that silently corrupts results if wrong):
%   gene k occupies bit k-1, so gene 1 is the LEAST significant bit and the
%   state index is  idx = sum_k b_k * 2^(k-1).  This matches the Python
%   package and the toolbox's basis-state strings -- see check_conventions.m.

    if nargin < 2, threshold = 0; end

    % ---- read (genes x cells) -------------------------------------------
    fid = fopen(fn, 'r');
    assert(fid > 0, 'qgrn:cannotOpen', 'cannot open %s', fn);
    cleanup = onCleanup(@() fclose(fid));
    fgetl(fid);                                       % discard header row
    names = {}; rows = {};
    while true
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        if isempty(strtrim(ln)), continue; end
        t = regexp(strtrim(ln), '\t', 'split');
        names{end+1,1} = t{1};                                    %#ok<AGROW>
        rows{end+1,1}  = str2double(t(2:end));                    %#ok<AGROW>
    end
    X = cell2mat(rows);                               % n x m
    [~, ncells] = size(X);

    % ---- binarize and order genes by activation count (descending) ------
    B = X > threshold;
    [~, ord] = sort(sum(B, 2), 'descend');            % MATLAB sort is stable
    B = B(ord, :);
    genes = names(ord)';
    act = sum(B, 2) / ncells;

    % ---- observed distribution over the 2^n activation states -----------
    pObs = qgrn_pobs(B);
end

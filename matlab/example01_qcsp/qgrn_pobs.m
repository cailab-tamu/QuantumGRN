function pObs = qgrn_pobs(B)
%QGRN_POBS  Observed distribution over the 2^n activation states.
%
%   POBS = QGRN_POBS(B) takes an n x ncells logical matrix (row k = gene k,
%   in whatever gene order you intend to use as the qubit mapping) and
%   returns the 2^n x 1 histogram of activation patterns, with |0...0>
%   zeroed and the rest rescaled.
%
%   Gene k occupies bit k-1, so gene 1 is the least significant bit.
%   Permuting the rows of B therefore changes the qubit mapping -- which is
%   precisely the degree of freedom probed by the ordering sweeps summarised
%   in README.md.

    [n, ncells] = size(B);
    lab = zeros(1, ncells);
    for k = 1:n
        lab = lab + B(k, :) * 2^(k-1);
    end
    pObs = accumarray(lab(:) + 1, 1, [2^n 1]) / ncells;
    pObs(1) = 0;
    pObs = pObs / sum(pObs);
end

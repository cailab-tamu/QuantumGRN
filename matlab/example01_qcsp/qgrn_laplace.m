function [ph, N] = qgrn_laplace(p, ncells, alpha)
%QGRN_LAPLACE  Laplace (add-alpha) smoothing of a probability vector.
%
%   [PH, N] = QGRN_LAPLACE(P, NCELLS, ALPHA) converts the distribution P back
%   to counts (NCELLS*P), adds ALPHA to every one of the 2^n bins, and
%   renormalises. N is the smoothed total count, needed by the chain rule in
%   qgrn_gradient.
%
%   Smoothing matters here because p_obs has exact zeros for activation
%   states never seen in the data, and KL divergence would otherwise diverge.

    if nargin < 3, alpha = 1; end
    d = ncells * p;
    N = sum(d);
    d = d + alpha;
    N = N + alpha * numel(p);
    ph = d / N;
end

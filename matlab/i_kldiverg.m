function [kl]=i_kldiverg(pt,po,rm0,m)
%I_KLDIVERG One-directional KL divergence KL(po || pt), as the paper defines it.
%
%   KL = I_KLDIVERG(PT,PO) returns sum(po.*log(po./pt)): the divergence OF
%   the model output PO FROM the target PT. Argument order is (target, model),
%   which is how the call sites in this repo pass them.
%
%   KL = I_KLDIVERG(PT,PO,TRUE) first drops the all-zero state and
%   renormalises both, as the paper does for p_obs and p_out.
%
%   DIRECTION NOW MATTERS. This used to return the symmetrised (Jeffreys)
%   average, 10*mean([KL(po||pt) KL(pt||po)]), which is invariant under
%   swapping the two arguments. The published loss (Eq. 5) is the single
%   direction KL(p_out || p_obs), so that is what this returns now.
%
%   Consequence for other callers: any site passing (model, target) rather
%   than (target, model) now computes the opposite divergence. There is one --
%   s5_asymmetri.m:122, i_kldiverg(bestpo,pt) -- which is reversed relative to
%   every other call in the repo. It was harmless while this was symmetric.
%
%   NO SCALE FACTOR. This returned 10*KL until 2026-09-08. The 10 is not in
%   the paper and it made every reported value ten times the divergence, so
%   comparing against a published KL (4.2504e-03) or against qgrn_objective.m,
%   which has no such factor, needed a silent division. Removing it is safe:
%   scaling an objective by a positive constant cannot move its argmin, and
%   measured on the 6-gene Model B fit the two agree to 2.7e-07 in theta and
%   7.8e-14 in KL, with the same 12 of 15 edges surviving |theta| > pi/180.
%   The only mechanism by which the scale mattered is that it also scaled the
%   gradient, so fminunc's absolute OptimalityTolerance bit a fraction of a
%   step earlier or later.
%
%   Two consequences of the removal, both handled: every displayed "kl=%f"
%   title across s4/s5/xs2/xs3/pkg_e_qgrn_* drops by 10x, and the
%   best-of-restarts tests in s5_asymmetri.m and y1_diffnet.m compared
%   against a hardcoded bestkl = 1, which was scale-dependent; those are now
%   bestkl = Inf.
%
%   KL = I_KLDIVERG(PT,PO,RM0,M) additionally applies the paper's Laplace
%   smoothing with alpha = 1 over M cells, phat = (p*M + 1)/(M + N), before
%   taking the divergence. Omit M and no smoothing is applied, which is what
%   every other script in this repo still gets.
%
%   Smoothing here is applied to the states that REMAIN after RM0, so N is 63
%   rather than the paper's 64. That is equivalent to the paper's version to
%   ~4e-5 relative: the (p*M + 1) numerators are identical, the denominator
%   cancels inside the log ratio, and the all-zero state the paper keeps
%   contributes exactly zero to the sum anyway (its smoothed value is
%   1/(M+64) in both distributions).
    if nargin<3, rm0=false; end
    if nargin<4, m=[]; end
    if rm0
        pt=pt(2:end);
        po=po(2:end);
        pt=pt./sum(pt);
        po=po./sum(po);
    end
    if ~isempty(m)
        N  = numel(po);                        % Laplace smoothing, alpha = 1
        po = (po*m + 1)/(m + N);
        pt = (pt*m + 1)/(m + N);
    end
    kl = sum(po .* (log(po)-log(pt)));         % KL(po || pt), Eq. 5
end
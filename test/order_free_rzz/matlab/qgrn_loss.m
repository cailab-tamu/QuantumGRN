function [L, Lkl, Lc] = qgrn_loss(theta, pout, pobs_s, m, lam)
%QGRN_LOSS  KL(smooth(pout) || smooth(pobs)) + lam * L_cons  (Eqs. 5-7 of the paper).
    ps  = (pout*m + 1) / (m + numel(pout));       % Laplace smoothing, alpha = 1
    Lkl = sum(ps .* log(ps ./ pobs_s));
    Lc  = sum(1 ./ (theta(:).^4 - (pi/2)^4).^2);
    L   = Lkl + lam * Lc;
end

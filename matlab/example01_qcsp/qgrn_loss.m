function [L, cons, err] = qgrn_loss(theta, pOut, pObs, ncells, ratio)
%QGRN_LOSS  KL term, boundary-constraint term and squared error.
%
%   [L, CONS, ERR] = QGRN_LOSS(THETA, POUT, POBS, NCELLS, RATIO)
%     L    : KL( p_out^ || p_obs^ ) on the Laplace-smoothed distributions
%     CONS : constraint keeping every angle inside (-pi/2, pi/2), scaled by RATIO
%     ERR  : plain sum of squared differences on the unsmoothed distributions
%
%   The constraint is  sum  ratio / (theta^4 - (pi/2)^4)^2 , which blows up as
%   any angle approaches +-pi/2. RATIO is recomputed every iteration as
%   L/constraint(theta,1) so the two terms stay comparable in magnitude.

    if nargin < 5, ratio = 1; end
    hOut = qgrn_laplace(pOut, ncells);
    hObs = qgrn_laplace(pObs, ncells);
    L    = sum(hOut .* log(hOut ./ hObs));
    cons = sum(sum( ratio ./ (theta.^4 - (pi/2)^4).^2 ));
    err  = sum((pOut - pObs).^2);
end

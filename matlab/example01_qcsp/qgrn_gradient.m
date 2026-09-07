function G = qgrn_gradient(theta, pOut, pObs, ncells, ratio, trainEncoder)
%QGRN_GRADIENT  Gradient of the total loss w.r.t. every circuit angle.
%
%   G = QGRN_GRADIENT(THETA, POUT, POBS, NCELLS, RATIO, TRAINENCODER)
%
%   Chain rule, term by term:
%
%     L = sum_i h_i log(h_i / g_i)      h = smoothed p_out, g = smoothed p_obs
%     dL/dh_i = 1 + log(h_i / g_i)
%     dh_i/dtheta = (ncells / N_hat) * dp_i/dtheta        (from qgrn_laplace)
%     dp_i/dtheta  comes from the parameter-shift rule in qgrn_dpraw
%
%   plus the derivative of the boundary constraint,
%     d/dtheta [ ratio / (theta^4 - (pi/2)^4)^2 ]
%       = -8 ratio theta^3 / (theta^4 - (pi/2)^4)^3
%
%   SYMMETRISATION: theta(a,b) and theta(b,a) are two distinct gates, but the
%   model treats an edge as undirected, so both receive the averaged gradient
%   (dp/dtheta_ab + dp/dtheta_ba)/2. Starting from a symmetric theta this keeps
%   the matrix symmetric for the whole run -- which is why the fitted network
%   is undirected.
%
%   NOTE ON FIDELITY TO THE PYTHON: the Python package forms dp/dtheta from the
%   RAW Born probabilities while evaluating the loss on the |0..0>-rescaled
%   ones, i.e. it omits the derivative of the rescaling step. This function
%   reproduces that choice deliberately so the two implementations agree
%   numerically; it is an approximation, not an identity.

    n = size(theta, 1);
    [hOut, Nhat] = qgrn_laplace(pOut,  ncells);
    hObs         = qgrn_laplace(pObs, ncells);

    coef = (1 + log(hOut ./ hObs)) * (ncells / Nhat);

    % derivative of the constraint term, for every angle
    derConst = (-8 * ratio * theta.^3) ./ (theta.^4 - (pi/2)^4).^3;
    if ~trainEncoder
        derConst(1:n+1:end) = 0;                       % clear the diagonal
    end

    % raw-probability derivatives, one parameter-shift evaluation per angle
    dp = cell(n, n);
    for k = 1:n
        for q = 1:n
            if k ~= q || trainEncoder
                dp{k,q} = qgrn_dpraw(theta, k, q);
            end
        end
    end

    G = zeros(n, n);
    if trainEncoder
        for k = 1:n
            G(k,k) = sum(coef .* dp{k,k}) + derConst(k,k);
        end
    end
    for a = 1:n
        for b = 1:n
            if a ~= b
                dsym = (dp{a,b} + dp{b,a}) / 2;
                G(a,b) = sum(coef .* dsym) + derConst(a,b);
            end
        end
    end
end

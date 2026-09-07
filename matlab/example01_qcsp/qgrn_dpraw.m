function d = qgrn_dpraw(theta, k, q)
%QGRN_DPRAW  d(raw probabilities)/d(theta(k,q)) by the parameter-shift rule.
%
%   D = QGRN_DPRAW(THETA, K, Q) returns a 2^n x 1 vector of exact derivatives,
%   computed only from circuit simulations at shifted parameter values -- no
%   hand-derived derivative circuits and no finite differences.
%
%   WHICH SHIFT RULE?  That depends on the gate's generator spectrum.
%
%     Ry(t) = exp(-i t/2 Y).  Y has eigenvalues +-1, one gap, so the familiar
%     TWO-term rule is exact:
%         df/dt = [ f(t+pi/2) - f(t-pi/2) ] / 2
%
%     c-Ry(t) = exp(-i t/2 |1><1| (x) Y).  The generator has eigenvalues
%     {0, 0, +1/2, -1/2}, giving gaps 1/2 AND 1, so f(t) contains two
%     frequencies and the two-term rule is WRONG. The FOUR-term rule is exact:
%         df/dt = c+ [f(t+pi/2)  - f(t-pi/2) ]
%               - c- [f(t+3pi/2) - f(t-3pi/2)]
%         c+- = (sqrt(2) +- 1) / (4 sqrt(2))
%
%   check_conventions.m demonstrates the difference: on this circuit the
%   two-term rule is off by ~6% while the four-term rule matches finite
%   differences to ~1e-12.
%
%   Reference: Wierichs, Izaac, Wang & Lin, "General parameter-shift rules for
%   quantum gradients", Quantum 6, 677 (2022).

    shift = @(s) setparam(theta, k, q, theta(k,q) + s);

    if k == q
        % encoder Ry: two-term rule
        [~, pp] = qgrn_pout(shift( pi/2));
        [~, pm] = qgrn_pout(shift(-pi/2));
        d = (pp - pm) / 2;
    else
        % controlled Ry: four-term rule
        cp = (sqrt(2) + 1) / (4*sqrt(2));
        cm = (sqrt(2) - 1) / (4*sqrt(2));
        [~, p1] = qgrn_pout(shift( pi/2));
        [~, p2] = qgrn_pout(shift(-pi/2));
        [~, p3] = qgrn_pout(shift( 3*pi/2));
        [~, p4] = qgrn_pout(shift(-3*pi/2));
        d = cp*(p1 - p2) - cm*(p3 - p4);
    end
end

function T = setparam(T, k, q, v)
    T(k,q) = v;
end

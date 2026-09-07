function [pOut, pRaw] = qgrn_pout(theta)
%QGRN_POUT  Simulate the QuantumGRN circuit and return its output distribution.
%
%   [pOut, pRaw] = QGRN_POUT(THETA)
%     pRaw : 2^n x 1 raw Born probabilities |<x|psi>|^2
%     pOut : same, with the |0...0> state zeroed and the rest rescaled
%            (the normalisation used throughout the paper)
%
%   Uses simulate() and querystates() from the Quantum Computing Toolbox.
%   Basis-state strings are indexed with qubit k as bit k-1 (see
%   qgrn_read_data for why).

    n = size(theta, 1);
    S = simulate(qgrn_circuit(theta));
    [states, probs] = querystates(S);
    states = string(states);

    pRaw = zeros(2^n, 1);
    for i = 1:numel(states)
        b = char(states(i)) == '1';                   % b(k) is qubit k
        idx = sum(b .* 2.^(0:n-1));
        pRaw(idx+1) = probs(i);
    end

    pOut = pRaw;
    pOut(1) = 0;
    pOut = pOut / sum(pOut);
end

function [p, C] = pout_qscgrn_qcsp(theta, bitorder)
%POUT_QSCGRN_QCSP  Model A via the MATLAB Quantum Computing Support Package.
%   Same circuit as pout_qscgrn: Ry encoder, then for each control k in
%   order, c-Ry to every target p ~= k.
    if nargin < 2, bitorder = 'q1lsb'; end
    n = size(theta,1);
    g = [];
    for k = 1:n, g = [g; ryGate(k, theta(k,k))]; end            %#ok<AGROW>
    for k = 1:n
        for q = 1:n
            if q ~= k, g = [g; cryGate(k, q, theta(k,q))]; end  %#ok<AGROW>
        end
    end
    C = quantumCircuit(g);
    p = qcsp_probs(C, n, bitorder);
end

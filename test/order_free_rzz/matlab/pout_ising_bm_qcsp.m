function [p, C] = pout_ising_bm_qcsp(theta, mix, bitorder)
%POUT_ISING_BM_QCSP  Model B via the MATLAB Quantum Computing Support Package.
%   Same circuit as pout_ising_bm: Ry encoder, all-pairs RZZ, Rx mixer.
%   Returns p_out and the quantumCircuit (use generateQASM(C) for hardware).
    if nargin < 3, bitorder = 'q1lsb'; end
    n = size(theta,1);
    g = [];
    for k = 1:n, g = [g; ryGate(k, theta(k,k))]; end            %#ok<AGROW>
    for k = 1:n-1
        for q = k+1:n, g = [g; rzzGate(k, q, theta(k,q))]; end   %#ok<AGROW>
    end
    for k = 1:n, g = [g; rxGate(k, mix)]; end                   %#ok<AGROW>
    C = quantumCircuit(g);
    p = qcsp_probs(C, n, bitorder);
end

function p = qcsp_probs(C, n, bitorder)
%QCSP_PROBS  Simulate a quantumCircuit and return p_out indexed like build_pobs.
%   bitorder = 'q1lsb' : character k of the basis-state string is qubit k,
%                        index = sum_k b_k 2^(k-1)   (this code's convention)
%   bitorder = 'q1msb' : character k is qubit n-k+1 (reversed).
%   |0>^n is set to 0 and the rest rescaled, as in the paper.
    S = simulate(C);
    [states, probs] = querystates(S);
    states = string(states);
    p = zeros(2^n, 1);
    for i = 1:numel(states)
        b = char(states(i)) == '1';
        if strcmp(bitorder, 'q1msb'), b = fliplr(b); end
        idx = sum(b .* 2.^(0:n-1));
        p(idx+1) = probs(i);
    end
    p(1) = 0; p = p / sum(p);
end

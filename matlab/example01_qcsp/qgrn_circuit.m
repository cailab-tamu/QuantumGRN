function C = qgrn_circuit(theta)
%QGRN_CIRCUIT  Build the QuantumGRN circuit with Quantum Computing Toolbox gates.
%
%   C = QGRN_CIRCUIT(THETA) returns a quantumCircuit for the n x n parameter
%   matrix THETA, where
%     THETA(k,k) is the encoder angle for gene k        -> ryGate
%     THETA(k,q) is the regulation angle for k -> q     -> cryGate (control k)
%
%   The circuit is the one in the paper:
%     L_enc : Ry(theta_kk) on every qubit          (encodes activation ratios)
%     L_reg : for each control k, a c-Ry to every target q ~= k
%
%   Gate order matters -- c-Ry gates with different controls do not commute,
%   which is exactly why this model is order-sensitive.

    n = size(theta, 1);
    gates = ryGate(1, theta(1,1));
    for k = 2:n
        gates = [gates; ryGate(k, theta(k,k))];                   %#ok<AGROW>
    end
    for k = 1:n
        for q = 1:n
            if q ~= k
                gates = [gates; cryGate(k, q, theta(k,q))];       %#ok<AGROW>
            end
        end
    end
    C = quantumCircuit(gates);
end

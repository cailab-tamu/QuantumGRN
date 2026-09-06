function psi = apply_1q(psi, U, k, n)
%APPLY_1Q  Apply 2x2 unitary U to qubit k of an n-qubit statevector.
    psi = reshape(psi, [2^(k-1), 2, 2^(n-k)]);
    psi = reshape(permute(psi, [2 1 3]), 2, []);
    psi = U * psi;
    psi = permute(reshape(psi, [2, 2^(k-1), 2^(n-k)]), [2 1 3]);
    psi = psi(:);
end

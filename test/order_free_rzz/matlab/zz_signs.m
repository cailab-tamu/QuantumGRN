function Z = zz_signs(n)
%ZZ_SIGNS  2^n x n matrix of Pauli-Z eigenvalues (+1 for bit 0, -1 for bit 1).
    idx = (0:2^n-1)';
    Z = zeros(2^n, n);
    for k = 1:n, Z(:,k) = 1 - 2*bitget(idx, k); end
end

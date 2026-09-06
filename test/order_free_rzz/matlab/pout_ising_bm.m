function p = pout_ising_bm(theta, mix, Z)
%POUT_ISING_BM  Model B: order-free Ising Born machine.
%   Ry(theta_kk) encoder -> prod_{k<p} RZZ(theta_kp) (commuting, diagonal)
%   -> Rx(mix) on every qubit -> measure. No gene ordering exists.
%   Z = zz_signs(n) (precomputed).
    n = size(theta,1); G = gate_defs();
    psi = zeros(2^n,1); psi(1) = 1;
    for k = 1:n, psi = apply_1q(psi, G.Ry(theta(k,k)), k, n); end
    phase = zeros(2^n,1);
    for k = 1:n-1
        for q = k+1:n, phase = phase + theta(k,q) * (Z(:,k) .* Z(:,q)); end
    end
    psi = psi .* exp(-0.5i * phase);              % RZZ(t) = exp(-i t/2 Z⊗Z)
    for k = 1:n, psi = apply_1q(psi, G.Rx(mix), k, n); end
    p = abs(psi).^2; p(1) = 0; p = p / sum(p);
end

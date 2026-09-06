import autograd.numpy as np
from autograd import grad
import numpy as onp, pandas as pd, itertools, sys

# ---------------- data ----------------
X = pd.read_csv('../../dataset/expr_matrix_pearsonresidual_7.txt', sep='\t', index_col=0)
Xb = (X.values > 0).astype(int)
genes_all = list(X.index)
m = Xb.shape[1]
# paper order: by activation ratio, decreasing
act = Xb.mean(1)
order = onp.argsort(-act)
def build_pobs(gene_idx):
    """gene_idx: list of row indices in the order qubit0..qubit(n-1). Bitstring = q0 q1 ... (q0 is MSB)."""
    n = len(gene_idx)
    B = Xb[gene_idx]
    lab = onp.zeros(m, dtype=int)
    for k in range(n):
        lab = lab*2 + B[k]
    p = onp.bincount(lab, minlength=2**n).astype(float)
    p[0] = 0; p /= p.sum()
    return p

def smooth(p, alpha=1.0, m=m):
    return (p*m + alpha) / (m + p.size*alpha)

# ---------------- gates (statevector, n qubits) ----------------
def apply_1q(psi, U, k, n):
    psi = psi.reshape([2]*n)
    psi = np.moveaxis(np.tensordot(U, psi, axes=([1],[k])), 0, k)
    return psi.reshape(-1)

def apply_ctrl_1q(psi, U, c, t, n):
    """controlled-U: apply U on target t where control c == 1"""
    psi = psi.reshape([2]*n)
    # split on control axis
    idx0 = [slice(None)]*n; idx0[c] = 0
    idx1 = [slice(None)]*n; idx1[c] = 1
    p1 = psi[tuple(idx1)]                       # n-1 dim tensor
    tt = t if t < c else t-1
    p1 = np.moveaxis(np.tensordot(U, p1, axes=([1],[tt])), 0, tt)
    p1 = np.expand_dims(p1, c)
    p0 = np.expand_dims(psi[tuple(idx0)], c)
    psi = np.concatenate([p0, p1], axis=c)
    return psi.reshape(-1)

def Ry(t):  return np.array([[np.cos(t/2), -np.sin(t/2)],[np.sin(t/2), np.cos(t/2)]])
def Rx(t):  return np.array([[np.cos(t/2), -1j*np.sin(t/2)],[-1j*np.sin(t/2), np.cos(t/2)]])

def zz_signs(n):
    bits = onp.array(list(itertools.product([0,1], repeat=n)))   # row index = bitstring int, q0 MSB
    z = 1 - 2*bits
    return z  # shape 2^n x n

def pout_from_psi(psi):
    p = np.abs(psi)**2
    p = np.concatenate([np.zeros(1), p[1:]])
    return p / np.sum(p)

# ---------------- Model A: original qscGRN (sequential c-Ry) ----------------
def modelA(theta, n):
    psi = np.zeros(2**n); psi = np.concatenate([np.ones(1), psi[1:]])
    for k in range(n):
        psi = apply_1q(psi, Ry(theta[k,k]), k, n)
    for k in range(n):               # regulation layers L_0..L_{n-1}, in order
        for p in range(n):
            if p != k:
                psi = apply_ctrl_1q(psi, Ry(theta[k,p]), k, p, n)
    return pout_from_psi(psi)

# ---------------- Model B: Ising Born machine (Ry enc -> RZZ all pairs -> Rx mixer) ----------------
def modelB(theta, n, mix, Z):
    psi = np.zeros(2**n); psi = np.concatenate([np.ones(1), psi[1:]])
    for k in range(n):
        psi = apply_1q(psi, Ry(theta[k,k]), k, n)
    # diagonal phase exp(-i * sum_{k<p} theta_kp z_k z_p / 2)  (RZZ convention)
    phase = np.zeros(2**n)
    for k in range(n):
        for p in range(k+1, n):
            phase = phase + theta[k,p]*Z[:,k]*Z[:,p]
    psi = psi * np.exp(-0.5j*phase)
    for k in range(n):
        psi = apply_1q(psi, Rx(mix), k, n)
    return pout_from_psi(psi)

# ---------------- loss ----------------
def kl(pout, pobs_s):
    ps = smooth(pout)
    return np.sum(ps*np.log(ps/pobs_s))
def lcons(theta):
    return np.sum(1.0/(theta**4 - (np.pi/2)**4)**2)

def fit(model_fn, theta0, pobs, train_diag, n_iter=1500, lr=0.05, lam_dyn=True, symmetric=True, verbose=False):
    pobs_s = smooth(pobs)
    n = theta0.shape[0]
    mask = onp.ones((n,n)) - (0 if train_diag else onp.eye(n))
    def loss(theta):
        return kl(model_fn(theta), pobs_s) + lam(theta)*lcons(theta)
    lam_val = [0.0]
    def lam(theta): return lam_val[0]
    g = grad(loss)
    theta = theta0.copy()
    for s in range(n_iter):
        pout = model_fn(theta)
        L_kl = kl(pout, pobs_s); L_c = lcons(theta)
        if lam_dyn: lam_val[0] = float(L_kl / L_c)     # rescale to same order (paper's dynamic lambda)
        G = g(theta)
        if symmetric: G = (G + G.T)/2
        theta = theta - lr*G*mask
        if float(L_kl) < 2**n*1e-4: break
        if verbose and s % 200 == 0: print(s, float(L_kl))
    return theta, float(L_kl), s

def init_theta(gene_idx):
    n = len(gene_idx)
    th = onp.zeros((n,n))
    for k in range(n):
        th[k,k] = 2*onp.arcsin(onp.sqrt(act[gene_idx[k]]))
    return th

def show(theta, names, thr=onp.pi/180*2):
    df = pd.DataFrame(onp.round(theta,3), index=names, columns=names)
    return df

if __name__ == '__main__':
    n = 6
    Z = zz_signs(n)
    base = list(order)   # paper order
    names = [genes_all[i] for i in base]
    print('activation order:', names, onp.round(act[base],3))
    pobs = build_pobs(base)

    # ---- A: paper order
    thA, L, s = fit(lambda t: modelA(t, n), init_theta(base), pobs, train_diag=False)
    print(f'\n[A] paper order  KL={L:.2e} iters={s}')
    print(show(thA, names))

    # ---- A: random permutations -> report in canonical gene order
    rng = onp.random.default_rng(0)
    stack = []
    for r in range(8):
        perm = list(rng.permutation(6))
        gidx = [base[i] for i in perm]
        th, L, s = fit(lambda t: modelA(t, n), init_theta(gidx), build_pobs(gidx), train_diag=False)
        # map back to base order
        inv = onp.argsort(perm)
        stack.append(th[onp.ix_(inv, inv)])
        print(f'[A] perm {r} order={[genes_all[i] for i in gidx]} KL={L:.2e}')
    stack = onp.array(stack)
    print('\n[A] mean over 8 orderings'); print(show(stack.mean(0), names))
    print('[A] std  over 8 orderings'); print(show(stack.std(0), names))
    iu = onp.triu_indices(n,1)
    signs = onp.sign(stack[:,iu[0],iu[1]])
    print('[A] sign-consistent edges (all 8 orderings agree):', int((onp.abs(signs.sum(0))==8).sum()), '/', len(iu[0]))
    onp.save('results/A_stack.npy', stack); onp.save('results/A_paper.npy', thA)

    # ---- B: Ising Born machine, several mixer angles, diag trained vs frozen
    for mix in [onp.pi/4, onp.pi/2]:
        for td in [False, True]:
            thB, L, s = fit(lambda t: modelB(t, n, mix, Z), init_theta(base), pobs, train_diag=td, n_iter=3000)
            print(f'\n[B] mix={mix:.3f} train_diag={td}  KL={L:.2e} iters={s}')
            print(show(thB, names))
            onp.save(f'results/B_mix{mix:.2f}_td{int(td)}.npy', thB)

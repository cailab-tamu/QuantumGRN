from qscgrn_vs_ising_bm import *
n=6; Z=zz_signs(n); base=list(order); names=[genes_all[i] for i in base]; pobs=build_pobs(base)
def fitfix(model_fn, th0, pobs, td, iters=1200):
    pobs_s=smooth(pobs); mask=onp.ones((n,n))-(0 if td else onp.eye(n))
    lamv=[0.0]
    def loss(t): return kl(model_fn(t),pobs_s)+lamv[0]*lcons(t)
    g=grad(loss); th=th0.copy()
    for s in range(iters):
        Lk=float(kl(model_fn(th),pobs_s)); lamv[0]=Lk/float(lcons(th))
        G=g(th); G=(G+G.T)/2; th=th-0.05*G*mask
    return th, float(kl(model_fn(th),pobs_s))
paper = onp.array([[0,-.070,0,0,.166,-.018],[-.070,0,.078,0,-.129,.084],[0,.078,0,.043,-.101,.067],[0,0,.043,0,0,0],[.166,-.129,-.101,0,0,-.078],[-.018,.084,.067,0,-.078,0]])
thA,L=fitfix(lambda t: modelA(t,n), init_theta(base), pobs, False); print(f'[A] paper order 1200 it KL={L:.2e}'); print(show(thA,names))
rng=onp.random.default_rng(1); stack=[]
for r in range(6):
    perm=list(rng.permutation(6)); gidx=[base[i] for i in perm]
    th,L=fitfix(lambda t: modelA(t,n), init_theta(gidx), build_pobs(gidx), False); inv=onp.argsort(perm); stack.append(th[onp.ix_(inv,inv)])
    print(f'[A] perm {r} KL={L:.2e}')
stack=onp.array(stack); print('[A] std over 6 orderings (max off-diag): %.4f'%stack.std(0)[onp.triu_indices(n,1)].max())
print('[A] max |theta_ij| range across orderings:'); print(show(stack.max(0)-stack.min(0),names))
for mix in [onp.pi/6, onp.pi/4, onp.pi/3]:
    thB,L=fitfix(lambda t: modelB(t,n,mix,Z), init_theta(base), pobs, True); print(f'\n[B] mix={mix:.3f} train_diag  KL={L:.2e}'); print(show(thB,names))
    onp.save(f'results/B2_{mix:.2f}.npy',thB)
onp.save('results/A2.npy',thA)
iu=onp.triu_indices(n,1)
print('\nSign concordance with Fig.3c (nonzero entries):')
for lab,th in [('A',thA)]+[(f'B{mix:.2f}',onp.load(f'results/B2_{mix:.2f}.npy')) for mix in [onp.pi/6,onp.pi/4,onp.pi/3]]:
    nz=paper[iu]!=0; agree=(onp.sign(th[iu][nz])==onp.sign(paper[iu][nz])).sum(); print(lab, agree,'/',nz.sum())
